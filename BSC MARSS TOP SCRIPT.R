
## NOTE: Code developed with assistance from AI (ChatGPT ~5, and Claude Sonnet ~4)
# ============================================================================
#                               TOP SCRIPT
# ============================================================================

## Two hidden states:
# Male stock abundance
# Female stock abundance

## Six observation time series:
# comm_male - observes male stock at time t
# comm_female - observes female stock at time t
# surv_male_recruit - observes male stock at time t (larger juveniles, nearly fishable)
# surv_female_recruit - observes female stock at time t
# surv_prerecruit_combined - predicts male stock at time t (observed at t-2)
# surv_prerecruit_combined - predicts female stock at time t (observed at t-4)
# ** The same prerecruit series predicts both male and female stocks, but with different lags

## The prerecruits are treated as another noisy measure of true abundance. They directly
# influence state estimation but add extra variance terms (R). The alternative is the
# prerecruits enter as covariates; these affect the dynamics of abundance or catchability
# i.e. they act as drivers (deterministically), not as noisy observations.BUT this is less
# ideal because there are numerous NAs in the survey which means those time points must be 
# excluded from the model fitting (or imputed, which will be inaccurate).

# B = Process model (population dynamics)
#     Current: "diagonal and unequal" - male and female dynamics independent
#     Describes: How states evolve over time (autocorrelation, growth)
#     Alternatives to test:
#       - Full 2x2 matrix: Allow male-female interactions
#       - "diagonal and equal": Same dynamics for both sexes
#       - "identity": Random walk (no mean reversion)

# C = Covariate effects on process (environmental drivers)
#     Current: 2x4 matrix - each state affected by temp, cond, cos, sin
#     Describes: How temperature, conductivity, seasonality affect population
#     Alternatives to test:
#       - "zero": Drop covariates entirely
#       - Constrained: Same effects for male/female (fewer parameters)
#       - Partial: Only seasonality (drop temp/cond)

# Z = Observation model (how data relate to states)
#     Current: 6x2 matrix with fixed 1's for commercial, estimated for surveys
#     Describes: Catchability/selectivity of each gear for each state
#     Alternatives to test:
#       - Different series have Z=1 (change which defines scale)
#       - Equal Z for recruit surveys (z3=z4)

# A = Observation offsets/intercepts
#     Current: "zero" - no systematic bias in observations
#     Describes: Fixed additive bias in each time series
#     Alternatives to test:
#       - Estimate offsets if observations seem systematically offset

# R = Observation error covariance
#     Current: "diagonal and unequal" - independent errors, different variance
#     Describes: Measurement error in each observation series
#     Alternatives to test:
#       - "equalvarcov": Allow correlations (e.g., male/female from same survey)
#       - "diagonal and equal": Same observation error for all

# Q = Process error covariance  
#     Current: "diagonal and unequal" - male/female process errors independent
#     Describes: Environmental stochasticity affecting populations
#     Alternatives to test:
#       - "equalvarcov": Correlated (both sexes hit by same events)
#       - "diagonal and equal": Same variability for both

# c = Covariate data matrix
#     Current: 4 rows (temp, cond, cos, sin) × n_time columns
#     Describes: The actual environmental measurements

## The structure of the MARSS models was influenced by the GLM models


library(dplyr)
library(ggplot2)
library(patchwork)
library(lubridate)
library(slider)
library(MARSS)

#setwd("C:/Users/smithj08/OneDrive - DPIE/10 BSC FIS/R Code")

# ============================================================================
#                               LOAD DATA
# ============================================================================

## Load survey data (saved from prep code); number per trap
dat.surv <- as.data.frame(readRDS("survey_dat_monthly.rds"))  # 423 rows

## Load commercial catch data (prepped in my separate script); kg per trap
dat.com <- as.data.frame(readRDS("comm_dat_monthly.rds"))

## Load Wallis conductivity data
dat.cond <- readRDS("monthly_observed_Conductivity.rds")  #use 'interpolated_Cond'

## Load VR4G temperature data
dat.temp <- readRDS("monthly_VR4G_Temp.rds")

## Load nearby ocean data
dat.cwq <- readRDS("copernicus_wq_500m.rds")

## Align time series; first month-year is 11-2018 for survey
dat.com <- subset(dat.com, 
                  !(CalYear < 2018 | (CalYear == 2018 & Month < 11) |
                      CalYear > 2024 | (CalYear == 2024 & Month > 5)))  # 67 rows

## Load functions
source("BSC MARSS FUNCTIONS.R")

# ============================================================================
#                             DATA PREPARATION
# ============================================================================

## *** NOTE: I could standardize comm CPUE using vessel RE

## Create standardized CPUE metrics for survey data
dat.surv1 <- dat.surv %>%
  mutate(
    # Date handling
    date = as.Date(paste(Year, Month, "01", sep="-")),
    yearmon = paste(Year, sprintf("%02d", Month), sep="-"),
    
    # CPUE calculations
    prerecruit_F_CPUE = total_prerecruit_F / EffortDays,
    prerecruit_M_CPUE = total_prerecruit_M / EffortDays,
    recruit_F_CPUE = total_recruit_F / EffortDays,
    recruit_M_CPUE = total_recruit_M / EffortDays,
    
    # Combined metrics
    total_prerecruit = total_prerecruit_F + total_prerecruit_M,
    total_recruit = total_recruit_F + total_recruit_M,
    prerecruit_CPUE = total_prerecruit / EffortDays,
    recruit_CPUE = total_recruit / EffortDays,
    total_CPUE = (total_prerecruit + total_recruit) / EffortDays,
    
    # Sex ratios
    prerecruit_sex_ratio = total_prerecruit_F / (total_prerecruit_F + total_prerecruit_M),
    recruit_sex_ratio = total_recruit_F / (total_recruit_F + total_recruit_M),
    
    # Replace NaN with NA
    prerecruit_sex_ratio = ifelse(is.nan(prerecruit_sex_ratio), NA, prerecruit_sex_ratio),
    recruit_sex_ratio = ifelse(is.nan(recruit_sex_ratio), NA, recruit_sex_ratio) )

## Create standardized CPUE metrics for commercial data
dat.com <- dat.com %>%
  mutate(
    # Date handling
    date = as.Date(paste(CalYear, Month, "01", sep="-")),
    yearmon = paste(CalYear, sprintf("%02d", Month), sep="-"),
    
    # CPUE calculations
    Male_CPUE = Total_Male_CatchWtKg / Total_Traps,
    Female_CPUE = Total_Female_CatchWtKg / Total_Traps,
    Total_CPUE = (Total_Male_CatchWtKg + Total_Female_CatchWtKg) / Total_Traps,
    
    # Total catch
    Total_CatchWtKg = Total_Male_CatchWtKg + Total_Female_CatchWtKg,
    
    # Sex ratio
    catch_sex_ratio = Total_Female_CatchWtKg / Total_CatchWtKg )

# ============================================================================
#                         NOMINAL TEMPORAL PATTERNS
# ============================================================================

dat.surv1.plot <- dat.surv1[dat.surv1$Site %in% c(1:4),]

## Plot survey time series by site
p1 <- ggplot(dat.surv1.plot, aes(x = date, y = prerecruit_CPUE, color = factor(Site))) +
  geom_line() +
  geom_point() +
  scale_color_viridis_d(name = "Site") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Prerecruit CPUE by Site",
       x = "Date", y = "CPUE (num/trap)") +
  theme_minimal()

p2 <- ggplot(dat.surv1.plot, aes(x = date, y = recruit_CPUE, color = factor(Site))) +
  geom_line() +
  geom_point() +
  scale_color_viridis_d(name = "Site") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Recruit CPUE by Site",
       x = "Date", y = "CPUE (num/trap)") +
  theme_minimal()
p2m <- ggplot(dat.surv1.plot, aes(x = date, y = recruit_M_CPUE, color = factor(Site))) +
  geom_line() +
  geom_point() +
  scale_color_viridis_d(name = "Site") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Male Recruit CPUE by Site ",
       x = "Date", y = "CPUE (num/trap)") +
  theme_minimal()
p2f <- ggplot(dat.surv1.plot, aes(x = date, y = recruit_F_CPUE, color = factor(Site))) +
  geom_line() +
  geom_point() +
  scale_color_viridis_d(name = "Site") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Female Recruit CPUE by Site",
       x = "Date", y = "CPUE (num/trap)") +
  theme_minimal()

## Plot commercial time series
p3 <- ggplot(dat.com, aes(x = date)) +
  geom_line(aes(y = Male_CPUE, color = "Male"), linewidth = 1) +
  geom_line(aes(y = Female_CPUE, color = "Female"), linewidth = 1) +
  geom_line(aes(y = Total_CPUE, color = "Total"), linewidth = 1) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_color_manual(values = c("Male" = "blue", "Female" = "red", "Total" = "black")) +
  labs(title = "Commercial CPUE",
       x = "Date", y = "CPUE (kg/trap)", color = "Category") +
  theme_minimal()

#print(p1)
#print(p2)
#print(p3)
#p1/p2/p3
p1/p2m/p2f/p3

# ============================================================================
#                       PREPARE MODELING DATASET
# ============================================================================

## Create modeling dataset (updated function call)
model_data <- create_modeling_dataset(
  survey_data = dat.surv1, 
  commercial_data = dat.com, 
  conductivity_data = dat.cond,
  temperature_data = dat.temp,
  cwq_data = dat.cwq )

## Create site-averaged (and some lagged) variables
model_data$prerecruit_CPUE_site1234_lag2_av <- rowMeans(
  model_data[c("prerecruit_CPUE_site1_lag2", 
               "prerecruit_CPUE_site2_lag2", 
               "prerecruit_CPUE_site3_lag2", 
               "prerecruit_CPUE_site4_lag2")])
model_data$prerecruit_CPUE_site1234_lag4_av <- rowMeans(
  model_data[c("prerecruit_CPUE_site1_lag4",
               "prerecruit_CPUE_site2_lag4",
               "prerecruit_CPUE_site3_lag4",
               "prerecruit_CPUE_site4_lag4")])
model_data$recruit_F_CPUE_site1234_lag0_av <- rowMeans(
  model_data[c("recruit_F_CPUE_site1_lag0", 
               "recruit_F_CPUE_site2_lag0", 
               "recruit_F_CPUE_site3_lag0", 
               "recruit_F_CPUE_site4_lag0")])
model_data$recruit_M_CPUE_site1234_lag0_av <- rowMeans(
  model_data[c("recruit_M_CPUE_site1_lag0", 
               "recruit_M_CPUE_site2_lag0", 
               "recruit_M_CPUE_site3_lag0", 
               "recruit_M_CPUE_site4_lag0")])

## Build a per-month 6-month prior running mean from dat.cond
cond_6mo <- dat.cond %>%
  arrange(YearMonth) %>%
  mutate(
    yearmon = format(YearMonth, "%Y-%m"),
    # mean over the 6 months BEFORE the current month (exclude current)
    Cond_mean_prev6m = slide_dbl(
      interpolated_Cond,
      ~ mean(.x, na.rm = TRUE),
      .before = 6,   # look back 6 positions
      .after  = -1,  # exclude the current month
      .complete = TRUE )
  ) %>%
  dplyr::select(yearmon, Cond_mean_prev6m)

model_data <- model_data %>%
  left_join(cond_6mo, by = "yearmon")

# plot(model_data$date, model_data$Temp_VR4G, type="l",
#      las=1, ylim=c(10,24))
# lines(model_data$date, model_data$bottom_temp_lag2, col="red")
# ^ these two temps are sufficient different they can both be included

# ============================================================================
#                             BUILD FULL MODEL
# ============================================================================

## Create observation matrix with 6 time series
y <- rbind(
  comm_male = model_data$Male_CPUE,                                   # Observes male[t]
  comm_female = model_data$Female_CPUE,                               # Observes female[t]
  surv_male_recruit = model_data$recruit_M_CPUE_site1234_lag0_av,     # Observes male[t]
  surv_female_recruit = model_data$recruit_F_CPUE_site1234_lag0_av,   # Observes female[t]
  surv_prerecruit_for_male = model_data$prerecruit_CPUE_site1234_lag2_av,    # Value from t-2
  surv_prerecruit_for_female = model_data$prerecruit_CPUE_site1234_lag4_av ) # Value from t-4

# y_sqrt <- rbind(
#   comm_male = sqrt(model_data$Male_CPUE),                                   # Observes male[t]
#   comm_female = sqrt(model_data$Female_CPUE),                               # Observes female[t]
#   surv_male_recruit = sqrt(model_data$recruit_M_CPUE_site1234_lag0_av),     # Observes male[t]
#   surv_female_recruit = sqrt(model_data$recruit_F_CPUE_site1234_lag0_av),   # Observes female[t]
#   surv_prerecruit_for_male = sqrt(model_data$prerecruit_CPUE_site1234_lag2_av),    # Value from t-2
#   surv_prerecruit_for_female = sqrt(model_data$prerecruit_CPUE_site1234_lag4_av) ) # Value from t-4

## Create covariate matrix: environmental variables AND seasonal indicators
covars <- rbind(
  temp = model_data$bottom_temp_lag2,
  cond = model_data$interpolated_Cond_lag0,  # model_data$Cond_mean_prev6m; model_data$interpolated_Cond_lag0
  cos1 = model_data$month_cos,
  sin1 = model_data$month_sin )
covars[1,1] <- mean(model_data$bottom_temp_lag2[model_data$month==11], na.rm=T)
# covars[2,1:5] <- mean(model_data$Cond_mean_prev6m, na.rm=T)
# ^ replace first Temp NA with November average
# ^ *** using Cond_prev6m lower AIC and smooths the latent states more BUT the errors don't work properly
# ^ Also, using Cond_prev6m makes the forecast performance worse, especially for males

## Create Base MARSS model
# B = process model
# C = covariate effects on process
# Z = observation model structure
# A = observation offsets
# R = observation error covariance
# Q = process error covariance
# c = covariate data

model_list_base <- list(
  B = matrix(list("b11", "b12",
                  "b21", "b22"), 2, 2, byrow=T),
  # B = "diagonal and unequal",  #worse
  
  C = matrix(list(
    "c11", "c12", "c13", "c14",
    "c21", "c22", "c23", "c24"), 2, 4, byrow=T),
  
  Z = matrix(list(
    1,     0,      # comm_male directly observes male stock (fix at 1 for identifiability)
    0,     1,      # comm_female directly observes female stock
    "z3",  0,      # surv_male_recruit observes male stock (estimate catchability)
    0,     "z4",   # surv_female_recruit observes female stock
    "z5",  0,      # prerecruits lag2 predict male stock
    0,     "z6"    # prerecruits lag4 predict female stock
  ), 6, 2, byrow = TRUE),
  
  A = "zero",        # No offset in observation equation
  #A = matrix(list("a1", "a2", "a3", "a4", "a5", "a6"), 6, 1),
  
  # Observation error: allow correlation between similar measurements
  # i.e. same survey of same pool of individuals
  R = matrix(list(
    "r11",  0,     0,     0,     0,     0,      # comm_male
    0,      "r22", 0,     0,     0,     0,      # comm_female
    0,      0,     "r33", 0,     "r35", 0,      # surv_male_recruit (correlate with prerecruit_male)
    0,      0,     0,     "r44", 0,     "r46",  # surv_female_recruit
    0,      0,     "r35", 0,     "r55", 0,      # prerecruit for male
    0,      0,     0,     "r46", 0,     "r66"   # prerecruit for female
  ), 6, 6, byrow = TRUE),
  # R = "diagonal and unequal",  # *** use this to compare with commercial series model
  
  # Process error: allow male-female correlation
  # i.e. a good recruitment year, or bad environment, will affect both sexes
  Q = matrix(c("q11", "q12",
               "q12", "q22"), 2, 2, byrow = TRUE),
  # Q = "diagonal and unequal",
  
  c = covars )

# ============================================================================
#                             FIT FULL MODEL
# ============================================================================

## Use first few observations of commercial CPUE as starting point
x0_init <- matrix(c(
  mean(model_data$Male_CPUE[1:3]),      # Initial male stock
  mean(model_data$Female_CPUE[1:3])     # Initial female stock
), 2, 1)

## Fit
fit_em <- MARSS(y = y,
                model = model_list_base, 
                inits = list(x0 = x0_init),
                control = list(maxit = 2000))  #EM algorithm doesn't converge but gets close

fit <- MARSS(y = y, 
             model = model_list_base, 
             inits = fit_em,
             method = "BFGS",
             control = list(maxit = 1000))  #converges and avoids local minima


# ============================================================================
#                           EVALUATE FULL MODEL
# ============================================================================

print(fit$AIC)  # 395.7 (complex R), 391.9 (complex R, Cond_prev6m), 408.4 (simple R); -16.76 (sqrt base)

## Residuals
#residuals(fit)
plot(fit)

## Estimates
tidy(fit)

# Some specific estimates
coef(fit, type="matrix")$Z
# z3 = 2.3, z4 = 11.0, z5 = 28.1, z6 = 77.4
# If z5 > z3: Prerecruits are better predictors than recruit survey
# If z5 ≈ z3: Both survey components equally informative
# If z5 >> z6: Lag-2 relationship is stronger than lag-4
# ... confirms male and female have different maturation/recruitment timing

coef(fit, type="matrix")$B
# b11, b22: Autoregression (stability vs. fluctuation)
# b12, b21: Cross-sex effects
# b11 = 0.571, b22 = 0.63; moderate persistence between months, i.e. moderate process noise
# b12 = -0.118: Higher female abundance -> decreases male abundance next month; could reflect sex-specific sampling bias or alternating movement patterns?
# b21 = +0.253: Higher male abundance -> increases female abundance next month; female stock tends to track male trends with some lag
# these effects improve fit and reduce residual autocorrelation (ACF1 drops from 0.20 to 0.19)

coef(fit, type="matrix")$Q  #process error
# q11, q22; small year-to-year random variability (< 0.003)
# q12: small shared process noise (-0.009)

coef(fit, type="matrix")$R
# diagonal elements are the residual variance; off-diagonal elements are the residual covariances
# diagonal values quantify how noisy or reliable each observation series is after accounting for the latent state trend
# ideally the off-diagonals would be close to zero - you don't want patterns in the residuals
# surv_male_recruit & surv_prerecruit_for_male = -0.094; strong negative covariance; the prerecruit and recruit indices for males fluctuate inversely after accounting for the state model
# surv_female_recruit & surv_prerecruit_for_female = -0.73; very strong neg covariance ***; uncertain cause
# diagonal elements show survey indices have much larger residual variance than commercial CPUE, consistent with their higher sampling noise and different scales
# the negative residual covariances could indicate the lags are imperfect or inconsistent

# "Because indices differ in units and scale, we interpret R on the correlation scale and report residual SDs relative to the series mean":
Rhat   <- coef(fit, type = "matrix")$R
sd_i   <- sqrt(diag(Rhat))
Dinv   <- diag(1/sd_i)
Rcorr  <- Dinv %*% Rhat %*% Dinv
rel_sd <- sd_i / rowMeans(y, na.rm = TRUE)  # residual SD relative to mean level
rel_sd

## Variance explained by state (per series)
fit_df <- fitted(fit)           # columns: .rownames, t, .fitted
y_vecs <- split(as.numeric(y), rownames(y))
f_vecs <- split(fit_df$.fitted, fit_df$.rownames)
VE <- sapply(names(y_vecs), function(k) {
  yy <- y_vecs[[k]]; ff <- f_vecs[[k]][seq_along(yy)]
  idx <- is.finite(yy) & is.finite(ff)
  1 - var(yy[idx] - ff[idx], na.rm = TRUE) / var(yy[idx], na.rm = TRUE)
})
VE

## Estimate the hidden states
states <- tsSmooth(fit)

## Compare with commercial CPUE (most likely indicator of stock abundance)
cor(states$.estimate[states$.rownames=="X1"], model_data$Male_CPUE, use="complete.obs")  # 0.78
cor(states$.estimate[states$.rownames=="X2"], model_data$Female_CPUE, use="complete.obs")  # 0.70
# ^^ For male CPUE: strong correlation. The hidden state closely tracks the commercial
# CPUE for males, meaning the MARSS model attributes most male CPUE variation to
# genuine stock fluctuations, not noise.
# ^^ For female CPUE: moderately strong correlation. The latent female stock signal
# tracks CPUE, but there’s more deviation — consistent with catchability or behaviour
# differences (e.g., seasonal availability, selectivity).
# The smoother state trajectories probably filters out much of the observation process noise.
# "The estimated latent states for males and females were strongly correlated with
# their respective commercial CPUE series (r = 0.78 and 0.70), indicating that the
# model successfully extracted an underlying abundance signal consistent with
# observed catches but filtered of short-term variability and observation noise."

## Plot latent states vs commercial catches

par(mfrow=c(2,1))
# Plot male stock trajectory
plot((states$.estimate[states$.rownames=="X1"]), type="l", lwd=2, 
     ylab="Male Stock (latent)", main="Estimated Male Stock State",
     ylim=c(0,0.6), las = 1, xaxt = "n", xlab="", las=1)
axis(1, at=1:length(model_data$yearmon), labels=model_data$yearmon, las=2)
lines((states$.estimate[states$.rownames=="X1"]) + 2*(states$.se[states$.rownames=="X1"]), lty=2)
lines((states$.estimate[states$.rownames=="X1"]) - 2*(states$.se[states$.rownames=="X1"]), lty=2)
points(model_data$Male_CPUE, col="red", pch=16)

# Plot female stock trajectory
plot((states$.estimate[states$.rownames=="X2"]), type="l", lwd=2, 
     ylab="Female Stock (latent)", main="Estimated Female Stock State",
     ylim=c(0,0.2), las = 1, xaxt = "n", xlab="", las=1)
axis(1, at=1:length(model_data$yearmon), labels=model_data$yearmon, las=2)
lines((states$.estimate[states$.rownames=="X2"]) + 2*(states$.se[states$.rownames=="X2"]), lty=2)
lines((states$.estimate[states$.rownames=="X2"]) - 2*(states$.se[states$.rownames=="X2"]), lty=2)
points(model_data$Female_CPUE, col="blue", pch=16)

## Plot fitted values for the n time series
fitted_vals <- fitted(fit)
obs_names <- c("comm_male", "comm_female", "surv_male_recruit", 
               "surv_female_recruit", "surv_prerecruit_for_male", 
               "surv_prerecruit_for_female")

par(mfrow = c(3,2), mar = c(4,4,2,1))
for (obs in obs_names) {
  subset_fit <- fitted_vals[fitted_vals$.rownames == obs, ]
  plot(subset_fit$t, (subset_fit$y), 
       pch = 16, col = "black", las = 1,
       xlab = "Time", ylab = "CPUE", xaxt = "n",
       main = obs, ylim=c(0,max(subset_fit$y, na.rm=T)))
  axis(1, at=1:length(model_data$yearmon), labels=model_data$yearmon, las=2)
  lines(subset_fit$t, (subset_fit$.fitted), col = "red", lwd = 2)
}

## Diagnostics of fit for six series
diagnostics <- compute_diagnostics(fit, "Base Model")
print(diagnostics)
# *** report MAPE, R2 and ACF1


# ============================================================================
#                         ANNUAL INDEX - RAW STATES
# ============================================================================

# --- RAW smoothed states ---
states <- tsSmooth(fit)  # data.frame: .rownames, t, .estimate, .se

# Keep X1/X2 and build tidy frame with labels
st_df <- subset(states, .rownames %in% c("X1","X2"))
st_df <- transform(
  st_df,
  sex      = ifelse(.rownames == "X1", "Male", "Female"),
  year     = model_data$year,     # assume aligned indexing
  estimate = .estimate,
  se       = .se )

# --- Annual summaries (mean and ~SE of annual mean) ---
# SE(mean) ≈ sqrt( sum(month_se^2) ) / n   (treating month estimates as approx. independent)
ann <- st_df %>%
  group_by(sex, year) %>%
  summarise(
    n      = sum(!is.na(estimate)),
    rho1   = if (n > 2) cor(estimate[-n], estimate[-1], use="complete.obs") else 0,
    rho1   = pmin(pmax(rho1, -0.9), 0.9),   # cap extreme rho values
    n_eff  = pmax(2, n * (1 - rho1) / (1 + rho1)),  # minimum of 2 effective obs
    mean   = mean(estimate, na.rm=TRUE),
    se_bar = sqrt(sum(se^2, na.rm=TRUE)) / n_eff,
    lo95   = mean - 1.96*se_bar,
    hi95   = mean + 1.96*se_bar,
    x_mid  = mean(t, na.rm=TRUE),
    .groups = "drop" )

# --- Plot: raw states + annual means/CIs overlaid ---
par(mfrow = c(2,1))

## MALE
ixM <- st_df$sex == "Male"
plot(st_df$estimate[ixM], type="l", lwd=2, col="grey40",
     ylab="Male Stock (latent)", main="Estimated Male Stock State — RAW (smoothed)",
     ylim=c(0, 0.6), xaxt="n", xlab="", las=1)
axis(1, at = 1:length(model_data$yearmon), labels = model_data$yearmon, las = 2, cex.axis=0.7)

lines(st_df$estimate[ixM] + 2*st_df$se[ixM], lty=2, col="grey70")
lines(st_df$estimate[ixM] - 2*st_df$se[ixM], lty=2, col="grey70")
points(model_data$Male_CPUE, col="grey70", pch=16)

annM <- subset(ann, sex == "Male" & year %in% 2019:2024)
with(annM, {
  arrows(x_mid, lo95, x_mid, hi95, angle=90, code=3, length=0.05, col="red")
  points(x_mid, mean, pch=19, cex=1.4, col="red")
  lines(x_mid, mean, col="red", lwd=2)
})

## FEMALE
ixF <- st_df$sex == "Female"
plot(st_df$estimate[ixF], type="l", lwd=2, col="grey40",
     ylab="Female Stock (latent)", main="Estimated Female Stock State — RAW (smoothed)",
     ylim=c(0, 0.2), xaxt="n", xlab="", las=1)
axis(1, at = 1:length(model_data$yearmon), labels = model_data$yearmon, las = 2, cex.axis=0.7)

lines(st_df$estimate[ixF] + 2*st_df$se[ixF], lty=2, col="grey70")
lines(st_df$estimate[ixF] - 2*st_df$se[ixF], lty=2, col="grey70")
points(model_data$Female_CPUE, col="grey70", pch=16)

annF <- subset(ann, sex == "Female" & year %in% 2019:2024)
with(annF, {
  arrows(x_mid, lo95, x_mid, hi95, angle=90, code=3, length=0.05, col="blue")
  points(x_mid, mean, pch=19, cex=1.4, col="blue")
  lines(x_mid, mean, col="blue", lwd=2)
})


# ============================================================================
#                      ANNUAL INDEX - SEASON CONTROLLED
# ============================================================================

## *** MARSS predict is not like normal predict
## This must be interpreted carefully. The monthly scale isn't useful
## BUT the decomposition does contain real signal about abundance

# covariate matrix with season turned off
c_cf <- covars
if ("cos1" %in% rownames(c_cf)) c_cf["cos1", ] <- 0
if ("sin1" %in% rownames(c_cf)) c_cf["sin1", ] <- 0

# produce season-neutral states (xtT) for the whole series
pr_cf <- predict(fit, type = "xtT", interval = "confidence",
                 newdata = list(c = c_cf))$pred
idx_states <- pr_cf$.rownames %in% c("X1", "X2")
state_cf <- pr_cf[idx_states, ]

##Now aggregate to annual
# ---- Build a tidy data frame with year labels ----
st_df <- state_cf %>%
  mutate(sex  = ifelse(.rownames=="X1", "Male", "Female"))

st_df$year <- rep(model_data$year,2)
  
# Annual summaries (mean and approx SE of the annual mean)
ann <- st_df %>%
  group_by(sex, year) %>%
  summarise(
    n      = sum(!is.na(estimate)),
    mean   = mean(estimate, na.rm=TRUE),
    se_bar = sqrt(sum(se^2, na.rm=TRUE)) / n,     # ≈ SE of annual mean
    lo95   = mean - 1.96*se_bar,
    hi95   = mean + 1.96*se_bar,
    x_mid  = mean(t, na.rm=TRUE),                 # x-position (mid of year)
    .groups = "drop" )

# ---- Plot, then overlay annual means + CI ----
par(mfrow=c(2,1))

## MALE
ixM <- st_df$sex=="Male"
plot(st_df$estimate[ixM], type="l", lwd=2,
     ylab="Male Stock (latent)", main="Estimated Male Stock State — season controlled",
     ylim=c(0,0.6), xaxt="n", xlab="", las=1, col="grey40")
axis(1, at=1:length(model_data$yearmon), labels=model_data$yearmon, las=2)

lines(st_df$estimate[ixM] + 2*st_df$se[ixM], lty=2, col="grey")
lines(st_df$estimate[ixM] - 2*st_df$se[ixM], lty=2, col="grey")
points(model_data$Male_CPUE, col="grey", pch=16)

annM <- subset(ann, sex=="Male")
annM <- annM[annM$year %in% 2019:2024,]
arrows(annM$x_mid, annM$lo95, annM$x_mid, annM$hi95,
       angle=90, code=3, length=0.05, col="red")
points(annM$x_mid, annM$mean, pch=19, cex=1.4, col="red")
lines(annM$x_mid, annM$mean, col="red", lwd=2)

## FEMALE
ixF <- st_df$sex=="Female"
plot(st_df$estimate[ixF], type="l", lwd=2,
     ylab="Female Stock (latent)", main="Estimated Female Stock State — season controlled",
     ylim=c(0,0.2), xaxt="n", xlab="", las=1, col="grey40")
axis(1, at=1:length(model_data$yearmon), labels=model_data$yearmon, las=2)

lines(st_df$estimate[ixF] + 2*st_df$se[ixF], lty=2, col="grey")
lines(st_df$estimate[ixF] - 2*st_df$se[ixF], lty=2, col="grey")
points(model_data$Female_CPUE, col="grey", pch=16)

annF <- subset(ann, sex=="Female")
annF <- annF[annF$year %in% 2019:2024,]
arrows(annF$x_mid, annF$lo95, annF$x_mid, annF$hi95,
       angle=90, code=3, length=0.05, col="blue")
points(annF$x_mid, annF$mean, pch=19, cex=1.4, col="blue")
lines(annF$x_mid, annF$mean, col="blue", lwd=2)


# ## How much variation is seasonal
# Cmat  <- coef(fit, type="matrix")$C
# C_sea <- Cmat[, match(c("cos1","sin1"), rownames(covars)), drop=FALSE]
# sea   <- t(C_sea %*% covars[c("cos1","sin1"), , drop=FALSE])  # T x 2
# sm <- tsSmooth(fit)
# states_mat <- matrix(NA, nrow = sum(sm$.rownames == "X1"), ncol = 2)
# states_mat[, 1] <- sm$.estimate[sm$.rownames == "X1"]  # male
# states_mat[, 2] <- sm$.estimate[sm$.rownames == "X2"]  # female
# colnames(states_mat) <- c("X1", "X2")
# # Compute proportion of variance explained by the seasonal component
# var_expl <- colMeans(sea^2) / colMeans(states_mat^2)
# var_expl  # 0.019, 0.042
# # ^ not sure if correct, but may show only 2-4% of latent variance is explained by the deterministic seasonal component in the state equation

## Correlation between states for full model and deseasoned model
x_rawM <- states$.estimate[states$.rownames=="X1"]
x_rawF <- states$.estimate[states$.rownames=="X2"]
cor(x_rawM, x_rawF)  # raw = 0.043
cor(x_no_season[,1], x_no_season[,2])  # deseasoned = 0.308
# massive increase in correlation
# -> most of the difference between the sexes in the raw states was phase-shifted seasonality
# -> consistent with both sexes experiencing the same underlying population trajectory
# and process noise (as the B and Q matrices allow)


# ============================================================================
#                           FORECAST FULL MODEL
# ============================================================================

## Forecast the latent states (M and F abundance)
# These states should track the true stock trajectory independent of observation
# noise; this is a measure of biological forecasting —> “stock status”. The
# latter's accuracy can't really be measured, but can be evaluated against:
# - correlation with commercial CPUE
# - precision of state estimates (SE); mean forecast variance

# t+1, t+2, t+3
# Compare to the GLM?

# -----------------------------
# Settings
# -----------------------------
H <- 6  # forecast horizon (months)
stopifnot(ncol(y) == ncol(covars))
Tn <- ncol(y)
Ttrain <- Tn - H

# -----------------------------
# Training data (fit to t = 1:(T-3))
# -----------------------------
y_tr <- y[, 1:Ttrain, drop = FALSE]
c_tr <- covars[, 1:Ttrain, drop = FALSE]
#d_tr <- d_obs[, 1:Ttrain, drop = FALSE]

# -----------------------------
# Fit Base model (EM -> BFGS)
# -----------------------------
x0_init <- matrix(c(
  mean(model_data$Male_CPUE[1:3],   na.rm = TRUE),
  mean(model_data$Female_CPUE[1:3], na.rm = TRUE)
), 2, 1)

model_base_tr <- model_list_base
#model_base_tr <- model_obsTemp
model_base_tr$c <- c_tr
#model_base_tr$d <- d_tr

fit_em <- MARSS(y = y_tr,
                model   = model_base_tr,
                inits   = list(x0 = x0_init),
                control = list(maxit = 2000),
                silent  = TRUE)

fit_tr <- MARSS(y = y_tr,
                model   = model_base_tr,
                inits   = fit_em,
                method  = "BFGS",
                control = list(maxit = 1000),
                silent  = TRUE)

# -----------------------------
# Build future covariates (three scenarios)
# -----------------------------
# 1) Observed: use covars[, (Ttrain+1):Tn]
c_future_obs <- covars[, (Ttrain+1):Tn, drop = FALSE]

# 2) Zero: all zeros
c_future_zero <- matrix(0, nrow = nrow(c_tr), ncol = H,
                        dimnames = list(rownames(c_tr), NULL))

# 3) Carry last: repeat last observed training vector H times
last_c <- c_tr[, ncol(c_tr), drop = FALSE]
c_future_carry <- matrix(rep(drop(last_c), H), nrow = nrow(c_tr), ncol = H,
                         dimnames = list(rownames(c_tr), NULL))

newdata_list <- list(
  obs   = list(c = c_future_obs),
  zero  = list(c = c_future_zero),
  carry = list(c = c_future_carry) )

# -----------------------------
# Forecast latent states (X1, X2) for h = 1..3 for each scenario
# -----------------------------
get_state_fc <- function(fit_obj, new_c) {
  pr <- predict(fit_obj, n.ahead = H, interval = "confidence",
                type = "xtT", newdata = new_c)$pred
  pr <- subset(pr, .rownames %in% c("X1","X2"))
  # Ensure horizon column h exists
  if (!"h" %in% names(pr)) {
    pr$h <- ave(pr$t, pr$.rownames, FUN = function(tt) seq_along(tt))
  }
  pr
}

fc_list <- lapply(newdata_list, function(nd) get_state_fc(fit_tr, nd))
# fc_list is a named list: $obs, $zero, $carry

# Keep ONLY the final forecast block (origin at end of training) ----
fc_final <- lapply(fc_list, function(df) {
  if (!all(c(".rownames","t","h","estimate") %in% names(df))) return(df[0, ])
  # Keep only X1/X2 and the requested horizons
  df <- subset(df, .rownames %in% c("X1","X2") & h %in% (Ttrain+1):(Ttrain+H))
  
  # Some MARSS builds return absolute t (Ttrain+1..Ttrain+H), others return relative t (1..H).
  # Normalise to absolute forecast index for plotting/summaries:
  if (length(df$t) && max(df$t, na.rm = TRUE) <= Ttrain) {
    df$t_abs <- Ttrain + df$h
  } else {
    df$t_abs <- df$t
  }
  df
})

# -----------------------------
# Smoothed states on training window (for plotting)
# -----------------------------
sm <- tsSmooth(fit_tr)
x1_tr <- sm$.estimate[sm$.rownames == "X1"]
x2_tr <- sm$.estimate[sm$.rownames == "X2"]
# ^ superfluous bc these are also in predict...

# -----------------------------
# Plot: 2 (sex) × 3 (scenario) plots
# -----------------------------
# scenarios <- c("obs","zero","carry")
scenarios <- c("obs","carry")  # exclude "zero" for paper
sexes <- c("X1","X2")
sex_lbl <- c(X1 = "Male latent state", X2 = "Female latent state")

# set the starting index for x-axis (e.g., 20 means start from the 20th month)
x_start <- 37
op <- par(mfrow = c(2, 2), mar = c(4,4,2,1))
for (sx in sexes) {
  for (sc in scenarios) {
    pf <- fc_list[[sc]]
    pf <- pf[pf$.rownames == sx, , drop = FALSE]
    pf <- pf[order(pf$t), , drop = FALSE]
    idx <- pf$t
    
    # training state
    x_tr <- if (sx == "X1") x1_tr else x2_tr
    
    # limit the time range for plotting
    idx_plot <- idx[idx >= x_start]
    x_tr_plot <- x_tr[x_start:length(x_tr)]
    
    # y-lims include smoothed training + forecast means and bands
    y_min <- min(x_tr_plot,
                 pf$estimate[pf$t >= x_start],
                 if ("Lo 95" %in% names(pf)) pf[["Lo 95"]][pf$t >= x_start] else pf$estimate[pf$t >= x_start],
                 na.rm = TRUE)
    y_max <- max(x_tr_plot,
                 pf$estimate[pf$t >= x_start],
                 if ("Hi 95" %in% names(pf)) pf[["Hi 95"]][pf$t >= x_start] else pf$estimate[pf$t >= x_start],
                 na.rm = TRUE) * 1.1
    if (sx == "X1") {y_max <- 0.45} else {y_max <- 0.17}
    
    plot(c(x_start, max(idx)), c(y_min, y_max), type = "n", las = 1,
         xlab = "Time", ylab = "Latent state", xaxt = "n",
         main = paste0(sex_lbl[[sx]], " — ", sc, " covars"))
    
    # custom x-axis labels from the subset
    axis(1, at = x_start:length(model_data$yearmon),
         labels = model_data$yearmon[x_start:length(model_data$yearmon)], las = 2)
    
    # 95% PI ribbon (if present)
    if (all(c("Lo 95","Hi 95") %in% names(pf))) {
      sel <- pf$t >= x_start
      polygon(c(pf$t[sel], rev(pf$t[sel])),
              c(pf[["Lo 95"]][sel], rev(pf[["Hi 95"]][sel])),
              border = NA, col = grDevices::adjustcolor("gray70", 0.5))
    }
    
    # 80% PI ribbon (if present)
    if (all(c("Lo 80","Hi 80") %in% names(pf))) {
      sel <- pf$t >= x_start
      polygon(c(pf$t[sel], rev(pf$t[sel])),
              c(pf[["Lo 80"]][sel], rev(pf[["Hi 80"]][sel])),
              border = NA, col = grDevices::adjustcolor("gray50", 0.5))
    }
    
    # smoothed training line (subset from x_start)
    lines(x_start:length(x_tr), x_tr_plot, lwd = 2)
    
    abline(v = Ttrain, lty = 3)
    
    # forecast mean
    lines(pf$t[pf$t >= Ttrain], pf$estimate[pf$t >= Ttrain], lwd = 2, lty = 2)
    points(pf$t[pf$t > Ttrain], pf$estimate[pf$t > Ttrain])
    
    # legend only on top-left panel
    if (sx == "X1" && sc == "obs") {
      legend("topleft",
             c("Smoothed (train)", "Forecast mean", "80% PI", "95% PI"),
             lty = c(1,2,NA,NA), lwd = c(2,2,NA,NA),
             pch = c(NA, NA, 15, 15), pt.cex = c(NA, NA, 2, 2),
             col = c("black","black","gray50","gray70"),
             bty = "n")
    }
  }
}
par(op)

# -----------------------------
# Precision summary: mean 95% width by horizon and scenario
# -----------------------------

# Build precision table across scenarios in fc_final
prec_tab <- do.call(rbind, lapply(names(fc_final), function(nm) {
  tmp <- prec_one(fc_final[[nm]], Ttrain = Ttrain, H = H, want80 = FALSE)
  if (nrow(tmp)) tmp$scenario <- nm
  tmp
}))

# Tidy labels & order
prec_tab$sex <- ifelse(prec_tab$state == "X1", "Male", "Female")
  prec_tab <- prec_tab[order(prec_tab$scenario, prec_tab$sex, prec_tab$h),
                       c("scenario","sex","h","mean_w95")]

print(prec_tab)
# The 95% interval widths (mean_w95) increase with forecast horizon (1 -> 3).
# Female: 0.079 -> 0.093 (3m) -> 0.11 (6m)
# Male: 0.229 -> 0.267 (3m) -> 0.270 (6m)
# The values are identical across obs, zero, and carry scenarios; this is
# expected: forecast precision is governed by the state-space dynamics,
# not by covariates
# The precision stabilises after ~3-6 months
# Credible intervals are similar when considering scale of predictions

meanM <- mean(sm$.estimate[sm$.rownames=="X1"])
meanF <- mean(sm$.estimate[sm$.rownames=="X2"])
meanM_prec <- mean(prec_tab$mean_w95[prec_tab$scenario=="obs" & prec_tab$sex=="Male"])
meanF_prec <- mean(prec_tab$mean_w95[prec_tab$scenario=="obs" & prec_tab$sex=="Female"])

meanM_prec/meanM  #3-m = 1.39; 6-m = 1.54
meanF_prec/meanF  #3-m = 1.21; 6-m = 1.38
# For 3-month window, the forecast uncertainty (95% interval width) is a little
# larger than the mean latent abundance itself. This means the model captures
# the trend direction of the latent abundance well, but absolute values are uncertain
# — expected when using short CPUE-based time series to estimate process states.
# The female series is slightly more precise (narrower relative intervals), which
# likely reflects lower process variance and/or stronger observation signals.
# The increase in width with forecast horizon (from h = 1 to h = 3) indicates
# modest uncertainty growth.
# "Although the state estimates are noisy in absolute terms, the relative growth
# in forecast uncertainty over three months was moderate, suggesting stable
# short-term predictions"


# ============================================================================
#                     VALUE OF THE FIS - FORECASTING TEST
# ============================================================================

## Forecast the commercial CPUE index, and compare to actual CPUE (M and F),
# AND compare these to a commercial-only model to examine the value of the 
# independent survey to predictability and stock status

## The forecasting of a noisy observation (the commercial CPUE series) tells
# us how well the model predicts what fishers will catch next month; this
# is a measure of operational forecasting —> “fishing quality”

# ---------------------------
# SETTINGS
# ---------------------------
SERIES_TO_PLOT <- c("comm_male","comm_female")  # data to forecast and plot
N_LAST <- 14           # forecasting start points (i.e. how many results to average)
HORIZON <- 3           # 1, 2, or 3 (with h=1 we can forecast to s=67, h=3 s=65)

## Commercial-only observation matrix and model (simple, matches your earlier setup)
y_comm <- rbind(
  comm_male   = model_data$Male_CPUE,
  comm_female = model_data$Female_CPUE )

model_comm_only <- list(
  B = model_list_base$B,
  Z = matrix(list(1,0,0,1), 2, 2, byrow=TRUE),
  A = "zero",
  R = "diagonal and unequal",
  Q = model_list_base$Q,
  C = model_list_base$C,
  c = covars )

## Edit the full model to allow a consistent comparison
model_list_base$R <- "diagonal and unequal"

# ---------------------------
# TARGETS
# ---------------------------
Tn <- ncol(y)  # full model y
N_LAST <- min(N_LAST, Tn - HORIZON)  # need space to forecast out to t+H-1
targets <- (Tn - N_LAST + 1):(Tn - (HORIZON - 1))

# ---------------------------
# STORAGE
# ---------------------------
pred_rows_full <- list()
pred_rows_comm <- list()

# --------- PLOTTING LAYOUT ---------
op <- par(mfrow=c(3,2),
          mar=c(3,3,2,1), oma=c(0,0,2,0))

for (s in targets) {
  print(paste("Time start = ",s))
  
  # ---------------------------
  # TRAINING DATA up to t-1
  # ---------------------------
  y_tr <- y[, 1:(s-1), drop=FALSE]
  c_tr <- covars[, 1:(s-1), drop=FALSE]
  
  # Future covariates for t..t+H-1 (needed for observation forecasts y_{t+h})
  c_fc <- covars[, s:(s + HORIZON - 1), drop=FALSE]
  
  # --------- FULL MODEL: EM -> BFGS ----------
  x0_init <- matrix(c(
    mean(model_data$Male_CPUE[1:3], na.rm=TRUE),
    mean(model_data$Female_CPUE[1:3], na.rm=TRUE)
  ), 2, 1)
  
  model_full_tr <- model_list_base
  model_full_tr$c <- c_tr
  
  print("fitting full")
  fit_em_full <- MARSS(y = y_tr,
                       model = model_full_tr,
                       inits = list(x0 = x0_init),
                       control = list(maxit = 2000), silent=TRUE)
  
  fit_full <- MARSS(y = y_tr,
                    model = model_full_tr,
                    inits = fit_em_full,
                    method = "BFGS",
                    control = list(maxit = 1000), silent=TRUE)
  
  # --------- COMM-ONLY: EM -> BFGS ----------
  model_comm_tr <- model_comm_only
  model_comm_tr$c <- c_tr
  
  print("fitting commercial")
  fit_em_comm <- MARSS(y = y_comm[, 1:(s-1), drop=FALSE],
                       model = model_comm_tr,
                       inits = list(x0 = x0_init),
                       control = list(maxit = 2000), silent=TRUE)
  
  # patch any zero/neg R diag in the EM par vector (2x1), to keep BFGS happy
  Rpar <- fit_em_comm$par$R  # dim = 2 x 1 with rows "(comm_male,comm_male)", "(comm_female,comm_female)"
  Rpar[1, 1] <- ifelse(!is.finite(Rpar[1,1]) || Rpar[1,1] <= 0, 1e-6, Rpar[1,1])  # bump zero/tiny
  Rpar[2, 1] <- ifelse(!is.finite(Rpar[2,1]) || Rpar[2,1] <= 0, 1e-6, Rpar[2,1])
  fit_em_comm$par$R <- Rpar
  
  fit_comm <- MARSS(y = y_comm[, 1:(s-1), drop=FALSE],
                    model = model_comm_tr,
                    inits = fit_em_comm,
                    method = "BFGS",
                    control = list(maxit = 1000), silent=TRUE)
  
  # ----------------------------
  # MULTI-STEP-AHEAD PREDICTIONS (t -> t+n)
  # ----------------------------
  new_c <- list(c = c_fc)
  
  pr_full <- predict(object = fit_full,
                     n.ahead = HORIZON,
                     interval = "prediction",
                     newdata = new_c)$pred
  pr_comm <- predict(object = fit_comm,
                     n.ahead = HORIZON,
                     interval = "prediction",
                     newdata = new_c)$pred
  
  # keep horizons & series of interest and add horizon column (h = 1..H)
  keep_f <- subset(pr_full, .rownames %in% SERIES_TO_PLOT & t >= s & t <= s + HORIZON - 1)
  keep_c <- subset(pr_comm, .rownames %in% SERIES_TO_PLOT & t >= s & t <= s + HORIZON - 1)
  keep_f$h <- keep_f$t - s + 1
  keep_c$h <- keep_c$t - s + 1
  keep_f$origin <- s
  keep_c$origin <- s
  pred_rows_full[[length(pred_rows_full)+1]] <- keep_f
  pred_rows_comm[[length(pred_rows_comm)+1]] <- keep_c
  
  # ---------------------------
  # FITTED LINES (to t-1) for plotting
  # ---------------------------
  fit_full_df <- fitted(fit_full)
  fit_comm_df <- fitted(fit_comm)
  
  # ---------- PLOT: full history up to s+H-1 (if available) ----------
  idx_plot_max <- min(Tn, s + HORIZON - 1)
  idx <- 1:idx_plot_max
  
  # ----- MALE panel -----
  {
    series <- "comm_male"
    y_all  <- model_data$Male_CPUE[idx]
    f_full <- subset(fit_full_df, .rownames == series & t <= (s-1))
    f_comm <- subset(fit_comm_df, .rownames == series & t <= (s-1))
    pf     <- subset(keep_f, .rownames == series) # h=1..H
    pc     <- subset(keep_c, .rownames == series)
    
    y_min <- min(y_all, f_full$.fitted, f_comm$.fitted, pf$estimate, pc$estimate, na.rm=TRUE)
    y_max <- max(y_all, f_full$.fitted, f_comm$.fitted, pf$estimate, pc$estimate, na.rm=TRUE) * 1.1
    
    plot(idx, y_all, pch=16, col="black", las=1,
         xlab="", ylab="CPUE", xaxt="n",
         ylim=range(y_min, y_max, finite=TRUE),
         xlim=c(1, Tn),
         main=paste0(series, " | forecasts @ t=", s, " (h=1..", HORIZON, ")"))
    #axis(1, at=idx, labels=idx)
    axis(1, at=1:length(model_data$yearmon), labels=model_data$yearmon, las=2)
    
    if (nrow(f_full)>0) lines(f_full$t, f_full$.fitted, lwd=2)
    if (nrow(f_comm)>0) lines(f_comm$t, f_comm$.fitted, lwd=2, lty=2)
    
    # forecast points for h=1..H (full=solid triangle; comm=triangle open)
    if (nrow(pf)>0) points(pf$t, pf$estimate, pch=17, col="red")
    if (nrow(pc)>0) points(pc$t, pc$estimate, pch=2,  col="red")
    
    # actual observations at those horizons
    have_obs <- pf$t[pf$t <= Tn]
    if (length(have_obs))
      points(have_obs, model_data$Male_CPUE[have_obs], pch=16, col="blue")
    
    # connect fitted @ s-1 to first forecast points
    if (nrow(f_full)>0 && any(pf$h==1))
      lines((s-1):(s + HORIZON - 1),
            c(f_full$.fitted[f_full$t==(s-1)], pf$estimate[pf$h %in% 1:HORIZON]),
            col="red", lwd=2)
    if (nrow(f_comm)>0 && any(pc$h==1))
      lines((s-1):(s + HORIZON - 1),
            c(f_comm$.fitted[f_comm$t==(s-1)], pc$estimate[pc$h %in% 1:HORIZON]),
            col="red", lwd=2, lty=2)
    
    abline(v = s - 1, lty = 3)
    
    legend("topleft",
           c("Observed", "Fitted full", "Fitted comm-only",
             "Pred full (h=1..H)", "Pred comm-only (h=1..H)", "Actual at h"),
           pch=c(16, NA, NA, 17, 2, 16), lty=c(NA,1,2,NA,NA,NA),
           col=c("black","black","black","red","red","blue"),
           lwd=c(NA,2,2,NA,NA,NA), bty="n", cex=0.8)
  }
  
  # ----- FEMALE panel -----
  {
    series <- "comm_female"
    y_all  <- model_data$Female_CPUE[idx]
    f_full <- subset(fit_full_df, .rownames == series & t <= (s-1))
    f_comm <- subset(fit_comm_df, .rownames == series & t <= (s-1))
    pf     <- subset(keep_f, .rownames == series)
    pc     <- subset(keep_c, .rownames == series)
    
    y_min <- min(y_all, f_full$.fitted, f_comm$.fitted, pf$estimate, pc$estimate, na.rm=TRUE)
    y_max <- max(y_all, f_full$.fitted, f_comm$.fitted, pf$estimate, pc$estimate, na.rm=TRUE) * 1.1
    
    plot(idx, y_all, pch=16, col="black", las=1,
         xlab="", ylab="CPUE", xaxt="n",
         ylim=range(y_min, y_max, finite=TRUE),
         xlim=c(1, Tn),
         main=paste0(series, " | forecasts @ t=", s, " (h=1..", HORIZON, ")"))
    #axis(1, at=idx, labels=idx)
    axis(1, at=1:length(model_data$yearmon), labels=model_data$yearmon, las=2)
    
    if (nrow(f_full)>0) lines(f_full$t, f_full$.fitted, lwd=2)
    if (nrow(f_comm)>0) lines(f_comm$t, f_comm$.fitted, lwd=2, lty=2)
    
    if (nrow(pf)>0) points(pf$t, pf$estimate, pch=17, col="red")
    if (nrow(pc)>0) points(pc$t, pc$estimate, pch=2,  col="red")
    
    have_obs <- pf$t[pf$t <= Tn]
    if (length(have_obs))
      points(have_obs, model_data$Female_CPUE[have_obs], pch=16, col="blue")
    
    # connect fitted @ s-1 to first forecast points
    if (nrow(f_full)>0 && any(pf$h==1))
      lines((s-1):(s + HORIZON - 1),
            c(f_full$.fitted[f_full$t==(s-1)], pf$estimate[pf$h %in% 1:HORIZON]),
            col="red", lwd=2)
    if (nrow(f_comm)>0 && any(pc$h==1))
      lines((s-1):(s + HORIZON - 1),
            c(f_comm$.fitted[f_comm$t==(s-1)], pc$estimate[pc$h %in% 1:HORIZON]),
            col="red", lwd=2, lty=2)
    
    abline(v = s - 1, lty = 3)
    
    legend("topleft",
           c("Observed", "Fitted full", "Fitted comm-only",
             "Pred full (h=1..H)", "Pred comm-only (h=1..H)", "Actual at h"),
           pch=c(16, NA, NA, 17, 2, 16), lty=c(NA,1,2,NA,NA,NA),
           col=c("black","black","black","red","red","blue"),
           lwd=c(NA,2,2,NA,NA,NA), bty="n", cex=0.8)
  }
  
}

mtext(sprintf("Rolling multi-step (h=1..%d): Full vs Commercial-only", HORIZON), outer=TRUE, line=0.5, cex=1)
par(op)

# ---------------------------
# ACCURACY BY HORIZON
# ---------------------------
pred_full_all <- do.call(rbind, pred_rows_full)
pred_comm_all <- do.call(rbind, pred_rows_comm)

# Build the evaluation frame df (adds obs, errors, etc.)
obs_vec <- ifelse(pred_full_all$.rownames == "comm_male",
                  model_data$Male_CPUE[pred_full_all$t],
                  model_data$Female_CPUE[pred_full_all$t])

df <- data.frame(
  series    = pred_full_all$.rownames,
  t         = pred_full_all$t,
  h         = pred_full_all$h,
  origin    = pred_full_all$origin,
  obs       = obs_vec,
  pred_full = pred_full_all$estimate,
  pred_comm = pred_comm_all$estimate )

df$err_full <- df$obs - df$pred_full
df$err_comm <- df$obs - df$pred_comm
df$abs_full <- abs(df$err_full)
df$abs_comm <- abs(df$err_comm)

# Diebold-Mariano test
dm_results <- dm_grid(df)
print(dm_results)
# Interpret: negative mean(d) (stat < 0) means FULL has lower loss than COMM
# Results align with the skill table (BELOW)
# significant for males for t+2 and t+3
# significant for females t+1 and t+2

# Accuracy measures and improvement
acc <- df %>%
  group_by(series, h) %>%
  summarise(
    RMSE_full = sqrt(mean(err_full^2, na.rm=TRUE)),
    RMSE_comm = sqrt(mean(err_comm^2, na.rm=TRUE)),
    MAE_full  = mean(abs_full, na.rm=TRUE),
    MAE_comm  = mean(abs_comm, na.rm=TRUE),
    skill_RMSE = (1 - RMSE_full/RMSE_comm)*100,
    skill_MAE  = (1 - MAE_full/MAE_comm)*100,
    .groups = "drop"
  ) %>%
  arrange(series, h)

cat("\n=== Accuracy by horizon (h=1..", HORIZON, ") ===\n", sep="")
print(acc)

# Prettier version joining improvment and D-M test 
pretty_skill <- acc |>
  dplyr::select(series, h, skill_RMSE, skill_MAE) |>
  dplyr::left_join(
    dm_results |>
      tidyr::pivot_wider(id_cols = c(series, h),
                         names_from = loss,
                         values_from = p_value,
                         names_prefix = "DMp_"),
    by = c("series","h")
  ) |>
  dplyr::mutate(
    skill_RMSE = sprintf("%.1f%%", skill_RMSE),
    skill_MAE  = sprintf("%.1f%%", skill_MAE),
    DMp_abs    = scales::pvalue(DMp_abs, accuracy = 0.001),
    DMp_squared= scales::pvalue(DMp_squared, accuracy = 0.001)
  ) |>
  dplyr::arrange(series, h)

print(pretty_skill)

## ^^ The full model (including survey data) improved the short-term forecast accuracy
## for female commercial CPUE by 28-32% (t+1), tapering with horizon 15-22% when t+3.
## This sows that the survey is adding short-term predictive signal for female CPUE,
## some of that persists to 2–3 months.

## Result for males is more nuanced; smaller gain at t+1 (3-21%), then higher gains
## at t+2 (21-34%) and t+3 (23-33%). The survey is helping more with the trajectory
## than the immediate step — perhaps picking up dynamics not visible in comm-only data.

## Summary: Adding fishery-independent data improved commercial CPUE forecasts by
## 15–32% across 1–3 months, with strongest gains for females at h=1 and for males at h=2–3.

## ** NOTE 1: This uses OBSERVED covariate values for Temp and Conductivity, which would need
## to be forecasted if this was implemented as a forecasting tool, which would reduce the
## accuracy of forecasts.

## ** NOTE 2: Because the MARSS framework estimates latent abundance states, these
## forecasts can be interpreted as predictions of underlying stock status. The commercial
## CPUE forecasts reflect short-term fishery performance (“catchability”), while the
## state forecasts reflect population abundance (“stock health”). Both are relevant,
## but for long-term management, the state forecasts are more robust indicators of
## stock status.


# ============================================================================
#                     VALUE OF THE FIS - REDUCES NOISE ETC
# ============================================================================

# model_list_base$R  #must be "diagonal and unequal"

## Fit full model with simplified R structure
fit_em_full <- MARSS(y = y,
                     model = model_list_base,
                     inits = list(x0 = x0_init),
                     control = list(maxit = 2000))

fit_full <- MARSS(y = y,
                  model = model_list_base,
                  inits = fit_em_full,
                  method = "BFGS",
                  control = list(maxit = 1000))

## Fit commercial only model
fit_em_comm <- MARSS(y = y_comm,
                     model = model_comm_only,
                     inits = list(x0 = x0_init),
                     control = list(maxit = 2000))

# patch any zero/neg R diag in the EM par vector (2x1), to keep BFGS happy
Rpar <- fit_em_comm$par$R  # dim = 2 x 1 with rows "(comm_male,comm_male)", "(comm_female,comm_female)"
Rpar[1, 1] <- ifelse(!is.finite(Rpar[1,1]) || Rpar[1,1] <= 0, 1e-6, Rpar[1,1])  # bump zero/tiny
Rpar[2, 1] <- ifelse(!is.finite(Rpar[2,1]) || Rpar[2,1] <= 0, 1e-6, Rpar[2,1])
fit_em_comm$par$R <- Rpar

fit_comm <- MARSS(y = y_comm,
                  model = model_comm_only,
                  inits = fit_em_comm,
                  method = "BFGS",
                  control = list(maxit = 1000))

tidy(fit_full)
tify(fit_comm)

