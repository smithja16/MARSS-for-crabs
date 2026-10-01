## ============================================================================
##  Generate a synthetic dataset with the structure of the Wallis Lake data
## ============================================================================
##  The real commercial catch data are confidential. This script simulates data
##  from the model itself, using the parameter estimates published in Table S2,
##  so the rest of the workflow can be run without them. The result has the same
##  dimensions, the same missing survey months and similar covariates, but it
##  is NOT the real data and will not reproduce the paper's numbers.
##
##  Because the data are simulated, the true latent states are known and are
##  saved too. That allows something the real analysis cannot do: check how
##  closely the model recovers the truth (see run_example.R).
##
##  Writes example_data/wallis_synthetic.csv and example_data/true_states.csv.
## ============================================================================

source("R/published_parameters.R")
set.seed(20181101)

n_months <- 67                                   # Nov 2018 - May 2024
n_extra  <- 4                                    # see "prerecruits" below
n_sim    <- n_months + n_extra
dates    <- seq(as.Date("2018-11-01"), by = "month", length.out = n_sim)
month    <- as.integer(format(dates, "%m"))
missing_survey_months <- 33:36                   # Jul-Oct 2021, as in the real data

ar1 <- function(n, phi, sd_innov) {
  e <- numeric(n); e[1] <- rnorm(1, 0, sd_innov / sqrt(1 - phi^2))
  for (i in 2:n) e[i] <- phi * e[i - 1] + rnorm(1, 0, sd_innov)
  e
}

## ---- covariates -------------------------------------------------------------
## Offshore bottom temperature, already lagged 2 months: a seasonal cycle that
## explains ~40% of its variance plus autocorrelated noise (as in the real
## series). Conductivity: mostly non-seasonal autocorrelated variation around
## its estuarine mean, with occasional freshwater drops.
cos_m <- cos(2 * pi * month / 12)
sin_m <- sin(2 * pi * month / 12)
bottom_temp_lag2 <- 17.08 - 0.44 * cos_m - 0.39 * sin_m + ar1(n_sim, 0.5, 0.44)
conductivity     <- 48.5 + 1.5 * cos_m + ar1(n_sim, 0.8, 3.2)

## The model uses standardised temperature and conductivity, centred and
## scaled over the modelling period
zs <- function(v) (v - mean(v[1:n_months])) / sd(v[1:n_months])
covars <- rbind(zs(bottom_temp_lag2), zs(conductivity), cos_m, sin_m)

## ---- latent states and observations, in model form --------------------------
##   x_t = B x_{t-1} + u + C c_t + w_t,   w_t ~ MVN(0, Q)
##   y_t = Z x_t + v_t,                   v_t ~ MVN(0, R)
p <- published
x <- matrix(NA_real_, 2, n_sim)
x_prev <- p$x0
for (t in seq_len(n_sim)) {
  x_prev <- p$B %*% x_prev + p$U + p$C %*% covars[, t, drop = FALSE] +
            t(chol(p$Q)) %*% rnorm(2)
  x[, t] <- x_prev
}
y <- p$Z %*% x + t(chol(p$R)) %*% matrix(rnorm(6 * n_sim), 6, n_sim)

## Gaussian errors permit negative values, which a catch rate cannot take; they
## are set to zero (survey catch rates are often zero anyway).
y[y < 0] <- 0

## ---- prerecruits: from model form back to the month they were collected ------
## In the model, row 5 at month t holds the male prerecruit catch rate collected
## at month t - 2, and row 6 the female catch rate collected at t - 4. A real
## dataset records each value against the month it was COLLECTED, so the rows
## are shifted back here. run_example.R shifts them forward again.
## This is why n_extra months were simulated: the prerecruits collected in the
## last four months inform states beyond the modelling period.
prerecruit_male   <- y[5, 3:(n_months + 2)]
prerecruit_female <- y[6, 5:(n_months + 4)]

## ---- assemble the as-collected table --------------------------------------
idx <- seq_len(n_months)
d <- data.frame(
  date                   = dates[idx],
  year                   = as.integer(format(dates[idx], "%Y")),
  month                  = month[idx],
  comm_cpue_male         = y[1, idx],
  comm_cpue_female       = y[2, idx],
  surv_recruit_male      = y[3, idx],
  surv_recruit_female    = y[4, idx],
  surv_prerecruit_male   = prerecruit_male,
  surv_prerecruit_female = prerecruit_female,
  bottom_temp_lag2       = bottom_temp_lag2[idx],
  conductivity           = conductivity[idx]
)
survey_cols <- grep("^surv_", names(d))
d[missing_survey_months, survey_cols] <- NA        # months the survey did not run

dir.create("example_data", showWarnings = FALSE)
num <- vapply(d, is.numeric, logical(1)) & !names(d) %in% c("year", "month")
d[num] <- lapply(d[num], round, 5)
write.csv(d, "example_data/wallis_synthetic.csv", row.names = FALSE)
write.csv(data.frame(date = dates[idx], true_state_male = x[1, idx],
                     true_state_female = x[2, idx]),
          "example_data/true_states.csv", row.names = FALSE)

message("Wrote example_data/wallis_synthetic.csv (", n_months, " months) ",
        "and example_data/true_states.csv")
