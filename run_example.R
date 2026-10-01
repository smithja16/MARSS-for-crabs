## ============================================================================
##  Worked example: integrating commercial and survey catch rates with MARSS,
##  and measuring what the survey adds
## ============================================================================
##  Runs on the synthetic data from simulate_data.R. Results are illustrative
##  and do not reproduce the paper's numbers.
##
##    source("simulate_data.R")   # once
##    source("run_example.R")
##
##  Four steps:
##    1. Fit the full model (commercial + survey).
##    2. Because the data are synthetic, check it against the true states.
##    3. Survey value by forecast accuracy (rolling-origin cross-validation).
##    4. Survey value by estimation precision (survey-removal test).
## ============================================================================

suppressMessages(library(MARSS))
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)

## ---- settings ----------------------------------------------------------------
## The paper used 12 origins (t = 54..65), both starts at every origin, and 30
## draws per coverage level. These defaults trade some precision for a run time
## of roughly 10-15 minutes; increase them for a closer replica of the design.
CV_ORIGINS <- seq(54, 65, by = 2)   # forecast origins (month index)
CV_HORIZON <- 3                     # months ahead
CV_STARTS  <- "z"                   # starts per fit inside the CV
SE_DRAWS   <- 20                    # random draws per coverage level
OUT_DIR    <- "output"
dir.create(OUT_DIR, showWarnings = FALSE)

if (!file.exists("example_data/wallis_synthetic.csv"))
  stop("Run simulate_data.R first.")
d     <- read.csv("example_data/wallis_synthetic.csv")
truth <- read.csv("example_data/true_states.csv")
prep  <- prepare_series(d)


## ============================================================================
##  1. Fit the full model
## ============================================================================
message("1. Fitting the full model (two starts)")
full <- fit_marss(prep$y, model_full, prep$covars)
cat(sprintf("\nFull model: logLik %.2f, %d parameters, best start: %s\n",
            full$logLik, full$fit$num.params, full$start))

## Parameter recovery. Only parameters that do not depend on how the covariates
## are scaled are compared (u and C change with scaling; B, Z, Q, R do not).
est <- coef(full$fit, type = "matrix")
p   <- published
rec <- data.frame(
  parameter = c("b11 male persistence", "b22 female persistence",
                "b21 male -> female", "b12 female -> male",
                "z3 male recruit loading", "z4 female recruit loading",
                "z5 male prerecruit loading", "z6 female prerecruit loading",
                "process error correlation",
                "r11 commercial male", "r22 commercial female"),
  true = c(p$B[1, 1], p$B[2, 2], p$B[2, 1], p$B[1, 2],
           p$Z[3, 1], p$Z[4, 2], p$Z[5, 1], p$Z[6, 2],
           p$Q[1, 2] / sqrt(p$Q[1, 1] * p$Q[2, 2]), p$R[1, 1], p$R[2, 2]),
  estimate = c(est$B[1, 1], est$B[2, 2], est$B[2, 1], est$B[1, 2],
               est$Z[3, 1], est$Z[4, 2], est$Z[5, 1], est$Z[6, 2],
               est$Q[1, 2] / sqrt(est$Q[1, 1] * est$Q[2, 2]), est$R[1, 1], est$R[2, 2]))
cat("\nParameter recovery (67 months of data, 31 parameters):\n")
print(transform(rec, true = signif(true, 3), estimate = signif(estimate, 3)), row.names = FALSE)
cat("Expect only rough agreement: 67 months is a short series for this many parameters.\n")


## ============================================================================
##  2. Against the truth: does the survey improve the abundance estimates?
## ============================================================================
## Something the real analysis cannot do. The commercial-only model is fitted
## with its observation variances fixed at the full model's values, as in step 3.
message("2. Comparing estimated and true states")
Rf <- est$R
comm <- fit_marss(prep$y[1:2, ], function(cv) model_comm(cv, c(Rf[1, 1], Rf[2, 2])),
                  prep$covars)

## Smoothed states and their standard errors, by POSITION (row 1 male, row 2
## female). MARSS names states differently depending on the model - "X1" here,
## "X.comm_male" in the commercial-only model - so names are not relied on.
states_of <- function(fit) {
  V <- MARSSkf(fit)$VtT
  list(est = fit$states, se = sqrt(rbind(V[1, 1, ], V[2, 2, ])))
}

state_check <- function(fit, label) {
  s <- states_of(fit)
  do.call(rbind, lapply(1:2, function(i) {
    tr <- truth[[c("true_state_male", "true_state_female")[i]]]
    data.frame(model = label, state = c("male", "female")[i],
               rmse_vs_truth = sqrt(mean((s$est[i, ] - tr)^2)),
               correlation   = cor(s$est[i, ], tr),
               coverage_95   = mean(abs(s$est[i, ] - tr) <= 1.96 * s$se[i, ]))
  }))
}
sc <- rbind(state_check(full$fit, "full"), state_check(comm$fit, "commercial only"))
cat("\nSmoothed states against the true (simulated) states:\n")
print(transform(sc, rmse_vs_truth = round(rmse_vs_truth, 4),
                correlation = round(correlation, 3), coverage_95 = round(coverage_95, 2)),
      row.names = FALSE)
cat("coverage_95 is the share of months whose true state lies inside the 95% interval.\n")


## ============================================================================
##  3. Survey value: forecast accuracy
## ============================================================================
message("3. Cross-validation (", length(CV_ORIGINS), " origins)")
cv <- run_cv(prep, CV_ORIGINS, CV_HORIZON, starts = CV_STARTS)
sc_cv <- score_cv(cv)

cat("\nMean absolute error of commercial CPUE forecasts:\n")
w <- reshape(sc_cv$mae, idvar = c("series", "h"), timevar = "model", direction = "wide")
names(w) <- sub("^MAE[.]", "", names(w))
print(w, row.names = FALSE, digits = 3)

cat("\nImprovement of the full model over the baseline (% reduction in MAE):\n")
print(transform(sc_cv$improvement, improvement_pct = round(improvement_pct, 1))[
  , c("series", "h", "MAE_full", "MAE_baseline", "improvement_pct")], row.names = FALSE, digits = 3)

cat("\nDiebold-Mariano tests, full model vs baseline (small-sample corrected):\n")
dm <- rbind(dm_test(cv, loss = "absolute"), dm_test(cv, loss = "squared"))
print(transform(dm, statistic = round(statistic, 2), p = signif(p, 2)), row.names = FALSE)
cat("With", length(CV_ORIGINS), "origins these tests have little power.\n")


## ============================================================================
##  4. Survey value: estimation precision
## ============================================================================
message("4. Survey-removal test")
se <- survey_removal_test(full$fit, prep, n_draws = SE_DRAWS)
cat("\nMean state standard error as survey coverage is reduced (parameters fixed):\n")
print(transform(se, se_male = signif(se_male, 3), se_female = signif(se_female, 3),
                male_increase_pct = round(male_increase_pct, 1),
                female_increase_pct = round(female_increase_pct, 1)), row.names = FALSE)


## ============================================================================
##  Figures
## ============================================================================
sf <- states_of(full$fit); scm <- states_of(comm$fit)
png(file.path(OUT_DIR, "states_vs_truth.png"), width = 9, height = 7, units = "in", res = 150)
op <- par(mfrow = c(2, 1), mar = c(3, 4.5, 2, 1))
for (i in 1:2) {
  est <- sf$est[i, ]; s_i <- sf$se[i, ]; estc <- scm$est[i, ]   # not `se`: that holds step 4's results
  tr <- truth[[c("true_state_male", "true_state_female")[i]]]
  yl <- range(est - 1.96 * s_i, est + 1.96 * s_i, tr, prep$y[i, ], na.rm = TRUE)
  plot(prep$date, tr, type = "n", ylim = yl, las = 1, xlab = "", ylab = "CPUE (kg per trap)")
  polygon(c(prep$date, rev(prep$date)),
          c(est - 1.96 * s_i, rev(est + 1.96 * s_i)),
          border = NA, col = adjustcolor("grey60", 0.4))
  points(prep$date, prep$y[i, ], pch = 16, col = "grey55", cex = 0.7)
  lines(prep$date, tr, lwd = 2, col = "#1B7837")
  lines(prep$date, est,  lwd = 2)
  lines(prep$date, estc, lwd = 1.5, lty = 2, col = "#B2182B")
  mtext(c("(a) Male", "(b) Female")[i], side = 3, adj = 0, line = 0.3, font = 2)
  if (i == 1) legend("topright", bty = "n", cex = 0.8,
                     c("True state", "Full model (95% CI)", "Commercial only", "Commercial CPUE"),
                     lty = c(1, 1, 2, NA), lwd = c(2, 2, 1.5, NA), pch = c(NA, NA, NA, 16),
                     col = c("#1B7837", "black", "#B2182B", "grey55"))
}
par(op); dev.off()

png(file.path(OUT_DIR, "survey_removal.png"), width = 6, height = 4.5, units = "in", res = 150)
plot(100 * se$coverage, se$male_increase_pct, type = "b", pch = 16, col = "#B2182B", las = 1,
     ylim = range(0, se$male_increase_pct, se$female_increase_pct),
     xlab = "Survey coverage (% of months retained)",
     ylab = "Increase in state standard error (%)", xlim = c(100, 0))
lines(100 * se$coverage, se$female_increase_pct, type = "b", pch = 16, col = "#2166AC")
legend("topleft", bty = "n", c("Male", "Female"), pch = 16, col = c("#B2182B", "#2166AC"))
dev.off()

message("Done. Figures written to ", OUT_DIR, "/")
