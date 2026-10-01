## ============================================================================
##  Survey value, test 2: precision of the abundance estimates
## ============================================================================
##  Forecast accuracy (survey_value_cv.R) is about predicting catch rates. This
##  test is about estimating abundance. Using ONE fitted model, survey
##  observations are removed from randomly chosen months and the Kalman smoother
##  is re-run with every parameter held at its fitted value. Nothing changes
##  except the information available, so the rise in state standard error is
##  attributable to losing the survey.
##
##  Removing a survey month means removing everything that month contributed:
##  the recruit rows in that column, and the prerecruit values that were shifted
##  forward 2 (male) and 4 (female) months. Blanking whole columns instead would
##  remove parts of three different survey months at once.
##
##  Only interior months are blanked and scored, because uncertainty near the
##  start and end of a series is dominated by position rather than data.
## ============================================================================

## State standard errors after removing the given survey months, parameters fixed
smooth_without <- function(fit, prep, months) {
  g <- fit
  for (tau in months) for (r in prep$survey_rows) {
    col <- tau + prep$row_lag[r]
    if (col >= 1 && col <= prep$n) {
      g$marss$data[r, col] <- NA
      g$model$data[r, col] <- NA
    }
  }
  stopifnot(identical(g$par, fit$par))          # the parameters must not move
  V <- MARSSkf(g)$VtT
  t(sqrt(apply(V, 3, diag)))                    # months x 2 (male, female)
}

survey_removal_test <- function(fit, prep, coverage = c(1, 0.8, 0.6, 0.4, 0.2, 0),
                                n_draws = 20, seed = 1) {
  set.seed(seed)
  interior <- 7:(prep$n - 2)
  surveyed <- intersect(interior, which(!is.na(prep$y[3, ])))

  rows <- list()
  for (cv in coverage) {
    n_remove <- round((1 - cv) * length(surveyed))
    draws <- if (n_remove %in% c(0, length(surveyed))) 1 else n_draws
    for (k in seq_len(draws)) {
      gone <- if (n_remove == 0) integer(0) else sample(surveyed, n_remove)
      se   <- smooth_without(fit, prep, gone)
      rows[[length(rows) + 1]] <- data.frame(
        coverage = cv, draw = k,
        se_male = mean(se[interior, 1]), se_female = mean(se[interior, 2]))
    }
  }
  res <- do.call(rbind, rows)

  s <- aggregate(cbind(se_male, se_female) ~ coverage, data = res, FUN = mean)
  full <- s[s$coverage == 1, ]
  s$male_increase_pct   <- 100 * (s$se_male   / full$se_male   - 1)
  s$female_increase_pct <- 100 * (s$se_female / full$se_female - 1)
  s[order(-s$coverage), ]
}
