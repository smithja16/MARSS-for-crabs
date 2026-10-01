## ============================================================================
##  Survey value, test 1: forecast accuracy by rolling-origin cross-validation
## ============================================================================
##  At each forecast origin s, both models are refitted to months 1..(s-1) and
##  used to forecast commercial CPUE at s, s+1, s+2.
##
##  WHY FORECASTS AND NOT AIC. The full and commercial-only models are fitted to
##  different observation matrices, so their likelihoods are on different
##  scales and cannot be compared. Forecast error on held-out data can. Only
##  the two commercial series are scored, so every model is judged on exactly
##  the same values however many series it was fitted to.
##
##  THE BASELINE. Fitted freely, the commercial-only model tends to drive its
##  observation variance to zero, which makes its forecasts chase every month's
##  sampling noise. It is therefore given the full model's estimates of the two
##  commercial observation variances, refitted at each origin. That helps the
##  baseline - a conservative choice. The unrepaired version can be included to
##  show the difference it makes.
##
##  FUTURE COVARIATES are the observed values, so absolute forecast skill is
##  optimistic for both fitted models. The comparison between them is fair
##  because both receive the same covariates; the seasonal naive does not.
## ============================================================================

run_cv <- function(prep, origins, horizon = 3, starts = "z", include_unrepaired = TRUE) {
  stopifnot(max(origins) + horizon - 1 <= prep$n, min(origins) > 13)
  comm <- c("comm_male", "comm_female")
  rows <- list()

  for (s in origins) {
    tr <- 1:(s - 1); fc <- s:(s + horizon - 1)
    cv_tr <- covars_at(prep$covars, tr)
    message(sprintf("  origin %d/%d (training to %s)", match(s, origins), length(origins),
                     format(prep$date[s - 1], "%b %Y")))

    full <- fit_marss(prep$y[, tr], model_full, cv_tr, starts = starts)
    Rf   <- coef(full$fit, type = "matrix")$R
    base <- fit_marss(prep$y[1:2, tr], function(cv) model_comm(cv, c(Rf[1, 1], Rf[2, 2])),
                      cv_tr, starts = starts)

    pred <- list(
      full     = forecast_recursion(full$fit, prep$covars[[full$start]][, fc, drop = FALSE])[comm, , drop = FALSE],
      baseline = forecast_recursion(base$fit, prep$covars[[base$start]][, fc, drop = FALSE]),
      seasonal_naive = prep$y[1:2, fc - 12, drop = FALSE])
    if (include_unrepaired) {
      free <- fit_marss(prep$y[1:2, tr], function(cv) model_comm(cv, NULL), cv_tr,
                        starts = starts, patch_R = TRUE)
      pred$baseline_unrepaired <- forecast_recursion(
        free$fit, prep$covars[[free$start]][, fc, drop = FALSE])
    }

    for (m in names(pred)) for (i in 1:2)
      rows[[length(rows) + 1]] <- data.frame(
        model = m, series = comm[i], origin = s, h = seq_len(horizon),
        actual = prep$y[i, fc], pred = as.numeric(pred[[m]][i, ]))
  }
  out <- do.call(rbind, rows)
  out$abs_err <- abs(out$pred - out$actual)
  out
}

## ----------------------------------------------------------------------------
## Summaries
## ----------------------------------------------------------------------------
## Mean absolute error by model, series and lead time, and the % reduction for
## the full model relative to the (repaired) baseline.
score_cv <- function(cv, baseline = "baseline") {
  mae <- aggregate(abs_err ~ model + series + h, data = cv, FUN = mean)
  names(mae)[4] <- "MAE"
  b <- mae[mae$model == baseline, c("series", "h", "MAE")]
  f <- mae[mae$model == "full",    c("series", "h", "MAE")]
  imp <- merge(f, b, by = c("series", "h"), suffixes = c("_full", "_baseline"))
  imp$improvement_pct <- 100 * (imp$MAE_baseline - imp$MAE_full) / imp$MAE_baseline
  list(mae = mae[order(mae$series, mae$h, mae$model), ], improvement = imp)
}

## Diebold-Mariano test of full vs baseline, per series and lead time.
## Forecasts from successive origins overlap at h > 1, so a Newey-West variance
## (lag h-1) is used; with few origins the Harvey-Leybourne-Newbold small-sample
## correction is applied and the statistic compared with a t distribution.
dm_test <- function(cv, baseline = "baseline", loss = c("absolute", "squared")) {
  loss <- match.arg(loss)
  out <- list()
  for (sr in unique(cv$series)) for (hh in sort(unique(cv$h))) {
    a <- cv[cv$model == "full"   & cv$series == sr & cv$h == hh, ]
    b <- cv[cv$model == baseline & cv$series == sr & cv$h == hh, ]
    m <- merge(a, b, by = "origin", suffixes = c("_f", "_b"))
    d <- if (loss == "absolute") m$abs_err_f - m$abs_err_b
         else (m$pred_f - m$actual_f)^2 - (m$pred_b - m$actual_b)^2
    n <- length(d)
    if (n < 3 || sd(d) == 0) next
    lmfit <- lm(d ~ 1)
    se    <- sqrt(sandwich::NeweyWest(lmfit, lag = hh - 1, prewhite = FALSE, adjust = TRUE)[1, 1])
    stat  <- coef(lmfit)[[1]] / se * sqrt((n + 1 - 2 * hh + hh * (hh - 1) / n) / n)
    out[[length(out) + 1]] <- data.frame(series = sr, h = hh, n = n, loss = loss,
                                         statistic = stat, p = 2 * pt(-abs(stat), n - 1))
  }
  do.call(rbind, out)
}
