## ============================================================================
##  Fitting and forecasting
## ============================================================================

## ----------------------------------------------------------------------------
## Two-stage fit from multiple starts
## ----------------------------------------------------------------------------
## Each fit runs EM first and then quasi-Newton (BFGS) initialised from it: EM
## gets close quickly but may stall, and BFGS finishes the job.
##
## MARSS can report convergence at a local optimum, so the fit is repeated from
## more than one start and the best likelihood kept. Here the starts are the
## same model with standardised and with raw covariates: these are exactly
## equivalent (the drift u absorbs the difference), so their likelihoods are
## directly comparable, but they send the optimiser down different paths.
##
##   y          observation matrix
##   build      function(covars) returning a MARSS model list
##   covars     named list of covariate matrices, one per start
##   patch_R    floor non-positive EM variances before BFGS (only needed for
##              the commercial-only model with free R; see model_spec.R)
##
## Returns the best fit and which start produced it. The start matters later:
## forecasts must use covariates on the same scale the model was fitted on.

fit_marss <- function(y, build, covars, starts = names(covars),
                      em_maxit = 2000, bfgs_maxit = 1000, patch_R = FALSE) {
  x0 <- matrix(c(mean(y[1, 1:3], na.rm = TRUE), mean(y[2, 1:3], na.rm = TRUE)), 2, 1)
  best <- list(fit = NULL, start = NA_character_, logLik = -Inf)

  for (st in starts) {
    model <- build(covars[[st]])
    em <- tryCatch(MARSS(y, model = model, inits = list(x0 = x0),
                         control = list(maxit = em_maxit), silent = TRUE),
                   error = function(e) NULL)
    if (is.null(em)) next
    if (patch_R && !is.null(em$par$R))
      em$par$R[em$par$R <= 0 | !is.finite(em$par$R)] <- 1e-6
    fit <- tryCatch(MARSS(y, model = model, inits = em, method = "BFGS",
                          control = list(maxit = bfgs_maxit), silent = TRUE),
                    error = function(e) NULL)
    if (!is.null(fit) && is.finite(fit$logLik) && fit$logLik > best$logLik)
      best <- list(fit = fit, start = st, logLik = as.numeric(fit$logLik))
  }
  if (is.null(best$fit)) stop("fit_marss(): no start converged")
  best
}

## ----------------------------------------------------------------------------
## Forecast by direct recursion
## ----------------------------------------------------------------------------
## Iterate the state equation from the last estimated state, then map through
## the observation equation:
##   x_{t+h} = B x_{t+h-1} + u + C c_{t+h},    y_{t+h} = Z x_{t+h} + a
## `c_future` must be on the scale the model was fitted on (the `start`
## returned by fit_marss). This gives the same point forecasts as
## predict.marssMLE(); it is used because predict() fails in MARSS 3.11.10
## under R >= 4.3 for models with covariates in the observation equation.

forecast_recursion <- function(fit, c_future) {
  p <- coef(fit, type = "matrix")
  x <- fit$states[, ncol(fit$states), drop = FALSE]
  out <- matrix(NA_real_, nrow(p$Z), ncol(c_future),
                dimnames = list(rownames(fit$model$data), NULL))
  for (h in seq_len(ncol(c_future))) {
    x <- p$B %*% x + p$U + p$C %*% c_future[, h, drop = FALSE]
    out[, h] <- p$Z %*% x + p$A
  }
  out
}

## Subset each covariate matrix in a named list to the given months

covars_at <- function(covars, months) lapply(covars, function(m) m[, months, drop = FALSE])
