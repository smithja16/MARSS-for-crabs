## ============================================================================
##  The model: two latent states observed by six monthly catch-rate series
## ============================================================================
##
##  State equation      x_t = B x_{t-1} + u + C c_t + w_t,   w_t ~ MVN(0, Q)
##  Observation eq.     y_t = Z x_t + a     + v_t,           v_t ~ MVN(0, R)
##
##  x_t  male and female fishery-available abundance, on the commercial CPUE
##       scale (because Z = 1 for the commercial series)
##  y_t  the six catch-rate series in month t
##  c_t  covariates: bottom temperature (2-month lag), conductivity, and a
##       two-term Fourier seasonal cycle
## ============================================================================

## ----------------------------------------------------------------------------
## Data preparation: from the as-collected table to the model's y and c
## ----------------------------------------------------------------------------
## The observation equation can only relate y at month t to x at month t, so
## prerecruit lags are applied by SHIFTING THE DATA: the male prerecruit value
## collected at month t-2 is placed in column t, so that it informs the male
## state two months later (when those crabs reach legal size); females are
## shifted four months. Each prerecruit row is therefore dated by the state it
## informs, not by when it was collected. The shifted-in leading months are NA,
## which MARSS handles natively.
prepare_series <- function(d) {
  n <- nrow(d)
  lag_by <- function(v, k) c(rep(NA, k), v[seq_len(n - k)])

  y <- rbind(
    comm_male              = d$comm_cpue_male,
    comm_female            = d$comm_cpue_female,
    surv_recruit_male      = d$surv_recruit_male,
    surv_recruit_female    = d$surv_recruit_female,
    surv_prerecruit_male   = lag_by(d$surv_prerecruit_male,   2),
    surv_prerecruit_female = lag_by(d$surv_prerecruit_female, 4))

  ## Temperature and conductivity are standardised so that their coefficients
  ## are comparable and u is not confounded with the mean covariate level. The
  ## Fourier terms are already bounded and centred and are left as they are.
  ## The unstandardised version is kept because fitting from both is used as
  ## a cheap multi-start (see fitting.R).
  fourier <- rbind(cos = cos(2 * pi * d$month / 12), sin = sin(2 * pi * d$month / 12))
  zs <- function(v) (v - mean(v)) / sd(v)
  covars <- list(
    z   = rbind(temp = zs(d$bottom_temp_lag2), cond = zs(d$conductivity), fourier),
    raw = rbind(temp = d$bottom_temp_lag2,     cond = d$conductivity,     fourier))

  list(y = y, covars = covars, n = n, date = as.Date(d$date), month = d$month,
       survey_rows = 3:6,
       row_lag = c(NA, NA, 0, 0, 2, 4))   # months between collection and column
}

## ----------------------------------------------------------------------------
## Full model: commercial + survey series
## ----------------------------------------------------------------------------
## `list()` matrices mix fixed numbers and parameter names; a character matrix
## would make MARSS estimate a parameter literally called "1".
model_full <- function(covars) {
  n_cov <- nrow(covars)
  C <- matrix(list(), 2, n_cov)
  for (j in seq_len(n_cov)) { C[[1, j]] <- paste0("c1", j); C[[2, j]] <- paste0("c2", j) }

  list(
    ## Persistence (b11, b22) and cross-sex effects (b12, b21)
    B = matrix(list("b11", "b12",
                    "b21", "b22"), 2, 2, byrow = TRUE),
    ## State-specific drift
    U = matrix(c("u1", "u2"), 2, 1),
    ## Covariate effects on each state
    C = C, c = covars,
    ## Loadings. Commercial series are fixed at 1, which puts both states on the
    ## commercial CPUE scale; each survey series has its own loading, which
    ## absorbs the difference in units and catchability.
    Z = matrix(list(1,    0,
                    0,    1,
                    "z3", 0,
                    0,    "z4",
                    "z5", 0,
                    0,    "z6"), 6, 2, byrow = TRUE),
    ## No offsets: each survey series' scale is carried entirely by Z
    A = "zero",
    ## Observation errors: independent, except that each sex's recruit and
    ## prerecruit series come from the same traps and may share error
    R = matrix(list("r11", 0,     0,     0,     0,     0,
                    0,     "r22", 0,     0,     0,     0,
                    0,     0,     "r33", 0,     "r35", 0,
                    0,     0,     0,     "r44", 0,     "r46",
                    0,     0,     "r35", 0,     "r55", 0,
                    0,     0,     0,     "r46", 0,     "r66"), 6, 6, byrow = TRUE),
    ## Process errors, correlated between the sexes
    Q = matrix(list("q11", "q12",
                    "q12", "q22"), 2, 2, byrow = TRUE)
  )
}

## ----------------------------------------------------------------------------
## Commercial-only model: the baseline the survey is compared against
## ----------------------------------------------------------------------------
## Identical dynamics, but only the two commercial series. With one series per
## state, its observation and process variances are only weakly identified: if
## left free, the observation variance tends to collapse towards zero. The
## comparison therefore fixes them (`r_fixed`) at values estimated by the full
## model - deliberately helping the baseline. `r_fixed = NULL` estimates them,
## which reproduces the problem.
model_comm <- function(covars, r_fixed = NULL) {
  m <- model_full(covars)
  m$Z <- matrix(list(1, 0, 0, 1), 2, 2, byrow = TRUE)
  m$R <- if (is.null(r_fixed)) {
    matrix(list("r11", 0, 0, "r22"), 2, 2, byrow = TRUE)
  } else {
    matrix(list(r_fixed[1], 0, 0, r_fixed[2]), 2, 2, byrow = TRUE)
  }
  m
}
