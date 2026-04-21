
## NOTE: Code developed with assistance from AI (ChatGPT ~5, and Claude Sonnet ~4)

# ============================================================================
#                               FUNCTIONS
# ============================================================================


#########################################################
###  Function to create comprehensive lagged dataset  ###
#########################################################

create_modeling_dataset <- function(
    survey_data, 
    commercial_data, 
    conductivity_data,
    temperature_data,
    cwq_data) {
  
  #key_lags <- c(0, 1, 2, 3, 6, 12, 18)
  key_lags <- c(0, 1, 2, 4, 18)
  
  survey_wide <- survey_data %>%
    dplyr::select(yearmon, Site, date,
                  prerecruit_CPUE, prerecruit_F_CPUE, prerecruit_M_CPUE,
                  recruit_CPUE, recruit_F_CPUE, recruit_M_CPUE,
                  prerecruit_sex_ratio, recruit_sex_ratio) %>%
    arrange(Site, date)
  
  modeling_data <- commercial_data %>%
    dplyr::select(yearmon, date, Total_CPUE, Male_CPUE, Female_CPUE, Total_CatchWtKg,
                  Total_Male_CatchWtKg, Total_Female_CatchWtKg, Total_Traps) %>%
    arrange(date)
  
  # Prepare conductivity data for merging
  conductivity_prepared <- conductivity_data %>%
    dplyr::select(YearMonth, interpolated_Cond) %>%
    mutate(yearmon = format(YearMonth, "%Y-%m")) %>%
    arrange(YearMonth)
  
  # Prepare temperature data for merging
  temp_prepared <- temperature_data %>%
    dplyr::select(Year_Month, Mean_Temp) %>%
    rename(Temp_VR4G = Mean_Temp) %>%
    mutate(yearmon = format(Year_Month, "%Y-%m")) %>%
    arrange(Year_Month)
  
  # Prepare CWQ data for merging
  cwq_prepared <- cwq_data %>%
    dplyr::select(year_month, mean_potential_temp, mean_bottom_temp, mean_current_magnitude) %>%
    mutate(
      yearmon = format(year_month, "%Y-%m"),
      # Rename for clarity
      potential_temp = mean_potential_temp,
      bottom_temp = mean_bottom_temp,
      current = mean_current_magnitude
    ) %>%
    dplyr::select(yearmon, potential_temp, bottom_temp, current) %>%
    arrange(yearmon)
  
  # Add lagged survey predictors  
  for (lag in key_lags) {
    for (site in c(1:4)) {
      
      # Create survey data with dates moved FORWARD by lag months
      # This allows matching survey data from (current_month - lag) 
      # with commercial data from current_month
      site_data <- survey_wide %>%
        filter(Site == site) %>%
        arrange(date) %>%
        mutate(
          date_future = date + months(lag),    # Survey date + lag = commercial date
          yearmon_future = format(date_future, "%Y-%m") )
      
      lag_vars <- site_data %>%
        dplyr::select(yearmon_future,
                      prerecruit_CPUE, prerecruit_F_CPUE, prerecruit_M_CPUE,
                      recruit_CPUE, recruit_F_CPUE, recruit_M_CPUE,
                      prerecruit_sex_ratio, recruit_sex_ratio) %>%
        rename_with(~paste0(.x, "_site", site, "_lag", lag), 
                    -yearmon_future)
      
      modeling_data <- modeling_data %>%
        left_join(lag_vars, by = c("yearmon" = "yearmon_future"))
    }
  }
  
  # Add lagged conductivity data
  for (lag in key_lags) {
    # Create conductivity data with dates moved FORWARD by lag months
    cond_lag_data <- conductivity_prepared %>%
      arrange(YearMonth) %>%
      mutate(
        YearMonth_future = YearMonth + months(lag),
        yearmon_future = format(YearMonth_future, "%Y-%m")
      ) %>%
      dplyr::select(yearmon_future, interpolated_Cond) %>%
      rename_with(~paste0(.x, "_lag", lag), -yearmon_future)
    
    modeling_data <- modeling_data %>%
      left_join(cond_lag_data, by = c("yearmon" = "yearmon_future"))
  }
  
  # Add temperature data
  modeling_data <- modeling_data %>%
    left_join(temp_prepared %>% select(yearmon, Temp_VR4G), by = "yearmon")
  
  # Add lagged CWQ data (potential temp, bottom temp, current)
  for (lag in key_lags) {
    # Create CWQ data with dates moved FORWARD by lag months
    cwq_lag_data <- cwq_prepared %>%
      arrange(yearmon) %>%
      mutate(
        # Convert yearmon back to date for arithmetic
        date_temp = as.Date(paste0(yearmon, "-01")),
        date_future = date_temp + months(lag),
        yearmon_future = format(date_future, "%Y-%m")
      ) %>%
      dplyr::select(yearmon_future, potential_temp, bottom_temp, current) %>%
      rename_with(~paste0(.x, "_lag", lag), -yearmon_future)
    
    modeling_data <- modeling_data %>%
      left_join(cwq_lag_data, by = c("yearmon" = "yearmon_future"))
  }
  
  # Add temporal predictors
  modeling_data <- modeling_data %>%
    mutate(
      year = year(date),
      month = month(date),
      quarter = quarter(date),
      season = case_when(
        month %in% c(12, 1, 2) ~ "Summer",
        month %in% c(3, 4, 5) ~ "Autumn", 
        month %in% c(6, 7, 8) ~ "Winter",
        month %in% c(9, 10, 11) ~ "Spring"
      ),
      # Cyclical time variables
      month_sin = sin(2 * pi * month / 12),  #Fourier series
      month_cos = cos(2 * pi * month / 12),
      # Effort standardization
      effort_scaled = scale(Total_Traps)[,1]
    ) %>%
    # Remove rows with missing predictors
    # filter(complete.cases(.))
    
    return(modeling_data)
}



#########################################################
###  Summarise diagnostics                            ###
#########################################################

compute_diagnostics <- function(fit, model_name) {
  fitted_vals <- fitted(fit)
  resids <- residuals(fit, type = "tt1")
  states <- tsSmooth(fit)
  
  obs_names <- unique(fitted_vals$.rownames)
  
  metrics_list <- lapply(obs_names, function(obs) {
    fit_sub <- fitted_vals[fitted_vals$.rownames == obs, ]
    resid_sub <- resids[resids$.rownames == obs, ]
    
    # Calculate normalized metrics
    obs_mean <- mean(fit_sub$y, na.rm = TRUE)
    obs_sd <- sd(fit_sub$y, na.rm = TRUE)
    
    data.frame(
      Model = model_name,
      Series = obs,
      RMSE = sqrt(mean((fit_sub$y - fit_sub$.fitted)^2, na.rm = TRUE)),
      NRMSE = sqrt(mean((fit_sub$y - fit_sub$.fitted)^2, na.rm = TRUE)) / obs_sd,  # Normalized by SD
      MAPE = mean(abs((fit_sub$y - fit_sub$.fitted) / fit_sub$y) * 100, na.rm = TRUE),  # %
      R2 = cor(fit_sub$y, fit_sub$.fitted, use = "complete.obs")^2,
      ACF1 = acf(resid_sub$.resids, plot = FALSE, na.action = na.pass)$acf[2],
      ACF1sr = acf(resid_sub$.std.resids, plot = FALSE, na.action = na.pass)$acf[2],
      N_obs = sum(!is.na(fit_sub$y))
    )
  })
  
  metrics_df <- do.call(rbind, metrics_list)
  metrics_df$AIC <- fit$AIC
  metrics_df$LogLik <- fit$logLik
  metrics_df$N_params <- fit$num.params
  metrics_df$Mean_State_SE <- mean(states$.se, na.rm = TRUE)
  
  return(metrics_df)
}



#########################################################
###  Diebold-Mariano test for forecast accuracy       ###
#########################################################

## Compares loss from full model vs comm-only model

# Diebold–Mariano test comparing FULL vs COMM losses at (series, h)
dm_test_abs <- function(dat, series, h, loss = c("abs","squared")) {
  loss   <- match.arg(loss)
  series <- as.character(series)
  h      <- as.integer(h)
  
  # base indexing (no tidy-eval)
  sub <- dat[dat$series == series & dat$h == h, , drop = FALSE]
  
  if (nrow(sub) < 3L) {
    return(data.frame(series = series, h = h, stat = NA_real_, p_value = NA_real_, loss = loss))
  }
  
  d <- if (loss == "abs") (sub$abs_full - sub$abs_comm) else (sub$err_full^2 - sub$err_comm^2)
  d <- d[is.finite(d)]
  if (length(d) < 3L || sd(d) == 0) {
    return(data.frame(series = series, h = h, stat = NA_real_, p_value = NA_real_, loss = loss))
  }
  
  if (requireNamespace("sandwich", quietly = TRUE) && requireNamespace("lmtest", quietly = TRUE)) {
    fit <- stats::lm(d ~ 1)
    L   <- max(0L, h - 1L)  # HAC lag for h-step overlaps
    vc  <- sandwich::NeweyWest(fit, lag = L, prewhite = FALSE, adjust = TRUE)
    ct  <- lmtest::coeftest(fit, vcov. = vc)
    stat <- unname(ct[1, "t value"])
    pval <- unname(ct[1, "Pr(>|t|)"])
  } else {
    # fallback: one-sample t-test on mean(d)
    stat <- mean(d) / (stats::sd(d) / sqrt(length(d)))
    pval <- 2 * stats::pt(-abs(stat), df = length(d) - 1L)
  }
  
  data.frame(series = series, h = h, stat = stat, p_value = pval, loss = loss)
}

# Grid over all series × h
dm_grid <- function(dat) {
  sh <- unique(dat[, c("series", "h")])
  out_list <- lapply(seq_len(nrow(sh)), function(i) {
    s  <- sh$series[i]
    hh <- as.integer(sh$h[i])
    rbind(
      dm_test_abs(dat, s, hh, "abs"),
      dm_test_abs(dat, s, hh, "squared")
    )
  })
  do.call(rbind, out_list)
}



#########################################################
###  Precision calculation of latent state forecast   ###
#########################################################

prec_one <- function(df, Ttrain, H, want80 = FALSE) {
  if (!is.data.frame(df) || nrow(df) == 0) {
    out <- data.frame(.rownames = character(0),
                      h = integer(0),
                      mean_w95 = numeric(0),
                      stringsAsFactors = FALSE)
    if (want80) out$mean_w80 <- numeric(0)
    return(out)
  }
  
  # Relative horizon = 1..H from absolute time
  df$rel_h <- df$t_abs - Ttrain
  
  # Keep only final H horizons and the two latent states
  keep <- (df$.rownames %in% c("X1","X2")) & 
    (df$t_abs %in% (Ttrain+1):(Ttrain+H))
  d <- df[keep, , drop = FALSE]
  if (nrow(d) == 0) {
    out <- data.frame(.rownames = character(0),
                      h = integer(0),
                      mean_w95 = numeric(0),
                      stringsAsFactors = FALSE)
    if (want80) out$mean_w80 <- numeric(0)
    return(out)
  }
  
  # Compute widths
  d$w95 <- d[["Hi 95"]] - d[["Lo 95"]]
  if (want80) {
    if (!all(c("Lo 80","Hi 80") %in% names(d))) {
      d$w80 <- NA_real_
    } else {
      d$w80 <- d[["Hi 80"]] - d[["Lo 80"]]
    }
  }
  
  # Aggregate by state and (relative) horizon
  agg95 <- aggregate(w95 ~ .rownames + rel_h, data = d, FUN = function(z) mean(z, na.rm = TRUE))
  names(agg95) <- c("state","h","mean_w95")
  
  if (want80) {
    agg80 <- aggregate(w80 ~ .rownames + rel_h, data = d, FUN = function(z) mean(z, na.rm = TRUE))
    names(agg80) <- c("state","h","mean_w80")
    out <- merge(agg95, agg80, by = c("state","h"), all = TRUE)
  } else {
    out <- agg95
  }
  
  out
}
