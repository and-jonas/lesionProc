
#=============================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2025  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

#=============================================================================== -

# extracts temperature course for specific site and time interval
extract_covar <- function(df, from, to) {
  df_ <- df %>% filter(timestamp > from, timestamp <= to)
  return(df_)
}

# extracts variables from list columns
extract_covars_from_nested <- function(tbl, from, vars)
{
  dat_vars <- list()
  for(i in vars){
    dat_vars[[i]] <- do.call("c", lapply(tbl[[from]], "[[", i))
  }
  if(length(vars)>1){
    vars <- do.call("cbind", dat_vars)
    out <- cbind(tbl, vars)
    out <- as_tibble(out)
  } else {
    out <- as_tibble(dat_vars)
    out <- bind_cols(tbl, out)
  }
  return(out)
}

# four-parameter beta function
temperature_response <- function(T, tmin, topt, tmax, scale) {
  if (T <= tmin || T >= tmax) {
    return(0)
  } else if (T <= topt) {
    return(((T - tmin) / (topt - tmin))^scale)
  } else {
    return(((tmax - T) / (tmax - topt))^scale)
  }
}

# effective thermal time duration
get_effective_time <- function(data, tmin, tmax, topt, scale, dt){
  if (nrow(data) == 0) return(0)
  temperature_series <- data$temp
  effective_time <- sum(sapply(temperature_series, function(X) temperature_response(X, tmin, topt, tmax, scale) * dt))
  return(effective_time)
}

# get all values at the lag time point (t0) and calculate the differences
# between t1 and t0
get_lags_diffs <- function(data, vars) {
  # get lag values for all lesion features
  df_lags <- data %>%
    mutate(across(
      .cols = vars, 
      .fns = list(lag_ = ~lag(.)), 
      .names = "lag_{.col}"
    ))
  # get relevant differences
  df_lags <- df_lags %>% mutate(lag_timestamp = dplyr::lag(timestamp)) %>% 
    relocate(lag_timestamp, .after=timestamp)
  df_lags_deltas <- df_lags %>% 
    mutate(
      diff_time = difftime(timestamp, lag_timestamp, units = "hours"),
      diff_gdd = lesion_age_gdd - lag_lesion_age_gdd,
      diff_area = area - lag_area,
      # perimeter-normalized differences 
      diff_area_pp_xy = diff_area / lag_xy_perimeter_f,
      diff_area_pp_x = diff_area / lag_x_perimeter_f,
      diff_area_pp_y = diff_area / lag_y_perimeter_f
    )
  return(df_lags_deltas)
}

# get the pseudo R2 metric using standard error estimates from bootstrapping
compute_pseudoR2 <- function(data, response, qr_model){
  qr_summary <- summary(qr_model, se = "boot", R = 50)
  pseudo_r2 <- 1 - sum(qr_summary$residuals^2) / sum((data[[response]] - median(data[[response]]))^2)
  coef <- qr_summary$coefficients
  p <- sprintf("%.20e", coef[, 4])
  return(list(pseudo_r2, p))
}

# get the pseudo R2 metric for non-linear models using standard error estimates from bootstrapping
compute_pseudoR2_nl <- function(obj, response){
  resid_fit <- residuals(obj)
  rho_tau <- function(u, tau) {
    u * (tau - (u < 0))
  }
  numerator <- sum(rho_tau(resid_fit, tau=.5))
  denominator <- sum(rho_tau(mdat[[response]] - median(mdat[[response]]), tau=.5))
  pseudo_R2 <- 1 - numerator / denominator
  return(pseudo_R2)
}

# get p values using standard error estimates from bootstrapping
get_pval <- function(data, qr_model){
  qr_summary <- summary(qr_model, se = "boot", R = 50) # Use bootstrapped standard errors
  coef <- qr_summary$coefficients
  p <- coef[2, 4]
}

# exponential decay function
exp_decay <- function(y0, yf, alpha, t) {
  y <- yf + (y0 - yf) * exp(-alpha * t)
  return(y)
}

compute_deviance <- function(fit, data, tau = 0.5) {
  if (is.null(fit)) return(NA)
  # Observed and predicted values
  y <- data$diff_area_norm_gdd_peri_y
  y_hat <- unlist(predict(fit, newdata = data))
  # Residuals
  residuals <- y - y_hat
  # Weighted absolute residuals
  weighted_residuals <- ifelse(residuals > 0, tau * abs(residuals), (1 - tau) * abs(residuals))
  # Sum of weighted residuals
  sum(weighted_residuals, na.rm = TRUE)
}


# non-linear quantile regression with grid of starting values
nls_quantile_exp <- function(data, x, y, tau = 0.5, n_samples = 300) {
  # get starting values
  start_samples <- data.frame(
    yf = runif(n_samples, min = 0, max = .01),
    alpha = runif(n_samples, min = 0, max = 1)
  ) %>% 
    mutate(y0 = runif(n_samples, min = yf, max = yf + .1))
  # specify model formula
  formula = as.formula(paste0(y, "~exp_decay(y0, yf, alpha, ", x, ")"))
  # fit using all starting values and formula
  results <- start_samples %>%
    rowwise() %>%
    mutate(
      fit = list(tryCatch(
        fit <- nlrq(formula, data = data, tau = tau,
                    start = list(y0 = y0, yf = yf, alpha = alpha)),
        error = function(e) NULL
      ))
    )
  # get parameters and filter 
  r <- results %>%
    filter(!is.null(fit)) %>% 
    mutate(deviance = compute_deviance(fit, data = data, tau = 0.5)) %>% 
    ungroup() %>% 
    mutate(pars = purrr::map(.x = fit, .f = get_qr_pars)) %>% unnest(pars) %>% 
    filter(alpha_fitted > 0) %>% 
    filter(y0_fitted >= 0) %>% 
    filter(alpha_fitted <= 1) %>% 
    filter(yf_fitted >= 0)
  # get best fit
  best_fit <- r %>% 
    slice_min(deviance) %>%
    pull(fit)
  # return
  if (length(best_fit) > 0){
    return(best_fit[[1]])
  } else{
    return(NA)
  }
}

# linear quantile regression 
linear_quantile <- function(data, x, y) {
  formula = as.formula(paste0(y, "~", x))
  fit <- rq(formula, data = data)
  return(fit)
}

# extract model parameters from model object
get_qr_pars <- function(obj){
  if (is.null(obj)) return(NA)
  coefs <- coef(obj)
  t <- tibble(!!!coef(obj))
  names(t) <- paste0(names(t), "_fitted")
  return(t)
}

