#' Compute ADF Test Statistics
#'
#' Internal function to compute ADF test statistics from data.
#'
#' @param y Numeric vector. Time series data.
#' @param lags Integer. Number of lagged differences to include.
#' @param detrend Character. Detrending method: "constant" or "none".
#'
#' @return List containing:
#'   \item{adf_alpha}{Normalized bias statistic T(rho-1)}
#'   \item{adf_t}{t-statistic for rho}
#'   \item{rho}{Estimated rho coefficient}
#'   \item{residuals}{Regression residuals}
#'   \item{sigma2}{Residual variance}
#'
#' @keywords internal
#' @noRd
compute_adf_stats <- function(y, lags = 0, detrend = "constant") {
  n <- length(y)
  
  # First difference
  dy <- diff(y)
  
  # Lagged level
  y_lag1 <- y[-n]
  
  # Build design matrix
  if (lags > 0) {
    # Embed lagged differences
    dy_embed <- stats::embed(dy, lags + 1)
    dy_dep <- dy_embed[, 1]
    dy_lags <- dy_embed[, -1, drop = FALSE]
    
    # Align y_lag1
    y_lag1 <- y_lag1[(lags + 1):length(y_lag1)]
    
    # Design matrix
    if (detrend == "constant") {
      X <- cbind(1, y_lag1, dy_lags)
    } else {
      X <- cbind(y_lag1, dy_lags)
    }
  } else {
    # For lags=0: dy has length n-1, y_lag1 has length n-1
    # We need dy[t] ~ y[t-1]
    dy_dep <- dy
    y_lag1 <- y[-n]
    
    if (detrend == "constant") {
      X <- cbind(1, y_lag1)
    } else {
      X <- matrix(y_lag1, ncol = 1)
    }
  }
  
  # Effective sample size
  T_eff <- length(dy_dep)
  
  # OLS regression
  fit <- stats::lm.fit(X, dy_dep)
  
  # Extract coefficient on y_lag1
  if (detrend == "constant") {
    rho_hat <- fit$coefficients[2]
    rho_idx <- 2
  } else {
    rho_hat <- fit$coefficients[1]
    rho_idx <- 1
  }
  
  # Residuals and variance
  resid <- fit$residuals
  sigma2 <- sum(resid^2) / T_eff
  
  # Standard error of rho_hat
  XtX_inv <- solve(crossprod(X))
  se_rho <- sqrt(sigma2 * XtX_inv[rho_idx, rho_idx])
  
  # Test statistics
  adf_alpha <- T_eff * rho_hat  # Normalized bias
  adf_t <- rho_hat / se_rho      # t-statistic
  
  list(
    adf_alpha = adf_alpha,
    adf_t = adf_t,
    rho = rho_hat,
    residuals = resid,
    sigma2 = sigma2,
    T_eff = T_eff
  )
}

#' Compute M-type Test Statistics
#'
#' Internal function to compute M-type test statistics (MZ-alpha, MZ-t, MSB)
#' following Ng and Perron (2001) and modified for bounded series.
#'
#' @param y Numeric vector. Time series data.
#' @param sigma2_lr Numeric. Long-run variance estimate.
#' @param detrend Character. Detrending method: "constant" or "none".
#'
#' @return List containing:
#'   \item{mz_alpha}{MZ-alpha statistic}
#'   \item{mz_t}{MZ-t statistic}
#'   \item{msb}{MSB statistic}
#'
#' @keywords internal
#' @noRd
compute_m_stats <- function(y, sigma2_lr, detrend = "constant") {
  n <- length(y)
  
  # Detrend the series
  if (detrend == "constant") {
    y_det <- y - mean(y)
  } else {
    y_det <- y
  }
  
  # Compute required quantities
  # s^2 = (1/T^2) * sum(y_{t-1}^2)
  s2 <- sum(y_det[-n]^2) / n^2
  
  # For M statistics, we need:
  # MZ_alpha = (y_T^2 / T - sigma^2_lr) / (2 * s^2)
  # MZ_t = (y_T^2 / T - sigma^2_lr) / (2 * s * sigma_lr)
  # MSB = s / sigma_lr
  
  y_T <- y_det[n]
  sigma_lr <- sqrt(sigma2_lr)
  s <- sqrt(s2)
  
  mz_alpha <- (y_T^2 / n - sigma2_lr) / (2 * s2)
  mz_t <- (y_T^2 / n - sigma2_lr) / (2 * s * sigma_lr)
  msb <- s / sigma_lr
  
  list(
    mz_alpha = mz_alpha,
    mz_t = mz_t,
    msb = msb
  )
}

#' Estimate Long-Run Variance
#'
#' Internal function to estimate long-run variance using the spectral AR
#' method following Ng and Perron (2001).
#'
#' @param y Numeric vector. Time series data.
#' @param lags Integer. Number of lags in AR approximation.
#' @param detrend Character. Detrending method.
#'
#' @return List containing:
#'   \item{sigma2_lr}{Long-run variance estimate}
#'   \item{alpha1}{Sum of AR coefficients alpha(1)}
#'   \item{sigma2}{Innovation variance}
#'
#' @keywords internal
#' @noRd
estimate_lrvar <- function(y, lags = 0, detrend = "constant") {
  n <- length(y)
  dy <- diff(y)
  
  if (lags > 0) {
    # Embed for lagged differences
    dy_embed <- stats::embed(dy, lags + 1)
    dy_dep <- dy_embed[, 1]
    dy_lags <- dy_embed[, -1, drop = FALSE]
    y_lag1 <- y[(lags + 1):(n - 1)]
    
    if (detrend == "constant") {
      X <- cbind(1, y_lag1, dy_lags)
    } else {
      X <- cbind(y_lag1, dy_lags)
    }
    
    fit <- stats::lm.fit(X, dy_dep)
    
    # alpha(1) = 1 - sum of lag coefficients
    if (detrend == "constant") {
      lag_coefs <- fit$coefficients[-(1:2)]
    } else {
      lag_coefs <- fit$coefficients[-1]
    }
    alpha1 <- 1 - sum(lag_coefs)
    
  } else {
    y_lag1 <- y[-n]
    
    if (detrend == "constant") {
      X <- cbind(1, y_lag1)
    } else {
      X <- cbind(y_lag1)
    }
    
    fit <- stats::lm.fit(X, dy)
    alpha1 <- 1
  }
  
  # Residual variance
  resid <- fit$residuals
  sigma2 <- sum(resid^2) / length(resid)
  
  # Long-run variance
  sigma2_lr <- sigma2 / (alpha1^2)
  
  list(
    sigma2_lr = sigma2_lr,
    alpha1 = alpha1,
    sigma2 = sigma2
  )
}
