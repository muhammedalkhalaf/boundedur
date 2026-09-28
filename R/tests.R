#' De-mean a Series for the Bounded Unit Root Tests
#'
#' @param y Numeric vector.
#' @param detrend Character: "constant" (OLS de-meaning) or "none".
#' @return The de-meaned (or original) series.
#' @keywords internal
#' @noRd
.bur_detrend <- function(y, detrend) {
  if (detrend == "constant") y - mean(y) else y
}


#' ADF Regression on the De-meaned Series
#'
#' Regresses \eqn{\Delta \hat X_t} on \eqn{\hat X_{t-1}} and \eqn{k} lagged
#' differences, without deterministic terms, as in equation (3.7) of
#' Cavaliere and Xu (2014).
#'
#' @param d De-meaned series.
#' @param k Number of lagged differences.
#' @return List with the ADF statistics, \eqn{\hat\alpha(1)}, the AR
#'   long-run variance \eqn{s^2_{AR}(k)} and the lag coefficients.
#' @keywords internal
#' @noRd
compute_adf_stats <- function(d, k) {
  TT <- length(d)
  tt <- (k + 2L):TT
  dep <- d[tt] - d[tt - 1L]
  X <- matrix(d[tt - 1L], ncol = 1L)
  if (k > 0L) {
    for (i in seq_len(k)) X <- cbind(X, d[tt - i] - d[tt - i - 1L])
  }
  XtX_inv <- solve(crossprod(X))
  b <- as.numeric(XtX_inv %*% crossprod(X, dep))
  e <- dep - as.numeric(X %*% b)
  s2 <- sum(e^2) / length(dep)
  se <- sqrt(XtX_inv[1, 1] * s2)
  pi_hat <- b[1]
  lagcoef <- if (k > 0L) b[-1] else numeric(0)
  alpha1 <- 1 - sum(lagcoef)
  list(adf_alpha = TT * pi_hat / alpha1,
       adf_t = pi_hat / se,
       rho = pi_hat,
       alpha1 = alpha1,
       sigma2 = s2,
       sigma2_lr = s2 / alpha1^2,
       lag_coefs = lagcoef)
}


#' M Statistics with the Initial-Value Term
#'
#' \eqn{MZ_\alpha = (T^{-1}\hat X_T^2 - T^{-1}\hat X_0^2 - s^2_{AR}) /
#' (2 T^{-2} \sum \hat X_{t-1}^2)}, \eqn{MSB = (T^{-2} \sum \hat X_{t-1}^2 /
#' s^2_{AR})^{1/2}} and \eqn{MZ_t = MZ_\alpha \times MSB}, as in
#' Cavaliere and Xu (2014, Section 3).
#'
#' @param d De-meaned series.
#' @param sigma2_lr Long-run variance \eqn{s^2_{AR}(k)}.
#' @keywords internal
#' @noRd
compute_m_stats <- function(d, sigma2_lr) {
  TT <- length(d)
  sslag <- sum(d[-TT]^2)
  mz_alpha <- (d[TT]^2 / TT - d[1]^2 / TT - sigma2_lr) / (2 * sslag / TT^2)
  msb <- sqrt(sslag / (TT^2 * sigma2_lr))
  list(mz_alpha = mz_alpha, mz_t = mz_alpha * msb, msb = msb)
}
