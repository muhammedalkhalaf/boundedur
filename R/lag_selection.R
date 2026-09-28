#' Select Optimal Lag using MAIC Criterion
#'
#' Selects the number of lagged differences in the ADF regression with the
#' Modified Akaike Information Criterion (MAIC) of Ng and Perron (2001).
#'
#' @param y Numeric vector. Time series data.
#' @param maxlag Integer or \code{NULL}. Maximum lag to consider. If
#'   \code{NULL}, uses the rule \code{floor(12 * (T/100)^0.25)}.
#' @param detrend Character. De-meaning method: "constant" (OLS de-meaning)
#'   or "none". Default is "constant".
#'
#' @return A list with class \code{"lag_selection"} containing:
#'   \item{selected_lag}{Optimal lag selected by MAIC}
#'   \item{maic}{MAIC value at optimal lag}
#'   \item{all_maic}{Vector of MAIC values for all lags}
#'   \item{maxlag}{Maximum lag considered}
#'   \item{n}{Sample size}
#'
#' @details
#' For \eqn{k = 0, \ldots, k_{max}}, the ADF regression of
#' \eqn{\Delta \hat X_t} on \eqn{\hat X_{t-1}} and \eqn{k} lagged differences
#' is estimated on the common sample \eqn{t = k_{max} + 2, \ldots, T}, and
#' \deqn{MAIC(k) = \ln(\hat{\sigma}^2_k) + 2(\tau_T(k) + k)/(T - k_{max} - 1),}
#' with \eqn{\tau_T(k) = \hat\sigma_k^{-2} \hat\beta_0^2 \sum \hat X_{t-1}^2}
#' and \eqn{\hat\beta_0} the coefficient on \eqn{\hat X_{t-1}}.
#'
#' @references
#' Ng, S., & Perron, P. (2001). Lag length selection and the construction of
#' unit root tests with good size and power. \emph{Econometrica}, 69(6),
#' 1519-1554. \doi{10.1111/1468-0262.00256}
#'
#' @examples
#' # Generate random walk
#' set.seed(123)
#' y <- cumsum(rnorm(200))
#'
#' # Select lag
#' lag_sel <- select_lag_maic(y)
#' print(lag_sel)
#'
#' @export
select_lag_maic <- function(y, maxlag = NULL, detrend = "constant") {
  if (!is.numeric(y) || !is.vector(y)) {
    stop("'y' must be a numeric vector")
  }
  n <- length(y)
  if (n < 10) {
    stop("Insufficient observations (minimum 10 required)")
  }
  detrend <- match.arg(detrend, c("constant", "none"))
  if (is.null(maxlag)) {
    maxlag <- floor(12 * (n / 100)^0.25)
  }
  maxlag <- max(0L, min(as.integer(maxlag), n - 6L))

  d <- .bur_detrend(y, detrend)
  tt <- (maxlag + 2L):n
  nef <- length(tt)
  dep <- d[tt] - d[tt - 1L]
  R <- matrix(d[tt - 1L], ncol = 1L)
  if (maxlag > 0L) {
    for (i in seq_len(maxlag)) R <- cbind(R, d[tt - i] - d[tt - i - 1L])
  }
  sumy <- sum(R[, 1]^2)
  maic_values <- numeric(maxlag + 1L)
  for (k in 0:maxlag) {
    Xk <- R[, seq_len(k + 1L), drop = FALSE]
    fit <- stats::lm.fit(Xk, dep)
    s2 <- sum(fit$residuals^2) / nef
    tau_k <- fit$coefficients[1]^2 * sumy / s2
    maic_values[k + 1L] <- log(s2) + 2 * (k + tau_k) / nef
  }
  best_lag <- which.min(maic_values) - 1L

  result <- list(
    selected_lag = best_lag,
    maic = maic_values[best_lag + 1L],
    all_maic = maic_values,
    maxlag = maxlag,
    n = n
  )
  class(result) <- "lag_selection"
  result
}

#' @export
print.lag_selection <- function(x, ...) {
  cat("MAIC Lag Selection\n")
  cat("==================\n")
  cat("Sample size:", x$n, "\n")
  cat("Max lag considered:", x$maxlag, "\n")
  cat("Selected lag:", x$selected_lag, "\n")
  cat("MAIC value:", format(x$maic, digits = 4), "\n")
  invisible(x)
}
