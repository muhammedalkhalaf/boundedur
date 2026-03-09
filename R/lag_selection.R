#' Select Optimal Lag using MAIC Criterion
#'
#' Selects the optimal number of lags for the ADF regression using the
#' Modified Akaike Information Criterion (MAIC) of Ng and Perron (2001).
#'
#' @param y Numeric vector. Time series data.
#' @param maxlag Integer or \code{NULL}. Maximum lag to consider. If
#'   \code{NULL}, uses the rule \code{floor(12 * (T/100)^0.25)}.
#' @param detrend Character. Detrending method: "constant" (demean) or
#'   "none". Default is "constant".
#'
#' @return A list with class \code{"lag_selection"} containing:
#'   \item{selected_lag}{Optimal lag selected by MAIC}
#'   \item{maic}{MAIC value at optimal lag}
#'   \item{all_maic}{Vector of MAIC values for all lags}
#'   \item{maxlag}{Maximum lag considered}
#'
#' @details
#' The MAIC criterion is defined as:
#' \deqn{MAIC(k) = \ln(\hat{\sigma}^2_k) + 2(k+1)/T}
#' where \eqn{\hat{\sigma}^2_k} is the residual variance from the ADF
#' regression with \eqn{k} lags.
#'
#' This criterion provides better size properties than standard AIC for
#' unit root testing.
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
  # Validate inputs
  if (!is.numeric(y) || !is.vector(y)) {
    stop("'y' must be a numeric vector")
  }
  
  n <- length(y)
  if (n < 10) {
    stop("Insufficient observations (minimum 10 required)")
  }
  
  detrend <- match.arg(detrend, c("constant", "none"))
  
  # Default maxlag (Ng & Perron 2001)
  if (is.null(maxlag)) {
    maxlag <- floor(12 * (n / 100)^0.25)
  }
  maxlag <- min(maxlag, floor(n / 3))  # Safety limit
  
  # First difference
  dy <- diff(y)
  
  # Storage for MAIC values
  maic_values <- numeric(maxlag + 1)
  
  for (k in 0:maxlag) {
    if (k > 0) {
      # Embed lagged differences
      dy_embed <- stats::embed(dy, k + 1)
      dy_dep <- dy_embed[, 1]
      dy_lags <- dy_embed[, -1, drop = FALSE]
      y_lag1 <- y[(k + 1):(n - 1)]
      
      if (detrend == "constant") {
        X <- cbind(1, y_lag1, dy_lags)
      } else {
        X <- cbind(y_lag1, dy_lags)
      }
    } else {
      dy_dep <- dy
      y_lag1 <- y[-n]
      
      if (detrend == "constant") {
        X <- cbind(1, y_lag1)
      } else {
        X <- cbind(y_lag1)
      }
    }
    
    N_k <- length(dy_dep)
    
    # OLS
    fit <- stats::lm.fit(X, dy_dep)
    sigma2_k <- sum(fit$residuals^2) / N_k
    
    # MAIC criterion
    tau_k <- 2 * (k + 1) / N_k
    maic_values[k + 1] <- log(sigma2_k) + tau_k
  }
  
  # Select minimum
  best_lag <- which.min(maic_values) - 1
  
  result <- list(
    selected_lag = best_lag,
    maic = maic_values[best_lag + 1],
    all_maic = maic_values,
    maxlag = maxlag,
    n = n
  )
  
  class(result) <- "lag_selection"
  return(result)
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
