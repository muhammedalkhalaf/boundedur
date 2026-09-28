#' Unit Root Tests for Bounded Time Series
#'
#' Performs unit root tests for time series constrained within known bounds,
#' following the methodology of Cavaliere and Xu (2014). Provides modified
#' ADF and M-type test statistics with p-values computed via Monte Carlo
#' simulation of bounded Brownian motion.
#'
#' @param y Numeric vector. The time series to test.
#' @param lbound Numeric. Lower bound for the series.
#' @param ubound Numeric or \code{Inf}. Upper bound for the series. Use
#'   \code{Inf} for one-sided (lower bound only) testing.
#' @param test Character. Which test(s) to perform. One of:
#'   \itemize{
#'     \item \code{"all"}: All available tests (default)
#'     \item \code{"adf"}: Both ADF-alpha and ADF-t
#'     \item \code{"adf_alpha"}: ADF normalized bias test only
#'     \item \code{"adf_t"}: ADF t-statistic test only
#'     \item \code{"mz_alpha"}: MZ-alpha test only
#'     \item \code{"mz_t"}: MZ-t test only
#'     \item \code{"msb"}: MSB test only
#'   }
#' @param lags Integer or \code{NULL}. Number of lagged differences to
#'   include. If \code{NULL} (default), selected automatically using MAIC.
#' @param maxlag Integer or \code{NULL}. Maximum lag for automatic selection.
#'   If \code{NULL}, uses \code{floor(12 * (T/100)^0.25)}.
#' @param detrend Character. Detrending method:
#'   \itemize{
#'     \item \code{"constant"}: Demean the series (default)
#'     \item \code{"none"}: No detrending
#'   }
#' @param nsim Integer. Number of Monte Carlo replications for p-value
#'   computation. Default is 499.
#' @param nstep Integer or \code{NULL}. Number of discretization steps for
#'   Brownian motion simulation. If \code{NULL}, uses sample size T.
#' @param seed Integer or \code{NULL}. Random seed for reproducibility.
#'
#' @return An object of class \code{"boundedur"} containing:
#'   \item{statistics}{Named vector of test statistics}
#'   \item{p_values}{Named vector of p-values}
#'   \item{results}{Data frame with statistics, p-values, and decisions}
#'   \item{n}{Sample size}
#'   \item{lags}{Number of lags used}
#'   \item{lbound}{Lower bound}
#'   \item{ubound}{Upper bound}
#'   \item{c_lower}{Standardized lower bound parameter}
#'   \item{c_upper}{Standardized upper bound parameter}
#'   \item{sigma2_lr}{Autoregressive spectral estimate \eqn{s^2_{AR}(k)}
#'     of the long-run variance}
#'   \item{alpha1}{\eqn{\hat\alpha(1)}, one minus the sum of the lag
#'     coefficients}
#'   \item{x0}{Initial observation \eqn{X_0} used in the bound parameters}
#'   \item{detrend}{Detrending method used}
#'   \item{nsim}{Number of Monte Carlo replications}
#'   \item{call}{The matched call}
#'
#' @details
#' Standard unit root tests assume the series is unbounded, leading to
#' non-standard limiting distributions when bounds are present. This
#' function implements the bounded unit root tests of Cavaliere and Xu
#' (2014), which account for the effect of bounds on the limiting
#' distribution.
#'
#' The null hypothesis is that the series has a unit root while respecting
#' the bounds. The alternative is stationarity.
#'
#' @section Test Statistics:
#' \describe{
#'   \item{ADF-alpha}{\eqn{T(\hat{\rho} - 1)/\hat\alpha(1)} from the ADF
#'     regression of the de-meaned series (no deterministic terms)}
#'   \item{ADF-t}{The t-statistic for \eqn{\rho - 1 = 0} in that regression}
#'   \item{MZ-alpha}{\eqn{(T^{-1}\hat X_T^2 - T^{-1}\hat X_0^2 - s^2_{AR}) /
#'     (2T^{-2}\sum \hat X_{t-1}^2)}}
#'   \item{MZ-t}{\eqn{MZ_\alpha \times MSB}}
#'   \item{MSB}{\eqn{(T^{-2}\sum \hat X_{t-1}^2 / s^2_{AR})^{1/2}}}
#' }
#' All five reject for small values; ADF-alpha and MZ-alpha, and ADF-t and
#' MZ-t, share the same limiting distribution (Theorem 1).
#'
#' @section P-value Computation:
#' P-values follow Algorithm 1 of Cavaliere and Xu (2014): the bound
#' parameters are estimated by \eqn{\hat c = (b - X_0)/(s_{AR} T^{1/2})},
#' and the null distribution is simulated from a random walk regulated at
#' \eqn{\hat c} and \eqn{\bar c}. The number of replications (\code{nsim}) controls accuracy;
#' larger values give more precise p-values but increase computation time.
#'
#' @references
#' Cavaliere, G., & Xu, F. (2014). Testing for unit roots in bounded time
#' series. \emph{Journal of Econometrics}, 178(2), 259-272.
#' \doi{10.1016/j.jeconom.2013.08.026}
#'
#' Ng, S., & Perron, P. (2001). Lag length selection and the construction of
#' unit root tests with good size and power. \emph{Econometrica}, 69(6),
#' 1519-1554. \doi{10.1111/1468-0262.00256}
#'
#' @examples
#' # Generate bounded random walk (interest rate between 0 and 10)
#' set.seed(123)
#' n <- 200
#' y <- numeric(n)
#' y[1] <- 5
#' for (i in 2:n) {
#'   y[i] <- y[i-1] + rnorm(1, 0, 0.5)
#'   y[i] <- max(0, min(10, y[i]))  # Reflect at bounds
#' }
#'
#' # Test for unit root with known bounds
#' result <- boundedur(y, lbound = 0, ubound = 10, nsim = 199)
#' print(result)
#' summary(result)
#'
#' # One-sided bound (e.g., price level, lower bound = 0)
#' result_lower <- boundedur(y, lbound = 0, ubound = Inf, nsim = 199)
#'
#' @export
boundedur <- function(y, lbound, ubound = Inf,
                      test = c("all", "adf", "adf_alpha", "adf_t",
                               "mz_alpha", "mz_t", "msb"),
                      lags = NULL, maxlag = NULL,
                      detrend = c("constant", "none"),
                      nsim = 499, nstep = NULL, seed = NULL) {
  
  # Capture call
  cl <- match.call()
  
  # Validate inputs
  if (!is.numeric(y) || !is.vector(y)) {
    stop("'y' must be a numeric vector")
  }
  if (any(is.na(y))) {
    stop("'y' contains missing values")
  }
  
  n <- length(y)
  if (n < 10) {
    stop("Insufficient observations (minimum 10 required)")
  }
  
  # Validate bounds
  if (!is.numeric(lbound) || length(lbound) != 1) {
    stop("'lbound' must be a single numeric value")
  }
  if (!is.numeric(ubound) || length(ubound) != 1) {
    stop("'ubound' must be a single numeric value or Inf")
  }
  if (is.finite(ubound) && lbound >= ubound) {
    stop("'lbound' must be less than 'ubound'")
  }
  
  # Check series respects bounds
  y_min <- min(y)
  y_max <- max(y)
  if (y_min < lbound) {
    stop(sprintf("Series minimum (%.4f) is below lower bound (%.4f)", y_min, lbound))
  }
  if (is.finite(ubound) && y_max > ubound) {
    stop(sprintf("Series maximum (%.4f) is above upper bound (%.4f)", y_max, ubound))
  }
  
  # Match arguments
  test <- match.arg(test)
  detrend <- match.arg(detrend)
  
  # Set seed if specified
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Determine which tests to run
  run_adf_alpha <- test %in% c("all", "adf", "adf_alpha")
  run_adf_t <- test %in% c("all", "adf", "adf_t")
  run_mz_alpha <- test %in% c("all", "mz_alpha")
  run_mz_t <- test %in% c("all", "mz_t")
  run_msb <- test %in% c("all", "msb")
  
  # De-mean (OLS) as in Cavaliere and Xu (2014, Section 3)
  d <- .bur_detrend(y, detrend)

  # Lag selection
  if (is.null(lags)) {
    lag_result <- select_lag_maic(y, maxlag = maxlag, detrend = detrend)
    lags <- lag_result$selected_lag
  }
  lags <- as.integer(min(max(lags, 0L), n - 4L))

  # Set nstep default
  if (is.null(nstep)) {
    nstep <- n
  }

  # ADF regression and AR spectral long-run variance s^2_AR(k), eq. (3.8)
  adf_result <- compute_adf_stats(d, lags)
  sigma2_lr <- adf_result$sigma2_lr
  sigma_lr <- sqrt(sigma2_lr)

  # Standardized bound parameters, eq. (4.10): the initial value X_0 is
  # used; replacing it by the sample mean makes the estimators inconsistent
  # (Cavaliere and Xu, 2014, Remark 4.1).
  X0 <- y[1]
  c_lower <- (lbound - X0) / (sigma_lr * sqrt(n))
  c_upper <- if (is.finite(ubound)) (ubound - X0) / (sigma_lr * sqrt(n)) else Inf

  m_result <- compute_m_stats(d, sigma2_lr)

  # Monte Carlo null distribution (Algorithm 1). ADF-alpha and MZ-alpha share
  # one limit, ADF-t and MZ-t another; all reject for small values.
  S <- simulate_null_distribution(c_lower, c_upper, nsim, nstep)

  all_stats <- c(adf_alpha = unname(adf_result$adf_alpha),
                 adf_t = unname(adf_result$adf_t),
                 mz_alpha = m_result$mz_alpha,
                 mz_t = m_result$mz_t,
                 msb = m_result$msb)
  sim_col <- c(adf_alpha = "alpha", adf_t = "t", mz_alpha = "alpha",
               mz_t = "t", msb = "msb")
  all_p <- vapply(names(all_stats), function(nm) {
    mean(S[, sim_col[[nm]]] <= all_stats[[nm]])
  }, numeric(1))

  keep <- c(adf_alpha = run_adf_alpha, adf_t = run_adf_t,
            mz_alpha = run_mz_alpha, mz_t = run_mz_t, msb = run_msb)
  statistics <- all_stats[keep]
  p_values <- all_p[keep]

  # Create results data frame
  results <- data.frame(
    statistic = statistics,
    p_value = p_values,
    decision = ifelse(p_values < 0.05, "Reject H0", "Fail to reject"),
    row.names = names(statistics)
  )
  
  # Create output object
  out <- list(
    statistics = statistics,
    p_values = p_values,
    results = results,
    n = n,
    lags = lags,
    lbound = lbound,
    ubound = ubound,
    c_lower = c_lower,
    c_upper = c_upper,
    sigma2_lr = sigma2_lr,
    alpha1 = adf_result$alpha1,
    x0 = X0,
    detrend = detrend,
    nsim = nsim,
    nstep = nstep,
    call = cl
  )
  
  class(out) <- "boundedur"
  return(out)
}

#' @export
print.boundedur <- function(x, ...) {
  cat("\n")
  cat("Bounded Unit Root Tests\n")
  cat("=======================\n")
  cat("Reference: Cavaliere & Xu (2014), Journal of Econometrics 178, 259-272\n")
  cat("\n")
  cat("Sample size:", x$n, "\n")
  cat("Lags:", x$lags, "\n")
  
  if (is.finite(x$ubound)) {
    cat(sprintf("Bounds: [%.4f, %.4f]\n", x$lbound, x$ubound))
  } else {
    cat(sprintf("Lower bound: %.4f (one-sided)\n", x$lbound))
  }
  
  cat(sprintf("Standardized c_lower: %.4f\n", x$c_lower))
  if (is.finite(x$c_upper)) {
    cat(sprintf("Standardized c_upper: %.4f\n", x$c_upper))
  }
  cat(sprintf("Long-run variance: %.6f\n", x$sigma2_lr))
  cat("\n")
  cat("Test Results:\n")
  cat("-------------\n")
  
  for (i in seq_len(nrow(x$results))) {
    test_name <- rownames(x$results)[i]
    stat <- x$results$statistic[i]
    pval <- x$results$p_value[i]
    dec <- x$results$decision[i]
    cat(sprintf("%-12s %10.4f    p = %.4f    %s\n",
                test_name, stat, pval, dec))
  }
  
  cat("\n")
  cat("H0: Unit root (series is I(1) with bounds)\n")
  cat(sprintf("p-values based on %d Monte Carlo replications\n", x$nsim))
  
  invisible(x)
}

#' @export
summary.boundedur <- function(object, ...) {
  cat("\n")
  cat("Summary: Bounded Unit Root Tests\n")
  cat("================================\n\n")
  
  cat("Call:\n")
  print(object$call)
  cat("\n")
  
  cat("Test Configuration:\n")
  cat(sprintf("  Sample size (T):      %d\n", object$n))
  cat(sprintf("  Lags (MAIC):          %d\n", object$lags))
  cat(sprintf("  Detrending:           %s\n", object$detrend))
  cat(sprintf("  MC replications:      %d\n", object$nsim))
  cat(sprintf("  Discretization steps: %d\n", object$nstep))
  cat("\n")
  
  cat("Bounds:\n")
  if (is.finite(object$ubound)) {
    cat(sprintf("  Original:     [%.4f, %.4f]\n", object$lbound, object$ubound))
    cat(sprintf("  Standardized: [%.4f, %.4f]\n", object$c_lower, object$c_upper))
  } else {
    cat(sprintf("  Lower bound:  %.4f (one-sided)\n", object$lbound))
    cat(sprintf("  Standardized: %.4f\n", object$c_lower))
  }
  cat("\n")
  
  cat("Variance Estimation:\n")
  cat(sprintf("  Long-run variance: %.6f\n", object$sigma2_lr))
  cat(sprintf("  alpha(1):          %.6f\n", object$alpha1))
  cat("\n")
  
  cat("Test Results:\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")
  cat(sprintf("%-12s %12s %12s %16s\n",
              "Test", "Statistic", "p-value", "Decision (5%)"))
  cat(paste(rep("-", 60), collapse = ""), "\n")
  
  for (i in seq_len(nrow(object$results))) {
    test_name <- rownames(object$results)[i]
    stat <- object$results$statistic[i]
    pval <- object$results$p_value[i]
    dec <- object$results$decision[i]
    cat(sprintf("%-12s %12.4f %12.4f %16s\n", test_name, stat, pval, dec))
  }
  
  cat(paste(rep("-", 60), collapse = ""), "\n")
  cat("\n")
  cat("H0: Unit root (series is I(1) with bounds)\n")
  cat("H1: Stationary (series is I(0))\n")
  
  invisible(object)
}
