#' Simulate Bounded Brownian Motion
#'
#' Simulates a discretized reflected Brownian motion constrained between
#' bounds, following the methodology of Cavaliere and Xu (2014).
#'
#' @param n Integer. Number of time steps for discretization.
#' @param c_lower Numeric. Standardized lower bound parameter.
#' @param c_upper Numeric or \code{Inf}. Standardized upper bound parameter.
#'   Use \code{Inf} for one-sided (lower) bound only.
#'
#' @return A numeric vector of length \code{n + 1} containing the simulated
#'   bounded Brownian motion path, starting at 0.
#'
#' @details
#' The function simulates a standard Brownian motion and applies reflection

#' at the boundaries. For two-sided bounds, both upper and lower reflections
#' are applied. For one-sided bounds (\code{c_upper = Inf}), only lower
#' reflection is used.
#'
#' The standardized bound parameters \code{c_lower} and \code{c_upper} are
#' computed from the original bounds as:
#' \deqn{c = (b - X_0) / (\sigma \sqrt{T})}
#' where \eqn{b} is the bound, \eqn{X_0} is the initial value, \eqn{\sigma}
#' is the long-run standard deviation, and \eqn{T} is the sample size.
#'
#' @references
#' Cavaliere, G., & Xu, F. (2014). Testing for unit roots in bounded time
#' series. \emph{Journal of Econometrics}, 178(2), 259-272.
#' \doi{10.1016/j.jeconom.2013.08.012}
#'
#' @examples
#' # Simulate bounded Brownian motion with two-sided bounds
#' set.seed(123)
#' bm <- simulate_bounded_bm(n = 1000, c_lower = -2, c_upper = 2)
#' plot(bm, type = "l", main = "Bounded Brownian Motion")
#' abline(h = c(-2, 2), col = "red", lty = 2)
#'
#' # One-sided bound (lower only)
#' bm_lower <- simulate_bounded_bm(n = 1000, c_lower = -1, c_upper = Inf)
#'
#' @export
simulate_bounded_bm <- function(n, c_lower, c_upper = Inf) {
  # Validate inputs
  if (!is.numeric(n) || length(n) != 1 || n < 1 || n != floor(n)) {
    stop("'n' must be a positive integer")
  }
  if (!is.numeric(c_lower) || length(c_lower) != 1) {
    stop("'c_lower' must be a single numeric value")
  }
  if (!is.numeric(c_upper) || length(c_upper) != 1) {
    stop("'c_upper' must be a single numeric value or Inf")
  }
  if (is.finite(c_upper) && c_lower >= c_upper) {
    stop("'c_lower' must be less than 'c_upper'")
  }
  
  # Generate standard Brownian motion increments
  dW <- stats::rnorm(n, mean = 0, sd = 1 / sqrt(n))
  
  # Initialize path
  W <- numeric(n + 1)
  W[1] <- 0
  
  # Simulate with reflection at boundaries
  for (i in seq_len(n)) {
    W[i + 1] <- W[i] + dW[i]
    
    # Apply reflection (Skorokhod reflection)
    if (is.finite(c_upper)) {
      # Two-sided reflection
      while (W[i + 1] < c_lower || W[i + 1] > c_upper) {
        if (W[i + 1] < c_lower) {
          W[i + 1] <- 2 * c_lower - W[i + 1]
        }
        if (W[i + 1] > c_upper) {
          W[i + 1] <- 2 * c_upper - W[i + 1]
        }
      }
    } else {
      # One-sided reflection (lower bound only)
      while (W[i + 1] < c_lower) {
        W[i + 1] <- 2 * c_lower - W[i + 1]
      }
    }
  }
  
  return(W)
}

#' Compute Functionals of Bounded Brownian Motion
#'
#' Internal function to compute test statistics from simulated bounded
#' Brownian motion paths.
#'
#' @param W Numeric vector. Simulated bounded Brownian motion path.
#' @param test Character. Type of test statistic: "adf_alpha", "adf_t",
#'   "mz_alpha", "mz_t", or "msb".
#'
#' @return Numeric value of the test statistic functional.
#'
#' @keywords internal
#' @noRd
compute_bm_functional <- function(W, test) {
  n <- length(W) - 1
  dt <- 1 / n
  
  # Grid points for integration
  t_grid <- seq(0, 1, length.out = n + 1)
  
  # W(1) - final value
  W1 <- W[n + 1]
  
  # Integral of W^2 dt (using trapezoidal rule)
  int_W2 <- sum((W[-1]^2 + W[-(n + 1)]^2) / 2 * dt)
  
  # Integral of W dW (Ito integral approximation)
  int_WdW <- sum(W[-(n + 1)] * diff(W))
  
  # Test statistics based on bounded Brownian motion functionals
  # Following Cavaliere & Xu (2014) Equations 2.6-2.10
  switch(test,
    "adf_alpha" = {
      # ADF normalized bias: T(rho - 1) -> int W dW / int W^2 dt
      int_WdW / int_W2
    },
    "adf_t" = {
      # ADF t-statistic
      int_WdW / sqrt(int_W2)
    },
    "mz_alpha" = {
      # MZ_alpha: (W(1)^2 - 1) / (2 * int W^2 dt)
      (W1^2 - 1) / (2 * int_W2)
    },
    "mz_t" = {
      # MZ_t: (W(1)^2 - 1) / (2 * sqrt(int W^2 dt))
      (W1^2 - 1) / (2 * sqrt(int_W2))
    },
    "msb" = {
      # MSB: sqrt(int W^2 dt)
      sqrt(int_W2)
    },
    stop("Unknown test type: ", test)
  )
}

#' Simulate Critical Values Distribution
#'
#' Internal function to generate the null distribution of test statistics
#' via Monte Carlo simulation of bounded Brownian motion.
#'
#' @param c_lower Numeric. Standardized lower bound.
#' @param c_upper Numeric or Inf. Standardized upper bound.
#' @param test Character. Test type.
#' @param nsim Integer. Number of Monte Carlo replications.
#' @param nstep Integer. Number of discretization steps.
#'
#' @return Numeric vector of simulated test statistics under the null.
#'
#' @keywords internal
#' @noRd
simulate_null_distribution <- function(c_lower, c_upper, test, nsim, nstep) {
  stats <- numeric(nsim)
  
  for (i in seq_len(nsim)) {
    W <- simulate_bounded_bm(nstep, c_lower, c_upper)
    stats[i] <- compute_bm_functional(W, test)
  }
  
  return(stats)
}
