#' Simulate Bounded Brownian Motion
#'
#' Simulates the discretized regulated Brownian motion of Algorithm 1 in
#' Cavaliere and Xu (2014).
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
#' The path follows the recursion (4.11) of Cavaliere and Xu (2014):
#' \eqn{X_t = X_{t-1} + n^{-1/2}\varepsilon_t}, set to the bound whenever it
#' would cross it, with \eqn{X_0 = 0} and i.i.d. standard normal
#' \eqn{\varepsilon_t}.
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
#' \doi{10.1016/j.jeconom.2013.08.026}
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
  .bur_path(stats::rnorm(n), c_lower, c_upper)
}


#' Regulated random walk (internal)
#' @keywords internal
#' @noRd
.bur_path <- function(eps, c_lower, c_upper) {
  n <- length(eps)
  sq <- sqrt(n)
  X <- numeric(n + 1L)
  for (t in seq_len(n)) {
    x <- X[t] + eps[t] / sq
    if (x > c_upper) x <- c_upper
    if (x < c_lower) x <- c_lower
    X[t + 1L] <- x
  }
  X
}


#' Monte Carlo Null Distribution (Algorithm 1, steps i to iii)
#'
#' @return Matrix with columns \code{alpha} (limit of ADF-alpha and
#'   MZ-alpha), \code{t} (limit of ADF-t and MZ-t) and \code{msb}.
#' @keywords internal
#' @noRd
simulate_null_distribution <- function(c_lower, c_upper, nsim, nstep) {
  S <- matrix(NA_real_, nsim, 3L, dimnames = list(NULL, c("alpha", "t", "msb")))
  for (b in seq_len(nsim)) {
    X <- .bur_path(stats::rnorm(nstep), c_lower, c_upper)
    Xt <- X - mean(X)
    A <- mean(Xt[-1L]^2)
    La <- (Xt[nstep + 1L]^2 - Xt[1L]^2 - 1) / (2 * A)
    S[b, ] <- c(La, La * sqrt(A), sqrt(A))
  }
  S
}
