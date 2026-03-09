# boundedur

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/boundedur)](https://CRAN.R-project.org/package=boundedur)
<!-- badges: end -->

Unit root tests for bounded time series following Cavaliere and Xu (2014).

## Overview

Standard unit root tests (ADF, Phillips-Perron) assume the time series is unbounded, but many economic variables are naturally bounded:

- Interest rates (non-negative, often with upper caps)
- Exchange rates (within currency bands)
- Capacity utilization (0-100%)
- Unemployment rates (0-100%)
- Price indices (non-negative)

When bounds are binding or nearly so, standard tests have non-standard limiting distributions, leading to incorrect inference. The **boundedur** package implements the modified tests of Cavaliere and Xu (2014), which properly account for bounds.

## Installation

```r
# From CRAN (when available)
install.packages("boundedur")

# Development version from GitHub
# install.packages("remotes")
```

## Usage

```r
library(boundedur)

# Generate bounded random walk (e.g., interest rate between 0 and 10%)
set.seed(123)
n <- 200
y <- numeric(n)
y[1] <- 5
for (i in 2:n) {
  y[i] <- y[i-1] + rnorm(1, 0, 0.3)
  y[i] <- max(0, min(10, y[i]))  # Reflect at bounds
}

# Test for unit root with known bounds
result <- boundedur(y, lbound = 0, ubound = 10)
print(result)
summary(result)
```

### One-sided bounds

For series with only a lower bound (e.g., prices):

```r
result <- boundedur(y, lbound = 0, ubound = Inf)
```

### Select specific tests

```r
# ADF tests only
result <- boundedur(y, lbound = 0, ubound = 10, test = "adf")

# Single test
result <- boundedur(y, lbound = 0, ubound = 10, test = "mz_t")
```

### Control lag selection

```r
# Automatic (MAIC criterion, default)
result <- boundedur(y, lbound = 0, ubound = 10)

# Manual lag specification
result <- boundedur(y, lbound = 0, ubound = 10, lags = 4)

# Custom maximum lag
result <- boundedur(y, lbound = 0, ubound = 10, maxlag = 12)
```

## Available Tests

| Test | Description |
|------|-------------|
| `adf_alpha` | ADF normalized bias: T(ρ̂ - 1) |
| `adf_t` | ADF t-statistic |
| `mz_alpha` | Modified Phillips-Perron normalized bias |
| `mz_t` | Modified Phillips-Perron t-statistic |
| `msb` | Modified Sargan-Bhargava |

## Methodology

The tests modify standard unit root statistics to account for the effect of bounds on the limiting distribution. Under the null hypothesis of a unit root, the process is constrained by bounds and follows a **reflected Brownian motion** rather than standard Brownian motion.

P-values are computed via Monte Carlo simulation of the bounded Brownian motion null distribution, with bound parameters estimated from the data.

## References

Cavaliere, G., & Xu, F. (2014). Testing for unit roots in bounded time series. *Journal of Econometrics*, 178(2), 259-272. [doi:10.1016/j.jeconom.2013.08.012](https://doi.org/10.1016/j.jeconom.2013.08.012)

Ng, S., & Perron, P. (2001). Lag length selection and the construction of unit root tests with good size and power. *Econometrica*, 69(6), 1519-1554. [doi:10.1111/1468-0262.00256](https://doi.org/10.1111/1468-0262.00256)

## License

GPL (>= 3)
