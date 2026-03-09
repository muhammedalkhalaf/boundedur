# boundedur 1.0.0

* Initial CRAN release.

## Features

* `boundedur()`: Main function for bounded unit root tests
  - ADF-alpha and ADF-t tests
  - M-type tests (MZ-alpha, MZ-t, MSB)
  - Monte Carlo p-value computation
  - Support for one-sided and two-sided bounds

* `select_lag_maic()`: MAIC lag selection (Ng & Perron, 2001)

* `simulate_bounded_bm()`: Bounded Brownian motion simulation

## References

* Cavaliere, G., & Xu, F. (2014). Testing for unit roots in bounded time series. Journal of Econometrics, 178(2), 259-272.
