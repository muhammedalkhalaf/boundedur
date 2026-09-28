# boundedur 1.0.3

* Bug fix: the bound parameters were computed with the sample mean in place of the initial observation X_0. Cavaliere and Xu (2014, eq. 4.10 and Remark 4.1) use X_0 and show that the mean makes the estimators inconsistent.
* Bug fix: the ADF statistics were taken from a regression with a constant on the raw series, and ADF-alpha omitted the division by alpha(1). They now come from the regression of the de-meaned series without deterministic terms, with ADF-alpha = T * pi / alpha(1), as in equation (3.7).
* Bug fix: MZ-alpha and MZ-t omitted the -X_0^2 / T term of the numerator, so they did not share the limiting distribution of ADF-alpha and ADF-t.
* Bug fix: the MSB p-value was taken in the wrong tail; MSB rejects for small values.
* Bug fix: the Monte Carlo null distribution now follows Algorithm 1 (a random walk regulated at the estimated bounds, de-meaned before the functionals are formed); the previous version used mirror reflection and did not de-mean.
* Bug fix: the MAIC lag selection omitted the tau_T(k) term of Ng and Perron (2001) and used varying samples; it now uses the full MAIC on a common sample.
* All five statistics and both bound parameters agree to four decimals with the Stata command boundedur (SSC) on the same data.

# boundedur 1.0.2

* Corrected the DOI of Cavaliere and Xu (2014) to 10.1016/j.jeconom.2013.08.026 in DESCRIPTION, README and all R and Rd files.

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
