test_that("boundedur basic functionality works", {
  set.seed(123)
  
  # Generate bounded random walk
  n <- 100
  y <- numeric(n)
  y[1] <- 5
  for (i in 2:n) {
    y[i] <- y[i-1] + rnorm(1, 0, 0.5)
    y[i] <- max(0, min(10, y[i]))
  }
  
  # Run test
  result <- boundedur(y, lbound = 0, ubound = 10, nsim = 99, seed = 42)
  
  expect_s3_class(result, "boundedur")
  expect_equal(result$n, n)
  expect_equal(result$lbound, 0)
  expect_equal(result$ubound, 10)
  expect_true(result$lags >= 0)
  expect_true(all(result$p_values >= 0 & result$p_values <= 1))
})

test_that("boundedur handles one-sided bounds", {
  set.seed(456)
  
  n <- 80
  y <- cumsum(rnorm(n))
  y <- pmax(y, 0)  # Reflect at lower bound
  
  result <- boundedur(y, lbound = 0, ubound = Inf, nsim = 99, seed = 42)
  
  expect_s3_class(result, "boundedur")
  expect_equal(result$ubound, Inf)
  expect_true(is.infinite(result$c_upper))
})

test_that("boundedur validates inputs correctly", {
  y <- cumsum(rnorm(50))
  
  # Missing lower bound

  expect_error(boundedur(y), "argument \"lbound\" is missing")
  
  # Invalid bounds
  expect_error(boundedur(y, lbound = 10, ubound = 5),
               "must be less than")
  
  # Series outside bounds
  expect_error(boundedur(y, lbound = max(y) + 1, ubound = max(y) + 10),
               "below lower bound")
  
  # Too few observations
  expect_error(boundedur(1:5, lbound = 0, ubound = 10),
               "Insufficient observations")
  
  # NA values
  y_na <- c(1:10, NA, 12:20)
  expect_error(boundedur(y_na, lbound = 0, ubound = 25),
               "missing values")
})

test_that("select_lag_maic works correctly", {
  set.seed(789)
  y <- cumsum(rnorm(150))
  
  result <- select_lag_maic(y)
  
  expect_s3_class(result, "lag_selection")
  expect_true(result$selected_lag >= 0)
  expect_true(result$selected_lag <= result$maxlag)
  expect_equal(length(result$all_maic), result$maxlag + 1)
})

test_that("simulate_bounded_bm respects bounds", {
  set.seed(101)
  
  # Two-sided
  bm1 <- simulate_bounded_bm(1000, c_lower = -2, c_upper = 2)
  expect_true(all(bm1 >= -2))
  expect_true(all(bm1 <= 2))
  expect_equal(length(bm1), 1001)
  
  # One-sided (lower only)
  bm2 <- simulate_bounded_bm(1000, c_lower = -1, c_upper = Inf)
  expect_true(all(bm2 >= -1))
})

test_that("specific test selection works", {
  set.seed(222)
  y <- cumsum(rnorm(80))
  y <- pmax(0, pmin(10, y + 5))
  
  # ADF only
  result_adf <- boundedur(y, lbound = 0, ubound = 10, test = "adf", nsim = 49)
  expect_true(all(c("adf_alpha", "adf_t") %in% names(result_adf$statistics)))
  expect_false("mz_alpha" %in% names(result_adf$statistics))
  
  # Single test
  result_msb <- boundedur(y, lbound = 0, ubound = 10, test = "msb", nsim = 49)
  expect_equal(names(result_msb$statistics), "msb")
})

test_that("detrending options work", {
  set.seed(333)
  y <- cumsum(rnorm(100)) + 50
  y <- pmax(0, pmin(100, y))
  
  result_const <- boundedur(y, lbound = 0, ubound = 100,
                            detrend = "constant", nsim = 49)
  result_none <- boundedur(y, lbound = 0, ubound = 100,
                           detrend = "none", nsim = 49)
  
  expect_equal(result_const$detrend, "constant")
  expect_equal(result_none$detrend, "none")
})

test_that("print and summary methods work", {
  set.seed(444)
  y <- cumsum(rnorm(60))
  y <- pmax(0, pmin(10, y + 5))
  
  result <- boundedur(y, lbound = 0, ubound = 10, nsim = 49)
  
  # Print should work without error
  expect_output(print(result), "Bounded Unit Root Tests")
  expect_output(summary(result), "Summary:")
})
