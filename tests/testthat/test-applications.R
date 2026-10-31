# Regression values from the applications in the paper, computed with
# version 1.1.0.

housing_data <- function() {
  n <- nrow(housing)
  data.frame(housing[-1, ],
             L1.dGDP = housing$dGDP[-n], L1.dC = housing$dC[-n],
             L1.dHP = housing$dHP[-n], L1.FFR = housing$FFR[-n])
}

test_that("Application 1: house prices, all coefficients change", {
  fit <- sysbreak(cbind(dGDP, dHP, dI) ~ L1.dGDP + L1.dHP + L1.FFR,
                  data = housing_data(), m = 3, trim = 0.20)
  expect_equal(fit$testing$statistic,
               c(56.34708277, 93.63874981, 103.89482343),
               tolerance = 1e-7, ignore_attr = TRUE)
  expect_equal(fit$br_selected$breaks, 1)
  expect_equal(fit$estimation$breaks, 92)
  expect_equal(unname(confint(fit)[1, ]), c(84, 92, 100))
})

test_that("Application 1: house prices, restricted model", {
  fit <- sysbreak(cbind(dGDP, dHP, dI) ~ L1.dGDP + L1.dHP + L1.FFR,
                  data = housing_data(), m = 3, trim = 0.20,
                  coef_breaks = list(dGDP ~ L1.dHP, dHP ~ L1.dHP,
                                     dI ~ L1.dHP))
  expect_equal(fit$testing$statistic,
               c(21.54716268, 29.56236776, 30.85634814),
               tolerance = 1e-7, ignore_attr = TRUE)
  expect_equal(fit$estimation$breaks, 101)
})

test_that("Application 2: spillover index without breaks", {
  intl <- as.matrix(inflation[, c("Canada", "UK", "EuroArea", "US")])
  expect_equal(dynobreak(intl, lags = 1, h = 12)$total, 47.38456045,
               tolerance = 1e-9)
})

test_that("Application 2: inflation VAR, coefficient and covariance breaks", {
  skip_on_cran()   # about 90 seconds
  intl <- as.matrix(inflation[, c("Canada", "UK", "EuroArea", "US")])
  fit <- suppressWarnings(varbreak(intl, lags = 1, m = 5, trim = 0.15))
  expect_equal(fit$estimation$breaks, c(125, 220))
  fit_cov <- suppressWarnings(varbreak(intl, lags = 1, m = 5, trim = 0.15,
                                       brv = TRUE))
  expect_equal(fit_cov$estimation$breaks, c(125, 220, 311))
  sp <- dyspill(fit_cov, h = 12)
  expect_equal(vapply(sp$regimes, `[[`, 0, "total"),
               c(17.82546320, 22.09298563, 19.34489331, 40.86977186),
               tolerance = 1e-7)
})
