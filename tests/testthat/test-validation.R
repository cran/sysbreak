set.seed(11)
d <- data.frame(y1 = rnorm(40), y2 = rnorm(40), x = rnorm(40))

test_that("invalid levels are rejected", {
  expect_error(sysbreak(cbind(y1, y2) ~ x, data = d, m = 1, ci_level = 95),
               "ci_level")
  expect_error(sysbreak(cbind(y1, y2) ~ x, data = d, m = 1, test_level = 0.07),
               "test_level")
})

test_that("infeasible trimming and break numbers are rejected", {
  expect_error(sysbreak(cbind(y1, y2) ~ x, data = d, m = 1, trim = 0.6))
  expect_error(sysbreak(cbind(y1, y2) ~ x, data = d, m = 10, trim = 0.2))
})

test_that("formula input requires a data frame", {
  expect_error(sysbreak(cbind(y1, y2) ~ x, m = 1), "data")
})

test_that("unknown coefficient selectors are rejected", {
  expect_error(sysbreak(cbind(y1, y2) ~ x, data = d, m = 1,
                        coef_breaks = list(y1 ~ z)))
})

test_that("brtype needs a fit with breaks", {
  set.seed(3)
  d0 <- data.frame(y = rnorm(60), x = rnorm(60))
  fit <- sysbreak(y ~ x, data = d0, m = 1, trim = 0.2)
  expect_error(brtype(fit))
})
