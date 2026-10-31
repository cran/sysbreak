sim_sys <- function() {
  set.seed(11)
  bigt <- 60
  x <- rnorm(bigt)
  y1 <- 1 + 0.4 * x + rnorm(bigt, sd = 0.3)
  y2 <- -0.5 + 0.2 * x + rnorm(bigt, sd = 0.3)
  y1[31:60] <- y1[31:60] + 1.2 - 0.6 * x[31:60]
  y2[31:60] <- y2[31:60] - 0.8 + 0.4 * x[31:60]
  data.frame(y1 = y1, y2 = y2, x = x)
}

sim_var <- function() {
  set.seed(12)
  bigt <- 75
  y <- matrix(0, bigt, 2)
  for (tt in 2:bigt) {
    shift <- if (tt > 39) c(0.8, -0.6) else c(0, 0)
    y[tt, ] <- shift + c(0.3 * y[tt - 1, 1], -0.2 * y[tt - 1, 2]) +
      rnorm(2, sd = 0.3)
  }
  y
}

test_that("methods agree with the fit: sysbreak", {
  fit <- sysbreak(cbind(y1, y2) ~ x, data = sim_sys(), m = 1, trim = 0.2)
  check_methods(fit)
  expect_equal(fit$estimation$breaks, 30)
})

test_that("methods agree with the fit: restricted sysbreak", {
  fit <- sysbreak(cbind(y1, y2) ~ x, data = sim_sys(), m = 1, trim = 0.2,
                  coef_breaks = list(y1 ~ x, y2 ~ x))
  check_methods(fit)
  # the intercepts are held constant across regimes
  cf <- coef(fit)
  expect_equal(cf[["regime1:y1:(Intercept)"]], cf[["regime2:y1:(Intercept)"]])
})

test_that("methods agree with the fit: varbreak", {
  fit <- varbreak(sim_var(), lags = 1, m = 1, trim = 0.2)
  check_methods(fit)
  expect_equal(nobs(fit), 74)
})

test_that("both break-date confidence levels are available", {
  fit <- sysbreak(cbind(y1, y2) ~ x, data = sim_sys(), m = 1, trim = 0.2)
  c95 <- confint(fit)
  c90 <- confint(fit, level = 0.90)
  expect_true(all(c90[, 1] >= c95[, 1] & c90[, 3] <= c95[, 3]))
  expect_error(confint(fit, level = 0.80), "0.95 or 0.90")
})

test_that("print, summary and plot run", {
  fit <- sysbreak(cbind(y1, y2) ~ x, data = sim_sys(), m = 1, trim = 0.2)
  expect_output(print(fit))
  expect_output(print(summary(fit)), "Regime 2")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(plot(fit))
  expect_invisible(plot(fit, dates = seq_len(60) + 1990))
  expect_error(plot(fit, dates = seq_len(30) + 1990), "one entry per row")
})

test_that("summary labels the sequential tests and handles two breaks", {
  d <- sim_sys()
  d$y1[46:60] <- d$y1[46:60] - 1.5
  d$y2[46:60] <- d$y2[46:60] + 1.5
  fit <- sysbreak(cbind(y1, y2) ~ x, data = d, m = 2, trim = 0.2)
  s <- summary(fit)
  expect_identical(rownames(s$tests)[1:2], c("WDmax", "SEQ(2|1)"))
  expect_output(print(s), "SEQ(2|1)", fixed = TRUE)
  expect_output(print(s), "Break 1")
  if (identical(fit$br_selected$breaks, 2L)) {
    expect_output(print(s), "Break 2")
    expect_output(print(s), "Regime 3")
  }
})

test_that("a fit without breaks reports no estimates", {
  set.seed(3)
  d <- data.frame(y = rnorm(60), x = rnorm(60))
  fit <- sysbreak(y ~ x, data = d, m = 1, trim = 0.2)
  expect_null(fit$estimation)
  expect_error(coef(fit), "No break was selected")
  expect_equal(nobs(fit), 60)
  expect_output(print(summary(fit)), "No break was selected")
})
