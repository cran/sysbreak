sim_spill <- function() {
  set.seed(1)
  bigt <- 120
  y <- matrix(0, bigt, 2)
  colnames(y) <- c("y1", "y2")
  for (tt in 2:bigt) {
    a <- if (tt > 60) matrix(c(0.2, 0.6, 0.6, 0.2), 2) else diag(0.2, 2)
    y[tt, ] <- a %*% y[tt - 1, ] + rnorm(2, sd = 0.3)
  }
  y
}

test_that("spillover measures have their defining properties", {
  fit <- varbreak(sim_spill(), lags = 1, m = 1, trim = 0.2)
  sp <- dyspill(fit, h = 10)
  for (r in sp$regimes) {
    expect_equal(unname(rowSums(r$fevd)), rep(1, 2))
    expect_equal(sum(r$net), 0)
    expect_true(r$total >= 0 && r$total <= 100)
    expect_equal(unname(r$net), unname(r$to - r$from))
  }
  # the spillover rises after the break at 60
  expect_gt(sp$regimes[[2]]$total, sp$regimes[[1]]$total)
})

test_that("the bootstrap is reproducible and restores the seed", {
  fit <- varbreak(sim_spill(), lags = 1, m = 1, trim = 0.2)
  set.seed(99)
  before <- .Random.seed
  a <- dyspill(fit, h = 10, boot = 19, seed = 1)
  expect_identical(.Random.seed, before)
  b <- dyspill(fit, h = 10, boot = 19, seed = 1)
  expect_identical(a$boot$total, b$boot$total)
})

test_that("dynobreak has the same properties", {
  set.seed(1)
  bigt <- 120
  a <- matrix(c(0.2, 0.4, 0.4, 0.2), 2)
  y <- matrix(0, bigt, 2)
  for (tt in 2:bigt) y[tt, ] <- a %*% y[tt - 1, ] + rnorm(2, sd = 0.3)
  sp <- dynobreak(y, lags = 1, h = 10)
  expect_equal(unname(rowSums(sp$fevd)), rep(1, 2))
  expect_equal(sum(sp$net), 0)
  expect_true(sp$total >= 0 && sp$total <= 100)
})

test_that("dyspill rejects fits it does not support", {
  set.seed(11)
  d <- data.frame(y1 = rnorm(60), y2 = rnorm(60), x = rnorm(60))
  d$y1[31:60] <- d$y1[31:60] + 2
  fit <- sysbreak(cbind(y1, y2) ~ x, data = d, m = 1, trim = 0.2)
  expect_error(dyspill(fit, h = 10))
})
