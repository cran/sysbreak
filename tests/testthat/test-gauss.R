# Reference: the original GAUSS program mbreak2.prg (Qu and Perron 2007) run on
# data1, output in mbr-out.out. GAUSS prints three decimals.

d1 <- data.frame(y1 = data1[, 1], y2 = data1[, 2],
                 x1 = data1[, 4], x2 = data1[, 5])

test_that("restricted model reproduces GAUSS", {
  # only the intercepts change, as in mbreak2.prg
  fit <- sysbreak(list(y1 ~ x1, y2 ~ x2), data = d1, m = 2, trim = 0.2,
                  coef_breaks = list(y1 ~ 1, y2 ~ 1))
  expect_equal(fit$testing$statistic[2], 48.382, tolerance = 1e-3,
               ignore_attr = TRUE)
  expect_equal(fit$estimation$breaks, c(50, 99))
  expect_equal(unname(confint(fit)[, c(1, 3)]),
               rbind(c(43, 57), c(92, 106)))
  expect_equal(unname(confint(fit, level = 0.90)[, c(1, 3)]),
               rbind(c(45, 55), c(94, 104)))
  expect_equal(unname(fit$estimation$covariance[1:2, ]),
               rbind(c(0.709, 0.377), c(0.377, 0.803)), tolerance = 1e-3)
})

test_that("unrestricted model: one-break test reproduces GAUSS", {
  fit <- sysbreak(list(y1 ~ x1, y2 ~ x2), data = d1, m = 2, trim = 0.2)
  expect_equal(fit$testing$statistic[1], 17.939, tolerance = 1e-3,
               ignore_attr = TRUE)
  # Two breaks: GAUSS reports 56, 100 and SupLR 55.718, a local optimum; the
  # partition below has the higher likelihood (det Sigma 0.4007 vs 0.4070).
  expect_equal(fit$estimation$breaks, c(50, 100))
  expect_equal(fit$testing$statistic[2], 58.03423595, tolerance = 1e-7,
               ignore_attr = TRUE)
})
