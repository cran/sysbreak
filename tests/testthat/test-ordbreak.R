# Reference: data2 from the original GAUSS ordered-break program.

d2 <- data.frame(y1 = data2[, 1], y2 = data2[, 2], x1 = data2[, 3],
                 x2 = data2[, 4], x3 = data2[, 5])

test_that("ordered two-break model on data2", {
  fit <- ordbreak(y1 ~ x2, y2 ~ x3, data = d2, trim = 0.2)
  expect_equal(fit$estimation$breaks, c(43, 53), ignore_attr = TRUE)
  expect_output(print(fit))
})

test_that("ordcheck describes the coefficient layout", {
  expect_output(ordcheck(y1 ~ x2, y2 ~ x3, data = d2), "set")
})
