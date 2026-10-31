check_methods <- function(fit) {
  est <- fit$estimation
  cf <- coef(fit)
  v <- vcov(fit)
  expect_equal(unname(cf), as.vector(est$coefficients))
  expect_identical(names(cf), rownames(v))
  expect_equal(dim(v), rep(length(cf), 2))
  expect_equal(rownames(v), rownames(confint(fit, type = "coef")))
  expect_equal(unname(fitted(fit) + residuals(fit)), unname(fit$y))
  expect_equal(nobs(fit), nrow(fit$y))
  ci <- confint(fit)
  expect_equal(unname(ci[, 2]), est$breaks)
  expect_equal(unname(ci[, c(1, 3), drop = FALSE]), unname(est$ci_breaks))
  z <- stats::qnorm(0.975)
  se <- sqrt(diag(v))
  expect_equal(unname(confint(fit, type = "coef")),
               unname(cbind(cf - z * se, cf + z * se)))
  zs <- unlist(lapply(summary(fit)$coefficients, function(t) t[, "z value"]))
  expect_equal(unname(cf / se), unname(zs))
}
