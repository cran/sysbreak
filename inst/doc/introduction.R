## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)

## -----------------------------------------------------------------------------
library(sysbreak)

set.seed(3)
bigt <- 60
x <- rnorm(bigt)

beta <- rbind(
  c(1.0,  0.5),
  c(2.0, -0.2),
  c(0.5,  0.8)
)

y <- numeric(bigt)
y[1:20]  <- beta[1, 1] + beta[1, 2] * x[1:20]  + rnorm(20, sd = 0.4)
y[21:40] <- beta[2, 1] + beta[2, 2] * x[21:40] + rnorm(20, sd = 0.4)
y[41:60] <- beta[3, 1] + beta[3, 2] * x[41:60] + rnorm(20, sd = 0.4)

df <- data.frame(y = y, x = x)

## -----------------------------------------------------------------------------
fit <- sysbreak(y ~ x, data = df, m = 2)

## -----------------------------------------------------------------------------
summary(fit)

## -----------------------------------------------------------------------------
confint(fit)
coef(fit)

## -----------------------------------------------------------------------------
plot(fit)

## -----------------------------------------------------------------------------
set.seed(13)
y_r <- 1 + rep(c(0.5, -0.2, 0.8), each = 20) * x + rnorm(bigt, sd = 0.4)
df_r <- data.frame(y = y_r, x = x)

## -----------------------------------------------------------------------------
restricted_fit <- sysbreak(y ~ x, data = df_r, m = 2, coef_breaks = list(y ~ x))
summary(restricted_fit)

## -----------------------------------------------------------------------------
fit_hac <- sysbreak(y ~ x, data = df, m = 2, robust = "hac", ci_level = 0.90)

## -----------------------------------------------------------------------------
set.seed(14)
bigt_var <- 200
A <- matrix(c(0.5, 0.1, 0.2, 0.4), 2, 2)
chol1 <- t(chol(matrix(c(1, 0.2, 0.2, 1), 2, 2)))
chol2 <- t(chol(matrix(c(4, 1.4, 1.4, 1), 2, 2)))
y_var <- matrix(0, bigt_var, 2, dimnames = list(NULL, c("y1", "y2")))
for (tt in 2:bigt_var) {
  chol_t <- if (tt <= 100) chol1 else chol2
  y_var[tt, ] <- A %*% y_var[tt - 1, ] + chol_t %*% rnorm(2)
}

## -----------------------------------------------------------------------------
var_fit <- varbreak(y_var, lags = 1, m = 2, brv = TRUE, hetq = TRUE)
confint(var_fit)

## -----------------------------------------------------------------------------
round(brtype(var_fit), 4)

## -----------------------------------------------------------------------------
set.seed(5)
bigt_nvar <- 80
W <- matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE)
y_nvar <- matrix(0, bigt_nvar, 2)
colnames(y_nvar) <- c("node1", "node2")

for (tt in 2:bigt_nvar) {
  alpha <- if (tt <= 40) 0.05 else 0.75
  beta <- 0.1
  y_nvar[tt, ] <- alpha * as.numeric(W %*% y_nvar[tt - 1, ]) +
    beta * y_nvar[tt - 1, ] + rnorm(2, sd = 0.2)
}

nvar_fit <- netbreak(
  y = y_nvar,
  W = W,
  lags = 1,
  m = 1,
  coef_breaks = c("network", "own"),
  trim = 0.2,
  h = 10
)

nvar_fit$estimation$breaks
nvar_fit$spillover$regimes[[1]]$total

## ----eval = FALSE-------------------------------------------------------------
# spill <- dyspill(nvar_fit, h = 10)
# spill$regimes[[1]]$fevd
# spill$regimes[[1]]$from
# spill$regimes[[1]]$to
# spill$regimes[[1]]$net

## -----------------------------------------------------------------------------
set.seed(13)
bigt2 <- 60
x1 <- rnorm(bigt2)
x2 <- rnorm(bigt2)
y1 <- 1 + 0.5 * x1 + rnorm(bigt2, sd = 0.3)
y2 <- 1 + 0.5 * x2 + rnorm(bigt2, sd = 0.3)
y1[21:bigt2] <- y1[21:bigt2] + 1.5 * x1[21:bigt2]
y2[36:bigt2] <- y2[36:bigt2] - 1.5 * x2[36:bigt2]
df2 <- data.frame(y1 = y1, y2 = y2, x1 = x1, x2 = x2)

ordered_fit <- ordbreak(y1 ~ x1, y2 ~ x2, data = df2, trim = 0.2)

ordered_fit$estimation$breaks

## -----------------------------------------------------------------------------
sessionInfo()

