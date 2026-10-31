# sysbreak
## Structural Breaks in a System of Equations

## Overview

Methods for estimating and testing multivariate linear structural break models with coefficient breaks, covariance breaks, or both. The package also includes VAR, Network VAR and ordered-break routines, tests of what changes at each break, and generalized Diebold-Yilmaz spillover measures.

## Example

```r
library(sysbreak)

set.seed(1)
bigt <- 60
x <- rnorm(bigt)
beta <- rbind(c(1, 0.5), c(2, -0.2), c(0.5, 0.8))
y <- numeric(bigt)
y[1:20] <- beta[1, 1] + beta[1, 2] * x[1:20] + rnorm(20, sd = 0.4)
y[21:40] <- beta[2, 1] + beta[2, 2] * x[21:40] + rnorm(20, sd = 0.4)
y[41:60] <- beta[3, 1] + beta[3, 2] * x[41:60] + rnorm(20, sd = 0.4)
dat <- data.frame(y = y, x = x)

fit <- sysbreak(y ~ x, data = dat, m = 2, trim = 0.15)
fit$estimation$breaks
fit$testing$statistic
```

## References

Qu, Z., and Perron, P. (2007). "Estimating and Testing Structural Changes in Multivariate Regressions".
*Econometrica*, 75(2), 459-502.

Bai, J., and Perron, P. (2003). "Critical Values for Multiple Structural Change Tests".
*The Econometrics Journal*, 6(1), 72-78.

Diebold, F. X., and Yilmaz, K. (2012).
"Better to Give than to Receive: Predictive Directional Measurement of Volatility Spillovers".
*International Journal of Forecasting*, 28(1), 57-66.
