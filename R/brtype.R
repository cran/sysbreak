#' Methods to Examine What Changes at Each Break Date
#'
#' For each estimated break date, tests whether the volatilities, the correlations
#' and the regression coefficients differ between the two regimes on either
#' side of it.
#'
#' @param fit A result from \code{\link{sysbreak}}, \code{\link{varbreak}} or
#'   \code{\link{netbreak}} with at least one break.
#'
#' @return A data frame with one row per break: \code{breaks}, the break date
#'   as a row of the data; the statistics \code{volatility}, \code{correlation}
#'   and \code{coefficients}; and their p-values \code{volatility_p},
#'   \code{correlation_p} and \code{coefficients_p}. The degrees of freedom are
#'   \eqn{n}, \eqn{n(n-1)/2} and the number of coefficients allowed to change.
#'
#' @details
#' The volatility statistic is a two-sample Hotelling (1931) statistic on the squared
#' residuals. The correlation statistic is the Jennrich (1970) test for equal
#' correlation matrices, which assumes normal errors. The coefficient statistic
#' is a Wald test for equal coefficients; it covers only the coefficients
#' allowed to change, and is \code{NA} when no coefficient is allowed to change. 
#' Each test compares the two regimes adjacent to a break and treats the break 
#' dates as known.
#'
#' The coefficient test holds the error covariance constant across regimes,
#' even when \code{fit} allows it to break (\code{brv = TRUE}). Estimating the
#' covariance separately in each regime would make the test oversized in samples of
#' typical length; the constant covariance tends to make it conservative when the
#' volatility does change.
#'
#' @references
#' Hotelling, H. (1931). The generalization of Student's ratio. \emph{The
#' Annals of Mathematical Statistics}, 2(3), 360-378.
#'
#' Jennrich, R. I. (1970). An asymptotic chi-square test for the equality of
#' two correlation matrices. \emph{Journal of the American Statistical
#' Association}, 65, 904-912.
#'
#'
#' @examples
#' set.seed(11)
#' bigt <- 120
#' x <- rnorm(bigt)
#' y1 <- 1 + 0.4 * x + rnorm(bigt, sd = 0.3)
#' y2 <- -0.5 + 0.2 * x + rnorm(bigt, sd = 0.3)
#' y1[61:120] <- y1[61:120] + 1.2 - 0.6 * x[61:120]
#' dat <- data.frame(y1 = y1, y2 = y2, x = x)
#'
#' fit <- sysbreak(cbind(y1, y2) ~ x, data = dat, m = 2, trim = 0.2)
#' brtype(fit)
#'
#' @export
brtype <- function(fit) {
  est <- fit$estimation
  if (is.null(est$residuals)) {
    stop("`fit` must come from sysbreak(), varbreak() or netbreak() and have at least one break.",
         call. = FALSE)
  }
  res <- est$residuals
  n <- ncol(res)
  br <- match(est$breaks, fit$spec$y_rows)   # breaks as rows of `res`
  m <- length(br)
  q <- length(est$coefficients) / (m + 1)
  chg <- which(fit$spec$coef_names %in% cbchg(fit$spec))
  # estimates with the covariance held constant
  beta <- if (is.null(est$const)) est$coefficients else est$const$coefficients
  vcov <- if (is.null(est$const)) est$vcov else est$const$vcov

  out <- data.frame(breaks = est$breaks,
                    volatility = NA_real_, volatility_p = NA_real_,
                    correlation = NA_real_, correlation_p = NA_real_,
                    coefficients = NA_real_, coefficients_p = NA_real_)
  for (j in seq_len(m)) {
    vol <- svolstat(res, br, j)
    out$volatility[j] <- vol
    out$volatility_p[j] <- stats::pchisq(vol, n, lower.tail = FALSE)

    # correlations do not depend on scale, therefore the residuals need no
    # standardization
    if (n > 1) {
      corr <- scorrstat(res, br, j)
      out$correlation[j] <- corr
      out$correlation_p[j] <- stats::pchisq(corr, n * (n - 1) / 2,
                                            lower.tail = FALSE)
    }

    if (length(chg)) {
      a <- (j - 1) * q + chg   # regime j
      b <- j * q + chg         # regime j + 1
      d <- beta[b] - beta[a]
      V <- vcov[a, a] + vcov[b, b] - vcov[a, b] - vcov[b, a]
      wald <- drop(t(d) %*% pinv(V) %*% d)
      out$coefficients[j] <- wald
      out$coefficients_p[j] <- stats::pchisq(wald, length(chg),
                                             lower.tail = FALSE)
    }
  }
  out
}

#' Hotelling-Type Volatility Break Statistic
#'
#' @noRd
svolstat <- function(residual, breaks, break_index) {
  endpoints <- c(0, breaks, nrow(residual))
  left <- (endpoints[[break_index]] + 1):endpoints[[break_index + 1]]
  right <- (endpoints[[break_index + 1]] + 1):
    endpoints[[break_index + 2]]
  left_sq <- residual[left, , drop = FALSE]^2
  right_sq <- residual[right, , drop = FALSE]^2
  mean_left <- colMeans(left_sq)
  mean_right <- colMeans(right_sq)
  centered_left <- sweep(left_sq, 2, mean_left)
  centered_right <- sweep(right_sq, 2, mean_right)
  pooled <- (
    crossprod(centered_left) + crossprod(centered_right)
  ) / (length(left) + length(right) - 2)
  difference <- mean_left - mean_right
  as.numeric(
    length(left) * length(right) / (length(left) + length(right)) *
      crossprod(difference, pinv(pooled) %*% difference)
  )
}

#' Jennrich Correlation Break Statistic
#'
#' @noRd
scorrstat <- function(standardized, breaks, break_index) {
  endpoints <- c(0, breaks, nrow(standardized))
  left <- (endpoints[[break_index]] + 1):endpoints[[break_index + 1]]
  right <- (endpoints[[break_index + 1]] + 1):
    endpoints[[break_index + 2]]
  p1 <- stats::cor(standardized[left, , drop = FALSE])
  p2 <- stats::cor(standardized[right, , drop = FALSE])
  n1 <- length(left)
  n2 <- length(right)
  pooled <- ((n1 - 1) * p1 + (n2 - 1) * p2) / (n1 + n2 - 2)
  pooled <- scovreg(pooled)
  pooled_inv <- pinv(pooled)
  z <- sqrt(n1 * n2 / (n1 + n2)) *
    pooled_inv %*% (p1 - p2)
  psi <- diag(ncol(standardized)) + pooled * pooled_inv
  diagonal <- diag(z)
  max(0, as.numeric(
    0.5 * sum(z * t(z)) -
      crossprod(diagonal, pinv(psi) %*% diagonal)
  ))
}

#' Regularize a Covariance Matrix
#'
#' @noRd
scovreg <- function(x, tol = 1e-8) {
  x <- (x + t(x)) / 2
  factor <- tryCatch(chol(x), error = function(e) NULL)
  if (!is.null(factor)) {
    diagonal <- abs(diag(factor))
    well_scaled <- length(diagonal) &&
      min(diagonal)^2 > tol * max(diagonal)^2
    if (well_scaled) return(x)
  }

  ev <- eigen(x, symmetric = TRUE)
  floor_value <- max(tol, max(abs(ev$values)) * tol)
  ev$values <- pmax(ev$values, floor_value)
  ev$vectors %*% (ev$values * t(ev$vectors))
}
