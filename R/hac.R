#' Long-Run Covariance Matrix of the Scores
#'
#' Estimates the long-run covariance matrix of the rows of \code{res}, a matrix
#' of score contributions with one row per observation and \eqn{d} columns,
#' using the quadratic spectral kernel with the automatic bandwidth of Andrews
#' (1991); see \code{jhatpr()}. The estimate is divided by the number of
#' observations less \eqn{d}. When \code{prewhit = TRUE}, the scores are first
#' prewhitened by a VAR(1) fitted by least squares, the long-run covariance
#' \eqn{\hat J} of the VAR residuals is estimated, and the result is recolored
#' as \eqn{(I - \hat A)^{-1} \hat J (I - \hat A)^{-1 \prime}}, where
#' \eqn{\hat A} is the estimated VAR(1) coefficient matrix (Andrews and
#' Monahan, 1992).
#' When \code{robust = "hc"}, the estimate is \eqn{\sum_t v_t v_t' / (T - d)},
#' the bandwidth-zero case without autocovariance terms, appropriate when the
#' scores are serially uncorrelated; \code{prewhit} is then ignored.
#'
#' @references
#' Andrews, D. W. K. (1991). Heteroskedasticity and autocorrelation consistent
#' covariance matrix estimation. \emph{Econometrica}, 59(3), 817--858.
#'
#' Andrews, D. W. K., and Monahan, J. C. (1992). An improved
#' heteroskedasticity and autocorrelation consistent covariance matrix
#' estimator. \emph{Econometrica}, 60(4), 953--966.
#'
#' @noRd
correct <- function(res, prewhit, robust) {
  d <- ncol(res)
  nt <- nrow(res)
  if (robust == "hc") {
    return(crossprod(res) / (nt - d))
  }

  vmat <- res
  if (prewhit) {
    bmat <- matrix(0, nrow = d, ncol = d)
    vstar <- matrix(0, nrow = nt - 1, ncol = d)

    for (i in seq_len(d)) {
      x_lag <- vmat[1:(nt - 1), , drop = FALSE]
      y_now <- vmat[2:nt, i, drop = FALSE]
      b <- pinv(x_lag) %*% y_now
      bmat[, i] <- as.vector(b)
      vstar[, i] <- vmat[2:nt, i] - as.vector(vmat[1:(nt - 1), , drop = FALSE] %*% b)
    }

    jh <- jhatpr(vstar)
    
    # bmat[, i] holds equation i's coefficients, so the fitted VAR is
    # v_t = t(bmat) %*% v_{t-1}, i.e. A = t(bmat). Andrews-Monahan
    # recoloring is (I - A)^-1 J (I - A)^-1'. The GAUSS original
    # (mbreak.src:1892) omits this transpose; corrected in 1.0.8.
    recolor <- pinv(diag(d) - t(bmat))
    
    #recolor <- pinv(diag(d) - bmat)
    return(recolor %*% jh %*% t(recolor))
  }

  jhatpr(vmat)
}

#' Andrews HAC Long-Run Covariance
#'
#' Computes the HAC covariance estimate using the automatic bandwidth and
#' quadratic spectral kernel used by the original mbreak routines.
#'
#' @noRd
jhatpr <- function(vmat) {
  nt <- nrow(vmat)
  d <- ncol(vmat)

  st <- bandw(vmat)
  jhat <- t(vmat) %*% vmat

  for (j in seq_len(nt - 1)) {
    weight <- kern(j / st)
    jhat <- jhat + weight * t(vmat[(j + 1):nt, , drop = FALSE]) %*%
      vmat[1:(nt - j), , drop = FALSE]
  }

  for (j in seq_len(nt - 1)) {
    weight <- kern(j / st)
    jhat <- jhat + weight * t(vmat[1:(nt - j), , drop = FALSE]) %*%
      vmat[(j + 1):nt, , drop = FALSE]
  }

  jhat / (nt - d)
}

#' Automatic HAC Bandwidth
#'
#' Computes the plug-in bandwidth used by the HAC covariance estimator.
#'
#' @noRd
bandw <- function(vhat) {
  nt <- nrow(vhat)
  d <- ncol(vhat)

  a2n <- 0
  a2d <- 0
  for (i in seq_len(d)) {
    x_lag <- vhat[1:(nt - 1), i, drop = FALSE]
    y_now <- vhat[2:nt, i, drop = FALSE]
    b <- as.numeric(pinv(x_lag) %*% y_now)
    sig <- sum((vhat[2:nt, i] - b * vhat[1:(nt - 1), i])^2) / (nt - 1)
    a2n <- a2n + 4 * b * b * sig * sig / (1 - b)^8
    a2d <- a2d + sig * sig / (1 - b)^4
  }

  a2 <- a2n / a2d
  # 1.3221 * (a2 * nt)^0.2
  # Cap added 2026-10-10: the effective degrees of freedom of the HAC estimate,
  # about nt / st for the quadratic spectral kernel, are kept at least 2d; a
  # larger bandwidth makes the estimate near singular when nt is small relative
  # to d (e.g. short segments in seqf()), and the Wald statistic then diverges.
  st <- 1.3221 * (a2 * nt)^0.2
  min(st, nt / (2 * d))
}

#' Quadratic Spectral Kernel
#'
#' Evaluates the quadratic spectral kernel at a scalar argument.
#'
#' @noRd
kern <- function(x) {
  if (abs(x) < .Machine$double.eps) {
    return(1)
  }

  del <- 6 * pi * x / 5
  3 * (sin(del) / del - cos(del)) / (del * del)
}
