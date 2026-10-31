#' Compute Coefficient Standard Errors Conditional on Break Dates
#'
#' Computes variance-covariance matrices and standard errors for the estimated
#' regime coefficient vector, treating estimated break dates as fixed.
#'
#' @noRd
coefse <- function(maty, matz, n, m, bigt, br, beta, vv,
                                        R = NULL, robust,
                                        brv = FALSE,
                                        prewhit) {
  q <- ncol(matz)
  xbar <- pzbar(matz, m, br, bigt)
  beta <- matrix(beta, ncol = 1)

  R <- valR(R, m = m, q = q)
  xtheta <- xbar %*% R

  residual <- maty - xbar %*% beta
  weighted <- wsys(
    x = xtheta,
    residual = residual,
    vv = vv,
    n = n,
    m = m,
    bigt = bigt,
    br = br,
    brv = brv
  )

  vcov_theta <- if (robust != "none") {
    vcovhac(
      xstar = weighted$x,
      ustar = weighted$residuals,
      n = n,
      m = m,
      bigt = bigt,
      br = br,
      brv = brv,
      prewhit = prewhit,
      robust = robust
    )
  } else {
    pinv(crossprod(weighted$x))
  }

  vcov_beta <- R %*% vcov_theta %*% t(R)
  vcov_beta <- (vcov_beta + t(vcov_beta)) / 2
  vcov_theta <- (vcov_theta + t(vcov_theta)) / 2

  beta_names <- bnames(m = m, q = q)
  dimnames(vcov_beta) <- list(beta_names, beta_names)
  se_beta <- sqdiag(vcov_beta)
  names(se_beta) <- beta_names

  list(
    vcov = vcov_beta,
    se = se_beta,
    se_by_regime = byregime(se_beta, m, q),
    vcov_theta = vcov_theta,
    se_theta = sqdiag(vcov_theta)
  )
}

#' Multiply Regressors and Residuals by Sigma^(-1/2)
#'
#' Whitens the system: it multiplies each period's regressors and residuals by
#' Sigma^(-1/2), so that the transformed errors have identity covariance.
#' `coefse()` uses the result to compute the coefficient standard errors.
#'
#' @noRd
wsys <- function(x, residual, vv, n, m, bigt, br,
                                        brv) {
  xstar <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  ustar <- matrix(0, nrow = bigt, ncol = n)
  endpoints <- c(0, br, bigt)

  for (seg in seq_len(m + 1)) {
    sigma_idx <- if (brv) seg else 1
    sigma <- vv[((sigma_idx - 1) * n + 1):(sigma_idx * n), , drop = FALSE]
    sigma_inv_sqrt <- pinv(sqrm(sigma))

    t_start <- endpoints[[seg]] + 1
    t_end <- endpoints[[seg + 1]]
    for (tt in t_start:t_end) {
      rows <- obs_rows(tt, tt, n)
      xstar[rows, ] <- sigma_inv_sqrt %*% x[rows, , drop = FALSE]
      ustar[tt, ] <- as.vector(sigma_inv_sqrt %*% residual[rows, , drop = FALSE])
    }
  }

  list(x = xstar, residuals = ustar)
}

#' Robust Covariance Matrix of the Coefficient Estimates
#'
#' Computes the sandwich covariance matrix of the coefficient estimates,
#' robust to heteroskedasticity and, when \code{robust = "hac"}, serial correlation in the errors. The
#' inputs are the regressors and residuals multiplied by Sigma^(-1/2)
#' (see `wsys()`). The long-run variance of the scores is estimated by
#' `correct()`, regime by regime when `brv = TRUE`.
#'
#' @noRd
vcovhac <- function(xstar, ustar, n, m, bigt, br,
                                 brv, prewhit, robust) {
  score_mat <- score_contrib(xstar, ustar, n = n, bigt = bigt)
  xx <- crossprod(xstar) / bigt
  bread <- pinv(xx)

  if (brv) {
    endpoints <- c(0, br, bigt)
    meat <- matrix(0, nrow = ncol(score_mat), ncol = ncol(score_mat))
    for (seg in seq_len(m + 1)) {
      t_start <- endpoints[[seg]] + 1
      t_end <- endpoints[[seg + 1]]
      seg_sc <- score_mat[t_start:t_end, , drop = FALSE]
      # In regime `seg`, the scores are nonzero only for the coefficients that
      # appear in that regime; the columns for the other regimes' coefficients
      # are all equal to zero. The code keeps only the nonzero columns.
      # The reason is that the degrees-of-freedom adjustment in correct()
      # divides by the number of observations minus the number of columns it
      # receives. That number should count only the coefficients estimated in
      # this regime. If the columns of zeros were also included, too much would
      # be subtracted and the variance would be overstated.
      active <- which(colSums(abs(seg_sc)) > 0)
      weight <- (t_end - t_start + 1) / bigt
      meat[active, active] <- meat[active, active] +
        weight * correct(seg_sc[, active, drop = FALSE], prewhit = prewhit, robust = robust)
    }
  } else {
    meat <- correct(score_mat, prewhit = prewhit, robust = robust)
  }

  bread %*% (meat / bigt) %*% bread
}

#' Scores of the Coefficients, One Row per Period
#'
#' Returns a matrix whose row `t` is X_t' Sigma^(-1) u_t, the contribution of
#' period `t` to the first-order condition for the coefficients. It is computed
#' from the regressors and residuals multiplied by Sigma^(-1/2) (see `wsys()`).
#' `vcovhac()` estimates the variance of these rows, or their long-run variance
#' when `robust = "hac"`.
#'
#' @noRd
score_contrib <- function(xstar, ustar, n, bigt) {
  out <- matrix(0, nrow = bigt, ncol = ncol(xstar))
  for (tt in seq_len(bigt)) {
    rows <- obs_rows(tt, tt, n)
    out[tt, ] <- as.vector(
      t(xstar[rows, , drop = FALSE]) %*% matrix(ustar[tt, ], ncol = 1)
    )
  }
  out
}

#' Standard Errors from a Covariance Matrix
#'
#' Returns the square roots of the diagonal of a covariance matrix. Diagonal
#' entries that are negative only because of rounding error (between
#' -sqrt(machine epsilon) and 0) are set to zero, so that they give 0
#' and not NaN.
#'
#' @noRd
sqdiag <- function(x) {
  variances <- diag(x)
  variances[variances < 0 & variances > -sqrt(.Machine$double.eps)] <- 0
  sqrt(variances)
}

#' Labels for the Coefficients of All Regimes
#'
#' Returns the labels "regime1.coef1", "regime1.coef2", and so on: one for each
#' of the `q` coefficients in each of the `m + 1` regimes, in the order that
#' the coefficients are stacked. `coefse()` uses them to label the covariance
#' matrix and the standard errors.
#'
#' @noRd
bnames <- function(m, q) {
  paste0(
    "regime", rep(seq_len(m + 1), each = q),
    ".coef", rep(seq_len(q), times = m + 1)
  )
}
