# Robust SupF tests allowing for heteroskedastic and, optionally, serially correlated errors.

#' Robust SupF Statistic
#'
#' Computes the robust SupF statistic for coefficient breaks at fixed break
#' dates, the unrestricted case.
#' The covariance of the scores is estimated by \code{correct()} according to \code{robust}.
#'
#' @noRd
supft <- function(maty, matx, m, n, br, rbeta, bigt,
                  prewhit, robust) {
  pmatx <- pzbar(matx, m, br, bigt)
  res <- maty - pmatx %*% rbeta
  umat <- matrix(as.vector(res), ncol = n, byrow = TRUE)
  sigma <- crossprod(umat) / bigt
  bigsig <- pinv(sqrm(kronecker(diag(bigt), sigma)))
  xbstar <- bigsig %*% pmatx
  ubstar <- matrix(as.vector(bigsig %*% res), ncol = n, byrow = TRUE)

  tempmat <- matrix(0, nrow = bigt, ncol = ncol(xbstar))
  for (k in seq_len(bigt)) {
    idx <- obs_rows(k, k, n)
    tempmat[k, ] <- as.vector(
      t(xbstar[idx, , drop = FALSE]) %*% matrix(ubstar[k, ], ncol = 1)
    )
  }

  vmat <- correct(tempmat, prewhit = prewhit, robust = robust)
  xx <- crossprod(xbstar) / bigt
  vbeta <- pinv(xx) %*% (vmat / bigt) %*% pinv(xx)
  rmat <- rdiff(m, ncol(matx))

  ftest <- t(rbeta) %*% t(rmat) %*%
    pinv(rmat %*% vbeta %*% t(rmat)) %*% rmat %*% rbeta
  as.numeric((bigt - (m + 1) * ncol(matx)) * ftest / bigt)
}

#' Restricted Robust SupF Statistic
#'
#' Computes the robust SupF statistic for coefficient breaks when the
#' coefficient vector is reparameterized by a restriction matrix.
#' The covariance of the scores is estimated by \code{correct()} according to \code{robust}.
#'
#' @noRd
suprft <- function(maty, matx, m, n, br, rbeta, R, bigt,
                   prewhit, robust) {
  pmatx <- pzbar(matx, m, br, bigt)
  res <- maty - pmatx %*% rbeta
  umat <- matrix(as.vector(res), ncol = n, byrow = TRUE)
  sigma <- crossprod(umat) / bigt
  bigsig <- pinv(sqrm(kronecker(diag(bigt), sigma)))
  pmatx_r <- pmatx %*% R
  xbstar <- bigsig %*% pmatx_r
  ubstar <- matrix(as.vector(bigsig %*% res), ncol = n, byrow = TRUE)

  tempmat <- matrix(0, nrow = bigt, ncol = ncol(xbstar))
  for (k in seq_len(bigt)) {
    idx <- obs_rows(k, k, n)
    tempmat[k, ] <- as.vector(
      t(xbstar[idx, , drop = FALSE]) %*% matrix(ubstar[k, ], ncol = 1)
    )
  }

  vmat <- correct(tempmat, prewhit = prewhit, robust = robust)
  xx <- crossprod(xbstar) / bigt
  vdelta <- pinv(xx) %*% (vmat / bigt) %*% pinv(xx)
  vbeta <- R %*% vdelta %*% t(R)
  rmat <- rdiff(m, ncol(matx))

  ftest <- t(rbeta) %*% t(rmat) %*%
    pinv(rmat %*% vbeta %*% t(rmat)) %*% rmat %*% rbeta
  as.numeric((bigt - (m + 1) * ncol(matx)) * ftest / bigt)
}

#' Sequential Robust F Statistic
#'
#' Computes the sequential `l + 1 | l` robust F statistic by inserting one
#' additional coefficient break into each admissible segment.
#' The robust covariance follows \code{robust}, as in \code{supft()}.
#'
#' @noRd
seqf <- function(maty, matz, n, br, bigt, h,
                 prewhit, robust, trm = NULL, R = NULL) {
  if (is.null(trm)) {
    trm <- h / bigt
  }

  nseg <- length(br) + 1
  dv <- c(0, br, bigt)
  lrtest <- numeric(nseg)
  no_room <- 0
  q <- ncol(matz)
  if (is.null(R)) {
    R <- diag(q * 2L)
  } else {
    R <- valR(R, m = 1L, q = q)
  }
  restricted <- !isTRUE(all.equal(
    unname(R), diag(q * 2L), check.attributes = FALSE
  ))

  for (is in seq_len(nseg)) {
    seg_length <- dv[[is + 1]] - dv[[is]]
    if (seg_length >= 2 * h) {
      starti <- dv[[is]] + 1
      endi <- dv[[is + 1]]
      idx <- obs_rows(starti, endi, n)
      segy <- maty[idx, , drop = FALSE]
      segz <- matz[idx, , drop = FALSE]
      seg_bigt <- endi - starti + 1
      # The statistic is a sup over the admissible dates.
      # est() alternates between coefficients and dates.
      #
      # The bound is `h` = round(trm * bigt), matching sup(l|0), WDmax and
      # seqlr(); see the fuller note in tests.R.  seqf() is brv = FALSE, so
      # the near-singular covariance that motivates it there cannot arise
      # here -- this is for consistency of the minimum segment length.
      best <- -Inf
      for (tau in h:(seg_bigt - h)) {
        rb <- r_estim(segy, segz, seg_bigt, n, 1, tau, R,
                      brbeta = TRUE, brv = FALSE)$nbeta
        f_tau <- if (restricted) {
          suprft(segy, segz, 1, n, tau, rb, R, seg_bigt, prewhit = prewhit,
                 robust = robust)
        } else {
          supft(segy, segz, 1, n, tau, rb, seg_bigt, prewhit = prewhit,
                robust = robust)
        }
        if (f_tau > best) best <- f_tau
      }
      lrtest[is] <- best

      # Replaced 2026-08-11.  Previous iterative search:
      # fit <- est(
      #   segy, segz, n, 1, seg_bigt, trm, R,
      #   brbeta = TRUE,
      #   brv = FALSE
      # )
      # lrtest[is] <- if (restricted) {
      #   suprft(
      #     segy, segz, 1, n, fit$dx, fit$rbeta, R, seg_bigt,
      #     prewhit = prewhit
      #   )
      # } else {
      #   supft(
      #     segy, segz, 1, n, fit$dx, fit$rbeta, seg_bigt,
      #     prewhit = prewhit
      #   )
      # }
    } else {
      no_room <- no_room + 1
      lrtest[is] <- 0
    }
  }

  if (no_room == nseg) {
    warning("No segment has enough observations to insert an additional break.")
  }

  max(lrtest)
}

#' Regime Difference Matrix
#'
#' Builds the matrix that takes adjacent differences of regime-specific
#' coefficient vectors.
#'
#' @noRd
rdiff <- function(m, q) {
  rsub <- matrix(0, nrow = m, ncol = m + 1)
  for (j in seq_len(m)) {
    rsub[j, j] <- -1
    rsub[j, j + 1] <- 1
  }
  kronecker(rsub, diag(q))
}
