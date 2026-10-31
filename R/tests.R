# Likelihood-ratio tests and asymptotic critical-values.

#' SupLR Statistic for Breaks at Fixed Break Dates
#'
#' Computes the likelihood-ratio statistic comparing a fitted break model 
#' with the no-break null model.
#'
#' @noRd
suplr <- function(maty, matx, m, n, br, rbeta, bigt,
                  brbeta = TRUE, brv = FALSE,
                  tol = 1e-6, max_iter = 1000) {
  br <- if (m == 0) integer() else br
  pmatx <- pzbar(matx, m, br, bigt)

  res <- maty - pmatx %*% rbeta
  res <- matrix(as.vector(res), ncol = n, byrow = TRUE)
  endpoints <- c(0, br, bigt)

  if (brv) {
    ulr <- 0
    for (k in seq_len(m + 1)) {
      i <- endpoints[[k]] + 1
      j <- endpoints[[k + 1]]
      res_k <- res[i:j, , drop = FALSE]
      vvar <- crossprod(res_k) / (j - i + 1)
      ulr <- ulr + nllk(j - i + 1, n, vvar)
    }
  } else if (brbeta) {
    vvar <- crossprod(res[seq_len(bigt), , drop = FALSE]) / bigt
    ulr <- nllk(bigt, n, vvar)
  } else {
    stop("At least one of `brbeta` or `brv` must be `TRUE` for `suplr()`.")
  }

  null_fit <- fgls(matx, maty, n = n, tol = tol, max_iter = max_iter)
  rlr <- nllk(bigt, n, null_fit$sigma)

  as.numeric(2 * (ulr - rlr))
}

#' Sequential LR Test Statistic
#'
#' Computes the sequential `l + 1 | l` likelihood-ratio statistic by inserting
#' an additional break into each admissible segment of the current partition.
#'
#' @noRd
seqlr <- function(maty, matz, n, br, bigt, h,
                  brbeta = TRUE, brv = FALSE,
                  trm = NULL, R = NULL,
                  tol = 1e-6, max_iter = 1000) {
  if (is.null(trm)) {
    trm <- h / bigt
  }

  nseg <- length(br) + 1
  dv <- c(0, br, bigt)
  lrtest <- numeric(nseg)
  no_room <- 0
  q <- ncol(matz)
  if (is.null(R)) {
    R <- defR(1L, q, brbeta = brbeta, brv = brv)
  } else {
    R <- valR(R, m = 1L, q = q)
  }
  # Bartlett correction, added 2026-09-29: free parameters per equation in the
  # one-break model, and changing parameters per equation.
  k <- ncol(R) / n
  p <- (ncol(R) - q) / n

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
      # The bound is `h` = round(trm * bigt), the same minimum segment length
      # sup(l|0) and WDmax use (mainp.R:198), not round(trm * seg_bigt) as
      # est() would apply on a slice.  Otherwise the sequential test searches
      # sub-segments shorter than any regime the estimation path can return,
      # and under brv those are too short to estimate a regime covariance:
      # logdet() then rewards the near-singular fit. The cost
      # is that supseq() is evaluated at trm while the trimming actually
      # applied is h / seg_bigt, which is larger -- the test is conservative.
      best <- -Inf
      for (tau in h:(seg_bigt - h)) {
        rfit <- r_estim(segy, segz, seg_bigt, n, 1, tau, R,
                        brbeta = brbeta, brv = brv,
                        tol = tol, max_iter = max_iter)
        lr_tau <- suplr(segy, segz, 1, n, tau, rfit$nbeta, seg_bigt,
                        brbeta = brbeta, brv = brv,
                        tol = tol, max_iter = max_iter)
        if (lr_tau > best) best <- lr_tau
      }
      # lrtest[is] <- best
      #
      # Bartlett correction, added 2026-09-29.  With common regressors, p of
      # them changing in every equation, and a fixed date, LR = -T log(Lambda)
      # with Lambda ~ Wilks(n, p, T - k), and -(T - b) log(Lambda) is chi^2(np)
      # up to O(T^-2), b = k + (n + 1 - p) / 2 (Bartlett 1938; Anderson 2003,
      # Theorem 8.5.2).  So LR is scaled by (T_s - b) / T_s on each segment.
      # The (n + 1 - p) / 2 term comes from estimating Sigma and dominates
      # when n is large.  With shared coefficients or equation-specific
      # regressors (netbreak, restricted R), k and p are effective counts and
      # the factor is an approximation.  Coefficient breaks only: with brv =
      # TRUE the covariance part is left uncorrected (its factor depends on
      # tau); covariance-only tests (brbeta = FALSE) are not scaled.
      if (brbeta) {
        bart <- (seg_bigt - k - (n - p + 1) / 2) / seg_bigt
        if (bart <= 0) {
          warning("Segment too short for the Bartlett correction; treated as having no room for a break.")
          bart <- 0
        }
        lrtest[is] <- best * bart
      } else {
        lrtest[is] <- best
      }

      # Replaced 2026-08-11.  Previous iterative search:
      # fit <- est(
      #   segy, segz, n, 1, seg_bigt, trm, R,
      #   brbeta = brbeta,
      #   brv = brv,
      #   tol = tol,
      #   max_iter = max_iter
      # )
      # lrtest[is] <- suplr(
      #   segy, segz, 1, n, fit$dx, fit$rbeta, seg_bigt,
      #   brbeta = brbeta,
      #   brv = brv,
      #   tol = tol,
      #   max_iter = max_iter
      # )
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

#' WDmax Test Statistic
#'
#' Computes the weighted double-maximum statistic using
#' response-surface critical values.
#'
#' @param lrtest Numeric vector of test statistics for \eqn{m=1,\dots,M}.
#' @param q Integer number of changing parameters used for the critical values.
#' @param trm Numeric trimming parameter.
#'
#' @return Numeric vector of weighted double-maximum statistics at the 90%,
#'   95%, 97.5%, and 99% levels.
#'
#' @noRd
wdmax <- function(lrtest, q, trm) {
  cvs <- t(vapply(seq_along(lrtest), function(k) cvlr(k, q, trm), numeric(4)))
  out <- numeric(4)
  for (i in seq_along(lrtest)) {
    out <- pmax(out, lrtest[[i]] * (cvs[1, ] / cvs[i, ]))
  }
  out
}

#' Critical Values for SupLR/SupF Tests
#'
#' Computes response-surface critical values for the sup test with a
#' fixed number of breaks.
#'
#' @param k Integer number of breaks under the alternative.
#' @param q Integer number of changing parameters.
#' @param trm Numeric trimming parameter.
#'
#' @return Numeric vector of critical values at the 90%, 95%, 97.5%, and 99%
#'   quantiles.
#'
#' @references
#' Bai, J., and Perron, P. (2003). Critical values for multiple structural
#' change tests. \emph{The Econometrics Journal}, 6(1), 72--78.
#'
#' @examples
#' cvlr(k = 2, q = 4, trm = 0.15)
#'
#' @noRd
cvlr <- function(k, q, trm) {
  s <- q / 100
  x1 <- c(
    1, s, s^2, k, trm, s / trm, sqrt(s), sqrt(s) * k,
    k^2, trm * k
  )
  x2 <- c(
    1 / k, 1 / (trm * k), 1 / sqrt(q),
    trm^2 / (k * sqrt(q)), 1 / (k * sqrt(q)),
    1 / (k^2 * sqrt(q)), trm / k, trm^2, trm / k^2
  )
  gamma <- rbind(
    c(11.196782, 104.88712, -4.9281477, -0.28498119, -8.5716909,
      0.40158261, 27.390867, -0.87199799, 0.020878203, -8.2428143),
    c(12.071082, 105.10911, -5.2622719, -0.4656648, -7.1676013,
      0.3997595, 28.774711, -0.94067789, 0.031926828, -8.3158755),
    c(13.494456, 106.34499, -6.0409618, -0.74248308, -6.1608266,
      0.38309265, 29.070138, -0.94566534, 0.04745933, -8.4308959),
    c(14.502622, 107.38504, -6.6871769, -0.92695539, -3.8823578,
      0.34916238, 30.223949, -0.99855164, 0.059769359, -8.7681435)
  )
  delta <- rbind(
    c(-0.028469926, 0.0025439001, -0.64887893, 1.2362379,
      0.11933421, -0.11407794, 1.9155175, -2.5112806, -0.90555119),
    c(-0.01564771, 0.0024439047, -0.60418473, 0.88703789,
      0.099245773, -0.024750329, 1.8865373, -2.3992114, -0.8987764),
    c(-0.0063640379, 0.0026798169, -0.55157844, 0.60550801,
      -0.041545425, 0.11773852, 1.9113885, -2.3368593, -0.91767072),
    c(0.0065757163, 0.0029624892, -0.52237261, 0.28120251,
      -0.0095765891, 0.14713479, 1.7996437, -2.1991939, -0.8416053)
  )

  k * as.numeric(gamma %*% x1) * exp(as.numeric(delta %*% x2))
}

#' Critical Values for the WDmax Test
#'
#' Computes response-surface critical values for the weighted
#' double-maximum test.
#'
#' @param q Integer number of changing parameters.
#' @param trm Numeric trimming parameter.
#' @param M Integer maximum number of breaks.
#'
#' @return Numeric vector of critical values at the 90%, 95%, 97.5%, and 99%
#'   quantiles.
#'
#' @references
#' Bai, J., and Perron, P. (2003). Critical values for multiple structural
#' change tests. \emph{The Econometrics Journal}, 6(1), 72--78.
#'
#' @examples
#' cvwd(q = 4, trm = 0.15, M = 5)
#'
#' @noRd
cvwd <- function(q, trm, M) {
  s <- q / 100
  inv_m <- 1 / M
  inv_q <- 1 / sqrt(q)
  x1 <- c(
    1, sqrt(s), s, s^2, trm, inv_m, inv_m^2,
    sqrt(s) * inv_m, trm * inv_m, inv_q, inv_q * inv_m, trm^2
  )
  x2 <- c(trm, inv_m, inv_q)
  gamma <- rbind(
    c(7.8928802, 54.926508, 99.176192, 1.7556585, -10.031546,
      0.21919933, 0.51560732, -1.5198532, -8.334575, 10.903617,
      -3.608048, 0.36988449),
    c(8.8974447, 56.455373, 100.05552, 1.7260493, -8.1088949,
      -0.16749937, 0.50620779, 0.57425313, -7.4276452, 10.813087,
      -3.0963096, -4.9004828),
    c(10.334167, 56.126876, 102.61121, 1.1313853, -7.1731583,
      -0.26675787, 0.32134578, 1.7364387, -6.3947836, 9.7661516,
      -2.5539402, -7.1142979),
    c(11.683051, 55.260852, 105.82025, 0.41543143, -6.6109861,
      -0.21299647, 0.088659723, 3.9162668, -5.7368093, 7.4874611,
      -1.9188067, -6.0109173)
  )
  delta <- rbind(
    c(-0.091920203, -0.028809008, -0.97091414),
    c(-0.077161524, -0.040406271, -0.84529562),
    c(-0.069996831, -0.045689675, -0.72690982),
    c(-0.061073936, -0.059190555, -0.53098037)
  )

  as.numeric(gamma %*% x1) * exp(as.numeric(delta %*% x2))
}

#' Critical Values for Sequential Break Tests
#'
#' Computes response-surface critical values for the sequential
#' \eqn{l+1 | l} test.
#'
#' @param l Integer number of breaks under the null.
#' @param q Integer number of changing parameters.
#' @param trm Numeric trimming parameter.
#'
#' @return Numeric vector of critical values at the 90%, 95%, 97.5%, and 99%
#'   quantiles.
#'
#' @references
#' Bai, J., and Perron, P. (2003). Critical values for multiple structural
#' change tests. \emph{The Econometrics Journal}, 6(1), 72--78.
#'
#' @examples
#' supseq(l = 1, q = 4, trm = 0.15)
#'
#' @noRd
supseq <- function(l, q, trm) {
  s <- q / 100
  lp1 <- l + 1
  x1 <- c(
    1, sqrt(s), s, s^2, lp1, 1 / lp1, trm, trm^2,
    trm / sqrt(q), 1 / (sqrt(q) * lp1), 1 / (q * lp1)
  )
  gamma <- rbind(
    c(3.3365189, 63.520755, 78.566433, 8.2220678, 0.3347565,
      -10.629771, -17.125942, 3.4257456, 8.3782433, 20.470265,
      -12.417636),
    c(4.7825403, 64.960542, 80.781825, 7.3834729, 0.32532089,
      -10.018118, -16.579631, 3.2370728, 7.7487827, 19.122954,
      -11.703497),
    c(6.1728497, 66.052063, 82.915169, 6.6742274, 0.32249463,
      -9.4546005, -15.73274, 2.4382694, 6.896105, 17.94025,
      -11.035772),
    c(7.8969219, 68.058397, 84.463656, 6.3615179, 0.32062121,
      -8.9655163, -15.192873, 2.4865186, 6.1439611, 17.097296,
      -10.634641)
  )

  as.numeric(gamma %*% x1)
}

#' Gaussian Log Likelihood at the Segment MLE
#'
#' Computes the normal log likelihood contribution for a segment using the
#' segment covariance estimate.
#'
#' @noRd
nllk <- function(seg_len, n, sigma) {
  -((seg_len * n) / 2) * (log(2 * pi) + 1) -
    (seg_len / 2) * logdet(sigma)
}
