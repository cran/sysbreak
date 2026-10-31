# routines for a model with ordered breaks.

#' Select which coefficients are allowed to change
#'
#' Maps a `coef_breaks` specification to column indices of the
#' coefficient vector. `"all"` allows every coefficient to change. Otherwise
#' supply a formula per equation; `y ~ 0` means the coefficients in that equation stay fixed.
#'
#' @noRd
ordsel <- function(coef_breaks, x_list, data, offsets, which_eq) {
  if (identical(coef_breaks, "all")) {
    return(unlist(lapply(which_eq, function(i) {
      offsets[[i]] + seq_len(ncol(x_list[[i]]))
    }), use.names = FALSE))
  }

  if (!is.list(coef_breaks) || length(coef_breaks) != length(which_eq)) {
    stop("`coef_breaks` must be \"all\" or one formula per equation.",
         call. = FALSE)
  }

  out <- integer(0)
  for (k in seq_along(which_eq)) {
    i <- which_eq[[k]]
    cb <- coef_breaks[[k]]
    if (!inherits(cb, "formula")) {
      stop("Each `coef_breaks` entry must be a formula.", call. = FALSE)
    }
    tt <- stats::terms(cb, data = data)
    if (length(cb) == 3L) {
      tt <- stats::delete.response(tt)
    }
    nm <- colnames(stats::model.matrix(tt, data = data))
    if (!length(nm)) {
      next
    }
    pos <- match(nm, colnames(x_list[[i]]))
    if (anyNA(pos)) {
      stop(
        sprintf(
          "Equation %d does not contain %s, so no coefficient is allowed to change.",
          i, paste(sQuote(nm[is.na(pos)]), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    out <- c(out, offsets[[i]] + pos)
  }
  sort(out)
}


#' Ordered-Breaks Design For Two Sets of Equations
#'
#' Builds `y`, `z`, `S`, `q1` that separate the two sets of equations, and the
#' indices of the coefficients that are allowed to change.
#'
#' @noRd
ordspec <- function(set1, set2, data,
                    coef_breaks1 = "all", coef_breaks2 = "all") {
  set1 <- valfor(set1)
  set2 <- valfor(set2)
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (!length(set1) || !length(set2)) {
    stop("A set of equations must contain at least one equation.", call. = FALSE)
  }

  formulas <- c(set1, set2)
  neq <- length(formulas)
  n1 <- length(set1)

  y_list <- vector("list", neq)
  x_list <- vector("list", neq)
  y_names <- character(neq)
  for (i in seq_len(neq)) {
    frame_i <- stats::model.frame(formulas[[i]], data = data,
                                  na.action = stats::na.pass)
    terms_i <- stats::terms(formulas[[i]], data = data)
    y_i <- stats::model.response(frame_i)
    if (is.null(y_i) || !is.numeric(y_i)) {
      stop("Each equation must have a numeric response.", call. = FALSE)
    }
    y_list[[i]] <- as.numeric(y_i)
    x_list[[i]] <- stats::model.matrix(stats::delete.response(terms_i),
                                       data = frame_i)
    y_names[[i]] <- rname(formulas[[i]])
  }

  keep <- Reduce(`&`, c(lapply(y_list, stats::complete.cases),
                        lapply(x_list, stats::complete.cases)))
  if (!any(keep)) {
    stop("No complete observations remain after evaluating the formulas.",
         call. = FALSE)
  }
  if (!all(keep)) {
    warning(
      sprintf(
        "%d observations dropped for missing values. Break dates are mapped to the original data indices.",
        sum(!keep)
      ),
      call. = FALSE
    )
  }
  y <- do.call(cbind, lapply(y_list, function(v) v[keep]))
  colnames(y) <- make.unique(y_names)
  x_list <- lapply(x_list, function(x) x[keep, , drop = FALSE])

  z <- unionz(x_list)
  S <- surS(x_list, union_names = colnames(z))

  widths <- vapply(x_list, ncol, integer(1))
  offsets <- c(0L, cumsum(widths)[-neq])
  q1 <- sum(widths[seq_len(n1)])
  q2 <- sum(widths[(n1 + 1L):neq])

  sel1 <- ordsel(coef_breaks1, x_list, data, offsets, seq_len(n1))
  sel2 <- ordsel(coef_breaks2, x_list, data, offsets, (n1 + 1L):neq)
  if (!length(sel1) || !length(sel2)) {
    stop("Each set needs to allow at least one changing coefficient.", call. = FALSE)
  }

  list(
    y = y, z = z, S = S,
    q1 = q1, q2 = q2, n1 = n1,
    sel1 = sel1, sel2 = sel2,
    coef_names = colnames(S),
    eqz = stats::setNames(lapply(x_list, colnames), colnames(y)),
    formulas = vapply(formulas, flabel, character(1)),
    y_rows = which(keep)
  )
}


#' Searching for Ordered Break Dates
#'
#' Searches over ordered break candidates. Columns in `sel1` shift at the
#' first date and those in `sel2` at the second; the rest are held constant
#' across both dates.
#'
#' @noRd
brorder <- function(maty, matz, n, sel1, sel2, h, bigt, trun,
                    tol = 1e-8, max_iter = 1000) {
  q <- ncol(matz)
  k1 <- length(sel1)
  k2 <- length(sel2)
  ncoef <- q + k1 + k2
  block1 <- q + seq_len(k1)
  block2 <- q + k1 + seq_len(k2)

  vecmle <- matrix(-1e12, nrow = bigt * bigt, ncol = 3 + ncoef)
  zz <- cbind(matz, matz[, sel1, drop = FALSE], matz[, sel2, drop = FALSE])
  prev_v <- NULL

  for (ii in (h + 1):(bigt - h)) {
    jj <- 0
    while ((ii + jj) <= bigt - h) {
      z1 <- zz
      z1[1:(ii * n), block1] <- 0
      z1[1:((ii + jj) * n), block2] <- 0

      fit <- ofgls(
        segx = z1, segy = maty, n = n,
        init_s = if (jj == 0) NULL else prev_v,
        tol = tol, max_iter = max_iter
      )
      prev_v <- fit$sigma

      row_i <- as.integer((ii - 1) * bigt + jj)
      vecmle[row_i, 1] <- ii
      vecmle[row_i, 2] <- ii + jj
      vecmle[row_i, 3] <- -bigt / 2 * logdet(fit$sigma)
      vecmle[row_i, 4:(3 + ncoef)] <- as.vector(fit$beta)

      jj <- jj + 1
    }
  }

  vecmle2 <- vecmle
  valid_trunc <- (vecmle2[, 2] - vecmle2[, 1]) <= trun
  vecmle2[, 3] <- ifelse(valid_trunc, vecmle2[, 3], -1e12)

  best_unr <- vecmle[which.max(vecmle[, 3]), , drop = FALSE]
  best_tr <- vecmle2[which.max(vecmle2[, 3]), , drop = FALSE]
  H1like <- best_tr[1, 3]

  null_fit <- ofgls(matz, maty, n = n, init_s = NULL,
                    tol = tol, max_iter = max_iter)
  H0like <- -bigt / 2 * logdet(null_fit$sigma)
  lr <- 2 * (H1like - H0like)

  beta1 <- matrix(best_unr[1, 4:(3 + q)], ncol = 1)
  shift <- matrix(0, nrow = q, ncol = 1)
  shift[sel1, 1] <- best_unr[1, 3 + block1]
  shift[sel2, 1] <- best_unr[1, 3 + block2]

  delt1 <- matrix(0, nrow = q, ncol = 1)
  delt2 <- matrix(0, nrow = q, ncol = 1)
  delt1[sel1, 1] <- shift[sel1, 1]
  delt2[sel2, 1] <- shift[sel2, 1]

  # residuals of the ordered-break model at the estimated dates (brdate), from
  # the coefficients stored for that pair; interval3() estimates Sigma from
  # them. The no-break fit above is used only for H0like.
  z_hat <- zz
  z_hat[seq_len(best_unr[1, 1] * n), block1] <- 0
  z_hat[seq_len(best_unr[1, 2] * n), block2] <- 0
  res_hat <- maty - z_hat %*% best_unr[1, 4:(3 + ncoef)]
  umat <- matrix(as.vector(res_hat), ncol = n, byrow = TRUE)

  list(
    brdate = as.integer(best_unr[1, 1:2]),
    beta1 = beta1,
    beta2 = beta1 + shift,
    delt1 = delt1,
    delt2 = delt2,
    umat = umat,
    lr = as.numeric(lr),
    truncated_brdate = as.integer(best_tr[1, 1:2]),
    H1like = as.numeric(H1like),
    H0like = as.numeric(H0like),
    candidates = vecmle
  )
}

#' Confidence Intervals for Breaks
#'
#' Simulates the limiting distribution used to construct confidence intervals for
#' the two ordered break dates.
#'
#' @noRd
interval3 <- function(matz, n, start, last, dx1, dx2, delt1, delt2, res,
                      bigM, rep, hetq = TRUE,
                      step = 1000, interv = 0.02, seed = NULL,
                      verbose = FALSE) {
  if (!is.numeric(interv) || length(interv) != 1 || !is.finite(interv) || interv <= 0) {
    stop("`interv` must be a positive finite number.")
  }
  if (!is.null(seed)) {
    old_seed <- get0(".Random.seed", envir = .GlobalEnv,
                     inherits = FALSE, ifnotfound = NULL)
    on.exit({
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  sigma <- t(res) %*% res / nrow(res)
  if (hetq) {
    left_idx <- obs_rows(start, dx1, n)
    right_idx <- obs_rows(dx2 + 1, last, n)
    Qj <- t(matz[left_idx, , drop = FALSE]) %*%
      kronecker(diag(dx1 - start + 1), pinv(sigma)) %*%
      matz[left_idx, , drop = FALSE] / (dx1 - start + 1)
    Qjp1 <- t(matz[right_idx, , drop = FALSE]) %*%
      kronecker(diag(last - dx2), pinv(sigma)) %*%
      matz[right_idx, , drop = FALSE] / (last - dx2)

    P1j <- as.numeric(t(delt1) %*% Qj %*% delt1)
    P2j <- as.numeric(t(delt2) %*% Qj %*% delt2)
    P12j <- as.numeric(t(delt1) %*% Qj %*% delt2)
    P1jp1 <- as.numeric(t(delt1) %*% Qjp1 %*% delt1)
    P2jp1 <- as.numeric(t(delt2) %*% Qjp1 %*% delt2)
    P12jp1 <- as.numeric(t(delt1) %*% Qjp1 %*% delt2)
  } else {
    idx <- obs_rows(start, last, n)
    Qj <- t(matz[idx, , drop = FALSE]) %*%
      kronecker(diag(last - start + 1), pinv(sigma)) %*%
      matz[idx, , drop = FALSE] / (last - start + 1)

    P1j <- as.numeric(t(delt1) %*% Qj %*% delt1)
    P2j <- as.numeric(t(delt2) %*% Qj %*% delt2)
    P12j <- as.numeric(t(delt1) %*% Qj %*% delt2)
    P1jp1 <- P1j
    P2jp1 <- P2j
    P12jp1 <- P12j
  }

  bigN <- step * bigM
  grids <- as.integer(round(1 / interv))
  gridN <- grids * bigM
  grid_step <- as.integer(round(step * interv))
  sv <- -bigM + (1 / grids) * seq.int(0, 2 * grids * bigM)
  sv_len <- length(sv)
  V1 <- rep(sv, each = sv_len)
  V2 <- rep(sv, times = sv_len)
  obd <- matrix(0, nrow = rep, ncol = 3)

  for (simu in seq_len(rep)) {
    if (verbose && simu %% 10 == 0) {
      message("rep: ", simu)
    }

    bigW <- matrix(stats::rnorm(bigN * 4), nrow = bigN, ncol = 4)
    bigW <- apply(bigW, 2, cumsum) / sqrt(step)
    W <- bigW[seq.int(1, by = grid_step, length.out = gridN), , drop = FALSE]
    sbv1 <- c(base::rev(W[, 1]), 0, W[, 4])
    sbv2 <- c(base::rev(W[, 2]), 0, W[, 3])
    Bv1 <- rep(sbv1, each = sv_len)
    Bv2 <- rep(sbv2, times = sv_len)

    if (hetq) {
      Z1 <- (V1 < 0) * (V2 <= 0) * (V1 <= V2) *
        ((Bv1 + V1 / 2) + sqrt(P2j / P1j) * Bv2 + V2 / 2 * (P2j + 2 * P12j) / P1j)
      Z2 <- (V1 < 0) * (V2 > 0) *
        (Bv1 + V1 / 2 + sqrt(P2jp1 / P1j) * Bv2 - V2 / 2 * (P2jp1 / P1j))
      Z3 <- (V1 >= 0) * (V2 > 0) * (V1 <= V2) *
        (sqrt(P1jp1 / P1j) * Bv1 + sqrt(P2jp1 / P1j) * Bv2 -
           (P2jp1 / P1j) * V2 / 2 - ((P1jp1 + 2 * P12jp1) / P1j) * V1 / 2)
    } else {
      Z1 <- (V1 < 0) * (V2 <= 0) * (V1 <= V2) *
        ((Bv1 + V1 / 2) + sqrt(P2j / P1j) * Bv2 + V2 / 2 * (P2j + 2 * P12j) / P1j)
      Z2 <- (V1 < 0) * (V2 > 0) *
        (Bv1 + V1 / 2 + sqrt(P2j / P1j) * Bv2 - V2 / 2 * (P2j / P1j))
      Z3 <- (V1 >= 0) * (V2 > 0) * (V1 <= V2) *
        (Bv1 + sqrt(P2j / P1j) * Bv2 -
           (P2j / P1j) * V2 / 2 - ((P1j + 2 * P12j) / P1j) * V1 / 2)
    }
    Z4 <- (V1 > V2) * (-1)
    Z <- Z1 + Z2 + Z3 + Z4
    best <- which.max(Z)
    obd[simu, ] <- c(V1[[best]], V2[[best]], Z[[best]])
  }

  obd[, 1] <- sort(obd[, 1])
  obd[, 2] <- sort(obd[, 2])
  count <- sum(abs(abs(obd[, 1]) - bigM) < 1e-4) + sum(abs(abs(obd[, 2]) - bigM) < 1e-4)
  if (count >= rep * 0.05) {
    warning("Some simulated ordered-break bounds reached `bigM`.")
  }

  q025 <- max(1, floor(rep * 0.025))
  q05 <- max(1, floor(rep * 0.05))
  q95 <- max(1, floor(rep * 0.95))
  q975 <- max(1, floor(rep * 0.975))

  bound <- matrix(0, nrow = 4, ncol = 2)
  bound[1, ] <- floor(c(dx1, dx2) - obd[q025, 1:2] / P1j) + c(1, 1)
  bound[2, ] <- floor(c(dx1, dx2) - obd[q05, 1:2] / P1j) + c(1, 1)
  bound[3, ] <- floor(c(dx1, dx2) - obd[q95, 1:2] / P1j)
  bound[4, ] <- floor(c(dx1, dx2) - obd[q975, 1:2] / P1j)
  bound
}


#' Estimates the Ordered Breaks Model
#'
#' Estimates the ordered breaks model built by `ordspec()`.
#'
#' @noRd
mainop <- function(spec, trim = 0.2, hetq = TRUE,
                   max_break_gap = NULL, test_level = 0.05, ci_level = 0.95,
                   confidence_intervals = FALSE,
                   ci_grid_radius = 20, ci_replications = 5000,
                   ci_step = 1000, ci_grid_step = 0.02,
                   ci_seed = NULL, ci_verbose = FALSE,
                   tol = 1e-8, max_iter = 1000) {
  y <- spec$y
  z <- spec$z
  n <- ncol(y)
  bigt <- nrow(y)
  trm <- as.numeric(trim)[[1]]
  
  if (is.na(trm) || trm <= 0 || trm >= 0.5) {
    stop("`trim` must be in (0, 0.5).", call. = FALSE)
  }
  h <- floor(trm * bigt)
  if (is.null(max_break_gap)) {
    max_break_gap <- floor(sqrt(bigt))
  }
  trun <- as.integer(max_break_gap)[[1]]

  if (!is.numeric(test_level) || length(test_level) != 1L ||
      !any(abs(test_level - c(0.10, 0.05, 0.025, 0.01)) < 1e-12)) {
    stop("`test_level` must be one of 0.10, 0.05, 0.025, or 0.01.",
         call. = FALSE)
  }
  if (!is.numeric(ci_level) || length(ci_level) != 1L ||
      !any(abs(ci_level - c(0.90, 0.95)) < 1e-12)) {
    stop("`ci_level` must be either 0.90 or 0.95 (the coverage probability).",
         call. = FALSE)
  }

  ci_grid_radius <- as.integer(ci_grid_radius)[[1]]
  ci_replications <- as.integer(ci_replications)[[1]]
  ci_step <- as.integer(ci_step)[[1]]
  if (!is.numeric(ci_grid_step) || length(ci_grid_step) != 1 ||
      !is.finite(ci_grid_step) || ci_grid_step <= 0) {
    stop("`ci_grid_step` must be a positive finite number.", call. = FALSE)
  }
  if (!is.logical(ci_verbose) || length(ci_verbose) != 1 || is.na(ci_verbose)) {
    stop("`ci_verbose` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  transformed <- transf(y, z, spec$S)
  ordered <- brorder(
    maty = transformed$maty, matz = transformed$matz, n = n,
    sel1 = spec$sel1, sel2 = spec$sel2,
    h = h, bigt = bigt, trun = trun, tol = tol, max_iter = max_iter
  )

  bound <- NULL
  if (confidence_intervals) {
    bound <- interval3(
      matz = transformed$matz, n = n, start = 1, last = bigt,
      dx1 = ordered$brdate[[1]], dx2 = ordered$brdate[[2]],
      delt1 = ordered$delt1, delt2 = ordered$delt2,
      res = ordered$umat, bigM = ci_grid_radius, rep = ci_replications,
      hetq = hetq, step = ci_step, interv = ci_grid_step,
      seed = ci_seed, verbose = ci_verbose
    )
    # interval3 returns rows 1,2 as the 95 and 90 percent upper bounds and rows
    # 4,3 as the matching lower bounds; keep the requested level only.
    rows <- if (isTRUE(all.equal(ci_level, 0.95))) c(4, 1) else c(3, 2)
    bound <- cbind(bound[rows[[1]], ], bound[rows[[2]], ])
    colnames(bound) <- sprintf("%g%% %s", 100 * ci_level, c("lower", "upper"))
  }

  # cvlr() returns the 10, 5, 2.5 and 1 percent values in that order
  cvi <- which.min(abs(test_level - c(0.10, 0.05, 0.025, 0.01)))

  est <- list(
    breaks = ordered$brdate,
    beta1 = ordered$beta1,
    beta2 = ordered$beta2,
    residuals = ordered$umat,
    delta1 = ordered$delt1,
    delta2 = ordered$delt2,
    ci_breaks = bound
  )
 
  attr(est, "coef_names") <- spec$coef_names
  attr(est, "changing1") <- spec$coef_names[spec$sel1]
  attr(est, "changing2") <- spec$coef_names[spec$sel2]
  attr(est, "ci_level") <- ci_level
  class(est) <- c("ordbreak_estimation", "list")

  list(
    spec = list(
      trim = trm, bigt = bigt, n = n, h = h,
      q1 = spec$q1, q2 = spec$q2, n1 = spec$n1,
      changing1 = spec$coef_names[spec$sel1],
      changing2 = spec$coef_names[spec$sel2],
      hetq = hetq, trun = trun,
      bigM = ci_grid_radius, rep = ci_replications,
      confidence_intervals = confidence_intervals,
      test_level = test_level, ci_level = ci_level
    ),
    testing = list(
      sup_lr_ordered = ordered$lr,
      # the critical value counts the coefficients that change
      critical_values = as.numeric(
        cvlr(1, length(spec$sel1) + length(spec$sel2), h / bigt)
      )[cvi],
      truncated_breaks = ordered$truncated_brdate,
      H1like = ordered$H1like,
      H0like = ordered$H0like
    ),
    estimation = est
  )
}


#' Ordered Structural Breaks in a System of Equations
#'
#' Estimation of a multiple-equation system in which the coefficients of one
#' group of equations change before those of another group. The two break dates
#' are estimated jointly and are not restricted to coincide.
#'
#' @usage
#' ordbreak(set1, set2, data,
#'   coef_breaks1 = "all", coef_breaks2 = "all",
#'   trim = 0.2, hetq = TRUE, max_break_gap = NULL,
#'   test_level = 0.05, ci_level = 0.95,
#'   confidence_intervals = FALSE,
#'   ci_grid_radius = 20, ci_replications = 5000, ci_step = 1000,
#'   ci_grid_step = 0.02, ci_seed = NULL, ci_verbose = FALSE,
#'   tol = 1e-8, max_iter = 1000)
#'
#' @param set1 The equations whose coefficients change at the first date. A
#'   formula such as \code{y1 ~ x1 + x2} for a single equation, or
#'   \code{list(y1 ~ x1 + x2, y2 ~ x2 + x3)} for several. The number of
#'   regressors may differ across equations. An intercept is added to each
#'   equation unless the formula indicates otherwise.
#' @param set2 The equations whose coefficients change at the second date,
#'   written in the same way as \code{set1}. Each equation of the system belongs
#'   to exactly one of the two sets.
#' @param data Data frame used to evaluate the formulas.
#' @param coef_breaks1 An expression for which coefficients of \code{set1} are
#'   allowed to change at the first date. The default \code{"all"} allows all of
#'   them to change. Otherwise supply a list of formulas, one per equation of
#'   \code{set1}, where the left-hand side should specify the equation and the
#'   right-hand side selects the coefficients of that equation that may change:
#'   \itemize{
#'     \item \code{y1 ~ x1}: the intercept and the coefficient on \code{x1}
#'     change.
#'     \item \code{y1 ~ x1 - 1}: the coefficient on \code{x1} only.
#'     \item \code{y1 ~ 0}: no coefficient of that equation changes.
#'   }
#'   Coefficients left out are held constant across both dates. At least one
#'   coefficient of \code{set1} must be allowed to change.
#' @param coef_breaks2 The same, for \code{set2} and the second date.
#' @param trim Trimming parameter, which must lie in \eqn{(0, 0.5)}. The minimum
#'   segment length is \code{floor(trim * nrow(data))}.
#' @param hetq Logical. If \code{TRUE}, allow regressor moments to vary across
#'   regimes when computing confidence intervals for break dates.
#' @param max_break_gap Truncation used when computing the test statistic. If
#'   \code{NULL}, \code{floor(sqrt(T))} is used, as in the original GAUSS
#'   program.
#' @param test_level Significance level at which the critical value is
#'   reported. One of \code{0.10}, \code{0.05}, \code{0.025} or \code{0.01}.
#' @param ci_level Coverage probability of the reported confidence intervals
#'   for the break dates. Either \code{0.95} or \code{0.90}.
#' @param confidence_intervals Logical. If \code{TRUE}, simulate confidence
#'   intervals for the two break dates. This is the expensive part of the
#'   procedure and is therefore off by default.
#' @param ci_grid_radius Half-width of the grid on which the limiting
#'   distribution is simulated.
#' @param ci_replications Number of simulation replications.
#' @param ci_step Number of steps used to discretize each Brownian motion.
#' @param ci_grid_step Spacing of the simulation grid.
#' @param ci_seed Optional seed for the simulation. The seed of the calling
#'   session is restored on exit.
#' @param ci_verbose Logical. If \code{TRUE}, print simulation progress.
#' @param tol Convergence tolerance of the feasible GLS iterations.
#' @param max_iter Maximum number of feasible GLS iterations.
#'
#' @return An object of class \code{"ordbreak"} containing:
#' \describe{
#'   \item{\code{spec}}{The setup: the sample size
#'     \code{bigt}, the number of equations \code{n}, the minimum segment length
#'     \code{h}, the number of coefficients in each set \code{q1} and \code{q2},
#'     the number of equations in the first set \code{n1}, the names of the
#'     coefficients allowed to change at each date, \code{changing1} and
#'     \code{changing2}, the arguments \code{trim}, \code{hetq} and
#'     \code{trun}, the input formulas split into \code{set1} and \code{set2},
#'     the equation and coefficient names, the retained data rows
#'     \code{y_rows} after discarding missing values, and the selection matrix
#'     \code{S}.}
#'   \item{\code{testing}}{The SupLR statistic for ordered breaks,
#'     \code{sup_lr_ordered}, its \code{critical_values} at \code{test_level},
#'     the pair \code{truncated_breaks} the statistic is computed at, and the
#'     log-likelihoods \code{H1like} and \code{H0like} it is built from. The
#'     critical value is indexed by the number of coefficients allowed to
#'     change, not by the total number of coefficients.}
#'   \item{\code{estimation}}{The estimated \code{breaks}, given as rows of
#'     \code{data}; the coefficients \code{beta1} and \code{beta2} of the two
#'     regimes; the \code{residuals} of the model at the estimated breaks, one
#'     column per equation; and the coefficient shifts \code{delta1}
#'     and \code{delta2} occurring at the two dates; and \code{ci_breaks}, a
#'     two by two matrix when \code{confidence_intervals = TRUE} and
#'     \code{NULL} otherwise, also given as rows of \code{data}. Its rows are
#'     the two break dates and its columns the lower and upper bounds at
#'     \code{ci_level}.}
#' }
#' Note that the ordering is imposed rather than tested. The equations in
#' \code{set1} are taken to break first, so exchanging the two sets specifies a
#' different model. Both specifications are tested against the same no-break
#' null, and their values of \code{sup_lr_ordered} can be compared.
#'
#' @details
#' The coefficients are stacked equation by equation in the order the formulas
#' are given, so those of \code{set1} occupy the leading \code{q1} positions and
#' those of \code{set2} the remaining \code{q2}. \code{\link{ordcheck}} reports
#' this layout without estimating anything, which is worth doing before a long
#' run.
#'
#' The test is that of Qu and Perron (2007, Section 6.4) for locally ordered
#' breaks. The test for globally ordered breaks in their Section 6.5, whose
#' limiting distribution differs, is not implemented.
#'
#' @references
#' Qu, Z., and Perron, P. (2007). Estimating and testing structural changes
#' in multivariate regressions. \emph{Econometrica}, 75(2), 459--502.
#'
#' @examples
#' set.seed(13)
#' bigt <- 60
#' x1 <- rnorm(bigt)
#' x2 <- rnorm(bigt)
#' y1 <- 1 + 0.5 * x1 + rnorm(bigt, sd = 0.3)
#' y2 <- 1 + 0.5 * x2 + rnorm(bigt, sd = 0.3)
#' y1[21:bigt] <- y1[21:bigt] + 1.5 * x1[21:bigt]
#' y2[36:bigt] <- y2[36:bigt] - 1.5 * x2[36:bigt]
#' dat <- data.frame(y1 = y1, y2 = y2, x1 = x1, x2 = x2)
#'
#' fit <- ordbreak(y1 ~ x1, y2 ~ x2, data = dat, trim = 0.2)
#' fit$estimation$breaks
#'
#' @export
ordbreak <- function(set1, set2, data,
                     coef_breaks1 = "all", coef_breaks2 = "all",
                     trim = 0.2, hetq = TRUE, max_break_gap = NULL,
                     test_level = 0.05, ci_level = 0.95,
                     confidence_intervals = FALSE,
                     ci_grid_radius = 20, ci_replications = 5000,
                     ci_step = 1000, ci_grid_step = 0.02,
                     ci_seed = NULL, ci_verbose = FALSE,
                     tol = 1e-8, max_iter = 1000) {
  if (inherits(set1, "formula")) set1 <- list(set1)
  if (inherits(set2, "formula")) set2 <- list(set2)

  spec <- ordspec(set1, set2, data,
                  coef_breaks1 = coef_breaks1, coef_breaks2 = coef_breaks2)

  fit <- mainop(
    spec, trim = trim, hetq = hetq, max_break_gap = max_break_gap,
    test_level = test_level, ci_level = ci_level,
    confidence_intervals = confidence_intervals,
    ci_grid_radius = ci_grid_radius, ci_replications = ci_replications,
    ci_step = ci_step, ci_grid_step = ci_grid_step,
    ci_seed = ci_seed, ci_verbose = ci_verbose,
    tol = tol, max_iter = max_iter
  )
  # report dates as rows of `data`; an identity map unless rows were dropped
  yrow <- spec$y_rows
  nt <- length(yrow)
  fit$estimation$breaks <- yrow[fit$estimation$breaks]
  fit$testing$truncated_breaks <- yrow[fit$testing$truncated_breaks]
  if (!is.null(fit$estimation$ci_breaks)) {
    bd <- fit$estimation$ci_breaks
    nlo <- sum(bd[, 1] < 1 | bd[, 1] > nt)
    nhi <- sum(bd[, 2] < 1 | bd[, 2] > nt)
    if (nlo + nhi > 0) {
      side <- c(if (nlo) "lower", if (nhi) "upper")
      warning(
        sprintf(
          "The %s confidence bound%s fell outside the sample and %s truncated. The break dates are weakly identified.",
          paste(side, collapse = " and "),
          if (length(side) > 1) "s" else "",
          if (length(side) > 1) "were" else "was"
        ),
        call. = FALSE
      )
    }
    fit$estimation$ci_breaks[] <- yrow[pmin(pmax(bd, 1L), nt)]
  }
  fit$spec <- c(fit$spec, list(
    model_type = "ordered_formula",
    formulas = spec$formulas,
    set1 = spec$formulas[seq_len(spec$n1)],
    set2 = spec$formulas[(spec$n1 + 1L):length(spec$formulas)],
    y_names = colnames(spec$y),
    z_names = colnames(spec$z),
    eqz = spec$eqz,
    coef_names = spec$coef_names,
    y_rows = spec$y_rows,
    S = spec$S
  ))
  class(fit) <- c("ordbreak", "list")
  fit
}


#' Print an Ordered-Break Result
#'
#' @export
#' @noRd
print.ordbreak <- function(x, ...) {
  s <- x$spec
  #cat("mbreak ordered-break result\n")
  cat("------------------------------------------------------------------------\n")
  cat("The specifications for testing and estimation:\n")
  cat("-------------------------------------------\n")
  cat(sprintf("(1) Trimming = %.3f\n", s$trim))
  cat(sprintf("(2) T = %d\n", s$bigt))
  cat(sprintf("(3) Equations breaking first: %d, second: %d\n",
              s$n1, s$n - s$n1))
  cat(sprintf("(4) The number of coefficients is: q1 = %d, q2 = %d\n",
              s$q1, s$q2))
  cat(sprintf("(5) Maximum break gap = %d\n", s$trun))
  cat("-------------------------------------------\n")
  cat(sprintf("hetq = %s\n", s$hetq))
  cat("------------------------------------------------------------------------\n\n")

  stat <- as.numeric(x$testing$sup_lr_ordered)
  cv <- as.numeric(x$testing$critical_values)
  any_sig <- is.finite(stat) && is.finite(cv) && stat > cv
  cat(sprintf(" ordered SupLR   (%g%% level)\n", 100 * s$test_level))
  cat(sprintf("%-20s %11s %11s\n", "", "statistic", "crit. value"))
  cat(sprintf("%-20s %11.3f %11.3f%s\n", "    0 vs 2 breaks", stat, cv,
              if (any_sig) " *" else ""))
  if (!identical(as.integer(x$estimation$breaks),
                 as.integer(x$testing$truncated_breaks))) {
    cat(sprintf("    computed at breaks %s\n",
                paste(x$testing$truncated_breaks, collapse = ", ")))
  }
  cat("\n")

  if (any_sig) {
    cat(sprintf("* significant at the %g%% level\n\n", 100 * s$test_level))
  }
  print(x$estimation)
  invisible(x)
}


#' Print the Estimation Block of an Ordered-Break Result
#'
#' @export
#' @noRd
print.ordbreak_estimation <- function(x, ...) {
  cat("Output from the estimation procedure\n")
  cat("__________________________________\n\n")
  pbr(x$breaks, x$ci_breaks, ci_level = attr(x, "ci_level"))

  nm <- attr(x, "coef_names")
  b1 <- as.numeric(x$beta1)
  b2 <- as.numeric(x$beta2)
  at <- ifelse(nm %in% attr(x, "changing1"), "break 1",
               ifelse(nm %in% attr(x, "changing2"), "break 2", "fixed"))
  cat("\nThe estimated coefficients are:\n")
  cat(sprintf("%-20s %11s %11s %11s  %s\n", "", "before", "after", "change",
              "changes at"))
  for (i in seq_along(nm)) {
    cat(sprintf("%-20s %11.5f %11.5f %11.5f  %s\n", nm[i], b1[i], b2[i],
                b2[i] - b1[i], at[i]))
  }
  invisible(x)
}


#' Coefficient Layout of an Ordered-Break Specification
#'
#' Reports the coefficients of an ordered-break model, which set of equations they
#' belongs to, and which of them are allowed to change, without estimating
#' anything.
#'
#' @usage
#' ordcheck(set1, set2, data, coef_breaks1 = "all", coef_breaks2 = "all")
#'
#' @param set1 The equations whose coefficients change at the first date,
#'   written as in \code{\link{ordbreak}}.
#' @param set2 The equations whose coefficients change at the second date.
#' @param data Data frame used to evaluate the formulas.
#' @param coef_breaks1 Which coefficients of \code{set1} are allowed to change
#'   at the first date. See \code{\link{ordbreak}}.
#' @param coef_breaks2 The same, for \code{set2} and the second date.
#'
#' @return The layout is printed: one
#'   row per coefficient, giving its position in the coefficient vector,
#'   its name, the set it belongs to, whether it is allowed to change, and the
#'   date at which it changes.
#'
#' @examples
#' set.seed(13)
#' dat <- data.frame(y1 = rnorm(40), y2 = rnorm(40),
#'                   x1 = rnorm(40), x2 = rnorm(40))
#'
#' ordcheck(y1 ~ x1 + x2, y2 ~ x1 + x2, data = dat,
#'         coef_breaks1 = list(y1 ~ x1))
#'
#' @export
ordcheck <- function(set1, set2, data,
                     coef_breaks1 = "all", coef_breaks2 = "all") {
  if (inherits(set1, "formula")) set1 <- list(set1)
  if (inherits(set2, "formula")) set2 <- list(set2)

  spec <- ordspec(set1, set2, data,
                  coef_breaks1 = coef_breaks1, coef_breaks2 = coef_breaks2)
  q <- length(spec$coef_names)
  tab <- data.frame(
    column = seq_len(q),
    coefficient = spec$coef_names,
    set = ifelse(seq_len(q) <= spec$q1, 1L, 2L),
    changes = ifelse(seq_len(q) %in% c(spec$sel1, spec$sel2), "yes", "no"),
    at = ifelse(seq_len(q) %in% spec$sel1, "break 1",
                ifelse(seq_len(q) %in% spec$sel2, "break 2", "-")),
    stringsAsFactors = FALSE
  )
  cat("equations in set 1:", spec$n1,
      " | set 2:", length(spec$formulas) - spec$n1, "\n")
  cat("q1 =", spec$q1, " q2 =", spec$q2, " q =", q, "\n\n")
  print(tab, row.names = FALSE)
  invisible(spec)
}


#' Iterated Feasible GLS for the Stacked System
#'
#' Estimates the stacked system by iterated feasible GLS, alternating the
#' GLS coefficients and Sigma = U'U / T until Sigma converges. brorder()
#' uses it for each candidate pair of break dates (encoded as zeroed
#' columns in segx) and for the model without breaks. The first step is
#' OLS, or GLS weighted by init_s^(-1) if init_s is supplied; brorder()
#' passes the Sigma of the previous pair to reduce the number of iterations.
#'
#' @noRd
ofgls <- function(segx, segy, n, init_s = NULL, tol = 1e-8, max_iter = 1000) {
  if (is.null(init_s)) {
    b <- pinv(t(segx) %*% segx) %*% t(segx) %*% segy
  } else {
    wx <- matrix(pinv(init_s) %*% matrix(segx, n), nrow(segx))
    b <- pinv(crossprod(segx, wx)) %*% crossprod(wx, segy)
  }

  res <- segy - segx %*% b
  umat <- matrix(as.vector(res), ncol = n, byrow = TRUE)
  vvar <- t(umat) %*% umat / nrow(umat)

  upd <- function(prev_sigma) {
    wx  <- matrix(pinv(prev_sigma) %*% matrix(segx, n), nrow(segx))
    b   <<- pinv(crossprod(segx, wx)) %*% crossprod(wx, segy)
    res <<- segy - segx %*% b
    umat <<- matrix(as.vector(res), ncol = n, byrow = TRUE)
    t(umat) %*% umat / nrow(umat)
  }
  fp <- fixed_point(vvar, upd, tol, max_iter, "ofgls()")

  list(beta = b, sigma = fp$value, residuals = res, umat = umat,
       iterations = fp$iterations, converged = fp$converged)
}
