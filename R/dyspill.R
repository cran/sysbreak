#' Diebold-Yilmaz Generalized Spillover Measures
#'
#' Computes horizon-specific spillover measures from a fitted structural-break
#' VAR model. The forecast-error variance decomposition uses generalized shocks
#' and is therefore invariant to variable ordering.
#'
#' @usage
#' dyspill(fit, h = 10, boot = 0, seed = NULL)
#'
#' @param fit A fitted object returned by \code{varbreak()} or
#'   \code{netbreak()}.
#' @param h Positive integer for forecast horizon.
#' @param boot Number of bootstrap replications for standard errors and
#'   confidence intervals; \code{0} for none. Available for \code{varbreak()}
#'   and \code{netbreak()} fits.
#' @param seed Optional seed for the bootstrapping. The seed of the calling session
#'   is restored on exit.
#'
#' @return An object of class \code{"dyspill"} containing:
#' \describe{
#'   \item{\code{horizon}}{The forecast horizon \code{h}.}
#'   \item{\code{model}}{The source model type.}
#'   \item{\code{variables}}{Names of the endogenous variables.}
#'   \item{\code{break_grid}}{A data frame with one row per regime: the
#'     \code{regime} number; its first and last observation, \code{start} and
#'     \code{end}, counted in the estimation sample (which omits the first
#'     \code{lags} observations); and the coefficient and covariance regimes
#'     it belongs to, \code{coef_regime} and \code{cov_regime}.}
#'   \item{\code{regimes}}{A list with one entry per regime. Each entry
#'     contains:
#'     \itemize{
#'       \item \code{input_rows}: the first and last row of the data in the
#'         regime.
#'       \item \code{spectral_radius}: the largest modulus of the eigenvalues
#'         of the companion matrix of the VAR; the VAR is stable when it is
#'         below one.
#'       \item \code{fevd}: the generalized forecast-error variance
#'         decomposition, with rows normalized to sum to one.
#'       \item \code{total}: the total spillover index.
#'       \item \code{from}, \code{to} and \code{net}: the directional
#'         spillovers of each variable.
#'       \item \code{net_pairwise}: the matrix of net pairwise spillovers.
#'       \item \code{benchmark_total} and \code{benchmark_net}: \code{total}
#'         and \code{net} recomputed with all lag coefficients set to zero, so
#'         that they reflect the error correlations only (see Details).
#'       \item With \code{boot > 0}: the bootstrap standard errors
#'         \code{total_se} and \code{net_se}, and the percentile intervals
#'         \code{total_ci} and \code{net_ci}.
#'     }}
#'   \item{\code{boot}}{Returned only when \code{boot > 0}. A list with:
#'     \itemize{
#'       \item \code{replications}: the number of bootstrap replications.
#'       \item \code{level}: the coverage probability of the intervals, taken
#'         from \code{ci_level} of the fit.
#'       \item \code{total}: the bootstrap draws of the total index, a matrix
#'         with one row per replication and one column per regime.
#'       \item \code{net}: the bootstrap draws of the net spillovers, an array
#'         of replications by variables by regimes.
#'       \item \code{differences}: a data frame with one row per pair of
#'         regimes, giving the \code{difference} in the total index (the later
#'         regime minus the earlier one), its bootstrap standard error
#'         \code{se}, and the percentile interval \code{lower} and
#'         \code{upper}.
#'       \item \code{dropped}: the number of draws dropped in each regime
#'         because the refitted VAR was not stable.
#'     }}
#' }
#'
#' @details
#' For a regime-specific VAR,
#' \deqn{y_t = c_r + \sum_{\ell=1}^p A_{\ell,r}y_{t-\ell}+u_t,\qquad
#' E(u_tu_t')=\Sigma_r,}
#' let \eqn{\Psi_{0,r}=I} and
#' \deqn{\Psi_{s,r}=\sum_{\ell=1}^{\min(p,s)}
#' A_{\ell,r}\Psi_{s-\ell,r}.}
#' The generalized variance share from shock \eqn{j} to variable \eqn{i} is
#' \deqn{\theta_{ij}^{g}(H)=
#' \frac{\sigma_{jj}^{-1}\sum_{s=0}^{H-1}
#' (e_i'\Psi_s\Sigma e_j)^2}
#' {\sum_{s=0}^{H-1}e_i'\Psi_s\Sigma\Psi_s'e_i}.}
#' Rows are normalized to sum to one before spillover statistics are formed.
#' For variable \eqn{i}, \code{from} is the contribution received from all
#' other shocks, \code{to} is the contribution transmitted to all other
#' variables, and \code{net = to - from}; a positive net value identifies a
#' net transmitter. The total spillover index is the average off-diagonal
#' contribution, expressed in percent.
#'
#' The benchmark puts every lag matrix to zero and keeps \eqn{\Sigma_r}. The
#' measures then depend on the error correlations only; comparing them with
#' the fitted ones shows how much of the spillover comes from the dynamics.
#' With two variables the benchmark total is
#' \eqn{100\rho^2/(1+\rho^2)}{100 rho^2 / (1 + rho^2)}.
#'
#' The bootstrap is a recursive residual bootstrap within each regime: the
#' residuals of the regime are resampled, the series is rebuilt from the
#' fitted VAR of the regime, the model is refitted at the same break dates in
#' the same way as \code{fit}, and the measures are recomputed. The break dates
#' are treated as known. The intervals are percentile intervals, since the
#' draws are skewed in small samples. Draws in which a refitted VAR is not
#' stable are dropped.
#'
#' @references
#' Diebold, F. X., and Yilmaz, K. (2012). Better to give than to receive:
#' Predictive directional measurement of volatility spillovers.
#' \emph{International Journal of Forecasting}, 28(1), 57--66.
#'
#' @examples
#' set.seed(1)
#' bigt <- 120
#' y <- matrix(0, bigt, 2)
#' colnames(y) <- c("y1", "y2")
#' for (tt in 2:bigt) {
#'   a <- if (tt > 60) matrix(c(0.2, 0.6, 0.6, 0.2), 2) else diag(0.2, 2)
#'   y[tt, ] <- a %*% y[tt - 1, ] + rnorm(2, sd = 0.3)
#' }
#'
#' fit <- varbreak(y, lags = 1, m = 1, trim = 0.2)
#' dyspill(fit, h = 10, boot = 99, seed = 1)
#'
#' @export
dyspill <- function(fit, h = 10, boot = 0, seed = NULL) {
  h <- dyvalid(h)
  if (!is.numeric(boot) || length(boot) != 1L || is.na(boot) ||
      boot < 0 || boot != floor(boot)) {
    stop("`boot` must be a non-negative integer.", call. = FALSE)
  }
  if (boot > 0 && is.null(fit$estimation$regressors)) {
    stop("The bootstrap is available for varbreak() and netbreak() fits.",
         call. = FALSE)
  }
  model <- dyreg(fit)
  unstable <- integer()

  regimes <- lapply(seq_along(model$regimes), function(k) {
    regime <- model$regimes[[k]]
    fevd <- dyfevd(regime$a, regime$sigma, h)
    stability <- dystab(regime$a)
    stat <- dystat(fevd$normalized)
    # the same measures with every lag matrix set to zero: shows what the error
    # correlations of the regime produce without any dynamics
    zero <- lapply(regime$a, function(a) a * 0)
    bench <- dystat(dyfevd(zero, regime$sigma, h)$normalized)
    if (!stability$stable) {
      unstable <<- c(unstable, k)
    }

    list(
      input_rows = c(
        start = model$y_rows[[regime$start]],
        end = model$y_rows[[regime$end]]
      ),
      spectral_radius = stability$radius,
      fevd = fevd$normalized,
      total = stat$total,
      from = stat$from,
      to = stat$to,
      net = stat$net,
      net_pairwise = stat$pairwise,
      benchmark_total = bench$total,
      benchmark_net = bench$net
    )
  })

  if (length(unstable)) {
    warning(
      sprintf(
        "VAR dynamics are unstable in spillover regime(s) %s; finite-horizon measures were still computed.",
        paste(unstable, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  result <- list(
    horizon = h,
    model = model$type,
    variables = model$names,
    break_grid = model$grid,
    regimes = regimes
  )

  if (boot > 0) {
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
    draws <- dyboot(fit, h, boot)
    level <- if (is.null(fit$spec$ci_level)) 0.95 else fit$spec$ci_level
    probs <- c(1 - level, 1 + level) / 2
    for (k in seq_along(regimes)) {
      net_k <- matrix(draws$net[, , k], nrow = boot,
                      dimnames = list(NULL, dimnames(draws$net)[[2]]))
      result$regimes[[k]]$total_se <- stats::sd(draws$total[, k], na.rm = TRUE)
      result$regimes[[k]]$total_ci <- stats::quantile(draws$total[, k], probs,
                                                      na.rm = TRUE)
      result$regimes[[k]]$net_se <- apply(net_k, 2, stats::sd, na.rm = TRUE)
      result$regimes[[k]]$net_ci <- t(apply(net_k, 2, stats::quantile, probs,
                                            na.rm = TRUE))
    }

    # every pair of regimes: the later one minus the earlier one
    total <- vapply(regimes, function(x) x$total, numeric(1))
    pairs <- utils::combn(length(regimes), 2)
    differences <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(j) {
      i <- pairs[1, j]; k <- pairs[2, j]
      d <- draws$total[, k] - draws$total[, i]
      ci <- stats::quantile(d, probs, na.rm = TRUE)
      data.frame(regimes = paste0(k, "-", i), difference = total[k] - total[i],
                 se = stats::sd(d, na.rm = TRUE), lower = ci[[1]],
                 upper = ci[[2]])
    }))

    dropped <- colSums(is.na(draws$total))
    if (any(dropped > 0)) {
      warning(
        sprintf("Bootstrap draws were dropped because a refitted VAR was not stable: %s.",
                paste(sprintf("regime %d, %d", which(dropped > 0),
                              dropped[dropped > 0]), collapse = "; ")),
        call. = FALSE
      )
    }
    result$boot <- list(replications = boot, level = level,
                        total = draws$total, net = draws$net,
                        differences = differences, dropped = dropped)
  }

  class(result) <- c("dyspill", class(result))
  result
}

#' Bootstrap of the Spillover Measures
#'
#' Recursive residual bootstrap within each regime. In each iteration, the draws
#' are used to build each
#' regime from its fitted VAR and resampled residuals, refits the model at the
#' same break dates with the estimates in `fit`, and recomputes the total and
#' net spillovers for each regime.
#'
#' @noRd
dyboot <- function(fit, h, B) {
  est  <- fit$estimation
  spec <- fit$spec
  res  <- est$residuals
  z    <- est$regressors
  n    <- ncol(res)
  bigt <- nrow(res)
  br   <- match(est$breaks, spec$y_rows)   # breaks as rows of `res`
  nb   <- length(br)
  ends <- c(0L, br, bigt)

  # regressors to coefficients: the identity for a VAR, the network structure
  # for a network VAR
  S <- if (identical(spec$model_type, "network_var")) {
    nvarS(spec$W, spec$lags, spec$intercept, spec$y_names, colnames(z))$S
  } else {
    diag(n * ncol(z))
  }
  q <- ncol(S)

  # the restriction the fit imposed across regimes
  changing <- spec$coef_names %in% cbchg(spec)
  ols <- identical(spec$model_type, "var") && isTRUE(spec$brbeta) &&
    all(changing)
  R <- if (!isTRUE(spec$brbeta)) {
    defR(nb, q, brbeta = FALSE, brv = TRUE)
  } else if (all(changing)) {
    defR(nb, q, brbeta = TRUE, brv = spec$brv)
  } else {
    buildR(changing, spec$coef_names, nb)
  }

  # the coefficients of each regime on the regressors, a column every equation
  theta <- matrix(est$coefficients, ncol = nb + 1)
  coef <- lapply(seq_len(nb + 1), function(r) {
    matrix(S %*% theta[, r], nrow = ncol(z))
  })
  lag_cols <- if (isTRUE(spec$intercept)) -1 else seq_len(ncol(z))

  # resampled rows, drawn regime by regime
  idx <- lapply(seq_len(nb + 1), function(r) {
    len <- ends[r + 1] - ends[r]
    replicate(B, sample.int(len, len, replace = TRUE))
  })

  total <- matrix(NA_real_, B, nb + 1)
  net <- array(NA_real_, c(B, n, nb + 1),
               dimnames = list(NULL, spec$y_names, NULL))
  for (b in seq_len(B)) {
    ys <- matrix(0, bigt, n)
    zs <- matrix(0, bigt, ncol(z))
    for (r in seq_len(nb + 1)) {
      rows <- (ends[r] + 1):ends[r + 1]
      e <- res[rows, , drop = FALSE][idx[[r]][, b], , drop = FALSE]
      zt <- z[rows[1], ]   # the observed lags at the start of the regime
      for (t in seq_along(rows)) {
        zs[rows[t], ] <- zt
        ys[rows[t], ] <- zt %*% coef[[r]] + e[t, ]
        zt[lag_cols] <- c(ys[rows[t], ], zt[lag_cols])[seq_along(zt[lag_cols])]
      }
    }

    fit_b <- fit
    if (ols) {
      # every coefficient changes and the regressors are common to all
      # equations, so the estimator reduces to OLS within each regime
      b_r <- lapply(seq_len(nb + 1), function(r) {
        rows <- (ends[r] + 1):ends[r + 1]
        solve(crossprod(zs[rows, ]), crossprod(zs[rows, ], ys[rows, ]))
      })
      u <- ys
      for (r in seq_len(nb + 1)) {
        rows <- (ends[r] + 1):ends[r + 1]
        u[rows, ] <- ys[rows, ] - zs[rows, ] %*% b_r[[r]]
      }
      sigma <- lapply(seq_len(nb + 1), function(r) {
        rows <- if (isTRUE(spec$brv)) (ends[r] + 1):ends[r + 1] else seq_len(bigt)
        crossprod(u[rows, ]) / length(rows)
      })
      fit_b$estimation$coefficients <- unlist(lapply(b_r, as.vector))
      fit_b$estimation$covariance <- do.call(rbind, sigma)
    } else {
      tr <- transf(ys, zs, S)
      refit <- tryCatch(
        r_estim(tr$maty, tr$matz, bigt, n, nb, br, R,
                brbeta = spec$brbeta, brv = spec$brv),
        error = function(e) NULL
      )
      if (is.null(refit)) next
      fit_b$estimation$coefficients <- as.numeric(refit$nbeta)
      fit_b$estimation$covariance <- refit$nvv
    }
    model_b <- dyreg(fit_b)
    for (k in seq_along(model_b$regimes)) {
      regime <- model_b$regimes[[k]]
      if (!dystab(regime$a)$stable) next
      stat <- dystat(dyfevd(regime$a, regime$sigma, h)$normalized)
      total[b, k] <- stat$total
      net[b, , k] <- stat$net
    }
  }
  list(total = total, net = net)
}

#' Print Diebold-Yilmaz Spillover Measures
#'
#' @param x An object returned by \code{dyspill()}.
#' @param ... Additional arguments, currently ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @noRd
#' @export
print.dyspill <- function(x, ...) {
  cat("Diebold-Yilmaz generalized spillover measures\n")
  cat(sprintf("Model: %s | Horizon: %d\n", x$model, x$horizon))
  overview <- do.call(rbind, lapply(seq_along(x$regimes), function(k) {
    regime <- x$regimes[[k]]
    data.frame(
      regime = k,
      sample = paste0(regime$input_rows[[1]], ":", regime$input_rows[[2]]),
      spectral_radius = regime$spectral_radius,
      stable = is.finite(regime$spectral_radius) &&
        regime$spectral_radius < 1,
      total = regime$total,
      benchmark = regime$benchmark_total,
      check.names = FALSE
    )
  }))
  overview$spectral_radius <- round(overview$spectral_radius, 4)
  overview$total <- round(overview$total, 2)
  overview$benchmark <- round(overview$benchmark, 2)
  if (!is.null(x$boot)) {
    overview$se <- round(vapply(x$regimes, function(r) r$total_se, 0), 2)
    overview$lower <- round(vapply(x$regimes, function(r) r$total_ci[[1]], 0), 2)
    overview$upper <- round(vapply(x$regimes, function(r) r$total_ci[[2]], 0), 2)
  }
  print(overview, row.names = FALSE)
  if (!is.null(x$boot)) {
    cat(sprintf("\nBootstrap: %d replications, %g%% percentile intervals\n",
                x$boot$replications, 100 * x$boot$level))
    cat("Differences in the total index between regimes:\n")
    d <- x$boot$differences
    d[-1] <- round(d[-1], 2)
    print(d, row.names = FALSE)
  }
  invisible(x)
}

#' Validate Diebold-Yilmaz Inputs
#'
#' @noRd
dyvalid <- function(h, a = NULL, sigma = NULL) {
  if (!is.numeric(h) || length(h) != 1L || is.na(h) || !is.finite(h) ||
      h != floor(h) || h < 1) {
    stop("`h` must be a positive integer.", call. = FALSE)
  }
  h <- as.integer(h)
  if (is.null(a) && is.null(sigma)) {
    return(h)
  }
  if (!is.list(a) || !length(a)) {
    stop("At least one VAR lag matrix is required.", call. = FALSE)
  }
  n <- nrow(a[[1]])
  valid_a <- vapply(a, function(mat) {
    is.matrix(mat) && is.numeric(mat) && identical(dim(mat), c(n, n)) &&
      all(is.finite(mat))
  }, logical(1))
  if (!length(n) || is.null(n) || n < 1L || !all(valid_a)) {
    stop("All VAR lag matrices must be finite numeric square matrices.",
         call. = FALSE)
  }
  sigma <- as.matrix(sigma)
  if (!is.numeric(sigma) || !identical(dim(sigma), c(n, n)) ||
      any(!is.finite(sigma))) {
    stop("The innovation covariance must be a finite numeric square matrix.",
         call. = FALSE)
  }
  scale <- max(1, max(abs(sigma)))
  if (max(abs(sigma - t(sigma))) > 1e-8 * scale) {
    stop("The innovation covariance must be symmetric.", call. = FALSE)
  }
  sigma <- (sigma + t(sigma)) / 2
  if (any(diag(sigma) <= sqrt(.Machine$double.eps) * scale)) {
    stop("The innovation covariance must have a positive diagonal.",
         call. = FALSE)
  }
  list(h = h, a = a, sigma = sigma, n = n)
}

#' VAR Moving-Average Coefficients
#'
#' @noRd
dyma <- function(a, h) {
  checked <- dyvalid(h, a, diag(nrow(a[[1]])))
  n <- checked$n
  p <- length(a)
  psi <- vector("list", h)
  psi[[1L]] <- diag(n)
  if (h > 1L) {
    for (s in 1:(h - 1L)) {
      value <- matrix(0, n, n)
      for (lag in seq_len(min(p, s))) {
        value <- value + a[[lag]] %*% psi[[s - lag + 1L]]
      }
      psi[[s + 1L]] <- value
    }
  }
  psi
}

#' VAR Stability
#'
#' @noRd
dystab <- function(a) {
  n <- nrow(a[[1]])
  p <- length(a)
  companion <- do.call(cbind, a)
  if (p > 1L) {
    companion <- rbind(
      companion,
      cbind(diag(n * (p - 1L)), matrix(0, n * (p - 1L), n))
    )
  }
  radius <- max(Mod(eigen(companion, only.values = TRUE)$values))
  list(radius = radius, stable = is.finite(radius) && radius < 1)
}

#' Generalized Forecast-Error Variance Decomposition
#'
#' @noRd
dyfevd <- function(a, sigma, h) {
  checked <- dyvalid(h, a, sigma)
  a <- checked$a
  sigma <- checked$sigma
  n <- checked$n
  psi <- dyma(a, checked$h)
  numerator <- matrix(0, n, n)
  denominator <- numeric(n)

  for (response in psi) {
    response_sigma <- response %*% sigma
    numerator <- numerator + sweep(
      response_sigma^2, 2L, diag(sigma), "/"
    )
    denominator <- denominator + diag(
      response_sigma %*% t(response)
    )
  }
  if (any(!is.finite(denominator)) ||
      any(denominator <= sqrt(.Machine$double.eps))) {
    stop("The forecast-error variance denominator is not positive.",
         call. = FALSE)
  }
  raw <- sweep(numerator, 1L, denominator, "/")
  row_total <- rowSums(raw)
  if (any(!is.finite(row_total)) ||
      any(row_total <= sqrt(.Machine$double.eps))) {
    stop("The generalized FEVD cannot be row-normalized.", call. = FALSE)
  }
  normalized <- sweep(raw, 1L, row_total, "/")
  names <- rownames(sigma)
  if (is.null(names)) {
    names <- paste0("y", seq_len(n))
  }
  dimnames(raw) <- list(names, names)
  dimnames(normalized) <- list(names, names)
  list(normalized = normalized)
}

#' Diebold-Yilmaz Spillover Statistics
#'
#' @noRd
dystat <- function(fevd) {
  n <- nrow(fevd)
  from <- 100 * (rowSums(fevd) - diag(fevd))
  to <- 100 * (colSums(fevd) - diag(fevd))
  net <- to - from
  pairwise <- 100 * (t(fevd) - fevd)
  diag(pairwise) <- 0
  names(from) <- rownames(fevd)
  names(to) <- colnames(fevd)
  names(net) <- rownames(fevd)
  dimnames(pairwise) <- dimnames(fevd)
  list(
    total = sum(from) / n,
    from = from,
    to = to,
    net = net,
    pairwise = pairwise
  )
}

#' Extract the VAR of Each Regime from a Fit
#'
#' Returns for each regime the lag coefficient matrices and the error
#' covariance matrix for`dyspill()`. `dyreg.default()` handles fits
#' from `varbreak()` and `netbreak()`.
#'
#' @noRd
dyreg <- function(fit) {
  UseMethod("dyreg")
}

#' @exportS3Method
#' @noRd
dyreg.default <- function(fit) {
  type <- fit$spec$model_type
  if (identical(type, "var")) {
    return(dyvar(fit))
  }
  if (identical(type, "network_var")) {
    return(dynet(fit))
  }
  stop(
    "`fit` must be returned by varbreak() or netbreak().",
    call. = FALSE
  )
}

#' Regimes from the Coefficient and Covariance Break Dates
#'
#' Partitions the sample at the union of the coefficient and covariance break
#' dates. Each row of the returned data frame corresponds to one segment of
#' this partition and saves its first and last observations together with
#' the indices of the coefficient regime and the covariance regime in effect.
#' For `varbreak()` and `netbreak()`, the two sets of break dates
#' coincide, and the partition reduces to the regimes of the fitted model.
#'
#' @noRd
dygrid <- function(coef_breaks, cov_breaks, bigt) {
  clean <- function(breaks) {
    breaks <- sort(unique(as.integer(breaks)))
    breaks[!is.na(breaks) & breaks > 0L & breaks < bigt]
  }
  coef_breaks <- clean(coef_breaks)
  cov_breaks <- clean(cov_breaks)
  breaks <- sort(unique(c(coef_breaks, cov_breaks)))
  endpoints <- c(0L, breaks, as.integer(bigt))
  start <- endpoints[-length(endpoints)] + 1L
  end <- endpoints[-1L]
  data.frame(
    regime = seq_along(start),
    start = start,
    end = end,
    coef_regime = findInterval(start - 1L, coef_breaks) + 1L,
    cov_regime = findInterval(start - 1L, cov_breaks) + 1L
  )
}

#' Lag Coefficient Matrices of a VAR Regime
#'
#' Converts the coefficient vector of one regime into the matrices
#' \eqn{A_1, \ldots, A_p} of the VAR. The vector is ordered by equation: the
#' coefficients of the first equation precede those of the second, and within
#' each equation the intercept, when present, precedes the coefficients on lag
#' 1, lag 2, and so on. Row \eqn{i} of \eqn{A_\ell} holds the coefficients of
#' equation \eqn{i} on the variables at lag \eqn{\ell}.
#'
#' @noRd
dyamat <- function(beta, n, p, intercept, names) {
  beta <- as.numeric(beta)
  if (length(beta) %% n != 0L) {
    stop("The fitted coefficient vector is incompatible with the VAR size.",
         call. = FALSE)
  }
  q <- length(beta) / n
  offset <- as.integer(isTRUE(intercept))
  if (q < offset + p * n) {
    stop("The fitted object does not contain all requested VAR lag matrices.",
         call. = FALSE)
  }
  coefficients <- matrix(beta, nrow = n, byrow = TRUE)
  lapply(seq_len(p), function(lag) {
    columns <- offset + (lag - 1L) * n + seq_len(n)
    matrix <- coefficients[, columns, drop = FALSE]
    dimnames(matrix) <- list(names, names)
    matrix
  })
}

#' Error Covariance Matrix of a Given Regime
#'
#' Returns the \eqn{n \times n} error covariance matrix of regime
#' \code{regime}. The fitted covariance matrices are stored either as a list
#' or as a matrix with the regime matrices stacked vertically, where regime 
#' \eqn{r} occupies rows \eqn{(r - 1)n + 1} to \eqn{rn}. A
#' single \eqn{n \times n} matrix is a covariance common to all regimes and is
#' returned for every regime.
#'
#' @noRd
dycov <- function(covariance, n, regime, names = NULL) {
  if (is.list(covariance)) {
    if (regime > length(covariance)) {
      stop("The covariance regime is missing from the fitted object.",
           call. = FALSE)
    }
    matrix <- as.matrix(covariance[[regime]])
  } else {
    covariance <- as.matrix(covariance)
    if (nrow(covariance) == n) {
      matrix <- covariance
    } else {
      rows <- (regime - 1L) * n + seq_len(n)
      if (max(rows) > nrow(covariance)) {
        stop("The covariance regime is missing from the fitted object.",
             call. = FALSE)
      }
      matrix <- covariance[rows, , drop = FALSE]
    }
  }
  if (!is.null(names)) {
    dimnames(matrix) <- list(names, names)
  }
  matrix
}

#' Regime-Specific VARs from a varbreak Fit
#'
#' Constructs, for each regime of a \code{varbreak()} fit, the lag coefficient
#' matrices and the error covariance matrix required by \code{dyspill()},
#' together with the first and last observation of the regime. Returns the
#' regimes, the partition produced by \code{dygrid()}, and the dimensions of
#' the system.
#'
#' @noRd
dyvar <- function(fit) {
  if (is.null(fit$estimation)) {
    stop("The varbreak fit has no selected estimation result.", call. = FALSE)
  }
  estimation <- fit$estimation
  spec <- fit$spec
  names <- spec$y_names
  n <- length(names)
  bigt <- spec$bigt
  # `breaks` are reported as rows of `y`; the regime grid runs on the
  # estimation sample, which drops the first `lags` observations
  brk <- match(estimation$breaks, spec$y_rows)
  coef_breaks <- if (isTRUE(spec$brbeta)) {
    brk
  } else {
    integer()
  }
  cov_breaks <- if (isTRUE(spec$brv)) {
    brk
  } else {
    integer()
  }
  grid <- dygrid(coef_breaks, cov_breaks, bigt)
  # Reshape by the number of breaks actually estimated, not by `spec$m`, which
  # is only the maximum allowed: the wrapper may have selected fewer.
  nb <- length(estimation$breaks)
  q <- length(estimation$coefficients) / (nb + 1L)
  coefficients <- byregime(estimation$coefficients, nb, q)

  regimes <- lapply(seq_len(nrow(grid)), function(k) {
    row <- grid[k, ]
    list(
      start = row$start,
      end = row$end,
      coef_regime = row$coef_regime,
      cov_regime = row$cov_regime,
      vol_regime = NA_integer_,
      corr_regime = NA_integer_,
      a = dyamat(
        coefficients[row$coef_regime, ], n, spec$lags,
        spec$intercept, names
      ),
      sigma = dycov(estimation$covariance, n, row$cov_regime, names)
    )
  })
  list(
    type = "varbreak", n = n, p = spec$lags, bigt = bigt,
    names = names, y_rows = spec$y_rows,
    grid = grid, regimes = regimes
  )
}

#' Regime-Specific VARs from a netbreak Fit
#'
#' Constructs, for each regime of a \code{netbreak()} fit, the lag coefficient
#' matrices and the error covariance matrix required by \code{dyspill()},
#' together with the first and last observation of the regime. The lag
#' matrices follow from the network VAR restriction
#' \eqn{A_\ell = \alpha_\ell W + \beta_\ell I}, where \eqn{\alpha_\ell} and
#' \eqn{\beta_\ell} are the estimated network and own-lag coefficients of the
#' regime. Returns the regimes, the partition produced by \code{dygrid()}, and
#' the dimensions of the system.
#'
#' @noRd
dynet <- function(fit) {
  if (is.null(fit$estimation)) {
    stop("The netbreak fit has no selected estimation result.", call. = FALSE)
  }
  estimation <- fit$estimation
  spec <- fit$spec
  names <- spec$y_names
  n <- length(names)
  bigt <- spec$bigt
  # As in `dyvar()`: `breaks` are rows of `y`, the grid runs on the estimation
  # sample, which drops the first `lags` observations
  brk <- match(estimation$breaks, spec$y_rows)
  coef_breaks <- if (isTRUE(spec$brbeta)) {
    brk
  } else {
    integer()
  }
  cov_breaks <- if (isTRUE(spec$brv)) {
    brk
  } else {
    integer()
  }
  grid <- dygrid(coef_breaks, cov_breaks, bigt)
  # As in `dyvar()`: reshape by the number of breaks estimated, not `spec$m`.
  nb <- length(estimation$breaks)
  q <- length(estimation$coefficients) / (nb + 1L)
  coefficients <- byregime(estimation$coefficients, nb, q)
  parameter_names <- spec$coef_names
  network_index <- match(
    paste0("network_lag", seq_len(spec$lags)), parameter_names
  )
  own_index <- match(
    paste0("own_lag", seq_len(spec$lags)), parameter_names
  )
  if (anyNA(c(network_index, own_index))) {
    stop("The netbreak lag parameters cannot be identified.", call. = FALSE)
  }

  regimes <- lapply(seq_len(nrow(grid)), function(k) {
    row <- grid[k, ]
    beta <- coefficients[row$coef_regime, ]
    a <- lapply(seq_len(spec$lags), function(lag) {
      matrix <- beta[[network_index[[lag]]]] * spec$W +
        beta[[own_index[[lag]]]] * diag(n)
      dimnames(matrix) <- list(names, names)
      matrix
    })
    list(
      start = row$start,
      end = row$end,
      coef_regime = row$coef_regime,
      cov_regime = row$cov_regime,
      vol_regime = NA_integer_,
      corr_regime = NA_integer_,
      a = a,
      sigma = dycov(estimation$covariance, n, row$cov_regime, names)
    )
  })
  list(
    type = "netbreak", n = n, p = spec$lags, bigt = bigt,
    names = names, y_rows = spec$y_rows,
    grid = grid, regimes = regimes
  )
}

#' Diebold-Yilmaz Spillover Measures Without Breaks
#'
#' Fits a VAR by OLS on the whole sample and computes the generalized
#' spillover measures of \code{dyspill()}, for comparison with the measures
#' of a model with breaks.
#'
#' @param y Numeric matrix, one column per variable.
#' @param lags Positive integer number of lags.
#' @param h Positive integer forecast horizon.
#' @param intercept Logical; include an intercept.
#'
#' @return A list with the row-normalized generalized \code{fevd}, the total
#'   spillover index \code{total}, and the directional vectors \code{from},
#'   \code{to} and \code{net}. The covariance divides by the number of
#'   observations.
#'
#' @examples
#' set.seed(1)
#' bigt <- 120
#' a <- matrix(c(0.2, 0.4, 0.4, 0.2), 2)
#' y <- matrix(0, bigt, 2)
#' colnames(y) <- c("y1", "y2")
#' for (tt in 2:bigt) {
#'   y[tt, ] <- a %*% y[tt - 1, ] + rnorm(2, sd = 0.3)
#' }
#'
#' sp <- dynobreak(y, lags = 1, h = 10)
#' sp$total
#' sp$net
#'
#' @export
dynobreak <- function(y, lags, h, intercept = TRUE) {
  y <- as.matrix(y)
  n <- ncol(y)
  bigt <- nrow(y)
  z <- do.call(cbind, lapply(seq_len(lags), function(l)
    y[(lags + 1 - l):(bigt - l), , drop = FALSE]))
  if (intercept) z <- cbind(1, z)
  yy <- y[(lags + 1):bigt, , drop = FALSE]
  b <- solve(crossprod(z), crossprod(z, yy))
  u <- yy - z %*% b
  sigma <- crossprod(u) / nrow(u)
  lag_rows <- if (intercept) b[-1, , drop = FALSE] else b
  a <- lapply(seq_len(lags), function(l) t(lag_rows[(l - 1) * n + 1:n, ]))
  fevd <- dyfevd(a, sigma, h)$normalized
  stat <- dystat(fevd)
  list(fevd = fevd, total = stat$total, from = stat$from, to = stat$to,
       net = stat$net)
}

