#' Methods for Displaying Structural-Break Fits
#'
#' Fits returned by \code{\link{sysbreak}},
#' \code{\link{varbreak}} and \code{\link{netbreak}}. All refer to the
#' estimates at the selected number of breaks.
#'
#' @param object A fit returned by \code{sysbreak()}, \code{varbreak()} or
#'   \code{netbreak()}.
#' @param parm Breaks (\code{type = "breaks"}) or coefficients
#'   (\code{type = "coef"}) to report, as indices or row names.
#' @param level Coverage probability. By default the \code{ci_level} of the
#'   fit, the level of the intervals shown by \code{print()} and
#'   \code{summary()}. For break dates either \code{0.95} or \code{0.90}, the
#'   two levels computed at estimation; for coefficients any value in
#'   \eqn{(0, 1)}.
#' @param type \code{"breaks"} for the break dates or \code{"coef"} for the
#'   regression coefficients.
#' @param ... Not used.
#'
#' @return
#' \describe{
#'   \item{\code{coef}}{A named vector of the regime coefficients: all
#'     coefficients of regime 1, then those of regime 2, and so on. Each name
#'     joins the regime and the coefficient label used by \code{summary()},
#'     as in \code{"regime1:y1:x"} (equation \code{y1}, regressor \code{x})
#'     or \code{"regime2:mu:node3"} for \code{netbreak}. Coefficients held
#'     constant by \code{coef_breaks} are repeated in every regime.}
#'   \item{\code{vcov}}{The covariance matrix of \code{coef()}, with rows
#'     and columns in the same order and with the same names. When
#'     \code{coef_breaks} holds coefficients constant, the matrix is singular,
#'     since the repeated entries of a constant coefficient are perfectly
#'     correlated.}
#'   \item{\code{confint}}{For \code{type = "breaks"}, a matrix with one row
#'     per break: the lower bound, the estimated break date and the upper
#'     bound. For \code{type = "coef"}, the lower and upper
#'     bounds of confidence intervals, estimate \eqn{\pm} normal quantile times
#'     standard error, treating the break dates as known.}
#'   \item{\code{residuals}, \code{fitted}}{Matrices with one column per
#'     equation and one row per observation; the row names are the rows
#'     of the data (\code{spec$y_rows}).}
#'   \item{\code{nobs}}{The number of observations used.}
#' }
#' When no break is selected the regime estimates are not computed, and all
#' but \code{nobs} stop with an error.
#'
#' @seealso \code{\link{summary.sysbreak}}, \code{\link{plot.sysbreak}}
#'
#' @examples
#' set.seed(11)
#' bigt <- 40
#' x <- rnorm(bigt)
#' y1 <- 1 + 0.4 * x + rnorm(bigt, sd = 0.3)
#' y2 <- -0.5 + 0.2 * x + rnorm(bigt, sd = 0.3)
#' y1[21:40] <- y1[21:40] + 1.2 - 0.6 * x[21:40]
#' y2[21:40] <- y2[21:40] - 0.8 + 0.4 * x[21:40]
#' dat <- data.frame(y1 = y1, y2 = y2, x = x)
#' fit <- sysbreak(cbind(y1, y2) ~ x, dat, m = 1, trim = 0.3)
#'
#' coef(fit)
#' confint(fit)
#' confint(fit, level = 0.90)
#' confint(fit, type = "coef")
#' head(residuals(fit))
#' nobs(fit)
#'
#' @importFrom stats coef confint fitted nobs residuals vcov
#' @name sysbreak-methods
NULL

#' @rdname sysbreak-methods
#' @export
coef.sysbreak <- function(object, ...) {
  est <- estget(object)
  nm <- cfnames(object)
  stats::setNames(as.numeric(est$coefficients), cflab(nm))
}

#' @rdname sysbreak-methods
#' @export
vcov.sysbreak <- function(object, ...) {
  est <- estget(object)
  if (is.null(est$vcov)) {
    stop("The coefficient covariance matrix is not available for this fit.",
         call. = FALSE)
  }
  nm <- cfnames(object)
  lab <- cflab(nm)
  out <- est$vcov
  dimnames(out) <- list(lab, lab)
  out
}

#' @rdname sysbreak-methods
#' @export
confint.sysbreak <- function(object, parm, level = NULL,
                             type = c("breaks", "coef"), ...) {
  type <- match.arg(type)
  est <- estget(object)
  if (is.null(level)) {
    level <- if (is.null(object$spec$ci_level)) 0.95 else object$spec$ci_level
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  }
  probs <- c(1 - level, 1 + level) / 2
  pct <- paste(format(100 * probs, trim = TRUE, scientific = FALSE,
                      digits = 3), "%")

  if (type == "breaks") {
    if (isTRUE(all.equal(level, 0.95))) {
      cols <- 1:2
    } else if (isTRUE(all.equal(level, 0.90))) {
      cols <- 3:4
    } else {
      stop("For break dates `level` must be 0.95 or 0.90.", call. = FALSE)
    }
    ci <- attr(est, "ci_all")
    if (is.null(ci)) {
      # fits built without both levels carry only the requested one
      if (!isTRUE(all.equal(level, object$spec$ci_level))) {
        stop("Intervals at this level are not stored; refit with `ci_level` set to it.",
             call. = FALSE)
      }
      ci <- est$ci_breaks
      cols <- 1:2
    }
    out <- cbind(ci[, cols[1]], est$breaks, ci[, cols[2]])
    dimnames(out) <- list(paste0("break", seq_along(est$breaks)),
                          c(pct[1], "break", pct[2]))
  } else {
    v <- stats::vcov(object)
    cf <- stats::coef(object)
    se <- sqdiag(v)
    z <- stats::qnorm(probs[2])
    out <- cbind(cf - z * se, cf + z * se)
    dimnames(out) <- list(rownames(v), pct)
  }
  if (!missing(parm)) {
    ok <- if (is.character(parm)) {
      all(parm %in% rownames(out))
    } else {
      is.numeric(parm) && !anyNA(parm) && all(parm >= 1 & parm <= nrow(out))
    }
    if (!ok) {
      stop(sprintf("`parm` must select rows of the result, by index (1 to %d) or by name, such as \"%s\".",
                   nrow(out), rownames(out)[1]), call. = FALSE)
    }
    out <- out[parm, , drop = FALSE]
  }
  out
}

#' @rdname sysbreak-methods
#' @export
residuals.sysbreak <- function(object, ...) {
  est <- estget(object)
  out <- est$residuals
  rownames(out) <- yrows(object)
  out
}

#' @rdname sysbreak-methods
#' @export
fitted.sysbreak <- function(object, ...) {
  if (is.null(object$y)) {
    stop("The fit object does not store the data; refit with sysbreak(), varbreak() or netbreak().",
         call. = FALSE)
  }
  out <- object$y - stats::residuals(object)
  rownames(out) <- yrows(object)
  out
}

#' @rdname sysbreak-methods
#' @export
nobs.sysbreak <- function(object, ...) {
  as.integer(object$spec$bigt)
}

#' Summarize a Structural-Break Fit
#'
#' A report of a fit returned by \code{\link{sysbreak}},
#' \code{\link{varbreak}} or \code{\link{netbreak}}: the break number,
#' the break dates with their confidence intervals, and a coefficient table for
#' each regime. \code{print()} of the fit gives the full output.
#'
#' @param object A fit returned by \code{sysbreak()}, \code{varbreak()} or
#'   \code{netbreak()}.
#' @param x An object returned by \code{summary()}.
#' @param digits Number of significant digits in the coefficient tables.
#' @param ... Not used.
#'
#' @return An object of class \code{"summary.sysbreak"}, a list with the
#'   \code{model} type, the numbers of equations \code{n} and observations
#'   \code{bigt}, the \code{tests} for the break-number decision, the number
#'   of breaks \code{selected}, and, when at least one break is selected, the
#'   break dates with their confidence intervals (\code{breaks}), the rows of
#'   the data in each regime (\code{regimes}) and the coefficient tables
#'   (\code{coefficients}): estimate, standard error, \eqn{z} value and
#'   two-sided \eqn{p}-value, one table per regime, treating the break dates as
#'   known.
#'
#' @seealso \code{\link{sysbreak-methods}}, \code{\link{plot.sysbreak}}
#'
#' @examples
#' set.seed(11)
#' bigt <- 40
#' x <- rnorm(bigt)
#' y1 <- 1 + 0.4 * x + rnorm(bigt, sd = 0.3)
#' y2 <- -0.5 + 0.2 * x + rnorm(bigt, sd = 0.3)
#' y1[21:40] <- y1[21:40] + 1.2 - 0.6 * x[21:40]
#' y2[21:40] <- y2[21:40] - 0.8 + 0.4 * x[21:40]
#' dat <- data.frame(y1 = y1, y2 = y2, x = x)
#' fit <- sysbreak(cbind(y1, y2) ~ x, dat, m = 1, trim = 0.3)
#' summary(fit)
#'
#' @export
summary.sysbreak <- function(object, ...) {
  sp <- object$spec
  bs <- object$br_selected
  seqs <- bs$sequential
  num <- function(v) if (length(v) == 1L) as.numeric(v) else NA_real_
  tests <- data.frame(
    statistic = c(num(bs$wdmax), vapply(seqs, function(s) num(s$statistic), 0)),
    critical_value = c(num(bs$wdmax_cv), vapply(seqs, function(s) num(s$cv), 0)),
    significant = c(isTRUE(bs$wdmax_sig),
                    vapply(seqs, function(s) isTRUE(s$sig), TRUE)),
    row.names = c("WDmax", vapply(seqs, function(s) {
      sprintf("SEQ(%d|%d)", s$l + 1L, s$l)
    }, ""))
  )
  out <- list(model = class(object)[1], n = sp$n, bigt = sp$bigt, m = sp$m,
              trim = sp$trim, brv = sp$brv, robust = sp$robust, level = bs$level, tests = tests,
              selected = bs$breaks, ci_level = sp$ci_level,
              changing = if (isTRUE(sp$restricted)) cbchg(sp) else NULL)

  if (!is.null(object$estimation)) {
    lev <- if (is.null(sp$ci_level)) 0.95 else sp$ci_level
    out$breaks <- stats::confint(object, level = lev)
    out$regimes <- regspan(object)
    nm <- cfnames(object)
    cf <- byregime(object$estimation$coefficients, nm$nreg - 1L, nm$q)
    dimnames(cf) <- list(nm$regimes, nm$coefs)
    v <- object$estimation$vcov
    se <- if (is.null(v)) {
      matrix(NA_real_, nrow(cf), ncol(cf))
    } else {
      byregime(sqdiag(v), nrow(cf) - 1L, ncol(cf))
    }
    out$coefficients <- lapply(seq_len(nrow(cf)), function(r) {
      z <- cf[r, ] / se[r, ]
      cbind(Estimate = cf[r, ], `Std. Error` = se[r, ], `z value` = z,
            `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
    })
    names(out$coefficients) <- rownames(cf)
  }
  class(out) <- "summary.sysbreak"
  out
}

#' @rdname summary.sysbreak
#' @export
print.summary.sysbreak <- function(x, digits = max(3L, getOption("digits") - 3L),
                                   ...) {
  cat(sprintf("%s: %d equation%s, T = %d, m = %d, trim = %g, robust = %s\n",
              x$model, as.integer(x$n), if (isTRUE(x$n > 1)) "s" else "",
              as.integer(x$bigt), as.integer(x$m), x$trim, x$robust))
  cat(sprintf("Covariance matrix: %s across regimes\n\n",
              if (isTRUE(x$brv)) "may change" else "constant"))
  cat(sprintf("Tests for the number of breaks (%g%% level)\n", 100 * x$level))
  tt <- x$tests
  tab <- data.frame(formatC(tt$statistic, format = "f", digits = 3),
                    formatC(tt$critical_value, format = "f", digits = 3),
                    ifelse(tt$significant, "*", ""),
                    row.names = paste0("  ", rownames(tt)))
  names(tab) <- c("Statistic", "Critical value", "")
  print(tab)
  if (length(x$selected) != 1L || is.na(x$selected)) {
    cat("Breaks selected: undetermined\n")
  } else {
    cat(sprintf("Breaks selected: %d\n", as.integer(x$selected)))
  }
  if (is.null(x$breaks)) {
    cat("No break was selected, so there are no regime estimates.\n")
    return(invisible(x))
  }

  cat(sprintf("\nBreak dates (rows of the data) with %g%% confidence intervals\n",
              100 * x$ci_level))
  br <- x$breaks
  tab <- data.frame(br[, 2], br[, 1], br[, 3],
                    row.names = paste0("  Break ", seq_len(nrow(br))))
  names(tab) <- c("Estimate", "Lower", "Upper")
  print(tab)
  if (length(x$changing)) {
    cat(sprintf("\nCoefficients allowed to change: %s\n",
                paste(x$changing, collapse = ", ")))
  }
  for (r in seq_along(x$coefficients)) {
    cat(sprintf("\nRegime %d (rows %d to %d)\n", r, x$regimes[r, 1],
                x$regimes[r, 2]))
    stats::printCoefmat(x$coefficients[[r]], digits = digits,
                        signif.legend = r == length(x$coefficients))
  }
  invisible(x)
}

#' Plot a Structural-Break Fit
#'
#' Plots each series with the estimated break dates (dashed lines), their
#' confidence intervals (shaded) and, optionally, the fitted values of the
#' regime-specific model.
#'
#' @param x A fit returned by \code{\link{sysbreak}}, \code{\link{varbreak}}
#'   or \code{\link{netbreak}}.
#' @param which Series to plot, as indices or names. By default all series,
#'   or the first four when there are more.
#' @param dates Optional vector with one entry per row of the data passed to
#'   the fitting function, including the initial lags of a VAR and any rows
#'   dropped for missing values; the observation in row \code{r} is placed at
#'   \code{dates[r]}. Numeric values or dates are used as the horizontal axis;
#'   character values, such as \code{"1983:Q1"}, are used as labels for the
#'   axis.
#' @param fitted Logical. If \code{TRUE}, add the fitted values.
#' @param ... Graphical parameters. \code{type}, \code{col}, \code{lwd},
#'   \code{lty} and \code{pch} apply to the plotted series; \code{main} and
#'   \code{ylab} may be vectors with one entry per panel; the others, such as
#'   \code{xlab}, \code{ylim} or \code{las}, are passed to \code{plot()} when
#'   each panel is set up.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{sysbreak-methods}}, \code{\link{summary.sysbreak}}
#'
#' @examples
#' set.seed(11)
#' bigt <- 40
#' x <- rnorm(bigt)
#' y1 <- 1 + 0.4 * x + rnorm(bigt, sd = 0.3)
#' y2 <- -0.5 + 0.2 * x + rnorm(bigt, sd = 0.3)
#' y1[21:40] <- y1[21:40] + 1.2 - 0.6 * x[21:40]
#' y2[21:40] <- y2[21:40] - 0.8 + 0.4 * x[21:40]
#' dat <- data.frame(y1 = y1, y2 = y2, x = x)
#' fit <- sysbreak(cbind(y1, y2) ~ x, dat, m = 1, trim = 0.3)
#' plot(fit)
#'
#' @export
plot.sysbreak <- function(x, which = NULL, dates = NULL, fitted = TRUE, ...) {
  y <- x$y
  if (is.null(y)) {
    stop("The fit object does not store the data; refit with sysbreak(), varbreak() or netbreak().",
         call. = FALSE)
  }
  k <- ncol(y)
  if (is.null(which)) {
    which <- seq_len(min(k, 4L))
    if (k > 4L) {
      message(sprintf("Plotting the first 4 of %d series; use `which` to choose.", k))
    }
  }
  if (is.character(which)) {
    which <- match(which, colnames(y))
  }
  if (anyNA(which) || any(which < 1 | which > k)) {
    stop("`which` must select columns of the data.", call. = FALSE)
  }

  rows <- yrows(x)
  if (!is.null(dates) && length(dates) < max(rows)) {
    stop(sprintf(
      "`dates` must have one entry per row of the data: the fit uses rows up to %d, but `dates` has %d entries.",
      max(rows), length(dates)), call. = FALSE)
  }
  lab <- NULL
  if (is.null(dates)) {
    tt <- rows
  } else if (is.numeric(dates) || inherits(dates, c("Date", "POSIXt"))) {
    tt <- dates[rows]
  } else {
    tt <- rows
    lab <- as.character(dates)
  }

  est <- x$estimation
  br <- if (is.null(est)) integer() else match(est$breaks, rows)
  ci <- if (is.null(est)) NULL else {
    matrix(match(est$ci_breaks, rows), ncol = 2)
  }
  fv <- if (isTRUE(fitted) && !is.null(est)) stats::fitted(x) else NULL

  ## widen the bottom margin when an xlab or sub is supplied
  dots <- list(...)
  given <- function(a) !is.null(a) && any(nzchar(a))
  mar1 <- if (given(dots[["sub"]])) 5 else if (given(dots[["xlab"]])) 4 else 3
  op <- graphics::par(mfrow = c(length(which), 1L), mar = c(mar1, 4, 2, 1))
  on.exit(graphics::par(op))
  ## user values in ... replace these defaults; type, col, lwd, lty and pch
  ## apply to the series, and main and ylab may give one entry per panel
  panel <- function(j, i, ..., type = "l", xlab = "", ylab = colnames(y)[j],
                    main = colnames(y)[j],
                    xaxt = if (is.null(lab)) "s" else "n",
                    col = graphics::par("col"), lwd = graphics::par("lwd"),
                    lty = graphics::par("lty"), pch = graphics::par("pch")) {
    if (length(main) > 1L && length(main) == length(which)) main <- main[i]
    if (length(ylab) > 1L && length(ylab) == length(which)) ylab <- ylab[i]
    graphics::plot(tt, y[, j], type = "n", xlab = xlab, ylab = ylab,
                   main = main, xaxt = xaxt, ...)
    if (!is.null(lab)) {
      at <- pretty(rows)
      at <- at[at >= min(rows) & at <= max(rows)]
      graphics::axis(1, at = at, labels = lab[at])
    }
    if (length(br)) {
      usr <- graphics::par("usr")
      graphics::rect(tt[ci[, 1]], usr[3], tt[ci[, 2]], usr[4], border = NA,
                     col = grDevices::adjustcolor("red", alpha.f = 0.15))
      graphics::abline(v = tt[br], col = "red", lty = 2)
    }
    graphics::lines(tt, y[, j], type = type, col = col, lwd = lwd, lty = lty,
                    pch = pch)
    if (!is.null(fv)) {
      graphics::lines(tt, fv[, j], col = "steelblue")
    }
  }
  for (i in seq_along(which)) {
    panel(which[i], i, ...)
  }
  invisible(x)
}

#' Regime Estimates of a Fit, or an Error When There Are None
#'
#' @noRd
estget <- function(object) {
  est <- object$estimation
  if (is.null(est)) {
    stop("No break was selected, so there are no regime estimates to report; see the `testing` component.",
         call. = FALSE)
  }
  est
}

#' Regime and Coefficient Labels of a Fit
#'
#' @noRd
cfnames <- function(object) {
  est <- object$estimation
  nreg <- length(est$breaks) + 1L
  q <- length(est$coefficients) / nreg
  cn <- object$spec$coef_names
  if (is.null(cn) || length(cn) != q) {
    cn <- paste0("coef", seq_len(q))
  }
  list(regimes = paste0("regime", seq_len(nreg)), coefs = cn, q = q,
       nreg = nreg)
}

#' Labels regime:coefficient, in the order the coefficients are stacked
#'
#' @noRd
cflab <- function(nm) {
  paste(rep(nm$regimes, each = nm$q), nm$coefs, sep = ":")
}

#' Rows of the Data Used in Estimation
#'
#' @noRd
yrows <- function(object) {
  rows <- object$spec$y_rows
  if (is.null(rows)) seq_len(object$spec$bigt) else rows
}

#' First and Last Row of the Data in Each Regime
#'
#' @noRd
regspan <- function(object) {
  rows <- yrows(object)
  br <- object$estimation$breaks
  pos <- match(br, rows)
  out <- cbind(start = rows[c(1L, pos + 1L)], end = c(br, rows[length(rows)]))
  rownames(out) <- paste0("regime", seq_len(nrow(out)))
  out
}
