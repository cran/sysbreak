# Main procedure that implements the estimation and testing of breaks for a
# linear multiple-equation system.

#' Main Procedure for Breaks in a Linear Multiple-Equation System
#'
#' Estimates multivariate structural break models with a maximum number of
#' breaks. The procedure computes break estimates, likelihood-ratio or
#' robust Wald test statistics, sequential tests for determining the number of breaks, 
#' and confidence intervals for break dates and coefficients.
#'
#' @usage
#' mainp(y, z, m, S, R = NULL, R_path = NULL,
#'   trim = 0.15, verbose = FALSE,
#'   brbeta = TRUE, brv = FALSE, robust = c("none", "hc", "hac"),
#'   hetq = FALSE, prewhit = FALSE,
#'   sumq = NULL,
#'   restart = c("none", "brbp", "brbp_residual"),
#'   restart_B = 10, restart_init = NULL,
#'   restart_verbose = FALSE)
#'
#' @param y Numeric matrix, vector, or data frame of dependent variables with
#'   \eqn{T} rows and \eqn{n} columns.
#' @param z Numeric matrix, vector, or data frame of regressors with \eqn{T}
#'   rows. 
#' @param m Integer, maximum number of breaks to estimate.
#' @param S Selection matrix, mapping regressors into a structure convenient for
#'   estimation. For observation \eqn{t}, the transformation is
#'   \eqn{(I_n \otimes z_t')S}.
#' @param R Optional restriction matrix for the regime coefficients to specify
#'   which coefficients are allowed to change. The restricted
#'   coefficient vector is estimated under \eqn{\beta = Rb}.
#' @param R_path Optional list of restriction matrices, one for each of 1 to
#'   \eqn{m} breaks. Built automatically from
#'   \code{coef_breaks}.
#' @param trim Trimming parameter. The minimum segment length is
#'   \code{h = round(trim * nrow(y))}. We recommend \eqn{h \ge 3(n + k)}, where
#'   \eqn{n} is the number of equations and \eqn{k} the number of regressors per
#'   equation, and a larger value for larger systems. A warning is issued if
#'   \eqn{h < 3(n + k)}.
#' @param verbose Logical. If \code{TRUE}, print an estimation summary.
#' @param brbeta Logical. If \code{TRUE}, allow regression
#'   coefficients to change across regimes.
#' @param brv Logical. If \code{TRUE}, allow the disturbance
#'   covariance matrix to change across regimes.
#' @param robust Character, one of \code{"none"}, \code{"hc"} or \code{"hac"}; see
#'   \code{sysbreak()}.
#' @param hetq Logical. If \code{TRUE}, the regressor distribution is allowed
#'   to change across regimes; if \code{FALSE}, pooled
#'   regressor moments are used in the break-date confidence intervals.
#'   This affects only the construction of confidence intervals, not testing.
#' @param prewhit Logical. If \code{TRUE}, prewhiten the scores by a VAR(1)
#'   before the HAC estimation. Used only when \code{robust = "hac"}.
#' @param sumq Optional integer giving the number of changing
#'   coefficient restrictions used for critical values. If \code{NULL}, all
#'   transformed regressors are treated as changing.
#' @param restart Character string, one of \code{"none"},
#'   \code{"brbp"}, or \code{"brbp_residual"}, controlling whether bootstrap
#'   restarting is used.
#' @param restart_B Integer number of bootstrap restart iterations when
#'   \code{restart} is not \code{"none"}.
#' @param restart_init Optional integer vector of starting break dates
#'   for bootstrap restarting. If \code{NULL}, the usual estimate is
#'   used as the first restart point.
#' @param restart_verbose Logical. If \code{TRUE}, print bootstrap restart
#'   progress.
#'
#' @return An object of class \code{"sysbreak"} with the following
#'   components:
#' \describe{
#'   \item{\code{spec}}{A list containing \code{m}, \code{trim}, the sample size
#'     \code{bigt}, the number of equations \code{n}, the number of coefficients
#'     per regime \code{q}, \code{brbeta}, \code{brv}, \code{robust}, \code{hetq},
#'     \code{prewhit}, and whether the fitted coefficient specification is
#'     restricted (i.e., partial structural change).}
#'   \item{\code{testing}}{A list containing:
#'     \code{statistic}, a numeric vector of length \code{m} for
#'     \eqn{1,\ldots,m} breaks;
#'     \code{statistic_type}, either \code{"SupLR"}, \code{"SupF"}, or missing
#'     when the requested robust test is unavailable;
#'     \code{critical_values}, an \eqn{m \times 4} matrix for the 90, 95,
#'     97.5, and 99 percent levels;
#'     \code{wdmax}, a four-element vector because its weights depend on the
#'     90, 95, 97.5, and 99 percent critical-value levels, together with its
#'     four-element critical-value vector;
#'     \code{sequential}, a list of the \eqn{l+1\mid l} tests for
#'     \eqn{l=1,\ldots,m-1}; and
#'     \code{changing_parameters}, the number of parameters added by one
#'     additional break. Each sequential-test entry contains integer
#'     \code{l}, scalar \code{statistic}, and a four-element
#'     \code{critical_values} vector ordered as 90, 95, 97.5, and 99 percent.}
#'   \item{\code{estimation}}{A list containing:
#'     \code{breaks}, the \code{m} final break dates;
#'     \code{coefficients}, the regime-stacked coefficient vector;
#'     \code{vcov}, its variance-covariance matrix;
#'     \code{covariance}, the vertically stacked regime covariance matrices;
#'     and \code{ci_breaks}, an \eqn{m \times 4} matrix whose columns are the
#'     95 percent lower and upper bounds followed by the 90 percent lower and
#'     upper bounds.}
#' }
#' The coefficient vector has length \eqn{(m+1)q} and ordering
#' \eqn{(\hat\beta_1',\ldots,\hat\beta_{m+1}')'}. Coefficients constrained to
#' remain constant are repeated across regime blocks; \code{vcov} uses the
#' same ordering. The covariance component has \eqn{(m+1)n} rows and
#' \eqn{n} columns. When covariance breaks are disabled, the common covariance
#' estimate is repeated in each regime block.
#'
#' @details
#' The selected branch estimates the optimal break dates for up to
#' \eqn{1,\dots,m} breaks, computes test statistics and critical values, and
#' then returns the estimates for up to \eqn{m} breaks. When \code{R} is
#' supplied, only the restricted branch is run. It iterates between restricted
#' coefficient estimation and break-date updating. Moreover, with
#' \code{restart = "brbp"}, the estimation applies the resampling bootstrap to
#' get restarting values. With \code{restart = "brbp_residual"}, it applies
#' residual-bootstrap to get starting values.
#'
#' When \code{robust = "none"}, the tests are likelihood-ratio tests and may
#' test coefficient breaks, covariance breaks, or both. When \code{robust} is
#' \code{"hc"} or \code{"hac"}, the reported statistic is the robust Wald
#' (SupF) statistic for coefficient breaks only; robust tests for covariance
#' breaks are not implemented, since their critical values are case specific.
#' With \code{brbeta = TRUE} and \code{brv = TRUE}, the robust statistic
#' therefore tests only coefficient stability.
#'
#' @references
#' Qu, Z., and Perron, P. (2007). Estimating and testing structural changes
#' in multivariate regressions. \emph{Econometrica}, 75(2), 459-502.
#'
#' @examples
#' set.seed(1)
#' bigt <- 60
#' z <- cbind(1, rnorm(bigt))
#' beta <- rbind(c(1, 0.5), c(2, -0.2), c(0.5, 0.8))
#' y <- matrix(NA_real_, nrow = bigt, ncol = 1)
#' y[1:20, ] <- z[1:20, ] %*% beta[1, ] + rnorm(20, sd = 0.4)
#' y[21:40, ] <- z[21:40, ] %*% beta[2, ] + rnorm(20, sd = 0.4)
#' y[41:60, ] <- z[41:60, ] %*% beta[3, ] + rnorm(20, sd = 0.4)
#' fit <- mainp(
#'   y = y,
#'   z = z,
#'   m = 2,
#'   S = diag(ncol(z)),
#'   trim = 0.15
#' )
#' fit$estimation$breaks
#'
#' @noRd
mainp <- function(y, z, m, S, R = NULL, R_path = NULL,
                  trim = 0.15, verbose = FALSE,
                  brbeta = TRUE, brv = FALSE,
                  robust = c("none", "hc", "hac"), hetq = FALSE,
                  prewhit = FALSE, sumq = NULL,
                  restart = c("none", "brbp", "brbp_residual"),
                  restart_B = 10, restart_init = NULL,
                  restart_verbose = FALSE) {
  y <- asmat(y, "y")
  z <- asmat(z, "z")
  S <- asmat(S, "S", allow_vector = FALSE)
  valyz(y, z)  #check rows match
  n <- ncol(y)
  bigt <- nrow(y)
  m <- as.integer(m)
  trm <- as.numeric(trim)
  restart <- match.arg(restart)
  if (is.logical(robust)) {
    stop("`robust` must be \"none\", \"hc\" or \"hac\"; it replaces the logical `hac`.",
         call. = FALSE)
  }
  robust <- match.arg(robust)
  prewhit <- isTRUE(prewhit) && robust == "hac"
  restart_B <- as.integer(restart_B)
  if (length(restart_B) != 1 || is.na(restart_B) || restart_B < 0) {
    stop("`restart_B` must be a non-negative integer.", call. = FALSE)
  }
  if (restart != "none" && is.null(R)) {
    warning("`restart` has no effect unless some coefficients are held constant, so it was ignored.",
            call. = FALSE)
  }

  transformed <- transf(y, z, S)
  maty <- transformed$maty
  matz <- transformed$matz
  q <- ncol(matz)
  if (!is.null(R)) {
    R <- valR(R, m = m, q = q)
  }
  if (!is.null(R_path)) {
    if (is.null(R)) {
      stop("`R_path` cannot be supplied without `R`.", call. = FALSE)
    }
    if (!is.list(R_path) || length(R_path) != m) {
      stop("`R_path` must contain one restriction matrix for each break count from 1 to `m`.",
           call. = FALSE)
    }
    R_path <- lapply(seq_len(m), function(k) valR(R_path[[k]], m = k, q = q))
    if (!isTRUE(all.equal(unname(R_path[[m]]), unname(R),
                          check.attributes = FALSE))) {
      stop("The final matrix in `R_path` must equal `R`.", call. = FALSE)
    }
  } else if (!is.null(R) && m == 1L) {
    R_path <- list(R)
  }
  if (is.null(sumq)) {
    sumq <- q
  }
  sumq <- as.integer(sumq)
  h <- as.integer(round(trm * bigt))
  # dating_MLE() estimates the coefficients and covariance in every regime,
  # therefore h >= n + k is needed (k = q / n regressors per equation) to have
  # enough observations. Short regimes also bias the log-determinant
  # objective; warn if h < 3 (n + k).
  kbar <- q / n
  if (h < n + kbar) {
    stop(sprintf(paste0(
      "`trim` is too small: h = %d is below n + k = %g, and the covariance of ",
      "a regime cannot be estimated. Increase `trim`."), h, n + kbar),
      call. = FALSE)
  }
  if (h < 3 * (n + kbar)) {
    warning(sprintf(paste0(
      "`trim` may be too small: h = %d is less than 3(n + k) = %g (n = %d ",
      "equations, k = %g regressors per equation). The estimated break dates ",
      "tend to produce regimes close to the minimum length; consider a larger ",
      "`trim`."), h, 3 * (n + kbar), n, kbar), call. = FALSE)
  }
  restricted_model <- !is.null(R) || !isTRUE(brbeta)

  if (is.null(R)) {
    selected <- unrm(
      maty = maty,
      matz = matz,
      n = n,
      m = m,
      bigt = bigt,
      trm = trm,
      h = h,
      sumq = sumq,
      brbeta = brbeta,
      brv = brv,
      robust = robust,
      hetq = hetq,
      prewhit = prewhit
    )
  } else {
    selected <- resm(
      maty = maty,
      matz = matz,
      n = n,
      m = m,
      bigt = bigt,
      trm = trm,
      h = h,
      R = R,
      R_path = R_path,
      brbeta = brbeta,
      brv = brv,
      robust = robust,
      hetq = hetq,
      prewhit = prewhit,
      restart = restart,
      restart_B = restart_B,
      restart_init = restart_init,
      restart_verbose = restart_verbose
    )
  }

  selected$estimation$coefficients <- as.numeric(
    selected$estimation$coefficients
  )
  dimnames(selected$estimation$vcov) <- NULL

  result <- c(
    list(spec = list(
      m = m,
      trim = trm,
      bigt = bigt,
      n = n,
      q = q,
      brbeta = brbeta,
      brv = brv,
      robust = robust,
      hetq = hetq,
      prewhit = prewhit,
      restricted = restricted_model
    )),
    selected
  )
  class(result) <- c("sysbreak", class(result))

  if (verbose) {
    print(result)
  }

  result
}

#' Print the results from \code{\link{mainp}}
#'
#' Prints the main test statistics, break-date estimates, confidence
#' intervals, coefficient estimates, and covariance estimates from an object
#' returned by \code{\link{mainp}}.
#'
#' @export
#' @noRd
print.sysbreak <- function(x, ...) {
  pheader(x$spec, class(x)[1])
  restricted <- isTRUE(x$spec$restricted)
  pbranch(
    branch = x,
    title = if (restricted) {
      "Output from the restricted structural change procedures"
    } else {
      "Output from the unrestricted structural change procedures"
    },
    n = x$spec$n,
    restricted = restricted,
    coef_names = x$spec$coef_names,
    changing = cbchg(x$spec),
    br_selected = x$br_selected,
    ci_level = x$spec$ci_level
  )

  invisible(x)
}

#' Print the Testing Result for a Selected Significance Level
#'
#' Used when the break-selection procedure has been run: each
#' statistic is shown against the single critical value for the chosen level,
#' followed by the number of breaks selected.
#'
#' @noRd
ptestsel <- function(br_selected) {
  lab <- sprintf("%g%%", 100 * br_selected$level)
  verdict <- function(sig) {
    if (is.na(sig)) {
      "The test is unavailable"
    } else if (sig) {
      "The test is significant"
    } else {
      "The test is insignificant"
    }
  }
  block <- function(name, stat, cv, sig) {
    cat(sprintf("The %s statistic is: %.3f\n", name, stat))
    cat(sprintf("%s critical value at the %s level: %.3f\n", name, lab, cv))
    cat(verdict(sig), "\n", sep = "")
    cat("\n")
  }

  cat("\n")
  block("WDmax", br_selected$wdmax, br_selected$wdmax_cv,
        br_selected$wdmax_sig)

  if (length(br_selected$sequential)) {
    cat("\nSequential tests:\n")
    for (s in br_selected$sequential) {
      cat(sprintf("The sequential test for %d versus %d break(s) is: %.3f\n",
                  s$l, s$l + 1L, s$statistic))
      cat(sprintf("Critical value at the %s level: %.3f\n", lab, s$cv))
      cat(verdict(s$sig), "\n", sep = "")
      cat("\n")
    }
  }

  if (is.na(br_selected$breaks)) {
    cat("The number of breaks selected is undetermined\n")
  } else {
    cat(sprintf("The number of breaks selected is %d\n", br_selected$breaks))
  }
  cat("\n")

  invisible(NULL)
}

#' Names of the Coefficients that are Allowed to Change
#'
#' Returns the names of the coefficients a restricted specification lets change
#' across regimes, or `NULL` when no selector was recorded (for example when
#' `mainp()` was called directly with a restriction matrix).
#'
#' @noRd
cbchg <- function(spec) {
  info <- spec$coef_breaks
  if (is.null(info)) {
    return(NULL)
  }
  as.character(info$selected)
}

#' Print Main-Procedure Result Header
#'
#' Prints the compact header used by the S3 print method for `mainp()`.
#'
#' @noRd
pheader <- function(spec, model = "sysbreak") {
  cat(sprintf("%s result\n", model))
  cat("------------------------------------------------------------------------\n")
  cat("The basic specifications for testing and estimation:\n")
  cat("-------------------------------------------\n")
  cat(sprintf("(1) M = %d\n", spec$m))
  cat(sprintf("(2) Trimming = %.3f\n", spec$trim))
  cat(sprintf("(3) T = %d\n", spec$bigt))
  cat(sprintf("(4) brbeta = %s, brv = %s\n",
              spec$brbeta, spec$brv))
  cat(sprintf("(5) The number of coefficients in each regime is: %d\n",
              spec$q))
  cat("-------------------------------------------\n")
  cat(sprintf("robust = %s, hetq = %s, prewhit = %s\n",
              spec$robust, spec$hetq, spec$prewhit))
  cat("------------------------------------------------------------------------\n")
  invisible(NULL)
}

#' Prints the testing and estimation output of a `mainp()` result
#'
#' Prints the testing and estimation output of a `mainp()` result.
#'
#' @noRd
pbranch <- function(branch, title, n, restricted,
                    coef_names = NULL, changing = NULL, br_selected = NULL,
                    ci_level = NULL) {
  cat("\n")
  cat(title, "\n")
  cat("***************************************************\n")
  cat("\n")

  ptest(branch$testing, restricted = restricted, br_selected = br_selected)
  pestim(branch$estimation, n = n,
         coef_names = coef_names,
         changing = changing,
         restricted = restricted,
         ci_level = ci_level)

  invisible(NULL)
}

#' Print Main-Procedure Testing Output
#'
#' Prints SupLR/SupF, sequential, and WDmax test statistics.
#'
#' @noRd
ptest <- function(testing, restricted, br_selected = NULL) {
  if (is.null(testing)) {
    return(invisible(NULL))
  }

  cat("Output from the testing procedure\n")
  cat("__________________________________\n")
  cat("\n")

  if (!is.null(br_selected)) {
    ptestsel(br_selected)
    return(invisible(NULL))
  }

  if (!is.null(testing$statistic)) {
    statistic_type <- if (is.null(testing$statistic_type)) {
      "test"
    } else {
      testing$statistic_type
    }
    for (i in seq_along(testing$statistic)) {
      if (!is.na(testing$statistic[i])) {
        label <- if (restricted) {
          paste("restricted", statistic_type)
        } else {
          statistic_type
        }
        cat(sprintf("The %s test for 0 versus %d break(s) is: %.3f\n",
                    label, i, testing$statistic[i]))
      }
    }
  }

  pcv(testing$critical_values)

  if (!is.null(testing$wdmax)) {
    if (length(testing$wdmax) == 1) {
      cat(sprintf("The WDmax statistic is: %.3f\n", testing$wdmax))
    } else {
      pcv(testing$wdmax, label = "WDmax statistics")
    }
    pcv(testing$wdmax_critical_values,
        label = "WDmax critical values")
  }

  if (!is.null(testing$sequential)) {
    cat("Sequential tests:\n")
    for (item in testing$sequential) {
      cat(sprintf("The sequential test for %d versus %d break(s) is: %.3f\n",
                  item$l, item$l + 1, item$statistic))
      pcv(item$critical_values)
    }
  }

  cat("\n")
  invisible(NULL)
}

#' Print Critical Values
#'
#' Prints a vector or matrix of critical values when available.
#'
#' @noRd
pcv <- function(values, label = "Critical values") {
  if (is.null(values)) {
    return(invisible(NULL))
  }

  cat(label, ":\n", sep = "")
  print(round(values, 3))
  cat("\n")

  invisible(NULL)
}

#' Print Main-Procedure Estimation Output
#'
#' Prints estimated break dates, confidence intervals, coefficients, and
#' covariance matrices.
#'
#' @noRd
pestim <- function(estimation, n, coef_names = NULL, changing = NULL,
                   restricted = FALSE, ci_level = NULL) {
  if (is.null(estimation)) {
    return(invisible(NULL))
  }

  cat("Output from the estimation procedure\n")
  cat("__________________________________\n")
  cat("\n")

  m <- length(estimation$breaks)
  q <- length(estimation$coefficients) / (m + 1L)
  coefficients <- byregime(estimation$coefficients, m, q)
  standard_errors <- if (is.null(estimation$vcov)) {
    NULL
  } else {
    byregime(sqdiag(estimation$vcov), m, q)
  }

  pbr(estimation$breaks, estimation$ci_breaks, ci_level = ci_level)
  pcoef(
    coefficients,
    standard_errors,
    coef_names = coef_names,
    changing = changing,
    restricted = restricted,
    ci_level = ci_level
  )
  pcov(estimation$covariance, n = n)

  invisible(NULL)
}

#' Print Break-Date Estimates and Confidence Intervals
#'
#' Prints the estimated break dates with their confidence intervals.  When a
#' `ci_level` is supplied the two-column interval for that coverage rate is
#' printed under each break date; otherwise the 95 and 90 percent intervals
#' are printed side by side, as in rbreak.
#'
#' @noRd
pbr <- function(br, ci, ci_level = NULL) {
  if (is.null(br) || length(br) == 0) {
    return(invisible(NULL))
  }

  cat("--------------------------------------------------------------\n")
  cat("The estimates and the confidence intervals for the break dates\n")
  cat("--------------------------------------------------------------\n")
  one <- !is.null(ci_level) && !is.null(ci) && ncol(ci) >= 2
  for (i in seq_along(br)) {
    cat(sprintf("The estimate of break %d is: %d\n", i, br[i]))
    if (one && nrow(ci) >= i) {
      cat(sprintf("The CI at the %g%% level is: [%d, %d]\n",
                  100 * ci_level, ci[i, 1], ci[i, 2]))
    } else if (!is.null(ci) && nrow(ci) >= i && ncol(ci) >= 4) {
      cat(sprintf("The 95%% C.I. for break %d is: [%d, %d]\n",
                  i, ci[i, 1], ci[i, 2]))
      cat(sprintf("The 90%% C.I. for break %d is: [%d, %d]\n",
                  i, ci[i, 3], ci[i, 4]))
    }
  }
  cat("********************************************************\n")
  cat("\n")

  invisible(NULL)
}

#' Print Coefficients by Regime
#'
#' Prints regime-specific coefficient estimates and, when available, standard
#' errors.  Rows are labelled with the model's coefficient names when the
#' calling function recorded them; otherwise they are `coef1`, `coef2`,
#' and so on.  For a restricted specification an extra column marks which
#' coefficients are free to change across regimes (`break`) and which are held
#' constant (`fixed`).
#'
#' @noRd
pcoef <- function(coefficients, standard_errors = NULL,
                  coef_names = NULL, changing = NULL, restricted = FALSE,
                  ci_level = NULL) {
  if (is.null(coefficients)) {
    return(invisible(NULL))
  }

  q <- ncol(coefficients)
  labels <- if (!is.null(coef_names) && length(coef_names) == q) {
    as.character(coef_names)
  } else {
    paste0("coef", seq_len(q))
  }
  marker <- isTRUE(restricted) && !is.null(changing) &&
    !is.null(coef_names) && length(coef_names) == q
  # Normal confidence interval, estimate +/- z * standard error, treating the
  # break dates as fixed (as the standard errors themselves do).
  show_ci <- !is.null(ci_level) && !is.null(standard_errors)
  zcrit <- if (show_ci) {
    stats::qnorm(1 - (1 - ci_level) / 2)
  } else {
    NA_real_
  }

  cat("The estimated coefficients are:\n")
  for (i in seq_len(nrow(coefficients))) {
    cat(sprintf("For regime %d:\n", i))
    est <- as.vector(coefficients[i, ])
    out <- data.frame(
      estimate = round(est, 5),
      row.names = labels,
      check.names = FALSE
    )
    if (!is.null(standard_errors)) {
      se <- as.vector(standard_errors[i, ])
      out[["std.error"]] <- round(se, 5)
      if (show_ci) {
        out[[sprintf("%g%% lower", 100 * ci_level)]] <- round(est - zcrit * se, 5)
        out[[sprintf("%g%% upper", 100 * ci_level)]] <- round(est + zcrit * se, 5)
      }
    }
    if (marker) {
      out[["  "]] <- ifelse(labels %in% changing, "break", "fixed")
    }
    print(out)
  }
  cat("\n")

  invisible(NULL)
}

#' Print Covariance Matrices by Regime
#'
#' Prints stacked covariance estimates as one block per regime.
#'
#' @noRd
pcov <- function(covariance, n) {
  if (is.null(covariance) || is.na(n) || n <= 0) {
    return(invisible(NULL))
  }

  nreg <- floor(nrow(covariance) / n)
  if (nreg == 0) {
    return(invisible(NULL))
  }

  cat("The estimated covariance matrices of the errors are:\n")
  for (i in seq_len(nreg)) {
    rows <- ((i - 1) * n + 1):(i * n)
    cat(sprintf("For regime %d:\n", i))
    out <- format(covariance[rows, , drop = FALSE], digits = 4)
    # R left-aligns the default labels of a character matrix
    dimnames(out) <- list(sprintf("[%d,]", seq_len(n)), sprintf("[,%d]", seq_len(n)))
    print(out, quote = FALSE, right = TRUE)
  }
  cat("\n")

  invisible(NULL)
}

#' Unrestricted Main-Procedure Function
#'
#' Runs unrestricted estimation, testing, critical-value calculations,
#' sequential tests, and confidence intervals for `mainp()`.
#'
#' @noRd
unrm <- function(maty, matz, n, m, bigt, trm, h, sumq,
                 brbeta, brv, robust,
                 hetq, prewhit) {
  q <- ncol(matz)
  q_test <- qtest(n, sumq, brbeta, brv, robust = robust)
  lrtest <- matrix(NA_real_, nrow = m, ncol = 1)
  breaks_by_m <- vector("list", m)
  for (i in seq_len(m)) {
    R_i <- defR(
      i, q,
      brbeta = brbeta,
      brv = brv
    )
    fit_i <- est(
      maty, matz, n, i, bigt, trm, R_i,
      brbeta = brbeta,
      brv = brv
    )
    final_i <- r_estim(
      maty, matz, bigt, n, i, fit_i$dx, R_i,
      brbeta = brbeta,
      brv = brv
    )
    breaks_by_m[[i]] <- fit_i$dx
    lrtest[i, 1] <- tstat(
      maty = maty,
      matz = matz,
      m = i,
      n = n,
      br = fit_i$dx,
      beta = final_i$nbeta,
      bigt = bigt,
      brbeta = brbeta,
      brv = brv,
      robust = robust,
      prewhit = prewhit,
      R = R_i
    )
  }

  cvs <- cvlrmat(m, q_test, trm)
  colnames(cvs) <- c("90%", "95%", "97.5%", "99%")
  rownames(cvs) <- paste0("m=", seq_len(m))

  sequential <- NULL
  if (m > 1) {
    sequential <- vector("list", m - 1)
    for (i in seq_len(m - 1)) {
      stat_i <- if (is.na(q_test)) {
        NA_real_
      } else if (robust == "none") {
        seqlr(
          maty, matz, n, breaks_by_m[[i]], bigt, h,
          brbeta = brbeta,
          brv = brv,
          trm = trm
        )
      } else {
        seqf(
          maty, matz, n, breaks_by_m[[i]], bigt, h,
          prewhit = prewhit,
          robust = robust,
          trm = trm
        )
      }
      sequential[[i]] <- list(
        l = i,
        statistic = stat_i,
        critical_values = seq_na(i, q_test, trm)
      )
    }
  }

  R_m <- defR(
    m, q,
    brbeta = brbeta,
    brv = brv
  )
  final <- r_estim(
    maty, matz, bigt, n, m, breaks_by_m[[m]], R_m,
    brbeta = brbeta,
    brv = brv
  )
  coef_se <- coefse(
    maty = maty,
    matz = matz,
    n = n,
    m = m,
    bigt = bigt,
    br = breaks_by_m[[m]],
    beta = final$nbeta,
    vv = final$nvv,
    R = R_m,
    robust = robust,
    brv = brv,
    prewhit = prewhit
  )
  ci <- interval2(
    maty, matz, m, n, bigt, breaks_by_m[[m]], final$nbeta, final$nvv,
    hetq = hetq,
    robust = robust,
    brv = brv,
    brbeta = brbeta,
    prewhit = prewhit
  )

  list(
    testing = list(
      statistic = as.numeric(lrtest[, 1]),
      statistic_type = if (robust != "none" && !brbeta) {
        NA_character_
      } else if (robust != "none") {
        "SupF"
      } else {
        "SupLR"
      },
      critical_values = cvs,
      wdmax = wd_na(lrtest[, 1], q_test, trm),
      wdmax_critical_values = cvwd_na(q_test, trm, m),
      sequential = sequential,
      changing_parameters = q_test
    ),
    # The optimal partition for every candidate break number.  Element k holds
    # the k-break dates, so a caller that selects fewer than `m` breaks can
    # re-estimate at fixed dates without repeating the dynamic program.
    breaks_by_m = breaks_by_m,
    estimation = list(
      breaks = breaks_by_m[[m]],
      coefficients = final$nbeta,
      vcov = coef_se$vcov,
      covariance = final$nvv,
      ci_breaks = ci
    )
  )
}

#' Restricted Main-Procedure Function
#'
#' Runs the restricted estimation path, double-maximum and sequential tests,
#' and confidence intervals for `mainp()` when restriction matrices are
#' supplied.
#'
#' @noRd
resm <- function(maty, matz, n, m, bigt, trm, h, R,
                 R_path = NULL, brbeta, brv, robust,
                 hetq, prewhit,
                 restart = "none", restart_B = 20,
                 restart_init = NULL,
                 restart_verbose = FALSE) {
  q <- ncol(matz)
  R <- valR(R, m = m, q = q)
  complete_path <- !is.null(R_path)
  if (is.null(R_path)) {
    R_path <- vector("list", m)
    R_path[[m]] <- R
  } else {
    if (!is.list(R_path) || length(R_path) != m) {
      stop("`R_path` must contain one restriction matrix for each break number from 1 to `m`.",
           call. = FALSE)
    }
    R_path <- lapply(seq_len(m), function(k) valR(R_path[[k]], m = k, q = q))
  }

  cq <- qres(
    n = n, m = m, R = R,
    brbeta = brbeta,
    brv = brv,
    robust = robust
  )
  lrtest <- rep(NA_real_, m)
  breaks_by_m <- vector("list", m)
  fit_m <- NULL
  final <- NULL

  for (i in seq_len(m)) {
    R_i <- R_path[[i]]
    if (is.null(R_i)) {
      next
    }
    cq_i <- qres(
      n = n, m = i, R = R_i,
      brbeta = brbeta,
      brv = brv,
      robust = robust
    )
    if (!(is.na(cq_i) && is.na(cq)) && !identical(cq_i, cq)) {
      stop(
        "The number of changing parameters must be constant across `R_path`.",
        call. = FALSE
      )
    }

    if (restart %in% c("brbp", "brbp_residual")) {
      r_res <- brbp(
        maty = maty,
        matz = matz,
        n = n,
        m = i,
        bigt = bigt,
        trm = trm,
        R = R_i,
        method = restart,
        initial_breaks = if (i == m) restart_init else NULL,
        B = restart_B,
        brbeta = brbeta,
        brv = brv,
        verbose = restart_verbose
      )
      fit_i <- r_res$fit
    } else {
      fit_i <- est(
        maty, matz, n, i, bigt, trm, R_i,
        brbeta = brbeta,
        brv = brv
      )
    }

    final_i <- r_estim(
      maty, matz, bigt, n, i, fit_i$dx, R_i,
      brbeta = brbeta,
      brv = brv
    )
    breaks_by_m[[i]] <- fit_i$dx
    lrtest[[i]] <- tstat(
      maty = maty,
      matz = matz,
      m = i,
      n = n,
      br = fit_i$dx,
      beta = final_i$nbeta,
      bigt = bigt,
      brbeta = brbeta,
      brv = brv,
      robust = robust,
      prewhit = prewhit,
      R = R_i
    )

    if (i == m) {
      fit_m <- fit_i
      final <- final_i
    }
  }

  if (is.null(fit_m) || is.null(final)) {
    stop(
      "The restricted model for `m` breaks was not estimated.",
      call. = FALSE
    )
  }
  if (!complete_path && m > 1L) {
    warning(
      "`R` works only for `m` breaks. To test 1 to `m` breaks, supply `R_path` with one `R` for each.",
      call. = FALSE
    )
  }

  cvs <- cvlrmat(m, cq, trm)
  if (!complete_path && m > 1L) {
    cvs[seq_len(m - 1L), ] <- NA_real_
  }
  colnames(cvs) <- c("90%", "95%", "97.5%", "99%")
  rownames(cvs) <- paste0("m=", seq_len(m))

  sequential <- NULL
  if (complete_path && m > 1L) {
    sequential <- vector("list", m - 1L)
    local_R <- R_path[[1L]]
    for (i in seq_len(m - 1L)) {
      stat_i <- if (is.na(cq)) {
        NA_real_
      } else if (robust == "none") {
        seqlr(
          maty, matz, n, breaks_by_m[[i]], bigt, h,
          brbeta = brbeta,
          brv = brv,
          trm = trm,
          R = local_R
        )
      } else {
        seqf(
          maty, matz, n, breaks_by_m[[i]], bigt, h,
          prewhit = prewhit,
          robust = robust,
          trm = trm,
          R = local_R
        )
      }
      sequential[[i]] <- list(
        l = i,
        statistic = stat_i,
        critical_values = seq_na(i, cq, trm)
      )
    }
  }

  coef_se <- coefse(
    maty = maty,
    matz = matz,
    n = n,
    m = m,
    bigt = bigt,
    br = fit_m$dx,
    beta = final$nbeta,
    vv = final$nvv,
    R = R,
    robust = robust,
    brv = brv,
    prewhit = prewhit
  )
  ci <- interval2(maty, matz, m, n, bigt, fit_m$dx, final$nbeta, final$nvv,
                  hetq = hetq,
                  robust = robust,
                  brv = brv,
                  brbeta = brbeta,
                  prewhit = prewhit)

  list(
    testing = list(
      statistic = lrtest,
      statistic_type = if (robust != "none" && !brbeta) {
        NA_character_
      } else if (robust != "none") {
        "SupF"
      } else {
        "SupLR"
      },
      critical_values = cvs,
      wdmax = wd_na(lrtest, cq, trm),
      wdmax_critical_values = if (complete_path) {
        cvwd_na(cq, trm, m)
      } else {
        rep(NA_real_, 4)
      },
      sequential = sequential,
      changing_parameters = cq
    ),
    # See the note in `unrm()`: the optimal partition for every candidate
    # break number.
    breaks_by_m = breaks_by_m,
    estimation = list(
      breaks = fit_m$dx,
      coefficients = final$nbeta,
      vcov = coef_se$vcov,
      covariance = final$nvv,
      ci_breaks = ci
    )
  )
}

#' Procedure for Testing
#'
#' Chooses and computes the appropriate LR or robust Wald test statistic for a
#' fitted `mainp()` break model.
#'
#' @noRd
tstat <- function(maty, matz, m, n, br, beta, bigt,
                  brbeta, brv, robust,
                  prewhit, R = NULL) {
  if (robust == "none") {
    return(suplr(
      maty, matz, m, n, br, beta, bigt,
      brbeta = brbeta,
      brv = brv
    ))
  }

  if (brv && !brbeta) {
    warning(
      "Robust tests for covariance breaks are not available; returning `NA` for the test statistic.",
      call. = FALSE
    )
    return(NA_real_)
  }

  if (brv && brbeta) {
    warning(
      "With `robust` set to \"hc\" or \"hac\" and `brv = TRUE`, the reported robust test is for coefficient breaks only and assumes a constant error covariance.",
      call. = FALSE
    )
  }

  if (is.null(R) || identical(R, diag(nrow(R)))) {
    return(supft(
      maty, matz, m, n, br, beta, bigt,
      prewhit = prewhit, robust = robust
    ))
  }

  suprft(
    maty, matz, m, n, br, beta, R, bigt,
    prewhit = prewhit, robust = robust
  )
}

#' Restriction Matrix
#'
#' Constructs the default reparameterization matrix implied by the requested
#' coefficient-break and covariance-break specification.
#'
#' @noRd
defR <- function(m, q, brbeta, brv) {
  if (brbeta) {
    return(diag((m + 1) * q))
  }
  if (brv) {
    return(kronecker(matrix(1, nrow = m + 1, ncol = 1), diag(q)))
  }
  diag((m + 1) * q)
}

#' Number of Changing Parameters for the Unrestricted Case
#'
#' Computes the number of changing parameters used to select critical
#' values for the unrestricted case.
#'
#' @noRd
qtest <- function(n, sumq, brbeta, brv, robust) {
  if (robust != "none") {
    if (brbeta) {
      return(as.integer(sumq))
    }
    return(NA_integer_)
  }

  if (brv && !brbeta) {
    return(as.integer(n * (n + 1) / 2))
  }
  if (brv && brbeta) {
    return(as.integer(sumq + n * (n + 1) / 2))
  }
  as.integer(sumq)
}

#' Number of Changing Parameters for the Unrestricted Case
#'
#' Computes the implied number of changing parameters for a restricted model
#'
#' @noRd
qres <- function(n, m, R, brbeta, brv, robust) {
  cq <- nrow(R) / (m + 1) - (nrow(R) - ncol(R)) / m
  if (robust != "none") {
    if (!brbeta) {
      return(NA_integer_)
    }
    if (abs(cq - round(cq)) > .Machine$double.eps^0.5) {
      stop(sprintf(
        "`R` must have q + m*k columns, where k is the number of changing parameters at each break. It has %d.",
        ncol(R)
      ))
    }
    return(as.integer(round(cq)))
  }

  if (brv && !brbeta) {
    return(as.integer(n * (n + 1) / 2))
  }

  cq <- n * (n + 1) / 2 * brv + cq
  if (abs(cq - round(cq)) > .Machine$double.eps^0.5) {
    stop(sprintf(
      "`R` must have q + m*k columns, where k is the number of changing parameters at each break. It has %d.",
      ncol(R)
    ))
  }
  as.integer(round(cq))
}

#' Critical-Value Matrix or Missing Values
#'
#' Builds a matrix of SupLR/SupF critical values for break numbers up to
#' `k_max`, or an `NA` matrix if unavailable.
#'
#' @noRd
cvlrmat <- function(k_max, q, trm) {
  if (length(q) != 1 || is.na(q)) {
    return(matrix(NA_real_, nrow = k_max, ncol = 4))
  }
  t(vapply(seq_len(k_max), function(i) cvlr(i, q, trm), numeric(4)))
}

#' Label Critical-Value Levels
#'
#' Adds the probability labels shared by test statistics and critical-value
#' vectors in the testing output.
#'
#' @noRd
cvnames <- function(values) {
  stats::setNames(
    as.numeric(values),
    c("90%", "95%", "97.5%", "99%")
  )
}

#' Sequential Critical Values or Missing Values
#'
#' Returns sequential-test critical values when the parameter number is
#' available, otherwise returns missing values.
#'
#' @noRd
seq_na <- function(l, q, trm) {
  if (length(q) != 1 || is.na(q)) {
    return(cvnames(rep(NA_real_, 4)))
  }
  cvnames(supseq(l, q, trm))
}

#' WDmax or Missing Values
#'
#' Computes WDmax when test statistics and relevant parameters are
#' available.
#'
#' @noRd
wd_na <- function(lrtest, q, trm) {
  if (length(q) != 1 || is.na(q) ||
      length(lrtest) == 0 || anyNA(lrtest) || !all(is.finite(lrtest))) {
    return(cvnames(rep(NA_real_, 4)))
  }
  cvnames(wdmax(lrtest, q, trm))
}

#' WDmax Critical Values or Missing Values
#'
#' Returns WDmax critical values when the parameter count is available.
#'
#' @noRd
cvwd_na <- function(q, trm, M) {
  if (length(q) != 1 || is.na(q)) {
    return(cvnames(rep(NA_real_, 4)))
  }
  cvnames(cvwd(q, trm, M))
}

#' Coefficients by Regime
#'
#' Coefficients by Regime
#'
#' Converts a vector of length (m + 1) q, which stacks the q coefficients of
#' regime 1, then those of regime 2, up to regime m + 1, into an (m + 1) x q
#' matrix whose jth row contains the coefficients of regime j. Also applied to
#' the corresponding standard errors.
#'
#' @noRd
byregime <- function(beta, m, q) {
  matrix(as.vector(beta), nrow = m + 1, ncol = q, byrow = TRUE)
}

#' Validate the Trimming Parameter and the Number of Breaks
#'
#' The minimum segment length is `h = round(trim * T)`, so `m` breaks need `m +
#' 1` segments and therefore `(m + 1) * h` observations. A request that does not
#' fit is rejected: the user needs to choose between a smaller `m` and a smaller
#' `trim`.
#'
#' Without it, `a:b` may lead to b<a, so the program counts backwards in the
#' break search, giving either meaningless statistics or a subscript error.
#'
#' @noRd
mcheck <- function(m, trim, bigt) {
  if (!is.numeric(trim) || length(trim) != 1L || !is.finite(trim) ||
      trim <= 0 || trim >= 0.5) {
    stop("`trim` must be strictly between 0 and 0.5.", call. = FALSE)
  }
  m <- as.integer(m)
  if (length(m) != 1L || is.na(m) || m < 1L) {
    stop("`m` must be a positive integer.", call. = FALSE)
  }

  h <- as.integer(round(trim * bigt))
  if (h < 2L) {
    stop(sprintf(
      "`trim = %g` gives a minimum segment length of %d observation(s) in a sample of %d; at least 2 are needed.",
      trim, h, bigt), call. = FALSE)
  }

  # round() gives
  # h <= trim * bigt + 0.5 < bigt / 2 + 0.5, so 2 * h <= bigt and hence
  # floor(bigt / h) >= 2.
  feasible <- as.integer(floor(bigt / h)) - 1L
  if (m > feasible) {
    stop(sprintf(
      paste0(
        "`m = %d` breaks need %d observations at `trim = %g` (%d segments ",
        "of at least %d), but only %d are available. The largest number ",
        "of breaks this sample allows at this trimming is %d; reduce `m` ",
        "or `trim`."
      ),
      m, (m + 1L) * h, trim, m + 1L, h, bigt, feasible), call. = FALSE)
  }
  m
}

#' Critical-Value Column
#'
#' @noRd
cindex <- function(level) {
  supported <- c(0.10, 0.05, 0.025, 0.01)
  index <- which(abs(level - supported) < 1e-12)
  if (length(index) != 1L) {
    stop(
      "`level` must be one of 0.10, 0.05, 0.025, or 0.01.",
      call. = FALSE
    )
  }
  index
}

#' Select the Number of Breaks from the Testing Output
#'
#' Applies the Bai-Perron sequential procedure to a `testing` component.  The
#' WDmax test decides whether there is any break: if it is not
#' significant at `level`, the procedure stops with zero breaks.  Otherwise it
#' starts from one break and applies the \eqn{l+1 \mid l} sequential tests,
#' stopping at the first one that is not significant or at `m`.
#'
#' WDmax reweights each candidate statistic by
#' `cvlr(1) / cvlr(k)` before taking the maximum, which puts every break count
#' on a common scale; an unweighted maximum has no such correction and cannot
#' be compared with a single critical value.
#'
#' Returns a list with the selected number of `breaks` (`NA_integer_` when the
#' test is unavailable), the `level` used, the WDmax statistic with its
#' critical value and significance, the sequential steps taken, and a `reason`
#' string for the outcome.
#'
#' @noRd
brsel <- function(testing, level, m) {
  idx <- cindex(level)
  m <- as.integer(m)

  gt <- function(stat, cv) {
    is.finite(stat) && is.finite(cv) && stat > cv
  }
  pick <- function(v) {
    if (is.null(v) || length(v) < idx) NA_real_ else as.numeric(v[[idx]])
  }

  wdmax <- pick(testing$wdmax)
  wdmax_cv <- pick(testing$wdmax_critical_values)

  if (!is.finite(wdmax)) {
    return(list(
      breaks = NA_integer_, level = level, m = m,
      wdmax = wdmax, wdmax_cv = wdmax_cv, wdmax_sig = NA,
      sequential = list(),
      reason = "the WDmax test is unavailable for this specification"
    ))
  }

  wd_sig <- gt(wdmax, wdmax_cv)

  if (!wd_sig) {
    return(list(
      breaks = 0L, level = level, m = m,
      wdmax = wdmax, wdmax_cv = wdmax_cv, wdmax_sig = FALSE,
      sequential = list(),
      reason = "WDmax is not significant: no break is detected"
    ))
  }

  selected <- 1L
  steps <- list()
  if (m > 1L && length(testing$sequential)) {
    for (l in seq_len(m - 1L)) {
      item <- testing$sequential[[l]]
      if (is.null(item)) break
      cv <- pick(item$critical_values)
      stat <- as.numeric(item$statistic)
      sig <- gt(stat, cv)
      steps[[length(steps) + 1L]] <- list(
        l = l, statistic = stat, cv = cv, sig = sig
      )
      if (!sig) break
      selected <- l + 1L
    }
  }

  reason <- if (selected == m) {
    sprintf(
      "sequential tests remained significant up to the maximum of %d break(s)",
      m
    )
  } else {
    sprintf(
      "the %d versus %d sequential test is not significant",
      selected, selected + 1L
    )
  }

  list(
    breaks = selected, level = level, m = m,
    wdmax = wdmax, wdmax_cv = wdmax_cv, wdmax_sig = wd_sig,
    sequential = steps,
    reason = reason
  )
}

#' Break-Date Confidence Interval Columns
#'
#' `interval2()` gives the 95 percent bounds in columns 1 and 2 and the 90
#' percent bounds in columns 3 and 4.
#'
#' @noRd
cicols <- function(ci, ci_level) {
  if (is.null(ci)) {
    return(NULL)
  }
  cols <- if (isTRUE(all.equal(ci_level, 0.95))) 1:2 else 3:4
  out <- ci[, cols, drop = FALSE]
  colnames(out) <- sprintf("%g%% %s", 100 * ci_level, c("lower", "upper"))
  out
}

#' Break-Date Confidence Intervals at Both Levels, as Rows of the Data
#'
#' Maps the four `interval2()` columns to rows of the data, clamped to the
#' sample as for `ci_breaks`.  Kept as an attribute of `estimation` so that
#' `confint()` can report either level without refitting.
#'
#' @noRd
cirows <- function(ci, yrow) {
  if (is.null(ci) || ncol(ci) < 4L) {
    return(NULL)
  }
  out <- ci[, 1:4, drop = FALSE]
  out[] <- yrow[pmin(pmax(out, 1L), length(yrow))]
  colnames(out) <- c("95% lower", "95% upper", "90% lower", "90% upper")
  out
}

#' Testing Output at a Single Level
#'
#' The tests are computed at all four levels, ordered 90, 95, 97.5, 99 percent.
#' This code is used to keep `level` and drops the other three. `wdmax` is
#' itself level-specific: the reweighting divides by the critical values at the
#' level used. Run after `brsel()`, which reads all four.
#'
#' @noRd
tcols <- function(testing, level) {
  if (is.null(testing)) {
    return(NULL)
  }
  j <- cindex(level)
  one <- function(v) if (is.null(v) || length(v) < j) v else unname(v[[j]])

  if (!is.null(testing$critical_values)) {
    cv <- testing$critical_values
    testing$critical_values <- stats::setNames(cv[, j], rownames(cv))
  }
  testing$wdmax <- one(testing$wdmax)
  testing$wdmax_critical_values <- one(testing$wdmax_critical_values)
  if (length(testing$sequential)) {
    testing$sequential <- lapply(testing$sequential, function(s) {
      s$critical_values <- one(s$critical_values)
      s
    })
  }
  attr(testing, "level") <- level
  class(testing) <- c("sysbreak_testing", class(testing))
  testing
}

#' Print the Break-Number Decision
#'
#' The statistics and critical values are in `testing`; this prints only the
#' outcome that `brsel()` adds: the WDmax verdict, each sequential step and its
#' outcome, the selected number of breaks, and the reason. The procedure stops
#' at the first insignificant test, so fewer steps are listed than in
#' `testing$sequential`.
#'
#' @export
#' @noRd
print.sysbreak_br_selected <- function(x, ...) {
  lab <- sprintf("%g%%", 100 * x$level)
  verdict <- function(sig) {
    if (length(sig) != 1L || is.na(sig)) {
      "unavailable"
    } else if (sig) {
      "significant"
    } else {
      "not significant"
    }
  }


  line <- function(lhs, v, note = "") {
    cat(
      sub("[[:space:]]+$", "", sprintf("  %-12s %-16s%s", lhs, v, note)),
      "\n", sep = ""
    )
  }

  cat(sprintf("Break selection   (%s level, m = %d)\n\n", lab, as.integer(x$m)))
  line("WDmax", verdict(x$wdmax_sig))

  n_steps <- length(x$sequential)
  for (i in seq_len(n_steps)) {
    s <- x$sequential[[i]]
    last <- i == n_steps && !isTRUE(s$sig)
    line(sprintf("%d vs %d", s$l, s$l + 1L), verdict(s$sig),
         if (last) "  -> stop" else "")
  }

  cat("\n")
  if (length(x$breaks) != 1L || is.na(x$breaks)) {
    cat("  Breaks selected: undetermined\n")
  } else {
    cat(sprintf("  Breaks selected: %d\n", as.integer(x$breaks)))
  }
  if (!is.null(x$reason)) {
    cat(sprintf("  %s\n", x$reason))
  }
  invisible(x)
}

#' Estimation Output That Can Print on Its Own
#'
#' `pestim()` needs the number of equations and the coefficient labels, which
#' are in `spec`. Attaches them as attributes and adds the `sysbreak_estimation`
#' class, so `fit$estimation` prints on its own. Contents are unchanged.
#'
#' @noRd
eclass <- function(estimation, spec) {
  if (is.null(estimation)) {
    return(NULL)
  }
  attr(estimation, "n") <- spec$n
  attr(estimation, "coef_names") <- spec$coef_names
  attr(estimation, "changing") <- cbchg(spec)
  attr(estimation, "restricted") <- isTRUE(spec$restricted)
  attr(estimation, "ci_level") <- spec$ci_level
  class(estimation) <- c("sysbreak_estimation", class(estimation))
  estimation
}

#' Print Estimation Output
#'
#' @export
#' @noRd
print.sysbreak_estimation <- function(x, ...) {
  pestim(
    unclass(x),
    n = attr(x, "n"),
    coef_names = attr(x, "coef_names"),
    changing = attr(x, "changing"),
    restricted = isTRUE(attr(x, "restricted")),
    ci_level = attr(x, "ci_level")
  )
  invisible(x)
}

#' Print Break Tests
#'
#' Prints the single-level output from `tcols()` as a table: the tests of no
#' break against `k` breaks, WDmax, and the sequential tests of `l` against
#' `l + 1` breaks, each beside its critical value. Contents are unchanged.
#'
#' @export
#' @noRd
print.sysbreak_testing <- function(x, ...) {
  lvl <- attr(x, "level")
  lab <- if (is.null(lvl)) "" else sprintf("%g%%", 100 * lvl)
  # `statistic_type` is NA, not NULL, when a requested robust test is
  # unavailable
  type <- if (length(x$statistic_type) != 1L || is.na(x$statistic_type)) {
    "test"
  } else {
    x$statistic_type
  }

  num <- function(v) {
    if (length(v) == 1L && is.finite(v)) sprintf("%12.3f", v) else sprintf("%12s", "NA")
  }
  any_sig <- FALSE
  mark <- function(s, c) {
    if (length(s) == 1L && length(c) == 1L && is.finite(s) &&
        is.finite(c) && s > c) {
      any_sig <<- TRUE
      " *"
    } else {
      ""
    }
  }
  # Each block carries its own column header, so the left column always means
  # "which row within this block".
  colhead <- function() {
    cat(sprintf("%-24s %12s %12s\n", "", "statistic", "crit. value"))
  }
  row <- function(lhs, s, c) {
    cat(sprintf("%-24s %s %s%s\n", sprintf("    %s", lhs),
                num(s), num(c), mark(s, c)))
  }

  cat(sprintf("Structural break tests   (%s level)\n\n", lab))

  if (!is.null(x$statistic)) {
    cv <- if (is.null(x$critical_values)) {
      rep(NA_real_, length(x$statistic))
    } else {
      as.numeric(x$critical_values)
    }
    cat(sprintf(" no break vs k breaks (%s)\n", type))
    colhead()
    for (i in seq_along(x$statistic)) {
      row(sprintf("k = %d", i), x$statistic[i], cv[i])
    }
    cat("\n")
  }

  if (!is.null(x$wdmax)) {
    cat(" double maximum (WDmax)\n")
    colhead()
    lhs <- if (length(x$statistic)) {
      sprintf("M = %d", length(x$statistic))
    } else {
      "WDmax"
    }
    row(lhs, x$wdmax, x$wdmax_critical_values)
    cat("\n")
  }

  if (length(x$sequential)) {
    cat(sprintf(" sequential (l vs l+1, %s)\n", type))
    colhead()
    for (s in x$sequential) {
      row(sprintf("l = %d", s$l), s$statistic, s$critical_values)
    }
    cat("\n")
  }

  if (!is.null(x$changing_parameters)) {
    cat(sprintf("Changing parameters: %d\n", as.integer(x$changing_parameters)))
  }
  if (any_sig) {
    cat(sprintf("* significant at the %s level\n", lab))
  }
  invisible(x)
}
