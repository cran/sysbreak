#' Structural Breaks in a System of Equations
#'
#' Estimation of a multiple-equation system affected by structural breaks in the
#' regression coefficients and/or the variance covariance matrix of the
#' residuals.
#'
#' @usage
#' sysbreak(model, data = NULL, m, coef_breaks = "all",
#'   test_level = 0.05, ci_level = 0.95,
#'   trim = 0.15, verbose = FALSE, brv = FALSE,
#'   robust = c("none", "hc", "hac"), hetq = FALSE, prewhit = FALSE,
#'   restart = c("none", "brbp", "brbp_residual"),
#'   restart_B = 10, restart_verbose = FALSE)
#'
#' @param model A formula such as \code{cbind(y1, y2) ~ x1 + x2} for a system
#'   with common regressors across equations, or
#'   \code{list(y1 ~ x1 + x2, y2 ~ x2 + x3)} for a system allowing different
#'   regressors across equations. The number of regressors may differ across
#'   equations, e.g., \code{list(y1 ~ x1 + x2, y2 ~ x3)} An intercept is added
#'   to each equation unless the formula indicates otherwise.
#' @param data Data frame used to evaluate \code{model}.
#' @param m Integer, the maximum number of breaks allowed when testing for
#'   breaks.
#' @param coef_breaks An expression for which regression coefficients are
#'   allowed to change across regimes. The default \code{"all"} allows all
#'   regression coefficients to change. \code{"none"} or \code{NULL}
#'   restrictions all coefficients to be constant. \ This case requires setting
#'   \code{brv = TRUE}, so that the disturbance covariance matrix\ can change;
#'   otherwise the code will return an error message. For in between, supply a
#'   list of formulas, one per equation, where left-hand side should specify the
#'   equation and the right-hand side selects the coefficients of that equation
#'   that may change:
#'   \itemize{
#'     \item \code{y1 ~ x1, y2 ~ x2}: the coefficient on \code{x1} in the first
#'     equation and the coefficient on \code{x2} in the second equation can
#'     change.
#'     \item \code{y1 ~ 1, y2 ~ x2+ x4}: the intercept in the first equation and
#'     the two slope coefficients in the second equation.
#'     \item \code{y1 ~ 1 + x1, y2~0}: the intercept and \code{x1} in the first
#'     equation, and no coefficient in the second equation.
#'   }
#' @param test_level Significance level used to select the number of breaks.
#'   One of \code{0.10}, \code{0.05}, \code{0.025} or \code{0.01}. The
#'   double-maximum tests decide whether there is any break; if
#'   WDmax is not significant the procedure stops and no break is
#'   reported. Otherwise the sequential tests of \eqn{l} against \eqn{l+1}
#'   breaks are applied from one break upwards, stopping at the first
#'   insignificant test or at \code{m}.
#' @param ci_level Coverage probability of the reported confidence intervals
#'   for the break dates. Either \code{0.95} or \code{0.90}.
#' @param trim Trimming parameter. The minimum segment length is
#'   \code{h = round(trim * nrow(y))}. We recommend \eqn{h \ge 3(n + k)}, where
#'   \eqn{n} is the number of equations and \eqn{k} the number of regressors per
#'   equation, and a larger value for larger systems. A warning is issued if
#'   \eqn{h < 3(n + k)}.
#' @param verbose Logical. If \code{TRUE}, print additional information from the
#'   estimation.
#' @param brv Logical. If \code{TRUE}, allow the covariance matrix of the
#'   residuals to change across regimes.
#' @param robust Inference for the coefficient-break tests, the coefficient
#'   standard errors and the break-date confidence intervals. \code{"none"}
#'   gives Gaussian quasi-likelihood-ratio tests, valid when the errors are
#'   serially uncorrelated and conditionally homoskedastic. \code{"hc"} gives
#'   the robust SupF statistic of Qu and Perron (2007, eq. 20) with a
#'   heteroskedasticity-consistent covariance estimate without autocovariance
#'   terms, valid when the errors are a martingale difference sequence.
#'   \code{"hac"} uses instead the quadratic spectral kernel with the automatic
#'   bandwidth of Andrews (1991), allowing for serial correlation.
#' @param hetq Logical. If \code{TRUE}, allow regressor moments to vary across
#'   regimes when computing confidence intervals for break dates.
#' @param prewhit Logical. If \code{TRUE}, prewhiten the scores by a VAR(1)
#'   before the HAC estimation (Andrews and Monahan, 1992). Used only when
#'   \code{robust = "hac"}.
#' @param restart A bootstrap restarting procedure to check and improve
#'   convergence. It applies only when some coefficients are held constant;
#'   otherwise it is ignored with a warning. Set to
#'   \code{"none"} to turn it off, \code{"brbp"} for resampling bootstrap
#'   restarts, or \code{"brbp_residual"} for residual-based restarts.
#' @param restart_B Number of restarting values used when
#'   \code{restart != "none"}; 10 is recommended.
#' @param restart_verbose Logical. If \code{TRUE}, print restart progress.
#'
#' @return An object of class \code{"sysbreak"} containing:
#' \describe{
#'   \item{\code{spec}}{A summary of the setup, including the sample size
#'     \code{bigt}, the number of equations \code{n}, the total number of
#'     regression coefficients per regime \code{q}, the arguments \code{m},
#'     \code{trim}, \code{brv}, \code{robust}, \code{hetq}, \code{prewhit},
#'     \code{test_level} and \code{ci_level}, the input formula or formulas, the
#'     equation and coefficient names, the retained data rows \code{y_rows}
#'     after discarding missing values, and a record of the \code{coef_breaks}
#'     selector.}
#'   \item{\code{testing}}{The tests, reported at \code{test_level}:
#'     \code{statistic}, testing no break against \eqn{k} breaks for every
#'     \eqn{k = 1,\ldots,m}, with its \code{critical_values};
#'     \code{statistic_type}, \code{"SupLR"} or \code{"SupF"}, and \code{NA}
#'     when a requested robust test is unavailable; \code{wdmax} and
#'     \code{wdmax_critical_values}; the \code{sequential} tests of \eqn{l}
#'     against \eqn{l+1} breaks, each containing \code{l}, a \code{statistic}
#'     and its \code{critical_values}; and the number of
#'     \code{changing_parameters}, which includes the changing covariance
#'     parameters when \code{robust = "none"}. }
#'   \item{\code{br_selected}}{The outcome of the break-number determination at
#'     \code{test_level}: the selected number of \code{breaks}, the WDmax
#'     statistic with its critical value and significance flag, and the
#'     \code{sequential} steps taken. \code{breaks} is \code{0} when WDmax is
#'     not significant, and \code{NA} when the tests are unavailable for the
#'     requested specification.}
#'   \item{\code{estimation}}{The estimates related to the selected number of
#'    breaks: estimated \code{breaks}; their confidence intervals
#'    \code{ci_breaks}; \code{coefficients} by regimes and their confidence
#'    intervals; disturbance \code{covariance} matrices. This component is
#'    \code{NULL} when no break is detected, in which case only \code{testing}
#'    and \code{br_selected} are informative.}
#'   \item{\code{y}}{The responses used in estimation, rows \code{spec$y_rows}
#'     of \code{data}.}
#' }
#' Note that \code{m} is the \emph{maximum} number of breaks allowed. Testing
#' output is reported for every value from 1 through \code{m}, while
#' \code{estimation} corresponds to the number of breaks selected by the
#' procedure under \code{test_level}. The fit has \code{print()},
#' \code{\link[=summary.sysbreak]{summary()}},
#' \code{\link[=plot.sysbreak]{plot()}} and the accessors in
#' \code{\link{sysbreak-methods}}.
#'
#'
#' @details
#' The interpretation of \code{testing$statistic} depends on \code{robust}.
#' \itemize{
#'   \item \code{robust = "none"} gives the SupLR statistic. The code can test
#'     the coefficients specified by \code{coef_breaks}, the disturbance
#'     covariance when \code{brv = TRUE}, or both.
#'   \item \code{robust = "hc"} or \code{"hac"} gives the robust SupF
#'     statistic, for the coefficients only, because the covariance-break
#'     critical values are case dependent under heteroskedasticity or serial
#'     correlation; it assumes a constant error covariance. With
#'     \code{brv = TRUE}, the estimates, standard errors and confidence
#'     intervals account for the covariance breaks; under \code{"hc"} the
#'     intervals further require errors without ARCH-type dependence.
#' }
#'
#' @references
#' Andrews, D. W. K. (1991). Heteroskedasticity and autocorrelation consistent
#' covariance matrix estimation. \emph{Econometrica}, 59(3), 817--858.
#'
#' Andrews, D. W. K., and Monahan, J. C. (1992). An improved
#' heteroskedasticity and autocorrelation consistent covariance matrix
#' estimator. \emph{Econometrica}, 60(4), 953--966.
#'
#' Bai, J., and Perron, P. (1998). Estimating and testing linear models with
#' multiple structural changes. \emph{Econometrica}, 66(1), 47--78.
#'
#' Qu, Z., and Perron, P. (2007). Estimating and testing structural changes
#' in multivariate regressions. \emph{Econometrica}, 75(2), 459--502.
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
#'
#' fit <- sysbreak(cbind(y1, y2) ~ x, data = dat, m = 1, trim = 0.2)
#' fit$estimation$breaks
#'
#' @export
sysbreak <- function(
  model, data = NULL, m, coef_breaks = "all", test_level = 0.05,
  ci_level = 0.95, trim = 0.15, verbose = FALSE, brv = FALSE,
  robust = c("none", "hc", "hac"), hetq = FALSE, prewhit = FALSE,
  restart = c("none", "brbp", "brbp_residual"), restart_B = 10,
  restart_verbose = FALSE
) {

  ## ----- checking function input or error: ------

  # collapse `restart` to one of "none", "brbp" or "brbp_residual", default
  # to the first
  restart <- match.arg(restart)
  if (is.logical(robust)) {
    stop("`robust` must be \"none\", \"hc\" or \"hac\"; it replaces the logical `hac`.",
         call. = FALSE)
  }
  robust <- match.arg(robust)
  prewhit <- isTRUE(prewhit) && robust == "hac"
  # `test_level` is a significance level for testing; `ci_level` is a
  # coverage level for confidence interval.
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
  # check if the model has common regressors and generate regressor matrices
  # accordingly
  eqlist <- is.list(model) && length(model) > 0 &&
    all(vapply(model, inherits, logical(1), what = "formula"))
  if (eqlist) { # a list of equations case
    modelparts <- sur_mats(formulas = model, data = data)
    y <- modelparts$y
    z <- modelparts$z
    S <- modelparts$S
    eqz <- modelparts$eqz
    mspec <- list(model_type = "sur", formulas = modelparts$formulas)
  } else if (inherits(model, "formula")) { # a single formula case
    modelparts <- mb_mats(formula = model, data = data)
    y <- modelparts$y
    z <- modelparts$z
    eqz <- rep(list(colnames(z)), ncol(y))
    names(eqz) <- colnames(y)
    S <- diag(ncol(y) * ncol(z))
    mspec <- list(model_type = "common", formula = modelparts$formula)
  } else {
    stop(
      "`model` must be a formula such as `cbind(y1, y2) ~ x1 + x2`, or a list of equation formulas.",
      call. = FALSE
    )
  }
  # get coef names
  coef_names <- eqnames(colnames(y), eqz)
  # check if `trim` and `m` are compatible with sample size; error if too
  # large or small.
  m <- mcheck(m, trim, nrow(y))
  # turn the formula coef_breaks into form that mainp() takes
  br_info <- coef_restrict(
    coef_breaks = coef_breaks,
    m = m,
    y_names = colnames(y),
    eqz = eqz
  )
  # check if any coefficient or covariance is allowed to change
  if (!br_info$brbeta && !isTRUE(brv)) {
    stop(
      "`coef_breaks` selects no changing coefficients; use `brv = TRUE` for a pure covariance-break model.",
      call. = FALSE
    )
  }
  # ----done checking input----

  # Call the main procedure to get testing and estimation results. Here
  # verbose is set to false on purpose; results are printed later
  fit <- mainp(
    y = y,
    z = z,
    m = m,
    S = S,
    R = br_info$R, # restrictions
    R_path = br_info$R_path,
    trim = trim,
    verbose = FALSE,
    brbeta = br_info$brbeta,
    brv = brv,
    robust = robust,
    hetq = hetq,
    prewhit = prewhit,
    restart = restart,
    restart_B = restart_B,
    restart_verbose = restart_verbose
  )

  fit$spec <- c(fit$spec, mspec, list(
    y_names = colnames(y),
    coef_names = coef_names,
    y_rows = modelparts$y_rows,
    coef_breaks = br_info$info,
    test_level = test_level,
    ci_level = ci_level
  ))

  # `m` is the maximum number of breaks allowed.
  # Select the number of breaks as follows: the double-maximum tests decide
  # whether there is any break, then the sequential tests decide how many.
  br_selected <- brsel(fit$testing, level = test_level, m = m)
  class(br_selected) <- c("sysbreak_br_selected", class(br_selected))
  fit$br_selected <- br_selected
  # reduces the four level-indexed fields to the requested level and marks
  # the object for printing
  fit$testing <- tcols(fit$testing, test_level)

  if (is.na(br_selected$breaks) || br_selected$breaks == 0L) {
    # No break detected (or no test available); the testing output is kept.
    fit$estimation <- NULL
  } else {
    if (br_selected$breaks < m) {
      # The optimal partition for this count was already found during the
      # search, so only the estimation stage is repeated, at fixed dates.
      k <- br_selected$breaks
      br_k <- fit$breaks_by_m[[k]]
      brbeta <- fit$spec$brbeta
      nn <- ncol(y)
      bigt <- nrow(y)
      tr <- transf(y, z, S)
      maty <- tr$maty
      matz <- tr$matz

      # R is sized (k + 1) * q, so take the one built for k breaks
      R_k <- if (is.null(br_info$R_path)) {
        defR(k, ncol(matz), brbeta = brbeta, brv = brv)
      } else {
        br_info$R_path[[k]]
      }

      final <- r_estim(
        maty, matz, bigt, nn, k, br_k, R_k, brbeta = brbeta, brv = brv
      )

      coef_se <- coefse(
        maty = maty, matz = matz, n = nn, m = k, bigt = bigt, br = br_k,
        beta = final$nbeta, vv = final$nvv, R = R_k, robust = robust, brv = brv,
        prewhit = prewhit
      )

      ci <- interval2(
        maty, matz, k, nn, bigt, br_k, final$nbeta, final$nvv,
        hetq = hetq, robust = robust, brv = brv, brbeta = brbeta, prewhit = prewhit
      )

      fit$estimation <- list(
        breaks = br_k,
        coefficients = as.numeric(final$nbeta),
        vcov = unname(coef_se$vcov),
        covariance = final$nvv,
        ci_breaks = ci
      )
    }
    # residuals, one row per observation used; brtype() needs them
    tr <- transf(y, z, S)
    nb <- length(fit$estimation$breaks)
    res <- tr$maty - pzbar(tr$matz, nb, fit$estimation$breaks, nrow(y)) %*%
      fit$estimation$coefficients
    fit$estimation$residuals <- matrix(res, ncol = ncol(y), byrow = TRUE,
                                       dimnames = list(NULL, colnames(y)))
    # brtype() tests the coefficients with the covariance held constant, so a
    # brv = TRUE fit also keeps those estimates, at the same break dates
    if (brv && isTRUE(fit$spec$brbeta)) {
      R_c <- if (is.null(br_info$R_path)) {
        defR(nb, ncol(tr$matz), brbeta = TRUE, brv = FALSE)
      } else {
        br_info$R_path[[nb]]
      }
      c_fit <- r_estim(tr$maty, tr$matz, nrow(y), ncol(y), nb,
                       fit$estimation$breaks, R_c, brbeta = TRUE, brv = FALSE)
      c_se <- coefse(maty = tr$maty, matz = tr$matz, n = ncol(y), m = nb,
                     bigt = nrow(y), br = fit$estimation$breaks,
                     beta = c_fit$nbeta, vv = c_fit$nvv, R = R_c, robust = robust,
                     brv = FALSE, prewhit = prewhit)
      fit$estimation$const <- list(coefficients = as.numeric(c_fit$nbeta),
                                   vcov = unname(c_se$vcov))
    }
    ci_all <- fit$estimation$ci_breaks  # 95 and 90 percent; kept for confint()
    fit$estimation$ci_breaks <- cicols(ci_all, ci_level)
    # report dates as rows of `data`; an identity map unless rows were dropped
    yrow <- fit$spec$y_rows
    fit$estimation$breaks <- yrow[fit$estimation$breaks]
    ci <- fit$estimation$ci_breaks
    nlo <- sum(ci[, 1] < 1 | ci[, 1] > length(yrow))
    nhi <- sum(ci[, 2] < 1 | ci[, 2] > length(yrow))
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
    fit$estimation$ci_breaks[] <- yrow[pmin(pmax(ci, 1L), length(yrow))]
    fit$estimation <- eclass(fit$estimation, fit$spec)
    attr(fit$estimation, "ci_all") <- cirows(ci_all, yrow)
  }
  fit$breaks_by_m <- NULL
  # the responses used in estimation, rows `spec$y_rows` of `data`
  fit$y <- y

  if (verbose) {
    print(fit)
  }
  fit
}

#' Structural Breaks in a VAR
#'
#' Estimation of a VAR affected by structural breaks in the regression
#' coefficients and/or the variance covariance matrix of the residuals.
#'
#' @usage
#' varbreak(y, lags, m, intercept = TRUE, coef_breaks = "all",
#'   test_level = 0.05, ci_level = 0.95,
#'   trim = 0.15, verbose = FALSE, brv = FALSE,
#'   robust = c("none", "hc", "hac"), hetq = FALSE, prewhit = FALSE,
#'   restart = c("none", "brbp", "brbp_residual"),
#'   restart_B = 10, restart_verbose = FALSE, h = NULL)
#'
#' @param y Numeric response matrix with one column per variable.
#' @param lags Positive integer, number of lags.
#' @param m Integer, maximum number of breaks allowed.
#' @param intercept Logical. If \code{TRUE}, include an intercept in the
#'   regressors matrix.
#' @param coef_breaks Coefficient-break selector, applied after the VAR lag
#'   matrix is built. The default \code{"all"} lets all VAR coefficients change.
#'   Use \code{"none"} or \code{NULL} for no coefficient
#'   changes, for which \code{brv = TRUE} is required so that the covariance
#'   matrix of the residuals can change. Use a list such as
#'   \code{list(y1 ~ ., y2 ~ 0)} to restrict breaks to selected regression
#'   coefficients.
#' @param test_level Significance level used for selecting the number of breaks.
#'   One of \code{0.10}, \code{0.05}, \code{0.025} or \code{0.01}. WDmax
#'   decides whether there is any break; if it is not significant the
#'   procedure stops and no break is reported. Otherwise the sequential tests
#'   of \eqn{l} against \eqn{l+1} breaks are applied from one break upwards,
#'   stopping at the first insignificant test or at \code{m}.
#' @param ci_level Coverage probability of the reported confidence intervals
#'   for the break dates. Either \code{0.95} or \code{0.90}.
#' @param trim Trimming parameter. The minimum segment length is
#'   \code{h = round(trim * nrow(y))}. We recommend \eqn{h \ge 3(n + k)}, where
#'   \eqn{n} is the number of equations and \eqn{k} the number of regressors per
#'   equation, and a larger value for larger systems. A warning is issued if
#'   \eqn{h < 3(n + k)}.
#' @param verbose Logical. If \code{TRUE}, print additional progress
#'   information from the testing and estimation.
#' @param brv Logical. If \code{TRUE}, allow the disturbance covariance matrix
#'   to change across regimes.
#' @param robust Inference for the coefficient-break tests, the coefficient
#'   standard errors and the break-date confidence intervals. \code{"none"}
#'   gives Gaussian quasi-likelihood-ratio tests, valid when the errors are
#'   serially uncorrelated and conditionally homoskedastic. \code{"hc"} gives
#'   the robust SupF statistic of Qu and Perron (2007, eq. 20) with a
#'   heteroskedasticity-consistent covariance estimate without autocovariance
#'   terms, valid when the errors are a martingale difference sequence.
#'   \code{"hac"} uses instead the quadratic spectral kernel with the automatic
#'   bandwidth of Andrews (1991), allowing for serial correlation.
#'   With enough lags the VAR errors are serially uncorrelated, and
#'   \code{"hc"} suffices for the tests and standard errors; with
#'   \code{brv = TRUE} and conditionally heteroskedastic errors of the ARCH type,
#'   the confidence intervals require \code{"hac"}.
#' @param hetq Logical. If \code{TRUE}, allow regressor moments to vary across
#'   regimes when computing confidence intervals for break dates.
#' @param prewhit Logical. If \code{TRUE}, prewhiten the scores by a VAR(1)
#'   before the HAC estimation (Andrews and Monahan, 1992). Used only when
#'   \code{robust = "hac"}.
#' @param restart A bootstrap restarting procedure to check and improve
#'   convergence. It applies only when some coefficients are held constant;
#'   otherwise it is ignored with a warning. Set to
#'   \code{"none"} to turn it off, \code{"brbp"} for resampling bootstrap
#'   restarts, or \code{"brbp_residual"} for residual-based restarts.
#' @param restart_B Integer, number of restart values used when
#'   \code{restart != "none"}.
#' @param restart_verbose Logical. If \code{TRUE}, print restart progress.
#' @param h Optional positive integer for forecast horizon. If supplied,
#'   attach Diebold--Yilmaz generalized spillover measures in \code{spillover}.
#'
#' @return An object of class \code{c("varbreak", "sysbreak")} containing:
#' \describe{
#'   \item{\code{spec}}{A summary of the setup, including the sample size
#'     \code{bigt}, the number of equations \code{n}, the total number of
#'     coefficients per regime \code{q}, the arguments \code{m}, \code{lags},
#'     \code{intercept}, \code{trim}, \code{brv}, \code{robust}, \code{hetq},
#'     \code{prewhit}, \code{test_level} and \code{ci_level}, the equation and
#'     coefficient names, the retained input rows \code{y_rows} after discarding
#'     the lags, and the \code{coef_breaks} selector.}
#'   \item{\code{testing}}{The tests, reported at \code{test_level}:
#'     \code{statistic}, testing no break against \eqn{k} breaks for every
#'     \eqn{k = 1,\ldots,m}, with its \code{critical_values};
#'     \code{statistic_type}, \code{"SupLR"} or \code{"SupF"}, and \code{NA}
#'     when a requested robust test is unavailable; \code{wdmax} and
#'     \code{wdmax_critical_values}; the \code{sequential} tests of \eqn{l}
#'     against \eqn{l+1} breaks, each containing \code{l}, a \code{statistic}
#'     and its \code{critical_values}; and the number of
#'     \code{changing_parameters}, which includes the changing covariance
#'     parameters when \code{robust = "none"}. }
#'   \item{\code{br_selected}}{The outcome of the break-number determination at
#'     \code{test_level}: the selected number of \code{breaks}, the WDmax
#'     statistic with its critical value and significance flag, and the
#'     \code{sequential} steps taken. \code{breaks} is \code{0} when WDmax is
#'     not significant, and \code{NA} when the tests are unavailable for the
#'     requested specification.}
#'   \item{\code{estimation}}{The estimates related to the selected number of
#'    breaks: estimated \code{breaks}; their confidence intervals
#'    \code{ci_breaks}; \code{coefficients} by regimes and their confidence
#'    intervals; disturbance \code{covariance} matrices. This component is
#'    \code{NULL} when no break is detected, in which case only \code{testing}
#'    and \code{br_selected} are informative.}
#'   \item{\code{spillover}}{When \code{h} is supplied and at least one break
#'     is selected, regime-specific Diebold--Yilmaz generalized spillover
#'     measures returned by \code{\link{dyspill}}.}
#'   \item{\code{y}}{The responses used in estimation, rows \code{spec$y_rows}
#'     of \code{y}.}
#' }
#' Note that \code{m} is the \emph{maximum} number of breaks allowed. Testing
#' output is reported for every value from 1 through \code{m}, while
#' \code{estimation} corresponds to the number of breaks selected by the
#' procedure under \code{test_level}. The optional \code{spillover} component
#' is omitted when \code{h = NULL}. The fit has \code{print()},
#' \code{\link[=summary.sysbreak]{summary()}},
#' \code{\link[=plot.sysbreak]{plot()}} and the accessors in
#' \code{\link{sysbreak-methods}}.
#'
#' @details
#' The interpretation of \code{testing$statistic} depends on \code{robust}.
#' \itemize{
#'   \item \code{robust = "none"} gives the SupLR statistic. The code can test
#'     the coefficients specified by \code{coef_breaks}, the disturbance
#'     covariance when \code{brv = TRUE}, or both.
#'   \item \code{robust = "hc"} or \code{"hac"} gives the robust SupF
#'     statistic, for the coefficients only, because the covariance-break
#'     critical values are case dependent under heteroskedasticity or serial
#'     correlation; it assumes a constant error covariance. With
#'     \code{brv = TRUE}, the estimates, standard errors and confidence
#'     intervals account for the covariance breaks; under \code{"hc"} the
#'     intervals further require errors without ARCH-type dependence.
#' }
#'
#' @references
#' Andrews, D. W. K. (1991). Heteroskedasticity and autocorrelation consistent
#' covariance matrix estimation. \emph{Econometrica}, 59(3), 817--858.
#'
#' Andrews, D. W. K., and Monahan, J. C. (1992). An improved
#' heteroskedasticity and autocorrelation consistent covariance matrix
#' estimator. \emph{Econometrica}, 60(4), 953--966.
#'
#' Bai, J., and Perron, P. (1998). Estimating and testing linear models with
#' multiple structural changes. \emph{Econometrica}, 66(1), 47--78.
#'
#' Qu, Z., and Perron, P. (2007). Estimating and testing structural changes
#' in multivariate regressions. \emph{Econometrica}, 75(2), 459--502.
#'
#' @examples
#' set.seed(12)
#' bigt <- 45
#' y <- matrix(0, bigt, 2)
#' for (tt in 2:bigt) {
#'   shift <- if (tt > 24) c(0.8, -0.6) else c(0, 0)
#'   y[tt, ] <- shift +
#'     c(0.3 * y[tt - 1, 1], -0.2 * y[tt - 1, 2]) +
#'     rnorm(2, sd = 0.3)
#' }
#'
#' fit <- varbreak(y, lags = 1, m = 1, trim = 0.2)
#' fit$estimation$breaks
#'
#' @export
varbreak <- function(
  y, lags, m, intercept = TRUE, coef_breaks = "all", test_level = 0.05,
  ci_level = 0.95, trim = 0.15, verbose = FALSE, brv = FALSE,
  robust = c("none", "hc", "hac"), hetq = FALSE, prewhit = FALSE,
  restart = c("none", "brbp", "brbp_residual"), restart_B = 10,
  restart_verbose = FALSE, h = NULL
) {
  ## ----- checking function input or error: ------

  # collapse `restart` to one of "none", "brbp" or "brbp_residual", default
  # to the first
  restart <- match.arg(restart)
  if (is.logical(robust)) {
    stop("`robust` must be \"none\", \"hc\" or \"hac\"; it replaces the logical `hac`.",
         call. = FALSE)
  }
  robust <- match.arg(robust)
  prewhit <- isTRUE(prewhit) && robust == "hac"
  # `test_level` is a significance level for testing; `ci_level` is a
  # coverage level for confidence interval.
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
  modelparts <- var_mats(y = y, lags = lags, intercept = intercept)
  eqz <- rep(list(colnames(modelparts$z)), ncol(modelparts$y))
  names(eqz) <- colnames(modelparts$y)
  coef_names <- eqnames(colnames(modelparts$y), eqz)
  S <- diag(ncol(modelparts$y) * ncol(modelparts$z))

  # Validate `trim`, and reject `m` if `m + 1` segments of round(trim * T)
  # observations cannot fit in the sample.  Done before the restriction
  # matrices are built, since those are sized by the number of breaks.
  m <- mcheck(m, trim, nrow(modelparts$y))

  br_info <- coef_restrict(
    coef_breaks = coef_breaks,
    m = m,
    y_names = colnames(modelparts$y),
    eqz = eqz
  )
  if (!br_info$brbeta && !isTRUE(brv)) {
    stop(
      "`coef_breaks` selects no changing coefficients; use `brv = TRUE` for a pure covariance-break model.",
      call. = FALSE
    )
  }

  # `verbose = FALSE` here on purpose; see the note in `sysbreak()`.  The
  # completed object is printed at the end of this function instead.
  fit <- mainp(
    y = modelparts$y,
    z = modelparts$z,
    m = m,
    S = S,
    R = br_info$R,
    R_path = br_info$R_path,
    trim = trim,
    verbose = FALSE,
    brbeta = br_info$brbeta,
    brv = brv,
    robust = robust,
    hetq = hetq,
    prewhit = prewhit,
    restart = restart,
    restart_B = restart_B,
    restart_verbose = restart_verbose
  )

  fit$spec <- c(fit$spec, list(
    model_type = "var",
    lags = as.integer(lags),
    intercept = intercept,
    y_names = colnames(modelparts$y),
    coef_names = coef_names,
    y_rows = modelparts$y_rows,
    coef_breaks = br_info$info,
    test_level = test_level,
    ci_level = ci_level
  ))

  # `m` is the maximum number of breaks allowed.
  # Select the number of breaks as follows: WDmax decides whether there is any
  # break, then the sequential tests decide how many.
  br_selected <- brsel(fit$testing, level = test_level, m = m)
  class(br_selected) <- c("sysbreak_br_selected", class(br_selected))
  fit$br_selected <- br_selected

  # reduces the four level-indexed fields to the requested level and marks
  # the object for printing
  fit$testing <- tcols(fit$testing, test_level)

  if (is.na(br_selected$breaks) || br_selected$breaks == 0L) {
    # No break detected (or no test available); the testing output is kept.
    fit$estimation <- NULL
  } else {
    if (br_selected$breaks < m) {
      # The optimal partition was already found during the
      # search; only the estimation stage is repeated, at fixed dates.
      k <- br_selected$breaks
      br_k <- fit$breaks_by_m[[k]]
      brbeta <- fit$spec$brbeta
      nn <- ncol(modelparts$y)
      bigt <- nrow(modelparts$y)
      tr <- transf(modelparts$y, modelparts$z, S)
      maty <- tr$maty
      matz <- tr$matz

      # R is sized (k + 1) * q, so take the one built for k breaks
      R_k <- if (is.null(br_info$R_path)) {
        defR(k, ncol(matz), brbeta = brbeta, brv = brv)
      } else {
        br_info$R_path[[k]]
      }

      final <- r_estim(
        maty, matz, bigt, nn, k, br_k, R_k, brbeta = brbeta, brv = brv
      )

      coef_se <- coefse(
        maty = maty, matz = matz, n = nn, m = k, bigt = bigt, br = br_k,
        beta = final$nbeta, vv = final$nvv, R = R_k, robust = robust, brv = brv,
        prewhit = prewhit
      )

      ci <- interval2(
        maty, matz, k, nn, bigt, br_k, final$nbeta, final$nvv,
        hetq = hetq, robust = robust, brv = brv, brbeta = brbeta, prewhit = prewhit
      )

      fit$estimation <- list(
        breaks = br_k,
        coefficients = as.numeric(final$nbeta),
        vcov = unname(coef_se$vcov),
        covariance = final$nvv,
        ci_breaks = ci
      )
    }
    # residuals, one row per observation used; brtype() needs them
    tr <- transf(modelparts$y, modelparts$z, S)
    nb <- length(fit$estimation$breaks)
    res <- tr$maty - pzbar(tr$matz, nb, fit$estimation$breaks,
                           nrow(modelparts$y)) %*% fit$estimation$coefficients
    fit$estimation$residuals <- matrix(res, ncol = ncol(modelparts$y), byrow = TRUE,
                                       dimnames = list(NULL, colnames(modelparts$y)))
    # the regressors, same rows; the bootstrap in dyspill() starts from them
    fit$estimation$regressors <- modelparts$z
    # brtype() tests the coefficients with the covariance held constant, so a
    # brv = TRUE fit also keeps those estimates, at the same break dates
    if (brv && isTRUE(fit$spec$brbeta)) {
      R_c <- if (is.null(br_info$R_path)) {
        defR(nb, ncol(tr$matz), brbeta = TRUE, brv = FALSE)
      } else {
        br_info$R_path[[nb]]
      }
      c_fit <- r_estim(tr$maty, tr$matz, nrow(modelparts$y), ncol(modelparts$y),
                       nb, fit$estimation$breaks, R_c, brbeta = TRUE, brv = FALSE)
      c_se <- coefse(maty = tr$maty, matz = tr$matz, n = ncol(modelparts$y),
                     m = nb, bigt = nrow(modelparts$y),
                     br = fit$estimation$breaks, beta = c_fit$nbeta,
                     vv = c_fit$nvv, R = R_c, robust = robust, brv = FALSE,
                     prewhit = prewhit)
      fit$estimation$const <- list(coefficients = as.numeric(c_fit$nbeta),
                                   vcov = unname(c_se$vcov))
    }
    ci_all <- fit$estimation$ci_breaks  # 95 and 90 percent; kept for confint()
    fit$estimation$ci_breaks <- cicols(ci_all, ci_level)
    # report dates as rows of `y`; `y_rows` drops the first `lags` observations
    yrow <- fit$spec$y_rows
    fit$estimation$breaks <- yrow[fit$estimation$breaks]
    ci <- fit$estimation$ci_breaks
    nlo <- sum(ci[, 1] < 1 | ci[, 1] > length(yrow))
    nhi <- sum(ci[, 2] < 1 | ci[, 2] > length(yrow))
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
    fit$estimation$ci_breaks[] <- yrow[pmin(pmax(ci, 1L), length(yrow))]
    fit$estimation <- eclass(fit$estimation, fit$spec)
    attr(fit$estimation, "ci_all") <- cirows(ci_all, yrow)
  }
  fit$breaks_by_m <- NULL
  # the responses used in estimation, rows `spec$y_rows` of `y`
  fit$y <- modelparts$y
  class(fit) <- c("varbreak", class(fit))

  # Spillovers. This is run after regime estimation.
  if (!is.null(h)) {
    if (is.null(fit$estimation)) {
      warning("No break was selected; `h` is ignored and no spillover measure is reported.",
              call. = FALSE)
    } else {
      fit$spillover <- dyspill(fit, h = h)
    }
  }
  if (verbose) {
    print(fit)
  }
  fit
}

#' Structural Breaks in a network VAR
#'
#' Estimates structural breaks in a VAR whose coefficients are restricted by a
#' fixed network matrix. With one lag, the model is
#' \deqn{y_t = \mu_j + (\alpha_{1,j} W + \beta_{1,j} I)y_{t-1} + u_t,}
#' for observations in regime \eqn{j}. With multiple lags, the model is
#' \deqn{y_t = \mu_j + \sum_{\ell=1}^p
#'   (\alpha_{\ell,j} W + \beta_{\ell,j} I)y_{t-\ell} + u_t.}
#'
#' @usage
#' netbreak(y, W, lags, m, intercept = TRUE,
#'   coef_breaks = "all", row_normalize = FALSE,
#'   test_level = 0.05, ci_level = 0.95,
#'   trim = 0.15, verbose = FALSE, brv = FALSE,
#'   robust = c("none", "hc", "hac"), hetq = FALSE, prewhit = FALSE,
#'   restart = c("none", "brbp", "brbp_residual"),
#'   restart_B = 10, restart_verbose = FALSE, h = NULL)
#'
#' @param y Numeric matrix/data frame with one column per network node.
#' @param W Numeric square network matrix. Row \eqn{i} gives the weights used to
#'   form the network lag for equation \eqn{i}.
#' @param lags Positive integer, number of lags.
#' @param m The maximum number of breaks allowed.
#' @param intercept Logical. If \code{TRUE}, include equation-specific
#'   intercepts \eqn{\mu_i}.
#' @param coef_breaks Parameters allowed to change across regimes. The
#'   parameters are one intercept per node, \code{mu:<node>}; the network
#'   coefficients \eqn{\alpha_\ell}, named \code{network_lag1},
#'   \code{network_lag2}, and so on; and the own-lag coefficients
#'   \eqn{\beta_\ell}, named \code{own_lag1}, \code{own_lag2}, and so on.
#'   The default \code{"all"} lets every parameter change. \code{"none"} or
#'   \code{NULL} holds them all constant, which requires \code{brv = TRUE} so
#'   that the disturbance covariance matrix can change; otherwise the code
#'   will return an error message. For anything in between, supply a character
#'   vector of parameter names or of these groups:
#'   \itemize{
#'     \item \code{"intercepts"}, or \code{"mu"}: all \eqn{\mu_i}.
#'     \item \code{"network"}: all \eqn{\alpha_\ell}.
#'     \item \code{"own"}: all \eqn{\beta_\ell}.
#'     \item \code{"dynamic"}: every lag coefficient, network and own.
#'     \item \code{"lag1"}, \code{"lag2"}, ...: \eqn{\alpha_\ell} and
#'       \eqn{\beta_\ell} at that lag.
#'     \item \code{"alpha1"}, \code{"beta1"}, ...: \eqn{\alpha_\ell} or
#'       \eqn{\beta_\ell} on its own.
#'   }
#'   Entries are case-insensitive and may be combined, as in
#'   \code{c("intercepts", "network_lag1")}.
#' @param row_normalize Logical. If \code{TRUE}, row-normalize \code{W}; rows
#'   with zero row sums are left as zero rows.
#' @param test_level Significance level used for selecting the number of breaks.
#'   One of \code{0.10}, \code{0.05}, \code{0.025} or \code{0.01}. WDmax
#'   decides whether there is any break; if it is not significant the
#'   procedure stops and no break is reported. Otherwise the sequential tests
#'   of \eqn{l} against \eqn{l+1} breaks are applied from one break upwards,
#'   stopping at the first insignificant test or at \code{m}.
#' @param ci_level Coverage probability of the reported confidence intervals
#'   for the break dates. Either \code{0.95} or \code{0.90}.
#' @param trim Trimming parameter. The minimum segment length is
#'   \code{h = round(trim * nrow(y))}. We recommend \eqn{h \ge 3(n + k)}, where
#'   \eqn{n} is the number of equations and \eqn{k} the number of regressors per
#'   equation, and a larger value for larger systems. A warning is issued if
#'   \eqn{h < 3(n + k)}.
#' @param verbose Logical. If \code{TRUE}, print additional progress
#'   information from the testing and estimation.
#' @param brv Logical. If \code{TRUE}, allow the disturbance covariance matrix
#'   to change across regimes.
#' @param robust Inference for the coefficient-break tests, the coefficient
#'   standard errors and the break-date confidence intervals. \code{"none"}
#'   gives Gaussian quasi-likelihood-ratio tests, valid when the errors are
#'   serially uncorrelated and conditionally homoskedastic. \code{"hc"} gives
#'   the robust SupF statistic of Qu and Perron (2007, eq. 20) with a
#'   heteroskedasticity-consistent covariance estimate without autocovariance
#'   terms, valid when the errors are a martingale difference sequence.
#'   \code{"hac"} uses instead the quadratic spectral kernel with the automatic
#'   bandwidth of Andrews (1991), allowing for serial correlation.
#'   With enough lags the VAR errors are serially uncorrelated, and
#'   \code{"hc"} suffices for the tests and standard errors; with
#'   \code{brv = TRUE} and conditionally heteroskedastic errors of the ARCH type,
#'   the confidence intervals require \code{"hac"}.
#' @param hetq Logical. If \code{TRUE}, allow regressor moments to vary across
#'   regimes when computing confidence intervals for break dates.
#' @param prewhit Logical. If \code{TRUE}, prewhiten the scores by a VAR(1)
#'   before the HAC estimation (Andrews and Monahan, 1992). Used only when
#'   \code{robust = "hac"}.
#' @param restart A bootstrap restarting procedure to check and improve
#'   convergence. It applies only when some coefficients are held constant;
#'   otherwise it is ignored with a warning. Set to
#'   \code{"none"} to turn it off, \code{"brbp"} for resampling bootstrap
#'   restarts, or \code{"brbp_residual"} for residual-based restarts.
#' @param restart_B Integer, number of restart values used when
#'   \code{restart != "none"}.
#' @param restart_verbose Logical. If \code{TRUE}, print restart progress.
#' @param h Optional positive integer for forecast horizon. If supplied,
#'   attach Diebold--Yilmaz generalized spillover measures in \code{spillover}.
#'
#' @return An object of class \code{c("netbreak", "sysbreak")} containing:
#' \describe{
#'   \item{\code{spec}}{A summary of the setup, including the sample size
#'     \code{bigt}, the number of equations \code{n}, the total number of NVAR
#'     parameters per regime \code{q}, the arguments \code{m}, \code{lags},
#'     \code{intercept}, \code{row_normalize}, \code{trim}, \code{brv},
#'     \code{robust}, \code{hetq}, \code{prewhit}, \code{test_level} and
#'     \code{ci_level}, the network matrix \code{W}, the equation and parameter
#'     names, the retained input rows \code{y_rows} after discarding the lags,
#'     and a record of the \code{coef_breaks} selector.}
#'   \item{\code{testing}}{The tests, reported at \code{test_level}:
#'     \code{statistic}, testing no break against \eqn{k} breaks for every
#'     \eqn{k = 1,\ldots,m}, with its \code{critical_values};
#'     \code{statistic_type}, \code{"SupLR"} or \code{"SupF"}, and \code{NA}
#'     when a requested robust test is unavailable; \code{wdmax} and
#'     \code{wdmax_critical_values}; the \code{sequential} tests of \eqn{l}
#'     against \eqn{l+1} breaks, each containing \code{l}, a \code{statistic}
#'     and its \code{critical_values}; and the number of
#'     \code{changing_parameters}, which includes the changing covariance
#'     parameters when \code{robust = "none"}. }
#'   \item{\code{br_selected}}{The outcome of the break-number determination at
#'     \code{test_level}: the selected number of \code{breaks}, the WDmax
#'     statistic with its critical value and significance flag, and the
#'     \code{sequential} steps taken. \code{breaks} is \code{0} when WDmax is
#'     not significant, and \code{NA} when the tests are unavailable for the
#'     requested specification.}
#'   \item{\code{estimation}}{The estimates related to the selected number of
#'    breaks: estimated \code{breaks}; their confidence intervals
#'    \code{ci_breaks}; \code{coefficients} by regimes and their confidence
#'    intervals; disturbance \code{covariance} matrices. This component is
#'    \code{NULL} when no break is detected, in which case only \code{testing}
#'    and \code{br_selected} are informative.}
#'   \item{\code{spillover}}{When \code{h} is supplied and at least one break
#'     is selected, regime-specific Diebold--Yilmaz generalized spillover
#'     measures returned by \code{\link{dyspill}}.}
#'   \item{\code{y}}{The responses used in estimation, rows \code{spec$y_rows}
#'     of \code{y}.}
#' }
#' Note that \code{m} is the \emph{maximum} number of breaks allowed. Testing
#' output is reported for every value from 1 through \code{m}, while
#' \code{estimation} corresponds to the number of breaks selected by the
#' procedure under \code{test_level}. Parameters excluded by
#' \code{coef_breaks} are repeated across regime blocks. The optional
#' \code{spillover} component is omitted when \code{h = NULL}. The fit has
#' \code{print()}, \code{\link[=summary.sysbreak]{summary()}},
#' \code{\link[=plot.sysbreak]{plot()}} and the accessors in
#' \code{\link{sysbreak-methods}}.
#'
#' @details
#' The interpretation of \code{testing$statistic} depends on \code{robust}.
#' \itemize{
#'   \item \code{robust = "none"} gives the SupLR statistic. The code can test
#'     the parameters specified by \code{coef_breaks}, the disturbance
#'     covariance when \code{brv = TRUE}, or both.
#'   \item \code{robust = "hc"} or \code{"hac"} gives the robust SupF
#'     statistic, for the parameters only, because the covariance-break
#'     critical values are case dependent under heteroskedasticity or serial
#'     correlation; it assumes a constant error covariance. With
#'     \code{brv = TRUE}, the estimates, standard errors and confidence
#'     intervals account for the covariance breaks; under \code{"hc"} the
#'     intervals further require errors without ARCH-type dependence.
#' }
#'
#' @references
#' Andrews, D. W. K. (1991). Heteroskedasticity and autocorrelation consistent
#' covariance matrix estimation. \emph{Econometrica}, 59(3), 817--858.
#'
#' Andrews, D. W. K., and Monahan, J. C. (1992). An improved
#' heteroskedasticity and autocorrelation consistent covariance matrix
#' estimator. \emph{Econometrica}, 60(4), 953--966.
#'
#' Zhu, X., Pan, R., Li, G., Liu, Y., and Wang, H. (2017). Network vector
#' autoregression. \emph{The Annals of Statistics}, 45(3), 1096--1123.
#'
#' Qu, Z., and Perron, P. (2007). Estimating and testing structural changes
#' in multivariate regressions. \emph{Econometrica}, 75(2), 459--502.
#'
#' @examples
#' set.seed(5)
#' bigt <- 80
#' W <- matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE)
#' y <- matrix(0, bigt, 2)
#' colnames(y) <- c("y1", "y2")
#' for (tt in 2:bigt) {
#'   alpha <- if (tt <= 40) 0.05 else 0.75
#'   beta <- 0.1
#'   y[tt, ] <- alpha * as.numeric(W %*% y[tt - 1, ]) +
#'     beta * y[tt - 1, ] + rnorm(2, sd = 0.2)
#' }
#'
#' fit <- netbreak(
#'   y = y, W = W, lags = 1, m = 1,
#'   coef_breaks = c("network", "own"),
#'   trim = 0.2
#' )
#' fit$estimation$breaks
#'
#' @export
netbreak <- function(
  y, W, lags, m, intercept = TRUE, coef_breaks = "all",
  row_normalize = FALSE, test_level = 0.05, ci_level = 0.95, trim = 0.15,
  verbose = FALSE, brv = FALSE, robust = c("none", "hc", "hac"), hetq = FALSE, prewhit = FALSE,
  restart = c("none", "brbp", "brbp_residual"), restart_B = 10,
  restart_verbose = FALSE, h = NULL
) {

  restart <- match.arg(restart)
  if (is.logical(robust)) {
    stop("`robust` must be \"none\", \"hc\" or \"hac\"; it replaces the logical `hac`.",
         call. = FALSE)
  }
  robust <- match.arg(robust)
  prewhit <- isTRUE(prewhit) && robust == "hac"
  # `test_level` is a significance level for testing; `ci_level` is a
  # coverage level for confidence interval.
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
  modelparts <- nvar_mats(
    y = y,
    W = W,
    lags = lags,
    intercept = intercept,
    row_normalize = row_normalize
  )
  # Validate `trim`, and reject `m` if `m + 1` segments of round(trim * T)
  # observations cannot fit in the sample.  Done before the restriction
  # matrices are built, since those are sized by the number of breaks.
  m <- mcheck(m, trim, nrow(modelparts$y))

  br_info <- nvar_restrict(
    coef_breaks = coef_breaks,
    m = m,
    param_names = modelparts$param_names,
    groups = modelparts$param_groups
  )
  if (!br_info$brbeta && !isTRUE(brv)) {
    stop(
      "`coef_breaks` selects no changing coefficients; use `brv = TRUE` for a pure covariance-break model.",
      call. = FALSE
    )
  }

  # `verbose = FALSE` here on purpose; see the note in `sysbreak()`.  The
  # completed object is printed at the end of this function instead.
  fit <- mainp(
    y = modelparts$y,
    z = modelparts$z,
    m = m,
    S = modelparts$S,
    R = br_info$R,
    R_path = br_info$R_path,
    trim = trim,
    verbose = FALSE,
    brbeta = br_info$brbeta,
    brv = brv,
    robust = robust,
    hetq = hetq,
    prewhit = prewhit,
    sumq = br_info$sumq,
    restart = restart,
    restart_B = restart_B,
    restart_verbose = restart_verbose
  )

  fit$spec <- c(fit$spec, list(
    model_type = "network_var",
    lags = as.integer(lags),
    intercept = intercept,
    row_normalize = row_normalize,
    y_names = colnames(modelparts$y),
    coef_names = modelparts$param_names,
    y_rows = modelparts$y_rows,
    W = modelparts$W,
    coef_breaks = br_info[c("selector", "selected")],
    test_level = test_level,
    ci_level = ci_level
  ))

  # `m` is the maximum number of breaks allowed.
  # Select the number of breaks as follows: WDmax decides whether there is any
  # break, then the sequential tests decide how many.
  br_selected <- brsel(fit$testing, level = test_level, m = m)
  class(br_selected) <- c("sysbreak_br_selected", class(br_selected))
  fit$br_selected <- br_selected
  # reduces the four level-indexed fields to the requested level and marks
  # the object for printing
  fit$testing <- tcols(fit$testing, test_level)

  if (is.na(br_selected$breaks) || br_selected$breaks == 0L) {
    # No break detected (or no test available); the testing output is kept.
    fit$estimation <- NULL
  } else {
    if (br_selected$breaks < m) {
      # The optimal partition for this count was already found during the
      # search, so only the estimation stage is repeated, at fixed dates.
      k <- br_selected$breaks
      br_k <- fit$breaks_by_m[[k]]
      brbeta <- fit$spec$brbeta
      nn <- ncol(modelparts$y)
      bigt <- nrow(modelparts$y)
      # `S` pools parameters across equations here, so it is not the identity
      tr <- transf(modelparts$y, modelparts$z, modelparts$S)
      maty <- tr$maty
      matz <- tr$matz

      # R is sized (k + 1) * q, so take the one built for k breaks
      R_k <- if (is.null(br_info$R_path)) {
        defR(k, ncol(matz), brbeta = brbeta, brv = brv)
      } else {
        br_info$R_path[[k]]
      }

      final <- r_estim(
        maty, matz, bigt, nn, k, br_k, R_k, brbeta = brbeta, brv = brv
      )

      coef_se <- coefse(
        maty = maty, matz = matz, n = nn, m = k, bigt = bigt, br = br_k,
        beta = final$nbeta, vv = final$nvv, R = R_k, robust = robust, brv = brv,
        prewhit = prewhit
      )

      ci <- interval2(
        maty, matz, k, nn, bigt, br_k, final$nbeta, final$nvv,
        hetq = hetq, robust = robust, brv = brv, brbeta = brbeta, prewhit = prewhit
      )

      fit$estimation <- list(
        breaks = br_k,
        coefficients = as.numeric(final$nbeta),
        vcov = unname(coef_se$vcov),
        covariance = final$nvv,
        ci_breaks = ci
      )
    }
    # residuals, one row per observation used; brtype() needs them
    tr <- transf(modelparts$y, modelparts$z, modelparts$S)
    nb <- length(fit$estimation$breaks)
    res <- tr$maty - pzbar(tr$matz, nb, fit$estimation$breaks,
                           nrow(modelparts$y)) %*% fit$estimation$coefficients
    fit$estimation$residuals <- matrix(res, ncol = ncol(modelparts$y), byrow = TRUE,
                                       dimnames = list(NULL, colnames(modelparts$y)))
    # the regressors, same rows; the bootstrap in dyspill() starts from them
    fit$estimation$regressors <- modelparts$z
    # brtype() tests the coefficients with the covariance held constant, so a
    # brv = TRUE fit also keeps those estimates, at the same break dates
    if (brv && isTRUE(fit$spec$brbeta)) {
      R_c <- if (is.null(br_info$R_path)) {
        defR(nb, ncol(tr$matz), brbeta = TRUE, brv = FALSE)
      } else {
        br_info$R_path[[nb]]
      }
      c_fit <- r_estim(tr$maty, tr$matz, nrow(modelparts$y), ncol(modelparts$y),
                       nb, fit$estimation$breaks, R_c, brbeta = TRUE, brv = FALSE)
      c_se <- coefse(maty = tr$maty, matz = tr$matz, n = ncol(modelparts$y),
                     m = nb, bigt = nrow(modelparts$y),
                     br = fit$estimation$breaks, beta = c_fit$nbeta,
                     vv = c_fit$nvv, R = R_c, robust = robust, brv = FALSE,
                     prewhit = prewhit)
      fit$estimation$const <- list(coefficients = as.numeric(c_fit$nbeta),
                                   vcov = unname(c_se$vcov))
    }
    ci_all <- fit$estimation$ci_breaks  # 95 and 90 percent; kept for confint()
    fit$estimation$ci_breaks <- cicols(ci_all, ci_level)
    # report dates as rows of `y`; `y_rows` drops the first `lags` observations
    yrow <- fit$spec$y_rows
    fit$estimation$breaks <- yrow[fit$estimation$breaks]
    ci <- fit$estimation$ci_breaks
    nlo <- sum(ci[, 1] < 1 | ci[, 1] > length(yrow))
    nhi <- sum(ci[, 2] < 1 | ci[, 2] > length(yrow))
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
    fit$estimation$ci_breaks[] <- yrow[pmin(pmax(ci, 1L), length(yrow))]
    fit$estimation <- eclass(fit$estimation, fit$spec)
    attr(fit$estimation, "ci_all") <- cirows(ci_all, yrow)
  }
  fit$breaks_by_m <- NULL
  # the responses used in estimation, rows `spec$y_rows` of `y`
  fit$y <- modelparts$y
  class(fit) <- c("netbreak", class(fit))

  # Spillovers describe the regimes, so this runs after selection.  With no
  # break selected there is no estimation to build them from.
  if (!is.null(h)) {
    if (is.null(fit$estimation)) {
      warning("No break was selected, so `h` is ignored and no spillover is reported.",
              call. = FALSE)
    } else {
      fit$spillover <- dyspill(fit, h = h)
    }
  }
  if (verbose) {
    print(fit)
  }
  fit
}

#' Build Restriction Matrices from Coefficient-Break Selectors
#'
#' Translates a `coef_breaks` selector into the arguments `mainp()` needs.
#' This is the formula-interface counterpart of `nvar_restrict()`, which does
#' the same job for the Network VAR parameter groups.
#'
#' Returns `R` and the full `R_path` (`NULL` for both when no restriction is
#' needed), the implied `brbeta` flag, and the selector metadata that the
#' wrappers record in `spec$coef_breaks`.
#'
#' @noRd
coef_restrict <- function(coef_breaks, m, y_names, eqz) {
  if (is.null(coef_breaks)) {
    return(list(R = NULL, R_path = NULL, brbeta = FALSE, info = NULL))
  }

  # Which coefficients change does not depend on the number of breaks, so the
  # coefficient map and the selector are computed once.
  cmap <- coefmap(y_names, eqz)
  selected <- breaks_sel(
    coef_breaks = coef_breaks, cmap = cmap, y_names = y_names
  )
  q <- nrow(cmap)
  nchg <- sum(selected)
  info <- list(
    selector = breaks_lab(coef_breaks),
    selected = cmap$coefficient[selected]
  )

  if (nchg == q) {
    # Every coefficient changes: no restriction.
    return(list(R = NULL, R_path = NULL, brbeta = TRUE, info = info))
  }
  if (nchg == 0L) {
    # Nothing changes: a pure covariance-break model, which needs `brv`.
    return(list(R = NULL, R_path = NULL, brbeta = FALSE, info = info))
  }

  R_path <- lapply(seq_len(m), function(k) {
    buildR(selected, cmap$coefficient, k)
  })
  list(R = R_path[[m]], R_path = R_path, brbeta = TRUE, info = info)
}

#' Reparameterization Matrix for a Given Selection of Changing Parameters
#'
#' Builds the matrix `R` in `beta = R theta` for `m` breaks.  Regime one always
#' loads the base block; a parameter marked in `selected` gets its own column
#' in each later regime, while an unselected one keeps reusing the base column.
#'
#' Which parameters change is regime-invariant, so `selected` is computed once
#' by the caller and this routine is called for each candidate break count.
#' It serves both the formula selectors and the Network VAR parameter groups,
#' which differ only in the labels used for the dimnames.
#'
#' @noRd
buildR <- function(selected, labels, m) {
  m <- as.integer(m)
  if (length(m) != 1 || is.na(m) || m < 0) {
    stop("`m` must be a non-negative integer.")
  }

  q <- length(selected)
  sel_idx <- which(selected)
  nchg <- length(sel_idx)

  R <- matrix(0, nrow = (m + 1) * q, ncol = q + m * nchg)
  rownames(R) <- unlist(lapply(seq_len(m + 1), function(regime) {
    paste0("regime", regime, ":", labels)
  }), use.names = FALSE)
  extra_names <- character()
  if (nchg > 0L) {
    extra_names <- unlist(lapply(seq_len(m), function(break_i) {
      paste0("regime", break_i + 1, ":", labels[sel_idx])
    }), use.names = FALSE)
  }
  colnames(R) <- c(paste0("base:", labels), extra_names)

  sel_pos <- match(seq_len(q), sel_idx)
  for (regime in seq_len(m + 1)) {
    for (par_i in seq_len(q)) {
      row_i <- (regime - 1) * q + par_i
      if (regime == 1 || !selected[[par_i]]) {
        R[row_i, par_i] <- 1
      } else {
        col_i <- q + (regime - 2) * nchg + sel_pos[[par_i]]
        R[row_i, col_i] <- 1
      }
    }
  }
  R
}

#' Select Coefficients allowed to change
#'
#' @noRd
breaks_sel <- function(coef_breaks, cmap, y_names) {
  if (is.character(coef_breaks) && length(coef_breaks) == 1L) {
    mode <- match.arg(tolower(coef_breaks), c("all", "none"))
    return(rep(mode == "all", nrow(cmap)))
  }

  if (!is.list(coef_breaks) || length(coef_breaks) < 1L) {
    stop(
      "`coef_breaks` must be \"all\", \"none\", or a list of formulas.",
      call. = FALSE
    )
  }

  selected <- rep(FALSE, nrow(cmap))
  list_names <- names(coef_breaks)
  if (is.null(list_names)) {
    list_names <- rep("", length(coef_breaks))
  }

  for (i in seq_along(coef_breaks)) {
    spec <- coef_breaks[[i]]
    if (!inherits(spec, "formula")) {
      stop("Every element of `coef_breaks` must be a formula.")
    }

    if (length(spec) == 3L) {
      equation <- rname(spec)
      rhs <- spec[[3]]
    } else if (length(spec) == 2L && nzchar(list_names[[i]])) {
      equation <- list_names[[i]]
      rhs <- spec[[2]]
    } else {
      stop(
        "Each `coef_breaks` formula must be two-sided, or named if one-sided.",
        call. = FALSE
      )
    }

    eq_i <- match(equation, y_names)
    if (is.na(eq_i)) {
      stop(
        sprintf("Equation `%s` in `coef_breaks` is not one of: %s.",
                equation, paste(y_names, collapse = ", ")),
        call. = FALSE
      )
    }

    available <- cmap$term[cmap$eq_index == eq_i]
    terms_i <- break_terms(rhs, available, equation)
    if (length(terms_i) == 0L) {
      next
    }

    idx <- which(cmap$eq_index == eq_i &
                 cmap$term %in% terms_i)
    selected[idx] <- TRUE
  }

  selected
}

#' Takes the right-hand side of one coef_breaks formula and returns the
#' coefficient names it picks out for that equation.
#'
#' @noRd
break_terms <- function(rhs, available, equation) {
  rhs_text <- paste(deparse(rhs), collapse = " ")
  rhs_compact <- gsub("\\s+", "", rhs_text)

  if (rhs_compact %in% c("0", "-1")) {
    return(character())
  }
  if (rhs_compact == ".") {
    return(available)
  }

  terms_i <- attr(
    stats::terms(stats::as.formula(paste("~", rhs_text))), "term.labels"
  )
  if ("." %in% terms_i) {
    stop(
      "In `coef_breaks`, `.` is only supported by itself, as in `y1 ~ .`.",
      call. = FALSE
    )
  }
  rhs_terms <- trimws(strsplit(rhs_text, "\\+", fixed = FALSE)[[1]])
  if (any(rhs_terms == "1")) {
    terms_i <- c("(Intercept)", terms_i)
  }
  terms_i <- unique(terms_i)
  if (length(terms_i) == 0L) {
    return(character())
  }

  terms_i <- ifelse(tolower(terms_i) == "intercept", "(Intercept)", terms_i)
  missing <- setdiff(terms_i, available)
  if (length(missing) > 0L) {
    stop(
      sprintf(
        "Coefficient(s) %s in `coef_breaks` are not available for equation `%s`. Available coefficients: %s.",
        paste(missing, collapse = ", "),
        equation,
        paste(available, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  terms_i
}

#' Builds one row per coefficient in the stacked system
#'
#' @noRd
coefmap <- function(y_names, eqz) {
  if (length(y_names) != length(eqz)) {
    stop("`eqz` must have one element per equation.")
  }

  eq_index <- rep(seq_along(eqz), lengths(eqz))
  term <- unlist(eqz, use.names = FALSE)
  equation <- y_names[eq_index]
  data.frame(
    index = seq_along(term),
    eq_index = eq_index,
    equation = equation,
    term = term,
    coefficient = paste0(equation, ":", term),
    stringsAsFactors = FALSE
  )
}

#' Build Coefficient Names by returning the last column of coefmap
#'
#' @noRd
eqnames <- function(y_names, eqz) {
  coefmap(y_names, eqz)$coefficient
}

#' Label Coefficient-Break Selectors
#'
#' @noRd
breaks_lab <- function(coef_breaks) {
  if (is.character(coef_breaks) && length(coef_breaks) == 1L) {
    return(coef_breaks)
  }
  if (is.list(coef_breaks)) {
    return(vapply(coef_breaks, flabel, character(1)))
  }
  NA_character_
}

#' Build SUR Matrices
#'
#' @noRd
sur_mats <- function(formulas, data) {
  formulas <- valfor(formulas)
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame when `model` is a list of equation formulas.")
  }

  neq <- length(formulas)
  frames <- vector("list", neq)
  x_list <- vector("list", neq)
  y_list <- vector("list", neq)
  y_names <- character(neq)

  for (i in seq_along(formulas)) {
    frames[[i]] <- stats::model.frame(
      formulas[[i]],
      data = data,
      na.action = stats::na.pass
    )
    terms_i <- stats::terms(formulas[[i]], data = data)
    y_i <- stats::model.response(frames[[i]])
    if (is.null(y_i) || !is.numeric(y_i)) {
      stop("Each SUR equation must have a numeric response.")
    }

    x_i <- stats::model.matrix(
      stats::delete.response(terms_i), data = frames[[i]]
    )
    y_list[[i]] <- as.numeric(y_i)
    x_list[[i]] <- x_i
    y_names[[i]] <- rname(formulas[[i]])
  }

  keep <- Reduce(
    `&`,
    c(
      lapply(y_list, function(x) stats::complete.cases(x)),
      lapply(x_list, stats::complete.cases)
    )
  )
  if (!any(keep)) {
    stop("No complete observations remain after evaluating the SUR formulas.")
  }
  if (!all(keep)) {
    warning(
      sprintf(
        "%d observations dropped for missing values. Break dates index the retained rows; map with `spec$y_rows`.",
        sum(!keep)
      ),
      call. = FALSE
    )
  }

  y <- do.call(cbind, lapply(y_list, function(x) x[keep]))
  colnames(y) <- make.unique(y_names)
  x_list <- lapply(x_list, function(x) x[keep, , drop = FALSE])

  z <- unionz(x_list)
  S <- surS(x_list, union_names = colnames(z))
  eqz <- lapply(x_list, colnames)
  names(eqz) <- colnames(y)

  list(
    y = y,
    z = z,
    S = S,
    formulas = vapply(formulas, flabel, character(1)),
    eqz = eqz,
    y_rows = which(keep)
  )
}

#' Build VAR Matrices
#'
#' @noRd
var_mats <- function(y, lags, intercept = TRUE) {
  y <- asmat(y, "y")
  if (!all(is.finite(y))) {
    stop("`y` must not contain missing or non-finite values.", call. = FALSE)
  }
  lags <- as.integer(lags)
  if (length(lags) != 1 || is.na(lags) || lags < 1) {
    stop("`lags` must be a positive integer.")
  }
  if (nrow(y) <= lags) {
    stop("`lags` must be smaller than the number of observations in `y`.")
  }
  if (is.null(colnames(y))) {
    colnames(y) <- paste0("y", seq_len(ncol(y)))
  }

  y_rows <- (lags + 1):nrow(y)
  y_response <- y[y_rows, , drop = FALSE]
  z_parts <- list()
  z_names <- character()

  if (intercept) {
    z_parts[[length(z_parts) + 1]] <- matrix(1, nrow = length(y_rows), ncol = 1)
    z_names <- c(z_names, "(Intercept)")
  }

  for (lag in seq_len(lags)) {
    lagged <- y[(lags + 1 - lag):(nrow(y) - lag), , drop = FALSE]
    z_parts[[length(z_parts) + 1]] <- lagged
    z_names <- c(z_names, paste0("L", lag, ".", colnames(y)))
  }

  z <- do.call(cbind, z_parts)
  colnames(z) <- z_names
  list(y = y_response, z = z, y_rows = y_rows)
}

#' Build Network VAR Matrices
#'
#' @noRd
nvar_mats <- function(y, W, lags, intercept = TRUE, row_normalize = FALSE) {
  y <- asmat(y, "y")
  if (is.null(colnames(y))) {
    colnames(y) <- paste0("y", seq_len(ncol(y)))
  }
  W <- valW(W, n = ncol(y), row_normalize = row_normalize)
  dimnames(W) <- list(colnames(y), colnames(y))

  modelparts <- var_mats(y = y, lags = lags, intercept = intercept)
  S_info <- nvarS(
    W = W,
    lags = lags,
    intercept = intercept,
    y_names = colnames(modelparts$y),
    z_names = colnames(modelparts$z)
  )

  list(
    y = modelparts$y,
    z = modelparts$z,
    y_rows = modelparts$y_rows,
    W = W,
    S = S_info$S,
    param_names = S_info$param_names,
    param_groups = S_info$param_groups
  )
}

#' Validate and Optionally Row-Normalize Network Matrix
#'
#' @noRd
valW <- function(W, n, row_normalize = FALSE) {
  W <- asmat(W, "W", allow_vector = FALSE)
  if (nrow(W) != n || ncol(W) != n) {
    stop("`W` must be a square matrix with dimension equal to `ncol(y)`.",
         call. = FALSE)
  }
  if (any(!is.finite(W))) {
    stop("`W` must contain only finite numeric values.", call. = FALSE)
  }
  if (isTRUE(row_normalize)) {
    rs <- rowSums(W)
    nz <- abs(rs) > .Machine$double.eps
    W[nz, ] <- W[nz, , drop = FALSE] / rs[nz]
  }
  W
}

#' Build the S(W) Matrix for a Network VAR
#'
#' @noRd
nvarS <- function(W, lags, intercept, y_names, z_names) {
  n <- nrow(W)
  pz <- length(z_names)
  lags <- as.integer(lags)
  nint <- if (isTRUE(intercept)) n else 0L
  network_names <- paste0("network_lag", seq_len(lags))
  own_names <- paste0("own_lag", seq_len(lags))
  intercept_names <- if (nint > 0L) paste0("mu:", y_names) else character()
  param_names <- c(intercept_names, network_names, own_names)

  S <- matrix(0, nrow = n * pz, ncol = length(param_names))
  rownames(S) <- eqnames(y_names, rep(list(z_names), n))
  colnames(S) <- param_names

  if (isTRUE(intercept)) {
    for (eq in seq_len(n)) {
      raw_row <- (eq - 1) * pz + 1
      S[raw_row, eq] <- 1
    }
  }

  network_offset <- nint
  own_offset <- nint + lags
  z_offset <- if (isTRUE(intercept)) 1L else 0L
  for (lag in seq_len(lags)) {
    network_col <- network_offset + lag
    own_col <- own_offset + lag
    for (eq in seq_len(n)) {
      for (var in seq_len(n)) {
        z_col <- z_offset + (lag - 1) * n + var
        raw_row <- (eq - 1) * pz + z_col
        S[raw_row, network_col] <- W[eq, var]
        S[raw_row, own_col] <- as.numeric(eq == var)
      }
    }
  }

  groups <- list(
    all = param_names,
    dynamic = c(network_names, own_names),
    network = network_names,
    own = own_names,
    intercepts = intercept_names,
    intercept = intercept_names,
    mu = intercept_names
  )
  for (lag in seq_len(lags)) {
    groups[[paste0("lag", lag)]] <- c(network_names[[lag]], own_names[[lag]])
    groups[[paste0("alpha", lag)]] <- network_names[[lag]]
    groups[[paste0("beta", lag)]] <- own_names[[lag]]
  }

  list(S = S, param_names = param_names, param_groups = groups)
}

#' Build Restriction Matrices for Network VAR
#'
#' The Network VAR counterpart of `coef_restrict()`.  Which parameters change
#' is regime-invariant, so the selection is made once only.
#'
#' @noRd
nvar_restrict <- function(coef_breaks, m, param_names, groups) {
  m <- as.integer(m)
  if (length(m) != 1 || is.na(m) || m < 0) {
    stop("`m` must be a non-negative integer.", call. = FALSE)
  }

  selected <- nvar_breaks_sel(coef_breaks, param_names, groups)
  q <- length(param_names)
  nchg <- sum(selected)
  info <- list(
    selector = coef_breaks,
    sumq = nchg,
    selected = param_names[selected]
  )

  if (nchg == q) {
    return(c(list(R = NULL, R_path = NULL, brbeta = TRUE), info))
  }
  if (nchg == 0L) {
    return(c(list(R = NULL, R_path = NULL, brbeta = FALSE), info))
  }

  R_path <- lapply(seq_len(m), function(k) buildR(selected, param_names, k))
  c(list(R = R_path[[m]], R_path = R_path, brbeta = TRUE), info)
}

#' Select Network VAR Parameters Allowed to Change
#'
#' @noRd
nvar_breaks_sel <- function(coef_breaks, param_names, groups) {
  if (is.null(coef_breaks)) {
    coef_breaks <- "none"
  }
  if (!is.character(coef_breaks) || length(coef_breaks) < 1L) {
    stop("`coef_breaks` must be \"all\", \"none\", or a character vector.",
         call. = FALSE)
  }

  keys <- tolower(coef_breaks)
  if (length(keys) == 1L && keys == "all") {
    return(rep(TRUE, length(param_names)))
  }
  if (length(keys) == 1L && keys == "none") {
    return(rep(FALSE, length(param_names)))
  }

  sel_names <- character()
  param_lower <- tolower(param_names)
  group_names <- names(groups)
  group_lower <- tolower(group_names)
  for (key in keys) {
    group_i <- match(key, group_lower)
    if (!is.na(group_i)) {
      sel_names <- c(sel_names, groups[[group_names[[group_i]]]])
      next
    }

    param_i <- match(key, param_lower)
    if (!is.na(param_i)) {
      sel_names <- c(sel_names, param_names[[param_i]])
      next
    }

    stop(
      sprintf(
        "`coef_breaks` entry `%s` is not recognized. Available parameters are: %s.",
        key,
        paste(param_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  param_names %in% unique(sel_names)
}

#' Build Common-Regressor Matrices from a Formula
#'
#' @noRd
mb_mats <- function(formula, data) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must be a two-sided formula.")
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame when using formula input.")
  }

  frame <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  terms_i <- stats::terms(formula, data = data)
  y <- stats::model.response(frame)
  y <- asmat(y, "formula response")
  if (is.null(colnames(y))) {
    colnames(y) <- rnames(formula, ncol(y))
  }

  z <- stats::model.matrix(stats::delete.response(terms_i), data = frame)
  keep <- stats::complete.cases(y) & stats::complete.cases(z)
  if (!any(keep)) {
    stop("No complete observations remain after evaluating the formula.")
  }
  if (!all(keep)) {
    warning(
      sprintf(
        "%d observations dropped for missing values. Break dates index the retained rows; map with `spec$y_rows`.",
        sum(!keep)
      ),
      call. = FALSE
    )
  }

  list(
    y = y[keep, , drop = FALSE],
    z = z[keep, , drop = FALSE],
    formula = flabel(formula),
    y_rows = which(keep)
  )
}

#' Validate Response and Regressor Matrices
#'
#' @noRd
valyz <- function(y, z) {
  if (!is.numeric(y) || !is.numeric(z)) {
    stop("`y` and `z` must be numeric.")
  }
  if (nrow(y) != nrow(z)) {
    stop("`y` and `z` must have the same number of rows.")
  }
  invisible(NULL)
}

#' Validate SUR Formulas
#'
#' @noRd
valfor <- function(formulas) {
  if (!is.list(formulas)) {
    stop("`formulas` must be a list of formulas.")
  }
  if (length(formulas) < 1) {
    stop("`formulas` must contain at least one equation.")
  }
  if (!all(vapply(formulas, inherits, logical(1), what = "formula"))) {
    stop("Every element of `formulas` must be a formula.")
  }
  if (!all(vapply(formulas, length, integer(1)) == 3L)) {
    stop("Every SUR formula must include a response.")
  }
  formulas
}

#' Build Union Regressor Matrix
#'
#' @noRd
unionz <- function(x_list) {
  union_names <- unique(unlist(lapply(x_list, colnames), use.names = FALSE))
  z <- matrix(0, nrow = nrow(x_list[[1]]), ncol = length(union_names))
  colnames(z) <- union_names

  for (name in union_names) {
    values <- NULL
    for (x in x_list) {
      idx <- match(name, colnames(x))
      if (!is.na(idx)) {
        if (is.null(values)) {
          values <- x[, idx]
        } else if (
          !isTRUE(all.equal(values, x[, idx], check.attributes = FALSE))
        ) {
          stop(
            sprintf(
              "Regressor column `%s` has conflicting values across equations.",
              name
            ),
            call. = FALSE
          )
        }
      }
    }
    z[, name] <- values
  }

  z
}

#' Build SUR Selection Matrix
#'
#' @noRd
surS <- function(x_list, union_names) {
  n <- length(x_list)
  p <- length(union_names)
  q_total <- sum(vapply(x_list, ncol, integer(1)))
  S <- matrix(0, nrow = n * p, ncol = q_total)
  rownames(S) <- paste0(
    "eq", rep(seq_len(n), each = p),
    ":", rep(union_names, times = n)
  )

  sel_names <- character(q_total)
  col_offset <- 0
  for (eq in seq_len(n)) {
    x_names <- colnames(x_list[[eq]])
    for (j in seq_along(x_names)) {
      raw_col <- (eq - 1) * p + match(x_names[[j]], union_names)
      sel_col <- col_offset + j
      S[raw_col, sel_col] <- 1
      sel_names[[sel_col]] <- paste0("eq", eq, ":", x_names[[j]])
    }
    col_offset <- col_offset + length(x_names)
  }
  colnames(S) <- sel_names
  S
}

#' Formula Label
#'
#' @noRd
flabel <- function(formula) {
  paste(deparse(formula), collapse = " ")
}

#' Formula Response Name
#'
#' @noRd
rname <- function(formula) {
  paste(deparse(formula[[2]]), collapse = " ")
}

#' Formula Response Names
#'
#' @noRd
rnames <- function(formula, n) {
  response <- formula[[2]]
  if (is.call(response) && identical(response[[1]], as.name("cbind"))) {
    names <- vapply(as.list(response[-1]), function(x) {
      paste(deparse(x), collapse = " ")
    }, character(1))
    return(make.unique(names))
  }
  if (n == 1L) {
    return(rname(formula))
  }
  paste0(rname(formula), seq_len(n))
}
