#' Structural Breaks in a System of Equations
#'
#' Estimation and testing for multivariate structural break models with
#' coefficient breaks, covariance breaks, or both.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{sysbreak}}: Breaks in a linear system of equations,
#'     including common-regressor and SUR-style structural break models
#'   \item \code{\link{varbreak}}: VAR structural break models
#'   \item \code{\link{netbreak}}: Network VAR structural break models
#'   \item \code{\link{ordbreak}}: Ordered two-break models
#'   \item \code{\link{ordcheck}}: Coefficient layout of an ordered two-break model
#'   \item \code{\link{brtype}}: Tests of what changes at each break
#'   \item \code{\link{dyspill}}: Generalized Diebold-Yilmaz spillover measures
#'   \item \code{\link{dynobreak}}: The same spillover measures without breaks
#' }
#'
#' @section Break Types:
#' The package supports the following structural break specifications:
#' \itemize{
#'   \item Coefficient breaks
#'   \item Covariance breaks
#'   \item Joint coefficient and covariance breaks
#' }
#'
#' @references
#' Qu, Z., and Perron, P. (2007). Estimating and testing structural changes
#' in multivariate regressions. \emph{Econometrica}, 75(2), 459-502.
#'
#' Bai, J., and Perron, P. (2003). Critical values for multiple structural
#' change tests. \emph{The Econometrics Journal}, 6(1), 72-78.
#'
#' @keywords internal
"_PACKAGE"
