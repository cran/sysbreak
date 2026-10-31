#' GAUSS/MATLAB Example Data 1
#'
#' Example data used by the original GAUSS/MATLAB `mbreak2` script. In this
#' package the corresponding procedure is \code{\link{sysbreak}}.
#'
#' @format A numeric matrix with 150 rows and 5 columns. Columns 1 and 2 are
#' dependent variables, and columns 3 through 5 are regressors.
#'
#' @source Qu, Z., and Perron, P. (2007). Estimating and testing structural
#' changes in multivariate regressions. \emph{Econometrica}, 75(2), 459-502.
"data1"

#' GAUSS/MATLAB Example Data 2
#'
#' Example data used by the original GAUSS/MATLAB `obreak2` script for the
#' ordered-break procedure.
#'
#' @format A numeric matrix with 100 rows and 5 columns. Columns 1 and 2 are
#' dependent variables, and columns 3 through 5 are regressors.
#'
#' @source Qu, Z., and Perron, P. (2007). Estimating and testing structural
#' changes in multivariate regressions. \emph{Econometrica}, 75(2), 459-502.
"data2"

#' U.S. House Prices and the Macroeconomy
#'
#' Quarterly data for the house price application of Ahamada and Diaz Sanchez
#' (2013), 1959:Q4 to 2009:Q3.
#'
#' @format A data frame with 200 rows and 6 columns:
#' \describe{
#'   \item{date}{Quarter, such as \code{"1960:Q1"}.}
#'   \item{dGDP}{Log difference of real GDP.}
#'   \item{dC}{Log difference of real personal consumption expenditures.}
#'   \item{dHP}{Log difference of the real home price index.}
#'   \item{FFR}{Federal funds rate, in percent, first month of the quarter.}
#'   \item{dI}{Change in \code{FFR}.}
#' }
#'
#' @source FRED series GDPC1, PCEC, CPIAUCSL and FEDFUNDS, and the Shiller
#' nominal home price index, seasonally adjusted by log-STL. Consumption and
#' home prices are deflated by CPIAUCSL.
#'
#' @references Ahamada, I., and Diaz Sanchez, J. L. (2013). A retrospective
#' analysis of the house price macro-relationship in the United States.
#' \emph{International Journal of Central Banking}, 9(4), 153-174.
"housing"

#' International Inflation
#'
#' Monthly consumer price inflation for the G7 countries and the euro area,
#' 1973:04 to 2007:12, from Bataa, Osborn, Sensier and van Dijk (2013).
#'
#' @format A data frame with 417 rows and 9 columns: \code{date}, such as
#'   \code{"1973:04"}, followed by inflation for \code{Canada},
#'   \code{EuroArea}, \code{France}, \code{Germany}, \code{Italy},
#'   \code{Japan}, \code{UK} and \code{US}. Inflation is
#'   \eqn{100 \times \Delta \log} of the consumer price index, adjusted for
#'   outliers and seasonality by the authors.
#'
#' @source The replication files of Bataa et al. (2013).
#'
#' @references Bataa, E., Osborn, D. R., Sensier, M., and van Dijk, D. (2013).
#' Structural breaks in the international dynamics of inflation.
#' \emph{Review of Economics and Statistics}, 95(2), 646-659.
"inflation"

#' Simulated Network VAR
#'
#' One draw from Scenario 1 of Han and Lee (2024): a network VAR for 20 nodes
#' over 300 periods, with breaks after periods 100 and 200.
#'
#' @format A list with two components:
#' \describe{
#'   \item{y}{A 300 by 20 matrix of observations, with columns \code{node1} to
#'     \code{node20}.}
#'   \item{W}{The 20 by 20 row-normalized network matrix. Entry \eqn{(i, l)}
#'     is \eqn{1/n_i} if node \eqn{i} follows node \eqn{l} and 0 otherwise,
#'     where \eqn{n_i} is the number of nodes that \eqn{i} follows.}
#' }
#'
#' @details In regime \eqn{j},
#' \deqn{y_{it} = V_i'\gamma_j + \sum_{\ell} \alpha_{\ell,j} (W y_{t-\ell})_i +
#'   \sum_{\ell} \beta_{\ell,j} y_{i,t-\ell} + e_{it}, \qquad
#'   e_{it} \sim N(0, 0.1^2).}
#' Regime 1 has one lag, with \eqn{\alpha = -0.10} and \eqn{\beta = 0.20}.
#' Regimes 2 and 3 have two lags, with \eqn{\alpha = (0.20, -0.22)},
#' \eqn{\beta = (-0.12, 0.40)} and \eqn{\alpha = (-0.12, 0.10)},
#' \eqn{\beta = (0.25, -0.40)}. The covariates \eqn{V_i} are fixed, so
#' \eqn{V_i'\gamma_j} is a node intercept that changes across regimes. The
#' in-degrees of the network follow a power law, as in Zhu et al. (2017).
#'
#' @source Simulated with \code{set.seed(1)} following Han and Lee (2024).
#'
#' @references Han, Y., and Lee, T. C. M. (2024). Structural break detection
#' in non-stationary network vector autoregression models. \emph{IEEE
#' Transactions on Network Science and Engineering}, 11(5), 4134-4145.
#'
#' Zhu, X., Pan, R., Li, G., Liu, Y., and Wang, H. (2017). Network vector
#' autoregression. \emph{The Annals of Statistics}, 45(3), 1096-1123.
"netsim"
