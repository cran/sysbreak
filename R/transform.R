#' Convert User Input to a Numeric Matrix
#'
#' Converts numeric vectors, matrices, and data frames to numeric matrices while
#' preserving column names where possible.
#'
#' @noRd
asmat <- function(x, arg, allow_vector = TRUE) {
  if (is.null(x)) {
    stop(sprintf("`%s` must not be NULL.", arg), call. = FALSE)
  }

  if (is.data.frame(x)) {
    numeric_cols <- vapply(x, is.numeric, logical(1))
    if (!all(numeric_cols)) {
      stop(
        sprintf("`%s` must contain only numeric columns.", arg),
        call. = FALSE
      )
    }
    out <- as.matrix(x)
  } else if (is.matrix(x)) {
    out <- x
  } else if (allow_vector && is.numeric(x)) {
    out <- matrix(x, ncol = 1)
  } else {
    stop(
      sprintf("`%s` must be a numeric matrix, vector, or data frame.", arg),
      call. = FALSE
    )
  }

  if (!is.numeric(out)) {
    stop(sprintf("`%s` must be numeric.", arg), call. = FALSE)
  }

  out
}

#' Converts responses and regressors into the equation
#' system used by the mbreak estimation and test routines.
#'
#' @noRd
transf <- function(y, z, S) {
  bigt <- nrow(y)
  n <- ncol(y)
  maty <- matrix(as.vector(t(y)), ncol = 1)
  matz <- matrix(0, nrow = bigt * n, ncol = ncol(S))
  eye_n <- diag(n)

  for (i in seq_len(bigt)) {
    rows_i <- ((i - 1) * n + 1):(i * n)
    matz[rows_i, ] <- kronecker(eye_n, z[i, , drop = FALSE]) %*% S
  }

  list(maty = maty, matz = matz)
}

#' Expands a regressors matrix into regime-specific blocks for a given set
#' of break dates.
#'
#' @noRd
pzbar <- function(matz, m, br, bigt = NULL) {
  if (m == 0) {
    return(matz)
  }

  n <- nrow(matz) / bigt
  q <- ncol(matz)
  endpoints <- c(0, br, bigt)

  matzb <- matrix(0, nrow = nrow(matz), ncol = (m + 1) * q)
  for (seg in seq_len(m + 1)) {
    t_start <- endpoints[[seg]] + 1
    t_end <- endpoints[[seg + 1]]
    row_idx <- obs_rows(t_start, t_end, n)
    col_idx <- ((seg - 1) * q + 1):(seg * q)
    matzb[row_idx, col_idx] <- matz[row_idx, , drop = FALSE]
  }

  matzb
}
