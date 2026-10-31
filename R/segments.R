# Index helpers for calculations.

#' Row Indices
#'
#' Returns the row indices for observations `from` to `to` in a system with
#' `n` equations.
#'
#' @noRd
obs_rows <- function(from, to, n) {
  ((from - 1) * n + 1):(to * n)
}

#' Upper-Triangular Segment Index
#'
#' Maps a segment `(start, end)` into the compact vector index used to store
#' all possible segment objective values.
#'
#' @noRd
bigvec_index <- function(start, end, bigt) {
  (start - 1) * bigt - (start - 1) * (start - 2) / 2 + (end - start + 1)
}
