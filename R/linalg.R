#' Symmetric Matrix Square Root
#'
#' Computes the symmetric square root of a positive semidefinite matrix after
#' truncating small negative eigenvalues caused by numerical rounding.
#'
#' @noRd
sqrm <- function(M, tol = sqrt(.Machine$double.eps)) {
  M <- (M + t(M)) / 2

  eig <- eigen(M, symmetric = TRUE)
  values <- eig$values

  if (any(values < -tol)) {
    warning("`M` has negative eigenvalues; setting them to zero in `sqrm()`.")
  }
  values[values < 0] <- 0

  eig$vectors %*% diag(sqrt(values), nrow = length(values)) %*% t(eig$vectors)
}

#' Moore-Penrose Pseudoinverse
#'
#' Computes a generalized inverse using `MASS::ginv()` when available, with an
#' SVD if not available.
#'
#' @noRd
# Tolerance changed 2026-10-11 from sqrt(.Machine$double.eps).  The relative
# cutoff sqrt(eps) discarded singular values that are small only because the
# variables are measured in different units.
pinv <- function(x, tol = 1e-14) {
  if (requireNamespace("MASS", quietly = TRUE)) {
    return(MASS::ginv(x, tol = tol))
  }

  sv <- svd(x)
  if (length(sv$d) == 0) {
    return(matrix(0, nrow = ncol(x), ncol = nrow(x)))
  }

  keep <- sv$d > tol * max(sv$d)
  if (!any(keep)) {
    return(matrix(0, nrow = ncol(x), ncol = nrow(x)))
  }

  sv$v[, keep, drop = FALSE] %*%
    (diag(1 / sv$d[keep], nrow = sum(keep)) %*% t(sv$u[, keep, drop = FALSE]))
}

#' Log Determinant
#'
#' Computes the log determinant of a positive definite matrix and stops if the
#' determinant is not positive.
#'
#' @noRd
logdet <- function(x) {
  det_x <- determinant(x, logarithm = TRUE)

  if (det_x$sign <= 0) {
    stop("A covariance matrix has a non-positive determinant.")
  }

  as.numeric(det_x$modulus)
}

#' Safe Log Determinant
#'
#' Returns log|x| if det(x) > 0 and Inf otherwise. Used in the break-date
#' searches, which minimise sum_j (T_j / 2) log|Sigma_j|. Inf is used to exclude a
#' partition whose regime covariance is singular, which would otherwise give
#' log|x| = -Inf and be selected.
#'
#' @noRd
slogdet <- function(x) {
  det_x <- determinant(x, logarithm = TRUE)

  if (det_x$sign <= 0 || !is.finite(det_x$modulus)) {
    return(Inf)
  }

  as.numeric(det_x$modulus)
}

#' Iterate a Function to a Fixed Point
#'
#' Repeatedly applies `update()` until successive iterates differ by at most
#' `tol` in max-abs norm, or `max_iter` iterations have been reached.  
#'
#' @noRd
fixed_point <- function(x, update, tol, max_iter, label) {
  for (itr in seq_len(max_iter)) {
    prev <- x
    x <- update(prev)
    if (max(abs(x - prev)) <= tol) {
      return(list(value = x, iterations = itr, converged = TRUE))
    }
  }

  warning(
    sprintf(
      "%s did not converge within max_iter = %d (tol = %g).",
      label, as.integer(max_iter), tol
    ),
    call. = FALSE
  )

  list(value = x, iterations = max_iter, converged = FALSE)
}
