#' @include internal-globals.R
#' @include internal-asserts.R
NULL

# Internal utility helpers shared across migration scaffolds.

# dim that returns length of vector for vector
.dim <- function(x) {
  if (is.null(dim(x))) return(length(x))
  return(dim(x))
}

#' Project matrix onto positive definite cone
#' @param sigma Square numeric covariance matrix.
#' @param tol Small positive number for eigenvalue clamping
#'   (default: MVBU_PROB_TOL).
#' @return Symmetric positive definite matrix.
#' @noRd
.make_pos_def <- function(sigma, tol = MVBU_PROB_TOL) {
  sigma <- (sigma + t(sigma)) / 2
  eig <- eigen(sigma, symmetric = TRUE)
  vals <- pmax(eig$values, tol)
  eig$vectors %*% diag(vals, nrow = length(vals)) %*% t(eig$vectors)
}
