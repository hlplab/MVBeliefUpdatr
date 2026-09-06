# Scalar type predicates ---------------------------------------------------

.is_scalar <- function(x) {
  length(x) == 1L && !is.na(x)
}

.is_scalar_numeric <- function(x) {
  is.numeric(x) && .is_scalar(x)
}

.is_scalar_integer <- function(x) {
  is.integer(x) && .is_scalar(x)
}

.is_scalar_double <- function(x) {
  is.double(x) && .is_scalar(x)
}

.is_scalar_character <- function(x) {
  is.character(x) && .is_scalar(x)
}

.is_scalar_factor <- function(x) {
  is.factor(x) && .is_scalar(x)
}

.is_scalar_logical <- function(x) {
  is.logical(x) && .is_scalar(x)
}

.is_scalar_count <- function(x) {
  .is_scalar_numeric(x) && x >= 0 && x == floor(x)
}

.is_non_empty_scalar_character <- function(x) {
  .is_scalar_character(x) && nzchar(x)
}

# Compound predicates ------------------------------------------------------

.is_try_error <- function(x) {
  inherits(x, "try-error")
}

.is_equal <- function(x, y, check.attributes = FALSE, ...) {
  isTRUE(all.equal(x, y, check.attributes = check.attributes, ...))
}

.is_like_factor <- function(x) {
  is.factor(x) || is.character(x) || is.logical(x)
}

.is_sigma <- function(x) {
  .assert_true(!is.null(x), msg = "Expected a covariance matrix, but got NULL.")
  if (is.matrix(x)) {
    if (all(x == 0) || is.positive.definite(x)) return(TRUE) else return(FALSE)
  } else {
    if (.is_scalar_double(x)) return(TRUE) else return(FALSE)
  }
}