#' @include internal-asserts.R
NULL

# Internal utility helpers shared across migration scaffolds.

# dim that returns length of vector for vector
.dim <- function(x) {
  if (is.null(dim(x))) return(length(x))
  return(dim(x))
}

.replace_na_in_array <- function(x, fill = 0) {
  .assert_true(is.array(x), msg = "x must be an array.")
  .assert_true(.is_scalar(fill), msg = "fill must be a scalar value.")

  x[is.na(x)] <- fill
  return(x)
}
