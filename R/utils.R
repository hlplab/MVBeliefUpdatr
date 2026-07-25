# Internal utility helpers shared across migration scaffolds.

#' @keywords internal
.is_numeric_vector <- function(x) {
  is.numeric(x) && is.null(dim(x))
}

#' @keywords internal
.as_numeric_vector <- function(x) {
  if (!.is_numeric_vector(x)) {
    stop("Expected a numeric vector.", call. = FALSE)
  }
  as.numeric(x)
}
