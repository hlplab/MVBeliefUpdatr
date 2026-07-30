# Internal utility helpers shared across migration scaffolds.

#' @keywords internal
.mvbu_is_s7_class <- function(x, class_name) {
  if (!inherits(x, "S7_object")) {
    return(FALSE)
  }

  class_names <- class(x)
  class_names <- gsub(".*::", "", class_names)
  any(class_names == class_name)
}

#' @keywords internal
.is_numeric_vector <- function(x) {
  is.numeric(x) && is.null(dim(x))
}

is_scalar_character <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x[1])
}

is_character <- function(x) {
  is.character(x)
}

#' @keywords internal
.as_numeric_vector <- function(x) {
  if (!.is_numeric_vector(x)) {
    stop("Expected a numeric vector.", call. = FALSE)
  }
  as.numeric(x)
}
