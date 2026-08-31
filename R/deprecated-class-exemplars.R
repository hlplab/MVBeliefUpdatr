get_expected_columns_for_exemplars <- function()
  c("category", "exemplars", "sim_function")

#' Deprecated: is.exemplars
#'
#' Check whether \code{x} is a set of exemplar categories.
#'
#' @param x Object to be checked.
#' @param group Name of one or more group variables, each unique combination of which describes a set of exemplars. (default: NULL)
#' @param category Name of the category variable. (default: "category")
#'
#' @return A logical.
#'
#' @seealso TBD
#' @description Deprecated. Use the S7-based validators and constructors for exemplar categories instead.
#' @keywords internal
#' @export
is.exemplars <- function(x, group = NULL, verbose = F) {
  lifecycle::deprecate_warn(
    when = "0.0.3",
    what = "is.exemplars()",
    details = "the S7-based validators and constructors for exemplar categories"
  )
  name_of_x <- deparse(substitute(x))

  if (!is.data.frame(x)) {
    if (verbose) message("Object is not a data frame-like object.")
    return(FALSE)
  }

  if (!is.null(group)) {
    if (verbose) message("Checking whether ", name_of_x, " is a collection of exemplars within each unique combination of group values.")
    # Grouping is not needed for these structural checks.
  }

  if (any(get_expected_columns_for_exemplars() %nin% names(x))) {
    if (verbose) message("x is missing a required column: ", paste(get_expected_columns_for_exemplars, collapse = ","))
    return(FALSE)
  }

  # Check that category is a factor only after everything else is checked.
  if (!is.factor(x$category)) return(FALSE)

  return(TRUE)
}


