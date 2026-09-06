get_expected_columns_for_exemplars <- function()
  c("category", "exemplars", "sim_function")

#' Deprecated: is.exemplars
#'
#' @description `r lifecycle::badge("deprecated")`
#' `is.exemplars()` was deprecated in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0.
#' Please use S7 validators and [Exemplar_Model] instead.
#'
#' @param x Object to be checked.
#' @param group Name of one or more group variables, each unique combination of which describes a set of exemplars. (default: NULL)
#' @param verbose Logical. If `TRUE`, emits diagnostics.
#'
#' @return A logical.
#'
#' @seealso [Exemplar_Model]
#' @keywords internal
#' @export
is.exemplars <- function(x, group = NULL, verbose = F) {
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "is.exemplars()",
    details = "Use S7 validators and Exemplar_Model instead."
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


