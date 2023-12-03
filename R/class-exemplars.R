get_expected_columns_for_exemplars <- function()
  c("category", "exemplars", "sim_function")

#' Is this a set of exemplar categories?
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
#' @keywords TBD
#' @export
is.exemplars <- function(x, group = NULL, verbose = F) {
  name_of_x <- deparse(substitute(x))

  if (!is_tibble(x)) {
    if (verbose) message("Object is not a tibble.")
    return(FALSE)
  }

  if (!is.null(group)) {
    if (verbose) message("Checking whether ", name_of_x, " is a collection of exemplars within each unique combination of group values.")
    x %<>% group_by(!!! syms(group))
  }

  if (any(get_expected_columns_for_exemplars() %nin% names(x))) {
    if (verbose) message("x is missing a required column: ", paste(get_expected_columns_for_exemplars, collapse = ","))
    return(FALSE)
  }

  # Check that category is a factor only after everything else is checked.
  if (!is.factor(x$category)) return(FALSE)

  return(TRUE)
}


