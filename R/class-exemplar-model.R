get_expected_columns_for_exemplar_model <- function() append(get_expected_columns_for_exemplars(), get_expected_columns_for_model())

#' Is this an exemplar model?
#'
#' Check whether \code{x} is an exemplar model. Optionally, one can also check whether a lapse rate
#' and lapse bias is part of the model.
#'
#' @param x Object to be checked.
#' @param group Name of one or more group variables, each unique combination of which describes an exemplar model. (default: NULL)
#' @param category Name of the category variable. (default: "category")
#' @param verbose Should verbose output be provided? (default: `TRUE`)
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#' @export
is.exemplar_model <- function(x, group = NULL, verbose = F, tolerance = 1e-5) {
  name_of_x <- deparse(substitute(x))

  # When no groups are specified, infer groups from object.
  if (is.null(group)) {
    group <- setdiff(names(x), get_expected_columns_for_exemplar_model())
    if (length(group) == 0) group <- NULL else {
      if (verbose) message(name_of_x, "has additional columns beyond those expected:", paste(group, collapse = ", "), "Interpreting those columns as group variables.")
    }
  }

  if (!is.null(group)) {
    if (verbose) message("Checking whether ", name_of_x, " is an exemplar model within each unique combination of group values.")
    x %<>% group_by(!!! syms(group))
  }

  if (!is.exemplars(x, group = group, verbose = verbose)) {
    if (verbose) message("x does not contain exemplars.")
    return(FALSE)
  }

  # Only need to test for exemplars columns here since is.model is called below.
  if (any(get_expected_columns_for_exemplars() %nin% names(x))) {
    if (verbose) message("x is missing a required column: ", paste(get_expected_columns_for_exemplars, collapse = ","))
    return(FALSE)
  }

  if (any(!is.factor(x$category))) {
    if (verbose) message(paste("category must be a factor."))
    return(FALSE)
  }

  if (!is.model(x, group = group, verbose = verbose, tolerance = tolerance)) {
    return(FALSE)
  }

  return(TRUE)
}


