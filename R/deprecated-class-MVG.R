get_expected_columns_for_MVG <- function() c("category", "mu", "Sigma")

#' Deprecated: is.MVG
#'
#' @description `r lifecycle::badge("deprecated")`
#' `is.MVG()` was deprecated in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0.
#' Please use S7 validators and [MVG_IdealObserver] instead.
#'
#' @param x Object to be checked.
#' @param group Name of one or more group variables, each unique combination of which describes an MVG. (default: NULL)
#' @param category DEPRECATED Name of the category variable. (default: "category")
#' @param is.long Is this check assessing whether the ideal observer is in long format (`TRUE`) or wide format (`FALSE`)?
#' (default: `TRUE`)
#' @param verbose Logical. If `TRUE`, emits diagnostics.
#'
#' @return A logical.
#'
#' @seealso [MVG_IdealObserver]
#' @keywords internal
#' @export
is.MVG <- function(x, group = NULL, category = "category", is.long = T, verbose = F) {
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "is.MVG()",
    details = "Use S7 validators and MVG_IdealObserver instead."
  )
  name_of_x <- deparse(substitute(x))
  .assert_that(.is_non_NA_scalar_logical(is.long))

  if (!is.data.frame(x)) {
    if (verbose) message("Object is not a data frame-like object.")
    return(FALSE)
  }

  if (!is.long) {
    if (verbose) message("Currently only MVGs in long format can be recognized.")
    return(FALSE)
  }

  # REMOVED until a better solution is found for category handling since this does lead to problems when working with data frames
  # that use a different category name.
  # if (category %nin% names(x)) {
  #   if (verbose) message("x is missing a category column. Did you use another name for this column? You can use the category
  #           argument to specify the name of that column.")
  #   return(FALSE)
  # }

  if (any(get_expected_columns_for_MVG() %nin% names(x))) {
    if (verbose) message(paste("x is missing a required column: ", paste(get_expected_columns_for_MVG, collapse = ",")))
    return(FALSE)
  }

  # Check that category is a factor only after everything else is checked.
  if (any(!is.factor(get(category, x)))) return(FALSE)

  # Check that mu and Sigma contain the cue names and that those cue names match.
  names_mu = names(x$mu[[1]])
  names_Sigma = dimnames(x$Sigma[[1]])
  if (!all(
    names_Sigma[[1]] == names_Sigma[[2]],
    names_Sigma[[1]] == names_mu)) {
    if (verbose) message("Names of cue dimensions do not match between mu and Sigma.")
    return(FALSE)
  }

  return(TRUE)
}


