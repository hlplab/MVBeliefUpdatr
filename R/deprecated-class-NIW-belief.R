get_expected_columns_for_NIW_belief <- function()
  c("category", "m", "S", "kappa", "nu")

#' Deprecated: is.NIW_belief
#'
#' @description `r lifecycle::badge("deprecated")`
#' `is.NIW_belief()` was deprecated in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0.
#' Please use S7 validators and [NIW_IdealAdaptor] instead.
#'
#' @param x Object to be checked.
#' @param group Name of one or more group variables, each unique combination of which describes an NIW_belief. (default: NULL)
#' @param category Name of the category variable. (default: "category")
#' @param is.long Is this check assessing whether the belief is in long format (`TRUE`) or wide format (`FALSE`)?
#' (default: `TRUE`)
#' @param verbose Logical. If `TRUE`, emits diagnostics.
#'
#' @return A logical.
#'
#' @seealso [NIW_IdealAdaptor]
#' @keywords internal
#' @export
is.NIW_belief <- function(x, group = NULL, category = "category", is.long = T, verbose = F) {
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "is.NIW_belief()",
    details = "Use S7 validators and NIW_IdealAdaptor instead."
  )
  name_of_x <- deparse(substitute(x))
  .assert_that(.is_scalar_logical(is.long))

  if (!is.data.frame(x)) {
    if (verbose) message("Object is not a data frame-like object.")
    return(FALSE)
  }

  if (!is.long) {
    if (verbose) message("Currently only NIW beliefs in long format can be recognized.")
    return(FALSE)
  }

  if (category %nin% names(x)) {
    if (verbose) message(paste(deparse(substitute(x)), "is missing a category column. Did you use another name for this column? You can use the category
            argument to specify the name of that column."))
    return(FALSE)
  }

  if (any(c("kappa", "nu", "m", "S") %nin% names(x))) {
    if (verbose) message(paste(deparse(substitute(x)), "is missing at least one of kappa, nu, m, or S."))
    return(FALSE)
  }

  # Check that category is a factor only after everything else is checked.
  if (any(!is.factor(get(category, x)))) {
    if (verbose) message("category is not a factor.")
    return(FALSE)
  }

  # Check that m and S contain the cue names and that those cue names match. This also serves as
  # check that m and S have appropriate dimensions.
  names_m = names(x$m[[1]])
  names_S = dimnames(x$S[[1]])
  if (!all(
    names_S[[1]] == names_S[[2]],
    names_S[[1]] == names_m)) {
    if (verbose) message("Names of cue dimensions do not match between m and S.")
    return(FALSE)
  }

  # Check nu vs. dimensionality of S
  D_S = if (is.null(dim(x$S[[1]]))) 1 else dim(x$S[[1]])[1]
  if (any(x$nu <= D_S + 1))
    warning(paste0("At least one category had nu smaller than allowed (is ", min(x$nu), "; should be >", D_S + 1, ").\n"))

  return(TRUE)
}

