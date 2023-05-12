get_expected_columns_for_MVG_ideal_observer <- function() append(get_expected_columns_for_MVG(), get_expected_columns_for_model())

#' Is this an ideal observer with multivariate Gaussian (MVG) categories?
#'
#' Check whether \code{x} is an ideal observer with \link[=is.MVG]{multivariate Gaussian (MVG) categories}. Optionally, one can also check whether a lapse rate
#' and lapse bias is part of the ideal observer.
#'
#' @param x Object to be checked.
#' @param group Name of one or more group variables, each unique combination of which describes an MVG_ideal_observer. (default: NULL)
#' @param category Name of the category variable. (default: "category")
#' @param is.long Is this check assessing whether the ideal observer is in long format (`TRUE`) or wide format (`FALSE`)?
#' (default: `TRUE`)
#' @param with.lapse Does this ideal observer have a lapse rate? (default: `FALSE`)
#' @param with.lapse_bias Does this ideal observer have a lapse bias? (default: `FALSE`)
#' @param verbose Should verbose output be provided? (default: `TRUE`)
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
is.MVG_ideal_observer <- function(x, group = NULL, category = "category", is.long = T, with.lapse = if (with.lapse_bias) T else F, with.lapse_bias = F, verbose = F, tolerance = 1e-5) {
  name_of_x <- deparse(substitute(x))
  assert_that(all(is.flag(with.lapse), is.flag(with.lapse_bias)))

  # When no groups are specified, infer groups from object.
  if (is.null(group)) {
    group <- setdiff(names(x), get_expected_columns_for_MVG_ideal_observer())
    if (length(group) == 0) group <- NULL else {
      if (verbose) message(paste(name_of_x, "has additional columns beyond those expected:", paste(group, collapse = ", "), "Interpreting those columns as group variables."))
    }
  }

  if (!is.null(group)) {
    if (verbose) message("Checking whether ", name_of_x, " is an MVG_ideal_observer within each unique combination of group values.")
    x %<>% group_by(!!! syms(group))
  }

  if (!is.MVG(x, category = category, group = group, verbose = verbose)) {
    if (verbose) message("x does not contain multivariate Gaussian categories.")
    return(FALSE)
  }

  # Only need to test for MVG columns here since is.model is called below.
  if (any(get_expected_columns_for_MVG() %nin% names(x))) {
    if (verbose) message(paste("x is missing a required column: ", paste(get_expected_columns_for_MVG, collapse = ",")))
    return(FALSE)
  }

  if (any(!is.factor(get(category, x)))) return(FALSE)

  if (!is.model(x, group = group, verbose = verbose, tolerance = tolerance)) {
    return(FALSE)
  }

  return(TRUE)
}


