#' Is this a Normal-Inverse-Wishart (NIW) belief?
#'
#' Check whether \code{x} is a Normal-Inverse-Wishard (NIW) belief/set of NIW beliefs. An NIW belief describes a distribution of
#' \link[=is.MVG]{multivariate Gaussian categories}.
#'
#' @param x Object to be checked.
#' @param category Name of the category variable. (default: "category")
#' @param is.long Is this check assessing whether the belief is in long format (`TRUE`) or wide format (`FALSE`)?
#' (default: `TRUE`)
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
is.NIW_belief = function(x, category = "category", is.long = T, verbose = F) {
  assert_that(is.flag(is.long))

  if (
    any(
      !is.long,
      all(!is_tibble(x), !is.data.frame(x))
    )
  ) {
    if (verbose) message("Currently only NIW beliefs in long format can be recognized.")
    return(FALSE)
  }

  if (category %nin% names(x)) {
    if (verbose) message("x is missing a category column. Did you use another name for this column? You can use the category
            argument to specify the name of that column.")
    return(FALSE)
  }

  if (any(c("kappa", "nu", "m", "S") %nin% names(x))) {
    if (verbose) message("x is missing at least one of kappa, nu, m, or S.")
    return(FALSE)
  }

  # Check that category is a factor only after everything else is checked.
  if (any(!is.factor(get(category, x)))) return(FALSE)

  # Check that m and S contain the cue names and that those cue names match.
  names_m = names(x$m[[1]])
  names_S = dimnames(x$S[[1]])
  if (!all(
    names_S[[1]] == names_S[[2]],
    names_S[[1]] == names_m)) {
    if (verbose) message("Names of cue dimensions do not match between m and S.")
    return(FALSE)
  }

  return(TRUE)
}

