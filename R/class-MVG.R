#' Is this a set of multivariate Gaussian (MVG) categories?
#'
#' Check whether \code{x} is a set of multivariate Gaussian (MVG) categories.
#'
#' @param x Object to be checked.
#' @param category Name of the category variable. (default: "category")
#' @param is.long Is this check assessing whether the ideal observer is in long format (`TRUE`) or wide format (`FALSE`)?
#' (default: `TRUE`)
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
is.MVG = function(x, category = "category", is.long = T, verbose = F) {
  assert_that(is.flag(is.long))

  if (
    any(
      !is.long,
      all(!is_tibble(x), !is.data.frame(x))
    )
  ) {
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

  if (any(c("mu", "Sigma") %nin% names(x))) {
    if (verbose) message("x is missing either mu or Sigma.")
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


