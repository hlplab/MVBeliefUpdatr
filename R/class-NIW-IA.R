get_expected_columns_for_NIW_ideal_adaptor <- function()
  c(get_expected_columns_for_NIW_belief(), "prior", "lapse_rate", "lapse_bias", "Sigma_noise")

#' Is this an ideal adaptor with Normal-Inverse-Wishart (NIW) beliefs?
#'
#' Check whether \code{x} is an ideal adaptor with \link[=is.NIW_belief]{Normal-Inverse-Wishard (NIW) beliefs}. An ideal adaptor
#' describes a distribution over \link[=is.MVG_ideal_observer]{ideal observers with multivariate Gaussian categories}. In this
#' sense, an ideal adaptor describes uncertainty about the true ideal observer. Optionally, one can also check whether a lapse
#' rate and lapse bias is part of the ideal adaptor.
#'
#' So far, the ideal adaptor is assumed to have perfect certainty about the prior, lapse rate and lapse bias (if present).
#' Future implementations might allow uncertainty over these parameters.
#'
#' @param x Object to be checked.
#' @param category Name of the category variable. (default: "category")
#' @param is.long Is this check assessing whether the ideal adaptor is in long format (`TRUE`) or wide format (`FALSE`)?
#' (default: `TRUE`)
#' @param with.lapse Does this ideal adaptor have a lapse rate? (default: `FALSE`)
#' @param with.lapse_bias Does this ideal adaptor have a lapse bias? (default: `FALSE`)
#' @param verbose Should verbose output be provided? (default: `TRUE`)
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
is.NIW_ideal_adaptor = function(x, category = "category", is.long = T, with.prior = T, with.lapse = if (with.lapse_bias) T else F, with.lapse_bias = F, verbose = F, tolerance = 1e-5) {
  name_of_x <- deparse(substitute(x))
  assert_that(all(is.flag(with.lapse), is.flag(with.lapse_bias)))

  if (!is.NIW_belief(x)) {
    if (verbose) message(paste(deparse(substitute(x)), "does not contain NIW beliefs."))
    return(FALSE)
  }

  if (
    any(
      !with.prior | "prior" %nin% names(x),
      with.lapse & "lapse_rate" %nin% names(x),
      with.lapse_bias & "lapse_bias" %nin% names(x)
    )
  ) {
    if (verbose) message(paste(name_of_x, " is missing prior, lapse rate, or lapse bias."))
    return(FALSE)
  }

  if (any(!is.factor(get(category, x)))) return(FALSE)

  groups <- setdiff(names(x), get_expected_columns_for_NIW_ideal_adaptor())
  if (length(groups) > 0) {
    if (verbose) message(paste(name_of_x, "has additional columns beyond those expected:", paste(groups, collapse = ", "),
                               "Checking whether", name_of_x, "is an NIW_ideal_adaptor within each unique combination of those additional variables."))
    x %<>%
      group_by(!!! syms(groups))
  }

  return(is.model(x, verbose = verbose, tolerance = tolerance))
}



