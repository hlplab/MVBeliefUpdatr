#' Is this an ideal adaptor with Normal-Inverse-Wishart (NIW) beliefs?
#'
#' Check whether \code{x} is an ideal adaptor with \link[=is.NIW_belief]{Normal-Inverse-Wishard (NIW) beliefs}. An ideal adaptor
#' describes a distribution over \link[=is.MVG_ideal_observer]{ideal observers with multivariate Gaussian categories}. In this
#' sense, an ideal adaptor describes uncertainty about the true ideal observer. Optionally, one can also check whether a lapse
#' rate and bias is part of the ideal adaptor.
#'
#' So far, the ideal adaptor is assumed to have perfect certainty about the prior, lapse rate and bias (if present).
#' Future implementations might allow uncertainty over these parameters.
#'
#' @param x Object to be checked.
#' @param category Name of the category variable. (default: "category")
#' @param is.long Is this check assessing whether the ideal adaptor is in long format (`TRUE`) or wide format (`FALSE`)?
#' (default: `TRUE`)
#' @param with.lapse Does this ideal adaptor have a lapse rate? (default: `FALSE`)
#' @param with.bias Does this ideal adaptor have a bias? (default: `FALSE`)
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
is.NIW_ideal_adaptor = function(x, category = "category", is.long = T, with.lapse = if (with.bias) T else F, with.bias = F, verbose = F) {
  assert_that(all(is.flag(with.lapse), is.flag(with.bias)))

  if (!is.NIW_belief(x)) {
    if (verbose) message("x does not contain NIW beliefs.")
    return(FALSE)
  }

  if (
    any(
      "prior" %nin% names(x),
      with.lapse & "lapse_rate" %nin% names(x),
      with.bias & "bias" %nin% names(x)
    )
  ) {
    if (verbose) message("x is missing prior, lapse rate, or bias.")
    return(FALSE)
  }

  # Check that category is a factor only after everything else is checked.
  if (any(!is.factor(get(category, x)))) return(FALSE)

  # Check that the prior probabilities add up to 1
  if (sum(x$prior) != 1) {
    if (verbose) message(paste("Prior probabilities in x do not add up to 1: ", sum(x$prior)))
    return(FALSE)
  }

  return(TRUE)
}

#' @describeIn is.NIW_ideal_adaptor Also checks whether the ideal adaptor has a lapse term.
#' @export
is.NIW_ideal_adaptor_w_lapse = function(x, category = "category", is.long = T, with.bias = F) {
  is.NIW_ideal_adaptor(x, category = category, is.long = is.long, with.lapse = T, with.bias = with.bias)
}

#' @describeIn is.NIW_ideal_adaptor Also checks whether the ideal adaptor has a bias term.
#' @export
is.NIW_ideal_adaptor_w_bias = function(x, category = "category", is.long = T) {
  is.NIW_ideal_adaptor(x, category = category, is.long = is.long, with.lapse = T, with.bias = T)
}


