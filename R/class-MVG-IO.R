#' Is this an ideal observer with multivariate Gaussian (MVG) categories?
#'
#' Check whether \code{x} is an ideal observer with \link[=is.MVG]{multivariate Gaussian (MVG) categories}. Optionally, one can also check whether a lapse rate
#' and lapse bias is part of the ideal observer.
#'
#' @param x Object to be checked.
#' @param category Name of the category variable. (default: "category")
#' @param is.long Is this check assessing whether the ideal observer is in long format (`TRUE`) or wide format (`FALSE`)?
#' (default: `TRUE`)
#' @param with.lapse Does this ideal observer have a lapse rate? (default: `FALSE`)
#' @param with.lapse_bias Does this ideal observer have a lapse bias? (default: `FALSE`)
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
is.MVG_ideal_observer = function(x, category = "category", is.long = T, with.lapse = if (with.bias) T else F, with.lapse_bias = F, verbose = F) {
  assert_that(all(is.flag(with.lapse), is.flag(with.lapse_bias)))

  if (!is.MVG(x)) {
    if (verbose) message("x does not contain multivariate Gaussian categories.")
    return(FALSE)
  }

  if (
    any(
      "prior" %nin% names(x),
      with.lapse & "lapse_rate" %nin% names(x),
      with.lapse_bias & "lapse_bias" %nin% names(x)
    )
  ) {
    if (verbose) message("x is missing prior, lapse rate, or lapse bias.")
    return(FALSE)
  }

  # Check that the prior probabilities add up to 1
  if (sum(x$prior) != 1) {
    if (verbose) message(paste("Prior probabilities in x do not add up to 1: ", sum(x$prior)))
    return(FALSE)
  }

  return(TRUE)
}

#' @describeIn is.MVG_ideal_observer Also checks whether the MVG has a lapse term.
#' @export
is.MVG_ideal_observer_w_lapse = function(x, category = "category", is.long = T, with.lapse_bias = F) {
  is.MVG_ideal_observer(x, category = category, is.long = is.long, with.lapse = T, with.lapse_bias = with.bias)
}

#' @describeIn is.MVG_ideal_observer Also checks whether the MVG has a lapse bias term.
#' @export
is.MVG_ideal_observer_w_lapse_bias = function(x, category = "category", is.long = T) {
  is.MVG_ideal_observer(x, category = category, is.long = is.long, with.lapse = T, with.lapse_bias = T)
}


