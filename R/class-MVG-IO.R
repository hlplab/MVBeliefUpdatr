get_expected_columns_for_MVG_ideal_observer <- function()
  c(get_expected_columns_for_MVG, "prior", "lapse_rate", "lapse_bias", "Sigma_noise")

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
is.MVG_ideal_observer = function(x, category = "category", is.long = T, with.lapse = if (with.lapse_bias) T else F, with.lapse_bias = F, verbose = F) {
  assert_that(all(is.flag(with.lapse), is.flag(with.lapse_bias)))

  if (!is.MVG(x, category = category, verbose = verbose)) {
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

  if (any(!is.factor(get(category, x)))) return(FALSE)

  groups <- setdiff(names(x), get_expected_columns_for_MVG_ideal_observer())
  if (length(groups) > 0) {
    if (verbose) message(paste(deparse(substitute(x)), "has additional columns beyond those expected. Checking whether",
                               deparse(substitute(x)), "is an MVG_ideal_observer within each unique combination of those additional variables."))
    x %<>%
      group_by(!!! syms(groups))
  }

  # Check that the prior probabilities add up to 1
  if (any(x %>% summarise(sum_prior = sum(prior)) %>% pull(sum_prior) != 1)) {
    if (verbose) message(paste("Prior probabilities in", deparse(substitute(x)), "do not add up to 1: ", sum(x$prior)))
    return(FALSE)
  }

  # Check that the lapse rate is constant across categories
  if (with.lapse &
      any(x %>% summarise(n_unique_lapse_rates = length(unique(lapse_rate))) %>% pull(n_unique_lapse_rates) != 1)) {
    if (verbose) message(paste("Lapse rates in", deparse(substitute(x)), "are not constant across categories: ", paste(x$lapse_rate, collapse = ", ")))
    return(FALSE)
  }

  # Check that the lapse bias probabilities add up to 1
  if (with.lapse_bias &
      any(x %>% summarise(sum_lapse_bias = sum(lapse_bias)) %>% pull(sum_lapse_bias) != 1)) {
    if (verbose) message(paste("Lapse bias probabilities in", deparse(substitute(x)), "do not add up to 1: ", sum(x$lapse_bias)))
    return(FALSE)
  }

  return(TRUE)
}


