#' @rdname get_categorization_function
#' @export
get_categorization_function_from_NIW_ideal_adaptor = function(x, ...) {
  get_categorization_function(
    ms = x$m,
    Ss = x$S,
    kappas = x$kappa,
    nus = x$nu,
    lapse_rate = x$lapse_rate,
    ...
  )
}

#' Get categorization from an NIW belief
#'
#' Categorize a single observation based on an NIW belief. The decision rule can be specified to be either the
#' criterion choice rule, proportional matching (Luce's choice rule), or the sampling-based interpretation of
#' Luce's choice rule.
#'
#' @param x An observation.
#' @param belief An \code{\link[=is.NIW_belief]{NIW_belief}} object.
#' @param decision_rule Must be one of "criterion", "proportional", or "sampling".
#' @param simplify Should the output be simplified, and just the label of the selected category be returned? This
#' option is only available for the criterion and sampling decision rules. (default: `FALSE`)
#'
#' @return Either a vector of posterior probabilities of the same length as the number of categories in the NIW
#' belief object, or a character indicating the chosen category (if simplify = `TRUE`).
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
#'
get_categorization_from_NIW_ideal_adaptor = function(
  x,
  belief,
  decision_rule,
  simplify = F
) {
  # TO DO: check dimensionality of x with regard to belief.
  # NOTE: Currently this function is intended for 1 observation only. Category priors are not yet incorporated.
  assert_NIW_belief(belief)
  assert_that(decision_rule  %in% c("criterion", "proportional", "sampling"),
              msg = paste("Decision rule must be one of:", c("criterion", "proportional", "sampling")))

  p = get_posterior_predictive_from_NIW_belief(x = x, belief = belief, log = F) %>%
    pull(pp) %>%
    divide_by(sum(.))

  if (decision_rule == "criterion")
    p = ifelse(p == max(p), 1, 0) else if (decision_rule == "sampling")
      # could change the first 1 to be based on number of observations
      p = as.vector(rmultinom(1, 1, p))

  if (simplify) {
    assert_that(decision_rule  %in% c("criterion", "sampling"),
                msg = "For simplify = T, decision rule must be either criterion or sampling.")
    return(belief$category[which(p == 1)])
  } else return(p)
}
