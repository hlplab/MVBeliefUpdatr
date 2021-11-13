#' @rdname get_categorization_function
#' @export
get_categorization_function_from_NIW_ideal_adaptor = function(x, ...) {
  get_categorization_function(
    ms = x$m,
    Ss = x$S,
    kappas = x$kappa,
    nus = x$nu,
    lapse_rate = x$lapse_rate,
    bias = x$bias,
    ...
  )
}

#' Get categorization from an NIW belief
#'
#' Categorize a single observation based on an NIW belief. The decision rule can be specified to be either the
#' criterion choice rule, proportional matching (Luce's choice rule), or the sampling-based interpretation of
#' Luce's choice rule.
#'
#' @param x A vector of observations.
#' @param belief An \code{\link[=is.NIW_belief]{NIW_belief}} object.
#' @param decision_rule Must be one of "criterion", "proportional", or "sampling".
#' @param add_noise Determines whether multivariate Gaussian noise is added to the input.
#' If `NULL`, no noise is added during the updating. If "sample" then a sample of
#' noise is added to the input. If "marginalize" then each observation is transformed into the marginal distribution
#' that results from convolving the input with noise. This latter option might be helpful, for example, if one is
#' interested in estimating the consequences of noise across individuals. If add_noise is not `NULL` a Sigma_noise
#' column must be present in the NIW_belief object specified as the priors argument. (default: `NULL`)
#' @param add_lapse Determines the proportion of trials on which the listener lapses, not categorizing based on the
#' percept. If add_lapse is "sample", the lapse rate and bias from the ideal adaptor will be used. If NULL, the lapse
#' is assumed to be zero. (default: `NULL`)
#' @param simplify Should the output be simplified, and just the label of the selected category be returned? This
#' option is only available for the criterion and sampling decision rules. (default: `FALSE`)
#'
#' @return Either a tibble of observations with posterior probabilities for each category (in long format), or a
#' character vector indicating the chosen category in the same order as the observations in x (if simplify = `TRUE`).
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
  add_noise = NULL,
  add_lapse = NULL,
  simplify = F
) {
  # TO DO: check dimensionality of x with regard to belief.
  assert_NIW_belief(belief)
  assert_that(decision_rule  %in% c("criterion", "proportional", "sampling"),
              msg = "Decision rule must be one of: 'criterion', 'proportional', or 'sampling'.")
  assert_that(any(is.null(add_noise), add_noise %in% c("sample", "marginalize")),
              msg = 'add_noise must be one of "sample" or "marginalize"')
  assert_that(any(is.null(add_noise), "Sigma_noise" %in% names(belief)),
              msg = "Can't add noise: argument belief does not have column Sigma_noise.")
  if (!is.null(add_noise)) {
    assert_that(any(add_noise %in% "marginalize",
                    all(add_noise %in% c("sample"), is_weakly_greater_than(length(x), 1))),
                msg = "For noise sampling, x must be of length 1 or longer.")

    # Handle noise
    if (add_noise == "sample") {
      x = x + rmvnorm(n = length(x), sigma = belief$Sigma_noise[[1]])
    } else if (add_noise == "marginalize") message("Method 'marginalize' not yet implemented for add_noise.")
  }
  assert_that(any(is.null(add_lapse), is_scalar_character(add_lapse)),
              msg = "add_lapse must be NULL or a character.")
  if (is_scalar_character(add_lapse)) {
    assert_that(add_lapse == "sample", msg = "If add_lapse is a character, it must be 'sample'.")
    assert_that(assert_NIW_ideal_adaptor(belief), msg = "If add_lapse = 'sample', belief must be an  ideal adaptor.")
    lapse_rate = unique(belief$lapse_rate)
    lapse_biases = belief$bias
  } else {
    lapse_rate = 0
    lapse_biases = rep(0, get_nlevels_of_category_labels_from_model(belief))
  }

  if (!is.list(x)) x <- list(x) # in case a single x is handed as argument
  p = get_posterior_predictive_from_NIW_belief(x = x, belief = belief, log = F) %>%
    group_by(category) %>%
    mutate(
      observation = 1:length(x),
      x = x,
      prior = belief$prior[category]) %>%
    group_by(observation) %>%
    mutate(pp = (pp * prior) / sum(pp * prior))

  # Consider lapses
  p %<>%
    mutate(pp = ifelse(
      rep(
        rbinom(1, 1, lapse_rate),
        get_nlevels_of_category_labels_from_model(belief)),
      lapse_biases, # substitute lapse probabilities for posterior
      pp))

  if (decision_rule == "criterion") {
    p %<>%
      mutate(
        # tie breaker in case of uniform probabilities
        pp = ifelse(
          rep(
            sum(pp == max(pp)) > 1,
            get_nlevels_of_category_labels_from_model(belief)),
          pp + runif(
            get_nlevels_of_category_labels_from_model(belief),
            min = 0,
            max = 0),
          pp),
        # select most probable category
        response = ifelse(pp == max(pp), 1, 0))
  } else if (decision_rule == "sampling") {
    p %<>%
      mutate(response = rmultinom(1, 1, pp) %>% as.vector())
  } else if (decision_rule == "proportional") {
    p %<>%
      mutate(response = pp)
  }

  p %<>%
    ungroup() %>%
    select(-pp) %>%
    select(observation, x, category, response)

  if (simplify) {
    assert_that(decision_rule  %in% c("criterion", "sampling"),
                msg = "For simplify = T, decision rule must be either criterion or sampling.")
    return(p %>%
             filter(response == 1) %>%
             select(observation, category) %>%
             arrange(observation) %>%
             rename(response = category) %>%
             ungroup() %>%
             pull(response))
  } else return(p)
}
