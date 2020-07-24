#' Get cue labels from NIW belief object
#'
#' Get the names for all cues from an NIW belief object.
#'
#' @param x An NIW belief object.
#'
#' @export
get_cue_labels_from_NIW_belief = function(x) {
  return(names(x$M[[1]]))
}



#' @rdname get_posterior_predictive
#' @export
get_posterior_predictive_from_NIW_belief = function(
  x,
  belief,
  log = T,
  category = "category",
  category.label = NULL,
  wide = FALSE
) {
  assert_that(is.NIW_belief(belief))
  assert_that(any(is.null(category.label) | is.character(category.label)))

  if (is.null(category.label)) {
    belief %<>%
      droplevels()
    category.label = belief %>%
      pull(!! sym(category)) %>%
      unique()
  }

  pp = foreach(c = category.label) %do% {
    b = belief %>%
      filter(!! sym(category) == c)

    get_posterior_predictive(
      x = x,
      M = b$M[[1]], S = b$S[[1]], kappa = b$kappa[[1]], nu = b$nu[[1]], log = log) %>%
      as_tibble() %>%
      rename_all(~ if (log) "lpp" else "pp") %>%
      mutate(!! sym(category) := c)
  }

  pp = reduce(pp, rbind)
  if (wide)
    pp %<>%
    pivot_wider(
      values_from = if (log) "lpp" else "pp",
      names_from = !! sym(category),
      names_prefix = if (log) "lpp." else "pp.") %>%
    unnest()

  return(pp)
}

# If there's a grouping variable extract the pp for each level of that grouping variable
get_posterior_predictives_from_NIW_beliefs = function(
  x,
  belief,
  log = T,
  category = "category",
  category.label = NULL,
  grouping.var,
  wide = FALSE
) {
  if (is.null(grouping.var)) {
    return(get_posterior_predictive_from_NIW_belief(
      x,
      belief,
      log = log,
      category = category,
      category.label = category.label,
      wide = wide))
  } else {
    assert_that(grouping.var %in% names(x),
                msg = "Grouping variable not found in the NIW belief object.")

    foreach (i = unique(x[[grouping.var]])) %do% {
      pp = get_posterior_predictive_from_NIW_belief(
        x,
        belief %>% filter(!! sym(grouping.var) == i),
        log = log,
        category = category,
        category.label = category.label,
        wide = wide
      ) %>%
        mutate(!! sym(grouping.var) := i)
    }

    pp = reduce(pp, rbind)
    return(pp)
  }
}


#' @rdname get_categorization_function
#' @export
get_categorization_function_from_NIW_belief = function(x, ...) {
  get_categorization_function(
    Ms = x$M,
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
#' option is only availale for the criterion and sampling decision rules. (default: `FALSE`)
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
get_categorization_from_NIW_belief = function(
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
    return(belief$category[which(p) == 1])
  } else return(p)
}
