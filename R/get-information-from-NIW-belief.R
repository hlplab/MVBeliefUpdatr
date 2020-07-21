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
  attest_that(is.NIW_belief(belief))
  attest_that(any(is.null(category.label) | is.character(category.label)))

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
    attest_that(grouping.var %in% names(x),
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
