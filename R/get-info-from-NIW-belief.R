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

  posterior_predictive <- foreach(c = category.label) %do% {
    b <-
      belief %>%
      filter(!! sym(category) == c)

    get_posterior_predictive(
      x = x,
      m = b$m[[1]], S = b$S[[1]], kappa = b$kappa[[1]], nu = b$nu[[1]], log = log) %>%
      as_tibble() %>%
      rename_all(~ if (log) "log_posterior_predictive" else "posterior_predictive") %>%
      mutate(!! sym(category) := c)
  }

  posterior_predictive %<>% reduce(rbind)
  if (wide)
    posterior_predictive %<>%
    pivot_wider(
      values_from = if (log) "log_posterior_predictive" else "posterior_predictive",
      names_from = !! sym(category),
      names_prefix = if (log) "log_posterior_predictive." else "posterior_predictive.") %>%
    unnest()

  return(posterior_predictive)
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
      posterior_predictive <-
        get_posterior_predictive_from_NIW_belief(
          x,
          belief %>% filter(!! sym(grouping.var) == i),
          log = log,
          category = category,
          category.label = category.label,
          wide = wide
        ) %>%
        mutate(!! sym(grouping.var) := i)
    }

    posterior_predictive %<>% reduce(rbind)
    return(posterior_predictive)
  }
}



