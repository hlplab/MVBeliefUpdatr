#' @rdname get_NIW_posterior_predictive
#' @export
get_posterior_predictive_from_NIW_belief = function(
  x,
  model,
  log = T,
  noise_treatment = if (is.NIW_ideal_adaptor(model)) "marginalize" else "no_noise",
  category = "category",
  category.label = NULL,
  wide = FALSE
) {
  assert_that(is.NIW_belief(model))
  assert_that(any(is.null(category.label) | is.character(category.label)))
  assert_that(any(noise_treatment == "no_noise", is.NIW_ideal_adaptor(model)),
              msg = 'No noise matrix Sigma_noise found. If noise_treatment is not "no_noise", then model must be an NIW_ideal_adaptor.')

  if (is.null(category.label)) {
    model %<>%
      droplevels()
    category.label = model %>%
      pull(!! sym(category)) %>%
      unique()
  }

  posterior_predictive <- foreach(c = category.label) %do% {
    b <-
      model %>%
      filter(!! sym(category) == c)

    get_NIW_posterior_predictive(
      x = x,
      m = b$m[[1]],
      S = b$S[[1]],
      kappa = b$kappa[[1]],
      nu = b$nu[[1]],
      log = log,
      noise_treatment = noise_treatment,
      Sigma_noise = if (noise_treatment == "no_noise") NULL else b$Sigma_noise[[1]]) %>%
      as_tibble() %>%
      rename_with(~ if (log) "log_posterior_predictive" else "posterior_predictive") %>%
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

# If there's a grouping variable extract the posterior predictive for each level of that grouping variable
get_posterior_predictives_from_NIW_beliefs = function(
  x,
  model,
  log = T,
  noise_treatment = if (is.NIW_ideal_adaptor(model)) "marginalize" else "no_noise",
  category = "category",
  category.label = NULL,
  grouping.var,
  wide = FALSE
) {
  if (is.null(grouping.var)) {
    return(get_posterior_predictive_from_NIW_belief(
      x,
      model,
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
          model %>% filter(!! sym(grouping.var) == i),
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



