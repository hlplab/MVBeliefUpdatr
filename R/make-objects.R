#' Make multivariate Gaussian ideal observer(s) from data.
#'
#' Constructs an \code{\link[=is.MVG]{MVG}} or \code{\link[=is.MVG_ideal_observer]{MVG_ideal_observer}} object with category
#' information for all categories found in the data.
#'
#'
#' @param data The tibble or data.frame from which to construct the MVG or MVG ideal observer object.
#' @param group Optionally, one or more grouping variables can be specified. If group is not NULL, one MVG or
#' ideal observers will be derived for each level of \code{group}. (default: NULL)
#' @param category Name of variable in \code{data} that contains the category information. (default: "category")
#' @param cues Name(s) of variables in \code{data} that contain the cue information.
#' @param prior Optionally specify a prior probability for each category (in each group). (default: a uniform
#' prior over all categories).
#' @param lapse_rate Optionally specify a lapse rate. (default: \code{NA})
#' @param lapse_bias Optionally specify a lapse bias. (default: \code{NA})
#' @param Sigma_noise Optionally specify a (multivariate Gaussian) covariance matrix of perceptual noise. This
#' argument will be ignored if NULL. (default: NULL)
#' @param add_Sigma_noise_to_category_representation Should the perceptual noise be added to the category
#' representation (category variability)? If FALSE, then noise will be considered during categorization decisions
#' but will not be added to the MVG categories. If TRUE, then noise will also be to the category covariance matrix.
#' This is typically the desired result since perceptual noise would have shaped the ideal observers representations.
#' Will be ignored if Sigma_noise is NULL. (default: TRUE)
#' @param verbose If true provides more information. (default: FALSE)
#'
#' @return A tibble that is an MVG or MVG ideal observer object.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
#'
make_MVG_from_data = function(
  data,
  group = NULL,
  category = "category",
  cues,
  verbose = F
) {
  assert_that(is.data.frame(data) | is_tibble(data))
  assert_that(all(is.null(group) | all(is.character(group) | is_symbol(group))))
  assert_that(all(is.character(category) | is_symbol(category), length(category) == 1))
  assert_that(all(is.character(cues) | is_symbol(cues), length(cues) > 0))

  if (is.character(group)) group <- syms(group)
  if (is.character(category)) category <- sym(category)
  if (is.character(cues)) {
    cue_names <- cues
    cues <- syms(cues)
  }

  assert_that(as_name(category) %in% names(data),
              msg = paste0("Category variable (", as_name(category), ") not found in data."))
  assert_that(all(cue_names %in% names(data)),
              msg = paste0("Some cues not found in data: ", paste(setdiff(cue_names, intersect(cue_names, names(data))), collapse = ", ")))

  model <- data %>%
    select(!! category, !!! cues, !!! group) %>%
    mutate(cues = pmap(list(!!! cues),
                       .f = function(...) {
                         x = c(...)
                         names(x) = as.character(cues)
                         return(x) })) %>%
    { if (is.null(group)) group_by(., !! category) else group_by(., !!! group, !! category) } %>%
    summarise(
      mu = list(colMeans(cbind(!!! cues))),
      Sigma = list(cov(cbind(!!! cues))))

  if (!is.null(group))
    while(length(group) > 1) {
      group = group[2:length(group)]
      model %<>%
        group_by(., !!! group, !! category) %>%
        summarise_at(vars(starts_with(c("mu", "Sigma"))),
                     ~ list(reduce(.x, `+`) / length(.x)))
    }

  model %<>%
    mutate(
      !! category := factor(!! category)) %>%
    ungroup()

  if (is.null(group) & !is.MVG(model, category = as_name(category), verbose = verbose)) {
    warning("Something went wrong. The returned object is not an MVG. Try again with verbose = T?")
  } else if (!is.null(group)) {
    warning("Currently, groups of ideal observers are not checked for internal consistency.")
  }

  return(model)
}

#' @export
#' @rdname make_MVG_from_data
make_MVG_ideal_observer_from_data = function(
  data,
  group = NULL,
  category = "category",
  cues,
  prior = NA_real_,
  lapse_rate = NA_real_,
  lapse_bias = NA_real_,
  Sigma_noise = NULL,
  add_Sigma_noise_to_category_representation = T,
  verbose = F
) {
  model <- data %>% make_MVG_from_data(group = group, category = category, cues = cues, verbose = verbose)
  model %<>% lift_MVG_to_MVG_ideal_observer(group = group, category = category, prior = prior, lapse_rate = lapse_rate, lapse_bias = lapse_bias,
                                       Sigma_noise = Sigma_noise, add_Sigma_noise_to_category_representation = add_Sigma_noise_to_category_representation)

  if (is.null(group) & !is.MVG_ideal_observer(model, category = as_name(category), verbose = verbose)) {
    warning("Something went wrong. The returned object is not an MVG ideal observer. Try again with verbose = T?")
  } else if (!is.null(group)) {
    warning("Currently, groups of ideal observers are not checked for internal consistency.")
  }

  return(model)
}


#' Make NIW belief from data.
#'
#' Constructs an \code{\link[=NIW_belief]{NIW_belief}} or \code{\link[=NIW_ideal_adaptor]{NIW_ideal_adaptor}} object,
#' representing Normal-Inverse Wishart (NIW) parameters for all categories found in the data. This object can be used
#' as a prior for functions like \code{\link{update_NIW_beliefs_incrementally}}.
#'
#' Currently, \code{make_NIW_prior_from_data()} does not infer kappa or nu, nor does it fit hierarchical data. Rather
#' the function simply estimates the category mean and covariance matrix from the sample (\code{data}), assumes them
#' to be the expected category mean (mu) and covariance (Sigma), and derives the m and S parameters
#' of the NIW from mu and Sigma based on the user-provided kappa and nu. That means m = mu and S = Sigma * (nu - D -1),
#' where D is the dimensionality of the data.
#'
#' @param data The tibble or data.frame from which to construct the prior.
#' @param group Optionally, one or more grouping variables can be specified. If group is not NULL, one NIW_belief or
#' ideal adaptor will be derived for each level of \code{group}. (default: NULL)
#' @param category Name of variable in \code{data} that contains the category information. (default: "category")
#' @param cues Name(s) of variables in \code{data} that contain the cue information.
#' @param kappa The strength of the beliefs over the category mean (pseudocounts). (default: same as nu)
#' @param nu The strength of the beliefs over the category covariance matrix (pseudocounts). (default: number of
#' cues + 2)
#' @param prior Optionally specify a prior probability for each category (in each group). (default: a uniform
#' prior over all categories).
#' @param lapse_rate Optionally specify a lapse rate. (default: \code{NA})
#' @param lapse_bias Optionally specify a lapse bias. (default: \code{NA})
#' @param Sigma_noise Optionally specify a (multivariate Gaussian) covariance matrix of perceptual noise. This argument will be
#' ignored if NULL. (default: NULL)
#' @param add_Sigma_noise_to_category_representation Should the perceptual noise be added to the category
#' representation (category variability)? If FALSE, then noise will be considered during categorization decisions
#' but will not be added to the MVG categories. If TRUE, then noise will also be to the category covariance matrix.
#' This is typically the desired result since perceptual noise would have shaped the ideal observers representations.
#' Will be ignored if Sigma_noise is NULL. (default: TRUE)
#' @param verbose If true provides more information. (default: FALSE)
#' @param keep.category_parameters Should categories' mu and Sigma be included in the output (in addition to m
#' and S of the prior)? (default: FALSE)
#'
#' @return A tibble that is an NIW_belief or NIW_ideal_adaptor object.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
#'
make_NIW_belief_from_data = function(
  data,
  group = NULL,
  category = "category",
  cues,
  kappa = nu,
  nu = length(cues) + 2,
  keep.category_parameters = F,
  verbose = F
) {
  assert_that(is.data.frame(data) | is_tibble(data))
  assert_that(all(is.null(group) | all(is.character(group) | is_symbol(group))))
  assert_that(all(is.character(category) | is_symbol(category), length(category) == 1))
  assert_that(all(is.character(cues) | is_symbol(cues), length(cues) > 0))
  assert_that(all(is.numeric(kappa), is.numeric(nu), !is.na(kappa), !is.na(nu)))
  assert_that(nu > length(cues) + 1,
              msg = paste0("nu must be larger than dimensionality of cues + 1 (>", length(cues) + 1, ")."))

  if (is.character(group)) group = syms(group)
  if (is.character(category)) category = sym(category)
  if (is.character(cues)) cues = syms(cues)

  model <- data %>%
    select(!! category, !!! cues, !!! group) %>%
    mutate(cues = pmap(list(!!! cues),
                       .f = function(...) {
                         x = c(...)
                         names(x) = as.character(cues)
                         return(x) })) %>%
    { if (is.null(group)) group_by(., !! category) else group_by(., !!! group, !! category) } %>%
    summarise(
      mu = list(colMeans(cbind(!!! cues))), # list(reduce(cues, `+`) / length(cues)),
      Sigma = list(cov(cbind(!!! cues))))

  message("S is set so that the expected category covariance matrix Sigma matches the category covariance in the sample (given nu). ",
          "It might be safer to fit an Inverse-Wishart distribution to the entire set of covariance matrices.")
  model %<>%
    mutate(
      !! category := factor(!! category),
      kappa = kappa,
      nu = nu,
      m = mu,
      S = get_S_from_expected_Sigma(Sigma, nu)) %>%
    ungroup()

  if (!keep.category_parameters) data %<>% select(-c(mu, Sigma))
  if (!is.NIW_belief(model, category = as_name(category), verbose = verbose))
    warning("Something went wrong. The returned object is not an NIW belief.")

  return(model)
}

#' @rdname make_NIW_belief_from_data
#' @export
make_NIW_prior_from_data <- make_NIW_belief_from_data

#' @export
#' @rdname make_NIW_belief_from_data
make_NIW_ideal_adaptor_from_data = function(
  data,
  group = NULL,
  category = "category",
  cues,
  kappa = nu,
  nu = length(cues) + 2,
  prior = NA_real_,
  lapse_rate = NA_real_,
  lapse_bias = NA_real_,
  Sigma_noise = NULL,
  add_Sigma_noise_to_category_representation = T,
  keep.category_parameters = F
) {
  model <- data %>% make_NIW_belief_from_data(group = group, category = category, cues = cues, verbose = verbose)
  model %<>% lift_NIW_belief_to_NIW_ideal_adaptor(group = group, category = category, prior = prior, lapse_rate = lapse_rate, lapse_bias = lapse_bias,
                                                 Sigma_noise = Sigma_noise, add_Sigma_noise_to_category_representation = add_Sigma_noise_to_category_representation)

  if (!keep.category_parameters) data %<>% select(-c(mu, Sigma))
  if (!is.NIW_ideal_adaptor(model, category = as_name(category), verbose = verbose))
    warning("Something went wrong. The returned object is not an NIW ideal adaptor.")

  return(model)
}


#' Turn an MVG/NIW_belief object into an ideal observer/adaptor
#'
#' Make an ideal observer or adaptor out of an MVG or NIW_belief object, respectively, by providing the missing
#' information about the prior, lapse rate, bias, and perceptual noise (if any).
#'
#' @param x Either an MVG or NIW_belief object.
#' @param group Optionally, a grouping structure can be specified. If group structure is not NULL, one
#' NIW belief or ideal adaptor will be derived for each level of \code{group_structure}. (default: NULL)
#' @param category Name of variable in \code{data} that contains the category information. (default: "category")
#' @param kappa The strength of the beliefs over the category mean (pseudocounts). (default: same as nu)
#' @param nu The strength of the beliefs over the category covariance matrix (pseudocounts). (default: number of
#' cues + 2)
#' @param prior Optionally specify a prior probability for each category (in each group). (default: a uniform
#' prior over all categories).
#' @param lapse_rate Optionally specify a lapse rate. (default: \code{NA})
#' @param lapse_bias Optionally specify a lapse bias. (default: \code{NA})
#' @param Sigma_noise Optionally specify a (multivariate Gaussian) covariance matrix of perceptual noise. This argument will be
#' ignored if NULL. (default: NULL)
#' @param add_Sigma_noise_to_category_representation Should the perceptual noise be added to the category
#' representation (category variability)? If FALSE, then noise will be considered during categorization decisions
#' but will not be added to the MVG categories. If TRUE, then noise will also be to the category covariance matrix.
#' This is typically the desired result since perceptual noise would have shaped the ideal observers representations.
#' Will be ignored if Sigma_noise is NULL. (default: TRUE)
#' @param verbose If true provides more information. (default: FALSE)
#' @param keep.category_parameters Should categories' mu and Sigma be included in the output (in addition to m
#' and S of the prior)? (default: FALSE)
#'
#' @return A tibble that is an NIW_belief or NIW_ideal_adaptor object.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
lift_likelihood_to_model = function(
  x,
  group = NULL,
  category = "category",
  prior = NA_real_,
  lapse_rate = NA_real_,
  lapse_bias = NA_real_,
  Sigma_noise = NULL
) {
  if (is.character(group)) group = syms(group)
  if (is.character(category)) category = sym(category)
  assert_that(all(is.numeric(lapse_rate), is.numeric(lapse_bias), is.numeric(prior)),
              msg = "The category prior, lapse rate, and lapse bias must be numeric.")

  n.cat = length(unique(x %>% pull(!! category)))
  if (all(is.na(prior) | is.null(prior))) {
    message(paste0("No prior specified. Defaulting to uniform prior over the ", n.cat, " categories found in x."))
    prior = rep(1 / n.cat, n.cat)
  }

  assert_that(all(
    is.na(lapse_rate) | between(lapse_rate, 0, 1),
    is.na(lapse_bias) | between(lapse_bias, 0, 1),
    between(prior, 0, 1)),
    msg = "If not NA, the category prior, lapse rate, and lapse bias must have values between 0 and 1.")
  assert_that(sum(prior) == 1,
              msg = paste0("Priors must add up to 1. (instead: ", sum(prior), ")."))

  assert_that(is.null(Sigma_noise) | is.matrix(Sigma_noise),
              msg = "If not NULL, Sigma_noise must be a matrix.")

  x %<>%
    group_by(!!! group) %>%
    mutate(
      prior = prior,
      lapse_rate = lapse_rate,
      lapse_bias = lapse_bias,
      Sigma_noise = list(Sigma_noise))

  return(x)
}

#' @export
#' @rdname lift_likelihood_to_model
lift_MVG_to_MVG_ideal_observer = function(
  x,
  group = NULL,
  category = "category",
  prior = NA_real_,
  lapse_rate = NA_real_,
  lapse_bias = NA_real_,
  Sigma_noise = NULL,
  add_Sigma_noise_to_category_representation = T
) {
  x %<>% lift_likelihood_to_model(group = group, category = category, prior = prior, lapse_rate = lapse_rate, lapse_bias = lapse_bias, Sigma_noise = Sigma_noise)
  if (!is.null(first(x$Sigma_noise))) {
    assert_that(all(dim(first(x$Sigma_noise)) == dim(first(x$Sigma))),
                msg = "If not NULL, Sigma_noise must be a matrix of the same dimensionality as Sigma.")
    assert_that(all(dimnames(first(x$Sigma_noise)) == dimnames(first(x$Sigma))),
                msg = "If Sigma_noise is not NULL, the dimnames of Sigma_noise and Sigma must match.")

    if (add_Sigma_noise_to_category_representation)
      x %<>%
        mutate(Sigma = map2(Sigma, Sigma_noise, ~ .x + .y))
  }

  return(x)
}

#' @export
#' @rdname lift_likelihood_to_model
lift_NIW_belief_to_NIW_ideal_adaptor = function(
  x,
  group = NULL,
  category = "category",
  prior = NA_real_,
  lapse_rate = NA_real_,
  lapse_bias = NA_real_,
  Sigma_noise = NULL,
  add_Sigma_noise_to_category_representation = T
) {
  x %<>% lift_likelihood_to_model(group = group, category = category, prior = prior, lapse_rate = lapse_rate, lapse_bias = lapse_bias, Sigma_noise = Sigma_noise)
  if (!is.null(first(x$Sigma_noise))) {
    assert_that(all(dim(first(x$Sigma_noise)) == dim(first(x$S))),
                msg = "If Sigma_noise is not NULL, Sigma_noise must be a matrix of the same dimensionality as S.")
    assert_that(all(dimnames(first(x$Sigma_noise)) == dimnames(first(x$S))),
                msg = "If Sigma_noise is not NULL, the dimnames of Sigma_noise and S must match.")

    if (add_Sigma_noise_to_category_representation)
      x %<>%
      mutate(S = map2(get_expected_Sigma_from_S(S, nu), Sigma_noise, ~ get_S_from_expected_Sigma(.x + .y, nu)))
  }

  return(x)
}

#' Aggregate multiple ideal observers/adaptors into an average ('typical') ideal observer/adaptor
#'
#' Takes a grouped collection of ideal observers or adaptors (or MVG / NIW_belief objects) and averages them into a single
#' object. Aggregation proceeds sequentially from
#' the first to the last element of \code{group_structure}: first, a mean and covariance matrix will be obtained for all combinations
#' of grouping variables; then the first grouping variable is collapsed over, creating an average for all combinations of the remaining
#' grouping variables; and so on. For this process to have a reasonable outcome, it is important that the grouping variables are hierarchically
#' organized from more specific to less specific: the second group variable should be a superset of the first group variable; the third group
#' variable should be a superset of the second group variable; etc.
#'
#' @param x An MVG, MVG_ideal_observer, NIW_belief, or NIW_ideal_adaptor object.
#' @param group_structure The group structure that will be used for aggregation.
#'
#' @return The aggregated object.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
aggregate_models_by_group_structure = function(
  x,
  group_structure = NULL
) {
  assert_that(all(is.character(group_structure) | is.symbol(group_structure)),
              msg = "Group structure must be a vector characters or symbols.")
  if (is.character(group_structure)) group_structure = syms(group_structure)
  assert_that(all(is.character(aggregate_to_group) | is.symbol(aggregate_to_group), length(aggregate_to_group) == 1),
              msg = "aggregate_to_group must be a single character or symbol.")
  if (is.character(aggregate_to_group)) aggregate_to_group = sym(aggregate_to_group)
  assert_that(as_name(aggregate_to_group) %in% as_name(group_structure),
              msg = "aggregate_to_group must be contained in group_structure.")

  aggregate_what = case_when(
    is.NIW_ideal_adaptor(x) ~ c("kappa", "nu", "m", "S", "prior", "lapse_rate", "lapse_bias", "Sigma_noise"),
    is.NIW_belief(x) ~ c("kappa", "nu", "m", "S"),
    is.MVG_ideal_observer(x) ~ c("mu", "Sigma", "prior", "lapse_rate", "lapse_bias", "Sigma_noise"),
    is.MVG(x) ~ c("mu", "Sigma"),
    T ~ NA_character_)

  while(length(group_structure) > 1) {
    group_structure = group_structure[2:length(group_structure)]
    x %<>%
      group_by(!!! group_structure, !! category) %>%
      summarise_at(aggregate_what, ~ list(reduce(.x, `+`) / length(.x)))
  }
}



#' Make multivariate Gaussian exposure data.
#'
#' Returns a tibble of observations drawn from multivariate Gaussians, with one observation per row. Each row
#' provides the category label and cue values. If \code{keep.parameters = T} then the parameters (\code{N, mean, sigma})
#' are also returned.
#'
#' The input is expected to be lists/vectors of parameters with the n-th element of each list/vector specifying the
#' category label, number of observations, \code{mu}, and \code{Sigma} of the n-th Gaussian.
#'
#' @param Ns Integer vector, with each number specifying the number of observations to be drawn from the corresponding
#' Gaussian.
#' @param mus List of mean vectors, each specifying the mean of a multivariate Gaussian.
#' @param sigmas List of covariance matrices, each specifying the covariance of a multivariate Gaussian.
#' @param category.labels Character vector of category names, each specifying the category label of a multivariate Gaussian. If \code{NULL}
#' (default) then Gaussians will be numbered from 1:N.
#' @param cue.labels Character vector of cue names. If \code{NULL} (default) then the cues will be numbered cue1, cue2, ...
#' @param x Alternatively to providing mus, sigma, category, and cue labels, one can also specify an \code{\link{MVG}},
#' \code{\link{MVG_ideal_observer}}, \code{\link{NIW_belief}}, or \code{\link{NIW_ideal_adaptor}}
#' object, which contains all information (except for the Ns).
#' @param randomize.order Should the order of the data be randomized? (default: FALSE) This won't affect the final outcome of
#' NIW belief updating, but it will change the incremental updates (and thus, for example, visualizations of the update process).
#' @param keep.input_parameters Should the parameters handed to this function be included in the output? (default: FALSE)
#'
#' @return A tibble.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
make_MVG_exposure_data = function(
  Ns, mus, Sigmas,
  category.labels = NULL,
  cue.labels = NULL,
  x = NULL,
  randomize.order = F,
  keep.input_parameters = F
) {
  if (!is.null(x)) {
    if (is.MVG(x) | is.MVG_ideal_observer(x)) {
      return(make_MVG_exposure_data(
        Ns = Ns,
        mus = x$mu,
        Sigmas = x$Sigma,
        category.labels = get_category_labels_from_MVG(x),
        cue.labels = get_cue_labels_from_MVG(x),
        randomize.order = randomize.order,
        keep.input_parameters = keep.input_parameters))
    } else if (is.NIW_belief(x) | is.NIW_ideal_adaptor(x)) {
      return(make_MVG_exposure_data(
        Ns = Ns,
        mus = x$m,
        Sigmas = get_expected_Sigma_from_S(x$S, x$nu),
        category.labels = get_category_labels_from_model(x),
        cue.labels = get_cue_labels_from_model(x),
        randomize.order = randomize.order,
        keep.input_parameters = keep.input_parameters))
    }
  } else {
    assert_that(!is.null(mus), !is.null(Sigmas))
    assert_that(is.null(category.labels) | length(mus) == length(category.labels),
                msg = "Number of category labels mismatch number of mus.")
    assert_that(is.null(cue.labels) | length(mus[[1]]) == length(cue.labels),
                msg = "Number of cue labels mismatches dimensionality of mus.")

    if (is.null(category.labels)) category.labels = 1:length(mus)
    if (is.null(cue.labels)) cue.labels = paste0("cue", 1:length(mus[[1]]))

    x = tibble(category = category.labels, n = Ns, mu = mus, Sigma = Sigmas) %>%
      mutate(data = pmap(.l = list(n, mu, Sigma), .f = rmvnorm)) %>%
      mutate(data = map(data, ~ .x %>% as.data.frame() %>% rename_all(~ cue.labels))) %>%
      unnest(data)

    if (randomize.order)
      x = sample_frac(x, 1)

    if (keep.input_parameters) return(x) else return(x %>% select(category, everything(), -c(n, mu, Sigma)))
  }
}
