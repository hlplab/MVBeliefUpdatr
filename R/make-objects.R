#' Make exemplar models from data.
#'
#' Constructs an exemplar model for all categories found in the data.
#'
#' @param data The tibble or data.frame from which to construct the exemplar model.
#' @param group Optionally, a vector of one or more grouping variables. If group is not NULL, one MVG or
#' ideal observers will be derived for each level of \code{group}. (default: NULL)
#' @param category Name of variable in \code{data} that contains the category information. (default: "category")
#' @param cues Name(s) of variables in \code{data} that contain the cue information.
#' @param sim_function Similarity function that is used to calculate exemplar-to-examplar similarity. Defaults
#' to \code{exp(-(x - y)^2)}.
#' @param prior Optionally specify a prior probability for each category (in each group). (default: a uniform
#' prior over all categories).
#' @param lapse_rate Optionally specify a lapse rate. (default: \code{NA})
#' @param lapse_bias Optionally specify a lapse bias. (default: \code{NA})
#' @param Sigma_noise Optionally specify a (multivariate Gaussian) covariance matrix of perceptual noise. This
#' argument will be ignored if NULL. (default: NULL)
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
make_exemplars_from_data = function(
    data,
    group = NULL,
    category = "category",
    cues,
    sim_function = function(x, y) {
      j <- 2
      k <- 1
      distance <- (x - y)^j
      similarity <- exp(-distance * k) },
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

  if (verbose) if (!is.null(group))
    message(
      paste("Group specified. Making one exemplar cloud for each unique value of group.",
            paste(unique(data$group), collapse = ",")))

  model <-
    data %>%
    select(!! category, !!! cues, !!! group) %>%
    mutate(cues = pmap(list(!!! cues),
                       .f = function(...) {
                         x = c(...)
                         names(x) = as.character(cues)
                         return(x) })) %>%
    { if (is.null(group)) group_by(., !! category) else group_by(., !!! group, !! category) } %>%
    select(-c(!!! cues)) %>%
    nest(exemplars = cues) %>%
    mutate(sim_function = list(.env$sim_function))

  model %<>%
    rename(category = !! category) %>%
    mutate(category = factor(category)) %>%
    ungroup()

  if (!is.exemplars(model, group = group, verbose = verbose))
    warning("Something went wrong. The returned object is not a set of exemplars. Try again with verbose = T?")

  return(model)
}


#' @export
#' @rdname make_exemplars_from_data
make_exemplar_model_from_data = function(
    data,
    group = NULL,
    category = "category",
    cues,
    sim_function = function(x, y) {
      j <- 2
      k <- 1
      distance <- (x - y)^j
      similarity <- exp(-distance * k) },
    ...,
    verbose = F
) {
  model <-
    data %>%
    make_exemplars_from_data(group = group, category = category, cues = cues, sim_function = sim_function, verbose = verbose)

  model %<>% lift_exemplars_to_exemplar_model(group = group, ..., verbose = verbose)

  return(model)
}


#' Make multivariate Gaussian ideal observer(s) from data.
#'
#' Constructs an \code{\link[=is.MVG]{MVG}} or \code{\link[=is.MVG_ideal_observer]{MVG_ideal_observer}} object with category
#' information for all categories found in the data.
#'
#' @param data The tibble or data.frame from which to construct the MVG or MVG ideal observer object.
#' @param group Optionally, a vector of one or more grouping variables. If group is not NULL, one MVG or
#' ideal observers will be derived for each level of \code{group}. (default: NULL)
#' @param category Name of variable in \code{data} that contains the category information. (default: "category")
#' @param cues Name(s) of variables in \code{data} that contain the cue information.
#' @param prior Optionally specify a prior probability for each category (in each group). (default: a uniform
#' prior over all categories).
#' @param lapse_rate Optionally specify a lapse rate. (default: \code{NA})
#' @param lapse_bias Optionally specify a lapse bias. (default: \code{NA})
#' @param Sigma_noise Optionally specify a (multivariate Gaussian) covariance matrix of perceptual noise. This
#' argument will be ignored if NULL. (default: NULL)
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

  if (verbose) if (!is.null(group))
    message(
      paste("Group specified. Making one MVG for each unique value of group.",
            paste(unique(data$group), collapse = ",")))

  model <-
    data %>%
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

  model %<>%
    rename(category = !! category) %>%
    mutate(category = factor(category)) %>%
    ungroup()

  if (!is.MVG(model, group = group, verbose = verbose))
    warning("Something went wrong. The returned object is not an MVG. Try again with verbose = T?")

  return(model)
}

#' @export
#' @rdname make_MVG_from_data
make_MVG_ideal_observer_from_data = function(
  data,
  group = NULL,
  category = "category",
  cues,
  ...,
  verbose = F
) {
  model <-
    data %>%
    make_MVG_from_data(group = group, category = category, cues = cues, verbose = verbose)

  model %<>% lift_MVG_to_MVG_ideal_observer(group = group, ..., verbose = verbose)

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
#' @param verbose If true provides more information. (default: FALSE)
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
  verbose = F
) {
  assert_that(all(is.numeric(kappa), is.numeric(nu), !is.na(kappa), !is.na(nu)))
  assert_that(nu > length(cues) + 1,
              msg = paste0("nu must be larger than dimensionality of cues + 1 (>", length(cues) + 1, ")."))

  # if (is.character(group)) group = syms(group)
  # if (is.character(category)) category = sym(category)
  # if (is.character(cues)) cues = syms(cues)

  model <-
    data %>%
    make_MVG_from_data(group = group, category = category, cues = cues, verbose = verbose) %>%
    rename(m = mu, S = Sigma)

  message("S is set so that the expected category covariance matrix Sigma matches the category covariance in the sample (given nu). ",
          "It might be safer to fit an Inverse-Wishart distribution to the entire set of covariance matrices.")
  model %<>%
    mutate(
      kappa = kappa,
      nu = nu,
      S = get_S_from_expected_Sigma(S, nu)) %>%
    ungroup()

  if (!is.NIW_belief(model, group = group, verbose = verbose))
    warning("Something went wrong. The returned object is not an NIW belief. Try again with verbose = T?")

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
  ...,
  verbose = F
) {
  model <-
    data %>%
    make_NIW_belief_from_data(group = group, category = category, cues = cues, kappa = kappa, nu = nu, verbose = verbose)

  model %<>% lift_NIW_belief_to_NIW_ideal_adaptor(group = group, ..., verbose = verbose)

  return(model)
}


#' Example NIW priors.
#'
#' @export
make_NIW_example = function(example = 1) {
  if (example == 1) {
    message("An example belief for two categories in a 2D cue continuum that differ in means and correlatation, but not standard deviations.")
    tibble(
      category = c("A", "B"),
      kappa = 10,
      nu = 30,
      m = list(c("cue1" = -2, "cue2" = -2), c("cue1" = 2, "cue2" = 2)),
      S = list(matrix(c(1, .3, .3, 1), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2"))),
               matrix(c(1, -.3, -.3, 1), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2")))),
      lapse_rate = .05
    ) %>%
      mutate(category = factor(category))
  } else if (example == 2) {
    message("An example belief for two categories in a 2D cue continuum that differ in means and correlatation, but not standard deviations.
            Same as Example 1, but with independent perceptual noise along both cue dimensions.")
    tibble(
      category = c("A", "B"),
      kappa = 10,
      nu = 30,
      m = list(c("cue1" = -2, "cue2" = -2), c("cue1" = 2, "cue2" = 2)),
      S = list(matrix(c(1, .3, .3, 1), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2"))),
               matrix(c(1, -.3, -.3, 1), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2")))),
      lapse_rate = .05,
      Sigma_noise = list(matrix(c(1, 0, 0, .25), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2"))))
    ) %>%
      mutate(category = factor(category))
  }
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
#' @param lapse_rate Optionally specify a lapse rate. (default: 0)
#' @param lapse_bias Optionally specify a lapse bias for each category. (default: uniform bias for all categories)
#' @param Sigma_noise Optionally specify a (multivariate Gaussian) covariance matrix of perceptual noise. This argument will be
#' ignored if NULL. (default: NULL)
#' @param keep.category_parameters Should categories' mu and Sigma be included in the output (in addition to m
#' and S of the prior)? (default: FALSE)
#' @param verbose If true provides more information. (default: FALSE)
#'
#' @return A tibble that is an NIW_belief or NIW_ideal_adaptor object.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
#' @rdname lift_likelihood_to_model
lift_likelihood_to_model <- function(
  x,
  group = NULL,
  category = "category",
  prior = rep(1 / get_nlevels_of_category_labels_from_model(x), get_nlevels_of_category_labels_from_model(x)),
  lapse_rate = 0,
  lapse_bias = rep(1 / get_nlevels_of_category_labels_from_model(x), get_nlevels_of_category_labels_from_model(x)),
  Sigma_noise = NULL,
  verbose = F
) {
  if (is.character(group)) group = syms(group)
  if (is.character(category)) category = sym(category)
  assert_that(all(is.numeric(lapse_rate), is.numeric(lapse_bias), is.numeric(prior)),
              msg = "Category prior, lapse rate, and lapse bias must be numeric.")

  category_levels <- get_category_labels_from_model(x)
  n.cat <- get_nlevels_of_category_labels_from_model(x)
  if (!is.null(prior)) {
    assert_that(length(prior) == n.cat,
              msg = paste("Category prior must have as many elements as there are categories. Has", length(prior), "instead of needed", n.cat))
    if (!is.null(names(prior))) {
        assert_that(all(names(prior) == category_levels),
                    msg = paste("Names of category priors must match levels of", category, "in x."))
    } else {
      message(paste("Category priors were not named. Assuming that priors are provided in alphabetic order of", category, "in x."))
      names(prior) <- category_levels
    }
  }
  if (!is.null(lapse_bias)) {
    assert_that(length(lapse_bias) == n.cat,
                msg = paste("Lapse_bias must have as many elements as there are categories. Has", length(lapse_bias), "instead of needed", n.cat))
    if (!is.null(names(lapse_bias))) {
      assert_that(all(names(lapse_bias) == category_levels),
                  msg = paste("Names of lapse biases must match levels of", category, "in x."))
    } else {
      message(paste("Lapse biases were not named. Assuming that lapse biases are provided in alphabetic order of", category, "in x."))
      names(lapse_bias) <- category_levels
    }
  }

  if (all(is.na(lapse_rate) | is.null(lapse_rate))) {
    lapse_rate <- 0
  }

  x %<>%
    group_by(!!! group) %>%
    mutate(
      prior = .env$prior[as.character(!! sym(category))],
      lapse_rate = .env$lapse_rate,
      lapse_bias = .env$lapse_bias[as.character(!! sym(category))],
      Sigma_noise = list(.env$Sigma_noise))

  if (!is.model(x, group = group, verbose = verbose))
    warning("Something went wrong. The returned object is not a model. Try again with verbose = T?")

  return(x)
}


#' @export
#' @rdname lift_likelihood_to_model
lift_exemplars_to_exemplar_model <- function(
    x,
    group = NULL,
    prior = rep(1 / (n.cat <- get_nlevels_of_category_labels_from_model(x)), n.cat),
    lapse_rate = 0,
    lapse_bias = rep(1 / (n.cat <- get_nlevels_of_category_labels_from_model(x)), n.cat),
    Sigma_noise = NULL,
    verbose = F
) {
  x %<>% lift_likelihood_to_model(group = group, prior = prior, lapse_rate = lapse_rate, lapse_bias = lapse_bias, Sigma_noise = Sigma_noise)
  if (!is.null(first(x$Sigma_noise))) {
    assert_that(!is.null(dimnames(first(x$Sigma_noise))),
                msg = "Sigma_noise = must have non-NULL dimnames.")
    assert_that(map2(dimnames(Sigma_noise), dimnames(first(x$Sigma)), ~ .x == .y) %>% reduce(c) %>% all(),
                msg = "The dimnames of Sigma_noise must match the cue names in the list of exemplars.")
  }

  if (!is.exemplar_model(x, group = group, verbose = verbose))
    warning("Something went wrong. The returned object is not an exemplar model. Try again with verbose = T?")

  return(x)
}


#' @export
#' @rdname lift_likelihood_to_model
lift_MVG_to_MVG_ideal_observer = function(
  x,
  group = NULL,
  category = "category",
  prior = rep(1 / (n.cat <- get_nlevels_of_category_labels_from_model(x)), n.cat),
  lapse_rate = 0,
  lapse_bias = rep(1 / (n.cat <- get_nlevels_of_category_labels_from_model(x)), n.cat),
  Sigma_noise = NULL,
  verbose = F
) {
  x %<>% lift_likelihood_to_model(group = group, category = category, prior = prior, lapse_rate = lapse_rate, lapse_bias = lapse_bias, Sigma_noise = Sigma_noise)
  if (!is.null(first(x$Sigma_noise))) {
    assert_that(all(dim(Sigma_noise) == dim(first(x$Sigma))),
                msg = "Sigma_noise must be a matrix of the same dimensionality as Sigma.")
    assert_that(!is.null(dimnames(first(x$Sigma_noise))),
                msg = "Sigma_noise = must have non-NULL dimnames.")
    assert_that(map2(dimnames(Sigma_noise), dimnames(first(x$Sigma)), ~ .x == .y) %>% reduce(c) %>% all(),
                msg = "The dimnames of Sigma_noise and Sigma must match.")
  }

  if (!is.MVG_ideal_observer(x, group = group, verbose = verbose))
    warning("Something went wrong. The returned object is not an MVG ideal observer. Try again with verbose = T?")

  return(x)
}

#' @export
#' @rdname lift_likelihood_to_model
lift_NIW_belief_to_NIW_ideal_adaptor = function(
  x,
  group = NULL,
  category = "category",
  prior = rep(1 / (n.cat <- get_nlevels_of_category_labels_from_model(x)), n.cat),
  lapse_rate = 0,
  lapse_bias = rep(1 / (n.cat <- get_nlevels_of_category_labels_from_model(x)), n.cat),
  Sigma_noise = NULL,
  verbose = F
) {
  x %<>% lift_likelihood_to_model(group = group, category = category, prior = prior, lapse_rate = lapse_rate, lapse_bias = lapse_bias, Sigma_noise = Sigma_noise)
  if (!is.null(first(x$Sigma_noise))) {
    assert_that(all(dim(Sigma_noise) == dim(first(x$S))),
                msg = "Sigma_noise must be a matrix of the same dimensionality as S.")
    assert_that(!is.null(dimnames(first(x$Sigma_noise))),
                msg = "Sigma_noise = must have non-NULL dimnames.")
    assert_that(map2(dimnames(Sigma_noise), dimnames(first(x$S)), ~ .x == .y) %>% reduce(c) %>% all(),
                msg = "The dimnames of Sigma_noise and S must match.")
  }

  if (!is.NIW_ideal_adaptor(x, group = group, verbose = verbose))
    warning("Something went wrong. The returned object is not an NIW ideal adaptor. Try again with verbose = T?")

  return(x)
}

#' Turn an MVG_ideal_observer into an NIW_ideal_adaptor
#'
#' Make an ideal adaptor out of an ideal observer, so that the *expected* category mean and category covariance
#' matrix of the ideal adaptor match the ideal observers category mean and category covariance matrix.
#'
#' @param x An MVG_ideal_observer object.
#' @param group Optionally, a grouping structure can be specified. If group structure is not NULL, one
#' NIW belief or ideal adaptor will be derived for each level of \code{group_structure}. (default: NULL)
#' @param category Name of variable in \code{data} that contains the category information. (default: "category")
#' @param kappa The strength of the beliefs over the category mean (pseudocounts). (default: same as nu)
#' @param nu The strength of the beliefs over the category covariance matrix (pseudocounts). (default: number of
#' cues + 2)
#' @param verbose If true provides more information. (default: FALSE)
#'
#' @return A tibble that is an NIW_ideal_adaptor object.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
lift_MVG_ideal_observer_to_NIW_ideal_adaptor = function(
  x,
  group = NULL,
  category = "category",
  kappa, nu,
  verbose = F
) {
  assert_that(!is.null(kappa), !is.null(nu),
              msg = "kappa and nu must be provided.")

  x %<>%
    rename(m = mu, S = Sigma) %>%
    mutate(
      kappa = kappa,
      nu = nu,
      S = get_S_from_expected_Sigma(S, nu))

  if (!is.NIW_ideal_adaptor(x, group = group, verbose = verbose))
    warning("Something went wrong. The returned object is not an NIW ideal adaptor. Try again with verbose = T?")

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
#' @param category Name of variable in \code{data} that contains the category information. (default: "category")
#' @param group_structure The group structure that will be used for aggregation.
#'
#' @return The aggregated object.
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
aggregate_models_by_group_structure = function(
  x,
  category = "category",
  group_structure = NULL
) {
  assert_that(all(is.character(group_structure) | is.symbol(group_structure)),
              msg = "Group structure must be a vector characters or symbols.")
  if (is.character(group_structure)) group_structure = syms(group_structure)
  # assert_that(all(is.character(group_structure) | is.symbol(group_structure), length(group_structure) == 1),
  #             msg = "aggregate_to_group must be a single character or symbol.")
  # if (is.character(aggregate_to_group)) aggregate_to_group = sym(aggregate_to_group)
  # assert_that(as_name(aggregate_to_group) %in% as_name(group_structure),
  #             msg = "aggregate_to_group must be contained in group_structure.")

  aggregate_what_into_means <-
    if (is.NIW_ideal_adaptor(x, category = category)) c("m", "prior", "lapse_rate", "lapse_bias", "kappa", "nu", "S", "Sigma_noise") else
      if (is.NIW_belief(x, category = category)) c("m", "kappa", "nu", "S") else
        if (is.MVG_ideal_observer(x, category = category)) c("mu", "prior", "lapse_rate", "lapse_bias", "Sigma", "Sigma_noise") else
          if (is.MVG(x, category = category)) c("mu", "Sigma") else NA_character_

  # Consider geommetric mean for some variables in the future:
  # (but for covs that would mean first getting the taus, then getting geommetric mean. then transforming
  # back into covs (since covs can be negative)). If implemented add the following to docu:
  # Note that geometric means will be used for count and variance variables (e.g., Sigma, S, Sigma_noise, kappa, nu),
  # whereas arithmetic means will be used for all other types of variables (e.g., mu, m, prior, lapse rate and bias).
  # aggregate_what_into_geometric_means <-
  #   if (is.NIW_ideal_adaptor(x, category = category)) c("kappa", "nu", "S", "Sigma_noise") else
  #     if (is.NIW_belief(x, category = category)) c("kappa", "nu", "S") else
  #       if (is.MVG_ideal_observer(x, category = category)) c("Sigma", "Sigma_noise") else
  #         if (is.MVG(x, category = category)) c("Sigma") else NA_character_

  while(length(group_structure) > 0) {
    group_structure = group_structure[-1]
    x %<>%
      group_by(!!! syms(group_structure), !! sym(category)) %>%
      summarise(
        across(aggregate_what_into_means, ~ list(reduce(.x, `+`) / length(.x))))
    # Consider geommetric mean for some variables in the future:
    # across(aggregate_what_into_geometric_means, ~ list(exp(reduce(log(.x), `+`) / length(.x))))
  }

  return(x)
}


#' Sample multivariate Gaussian exposure data.
#'
#' Returns a tibble of observations drawn from multivariate Gaussians, with one observation per row. Each row
#' provides the category label and cue values. If \code{keep.input_parameters = T} then the parameters (\code{N, mean, sigma})
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
#' @param model \code{\link{MVG}}, \code{\link{MVG_ideal_observer}}, \code{\link{NIW_belief}}, or \code{\link{NIW_ideal_adaptor}} object.
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
sample_MVG_data = function(
  Ns, mus, Sigmas,
  category.labels = NULL,
  cue.labels = NULL,
  randomize.order = F,
  keep.input_parameters = F
) {
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


#' @rdname sample_MVG_data
#' @export
sample_MVG_data_from_model = function(
  Ns,
  model = NULL,
  randomize.order = F,
  keep.input_parameters = F
) {
  if (is.MVG(model) | is.MVG_ideal_observer(model)) {
    return(sample_MVG_data(
      Ns = Ns,
      mus = model$mu,
      Sigmas = model$Sigma,
      category.labels = get_category_labels_from_model(model),
      cue.labels = get_cue_labels_from_model(model),
      randomize.order = randomize.order,
      keep.input_parameters = keep.input_parameters))
  } else if (is.NIW_belief(model) | is.NIW_ideal_adaptor(model)) {
    return(make_MVG_data(
      Ns = Ns,
      mus = model$m,
      Sigmas = get_expected_Sigma_from_S(model$S, model$nu),
      category.labels = get_category_labels_from_model(model),
      cue.labels = get_cue_labels_from_model(model),
      randomize.order = randomize.order,
      keep.input_parameters = keep.input_parameters))
  } else
    message(paste("Unrecognized model class:", class(model)))
}
