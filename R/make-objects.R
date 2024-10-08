#' @importFrom dplyr rename_all
NULL

#' Example exemplar model.
#'
#' @rdname example-model
#' @export
example_exemplar_model <- function(
    example = 1,
    categories = c("/b/", "/p/"),
    cues = c("VOT", "f0_semitones", "vowel_duration")[1:example],
    prior = rep(1 / length(categories), length(categories)),
    lapse_rate = .05,
    lapse_bias = prior,
    Sigma_noise = matrix(rep(0, length(cues)^2), nrow = length(cues), dimnames = list(cues, cues)),
    verbose = T
) {
  if (example %in% 1:3) {
    if (verbose)
      message("An example exemplar model based on the data from Chodroff & Wilson (2018).\n",
            "Using ", length(categories), " categories and ", length(cues), " cue(s). ",
            "Lapse is set to ", lapse_rate, " with lapse biases ", paste(lapse_bias, collapse = ", "), ".")

    suppressMessages(
      make_MVG_ideal_observer_from_data(
        data = ChodroffWilson2018 %>% filter(category %in% categories),
        category = "category",
        cues = cues,
        prior = prior,
        lapse_rate = lapse_rate,
        lapse_bias = lapse_bias,
        Sigma_noise = Sigma_noise)) %>%
      sample_MVG_data_from_model(Ns = 50) %>%
      make_exemplar_model_from_data(
        cues = cues,
        prior = prior,
        lapse_rate = lapse_rate,
        lapse_bias = lapse_bias,
        Sigma_noise = NULL) %>%
      mutate(across(category, factor))

  } else if (example == 5) {
    if (verbose)
      message(paste("An example exemplar model with two categories and two cues.\n",
            "Lapse is set to ", lapse_rate, " with lapse biases ", paste(lapse_bias, collapse = ", "), "."))
    suppressMessages(example_MVG_ideal_observer(4)) %>%
      sample_MVG_data_from_model(Ns = 50) %>%
      make_exemplar_model_from_data(
        cues = c("cue1", "cue2"),
        prior = prior,
        lapse_rate = lapse_rate,
        lapse_bias = lapse_bias,
        Sigma_noise = NULL) %>%
      mutate(across(category, factor))
  }
}


#' Example MVG ideal observers.
#'
#' @importFrom dplyr %>%
#' @importFrom purrr map
#'
#' @rdname example-model
#' @export
example_MVG_ideal_observer <- function(
    example = 1,
    categories = c("/b/", "/p/"),
    cues = c("VOT", "f0_semitones", "vowel_duration")[1:example],
    prior = rep(1 / length(categories), length(categories)),
    lapse_rate = .05,
    lapse_bias = prior,
    Sigma_noise = matrix(rep(0, length(cues)^2), nrow = length(cues), dimnames = list(cues, cues)),
    verbose = T
) {
  if (example %in% 1:3) {
    if (verbose)
      message("An example MVG ideal observer based on the data from Chodroff & Wilson (2018).\n",
            "Using ", length(categories), " categories and ", length(cues), " cue(s). ",
            "Lapse is set to ", lapse_rate, " with lapse biases ", paste(lapse_bias, collapse = ", "), ".")

    make_MVG_ideal_observer_from_data(
      data = ChodroffWilson2018 %>% filter(category %in% categories),
      category = "category",
      cues = cues,
      prior = prior,
      lapse_rate = lapse_rate,
      lapse_bias = lapse_bias,
      Sigma_noise = Sigma_noise)
  } else if (example == 5) {
    if (verbose)
      message("An example MVG ideal observer with two categories and two cues.\n",
            "Lapse is set to ", lapse_rate, " with lapse biases ", paste(lapse_bias, collapse = ", "), ".")

    tibble(
      category = c("A", "B"),
      mu = list(c("cue1" = -2, "cue2" = -2), c("cue1" = 2, "cue2" = 2)),
      Sigma = list(
        matrix(c(3, 2.4, 2.4, 3), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2"))),
        matrix(c(3, -2.4, -2.4, 3), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2")))),
      prior = prior,
      lapse_rate = lapse_rate,
      lapse_bias = lapse_bias,
      Sigma_noise = list(
        matrix(rep(0, 4), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2"))),
        matrix(rep(0, 4), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2"))))) %>%
      mutate(across(category, factor))
  }
}


#' Example NIW ideal adaptors.
#'
#' @rdname example-model
#' @export
example_NIW_ideal_adaptor <- function(
    example = 1,
    categories = c("/b/", "/p/"),
    cues = c("VOT", "f0_semitones", "vowel_duration")[1:example],
    prior = rep(1 / length(categories), length(categories)),
    lapse_rate = .05,
    lapse_bias = prior,
    Sigma_noise = matrix(rep(0, length(cues)^2), nrow = length(cues), dimnames = list(cues, cues)),
    verbose = T
) {
  x <-
    example_MVG_ideal_observer(
      example = example,
      categories = categories,
      cues = cues,
      prior = prior,
      lapse_rate = lapse_rate,
      lapse_bias = lapse_bias,
      Sigma_noise = Sigma_noise,
      verbose = verbose)

  if (example %in% 1:3) {
    if (verbose)
      message("Lifted to an ideal adaptor with kappa = 10 and nu = 30.")

     x %>%
      lift_MVG_ideal_observer_to_NIW_ideal_adaptor(
        kappa = 10,
        nu = 30)
  } else if (example == 5) {
    message("Lifted to an ideal adaptor with kappa = 30 and nu = 10.")

    x %>%
      lift_MVG_ideal_observer_to_NIW_ideal_adaptor(
        kappa = 30,
        nu = 10)
  }
}


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
#' @export
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
    warning("NOTE: The returned object does not appear to be a set of exemplars. For more information, try again with verbose = T.")

  return(model)
}


#' @export
#' @rdname make_exemplars_from_data
make_exemplar_model_from_data = function(
    data,
    group = NULL,
    category = "category",
    cues,
    sim_function = function(x, y, j = 2, k = 1) {
      # If j = 1, this uses Manhattan city block distances
      # If j = 2, this uses Euclidean distances
      distance <- sum((x - y)^j)^(1/j)
      similarity <- exp(-distance * k) },
    # Alternatively, users might specify Mahanalobis distance
    # function(x, y, center, cov, k = 1) {
    #
    #   distance <- mahanalobis(x, center, cov)^-2 ...
    #   similarity <- exp(-distance * k) }
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
#' information for all categories found in the data. Currently, this functions does nothing fancy. It simply gets the mean
#' and covariance matrix of the cues for each category and group from the data. No cross-validation or other measures against
#' overfitting are implemented, though it is recommended that such methods are applied.
#'
#' Alternative approaches include the use of `brms::brm()` to fit a multivariate Normal model to the data. While this approach
#' allows fitting of both the means and variances of each category (including for hierarchically organized grouped data), it
#' currently does not provide a way to model category-specific correlations (or covariances) between cues. Instead, the
#' approach implemented in `brms` only models correlation at the population-level (residual correlations).
#'
#' Yet another alternative would be to write a separate `Stan` program specifically for this purpose. However, while this is
#' relatively straightforward for data from a single talker, a hierarchical model for grouped data essentially requires an
#' extension of the multivariate model approach implemented in `brms` and described in the preceding paragraph.
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
#'
#' @export
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
    warning("NOTE: The returned object does not appear to be an MVG. For more information, try again with verbose = T.")

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
#' @export
make_NIW_belief_from_data <- function(
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
    rename(all_of(c(m = "mu", S = "Sigma")))

  message("S is set so that the expected category covariance matrix Sigma matches the category covariance in the sample (given nu). ",
          "It might be safer to fit an Inverse-Wishart distribution to the entire set of covariance matrices.")
  model %<>%
    mutate(
      kappa = .env$kappa,
      nu = .env$nu,
      S = get_S_from_expected_Sigma(.data$S, .env$nu)) %>%
    ungroup()

  if (!is.NIW_belief(model, group = group, verbose = verbose))
    warning("NOTE: The returned object does not appear to be an NIW belief. For more information, try again with verbose = T.")

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


#' Turn an MVG/NIW_belief object into an ideal observer/adaptor
#'
#' Make an ideal observer or adaptor out of an MVG or NIW_belief object, respectively, by providing the missing
#' information about the prior, lapse rate, bias, and perceptual noise (if any).
#'
#' @param x Either an MVG or NIW_belief object.
#' @param group Optionally, a grouping structure can be specified. If group structure is not NULL, one
#' NIW belief or ideal adaptor will be derived for each level of \code{group_structure}. (default: NULL)
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
#' @export
#' @rdname lift_likelihood_to_model
lift_likelihood_to_model <- function(
  x,
  group = NULL,
  prior = rep(1 / get_nlevels_of_category_labels_from_model(x), get_nlevels_of_category_labels_from_model(x)),
  lapse_rate = 0,
  lapse_bias = rep(1 / get_nlevels_of_category_labels_from_model(x), get_nlevels_of_category_labels_from_model(x)),
  Sigma_noise = NULL,
  verbose = F
) {
  if (is.character(group)) group = syms(group)
  assert_that(all(is.numeric(lapse_rate), is.numeric(lapse_bias), is.numeric(prior)),
              msg = "Category prior, lapse rate, and lapse bias must be numeric.")

  category_levels <- get_category_labels_from_model(x)
  n.cat <- get_nlevels_of_category_labels_from_model(x)
  if (!is.null(prior)) {
    assert_that(length(prior) == n.cat,
              msg = paste("Category prior must have as many elements as there are categories. Has", length(prior), "instead of needed", n.cat))
    if (!is.null(names(prior))) {
        assert_that(all(names(prior) == category_levels),
                    msg = paste("Names of category priors must match levels of category in x."))
    } else if (!all(prior == first(prior))) {
      # If priors are the same there's no need for this message. This also prevents that the message is
      # displayed when the default uniform prior is used.
      message(paste("Category priors were not named. Assuming that priors are provided in alphabetic order of category in x."))
    }
    names(prior) <- category_levels
  }
  if (!is.null(lapse_bias)) {
    assert_that(length(lapse_bias) == n.cat,
                msg = paste("Lapse_bias must have as many elements as there are categories. Has", length(lapse_bias), "instead of needed", n.cat))
    if (!is.null(names(lapse_bias))) {
      assert_that(all(names(lapse_bias) == category_levels),
                  msg = paste("Names of lapse biases must match levels of category in x."))
    } else if (!all(lapse_bias == first(lapse_bias))) {
      # If biases are the same there's no need for this message. This also prevents that the message is
      # displayed when the default uniform biases are used.
      message(paste("Lapse biases were not named. Assuming that lapse biases are provided in alphabetic order of category in x."))
    }
    names(lapse_bias) <- category_levels
  }

  if (all(is.na(lapse_rate) | is.null(lapse_rate))) {
    lapse_rate <- 0
  }

  assert_that(is.null(Sigma_noise) | is.matrix(Sigma_noise),
              msg = "Sigma_noise must be NULL or a matrix (noise is independent of, and thus constant across, categories).")
  if (!is.null(Sigma_noise)) {
    assert_that(all(dim(Sigma_noise) == rep(get_cue_dimensionality_from_model(x), 2)),
                msg = paste("If not NULL, Sigma_noise must match the dimensionality of other parameters in the model (here: a",
                            get_cue_dimensionality_from_model(x), "x", get_cue_dimensionality_from_model(x), " matrix)."))
    assert_that(!is.null(first(dimnames(Sigma_noise))),
                msg = "If not NULL, Sigma_noise = must have non-NULL dimnames.")
    assert_that(map(dimnames(Sigma_noise), ~ .x == get_cue_labels_from_model(x)) %>% reduce(all),
                msg = "The dimnames of Sigma_noise must match those of the model to lift.")
  }

  x %<>%
    group_by(!!! group) %>%
    mutate(
      prior = .env$prior[as.character(.data$category)],
      lapse_rate = .env$lapse_rate,
      lapse_bias = .env$lapse_bias[as.character(.data$category)],
      Sigma_noise = list(.env$Sigma_noise))

  if (!is.MVBU_model(x, group = group, verbose = verbose))
    warning("NOTE: The returned object does not appear to be an MVBeliefUpdatr model. For more information, try again with verbose = T.")

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
  assert_that(is.exemplars(x, group = group, verbose = verbose))

  x %<>% lift_likelihood_to_model(group = group, prior = prior, lapse_rate = lapse_rate, lapse_bias = lapse_bias, Sigma_noise = Sigma_noise)

  if (!is.exemplar_model(x, group = group, verbose = verbose))
    warning("NOTE: The returned object does not appear to be an exemplar model. For more information, try again with verbose = T.")

  return(x)
}


#' @export
#' @rdname lift_likelihood_to_model
lift_MVG_to_MVG_ideal_observer = function(
  x,
  group = NULL,
  prior = rep(1 / (n.cat <- get_nlevels_of_category_labels_from_model(x)), n.cat),
  lapse_rate = 0,
  lapse_bias = rep(1 / (n.cat <- get_nlevels_of_category_labels_from_model(x)), n.cat),
  Sigma_noise = NULL,
  verbose = F
) {
  assert_that(is.MVG(x, group = group, verbose = verbose))

  x %<>% lift_likelihood_to_model(group = group, prior = prior, lapse_rate = lapse_rate, lapse_bias = lapse_bias, Sigma_noise = Sigma_noise)

  if (!is.MVG_ideal_observer(x, group = group, verbose = verbose))
    warning("NOTE: The returned object does not appear to be an MVG ideal observer. For more information, try again with verbose = T.")

  return(x)
}

#' @export
#' @rdname lift_likelihood_to_model
lift_NIW_belief_to_NIW_ideal_adaptor <- function(
  x,
  group = NULL,
  prior = rep(1 / (n.cat <- get_nlevels_of_category_labels_from_model(x)), n.cat),
  lapse_rate = 0,
  lapse_bias = rep(1 / (n.cat <- get_nlevels_of_category_labels_from_model(x)), n.cat),
  Sigma_noise = NULL,
  verbose = F
) {
  assert_that(is.NIW_belief(x, group = group, verbose = verbose))

  x %<>% lift_likelihood_to_model(group = group, prior = prior, lapse_rate = lapse_rate, lapse_bias = lapse_bias, Sigma_noise = Sigma_noise)

  if (!is.NIW_ideal_adaptor(x, group = group, verbose = verbose))
    warning("NOTE: The returned object does not appear to be an NIW ideal adaptor. For more information, try again with verbose = T.")

  return(x)
}

#' @export
#' @rdname lift_likelihood_to_model
lift_MVG_ideal_observer_to_NIW_ideal_adaptor <- function(
  x,
  group = NULL,
  kappa, nu,
  verbose = F
) {
  assert_that(is.MVG_ideal_observer(x, group = group, verbose = verbose))
  assert_that(!is.null(kappa), !is.null(nu),
              msg = "kappa and nu must be provided.")

  x %<>%
    rename(all_of(c(m = "mu", S = "Sigma"))) %>%
    mutate(
      kappa = .env$kappa,
      nu = .env$nu,
      S = get_S_from_expected_Sigma(.data$S, .env$nu))

  if (!is.NIW_ideal_adaptor(x, group = group, verbose = verbose))
    warning("NOTE: The returned object does not appear to be an NIW ideal adaptor. For more information, try again with verbose = T.")

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
#' @seealso TBD
#' @keywords TBD
#' @export
aggregate_models_by_group_structure = function(
  x,
  group_structure = NULL
) {
  assert_that(all(is.character(group_structure)),
              msg = "Group structure must be a vector of characters.")
  assert_that(all(group_structure %in% names(x)),
              msg = "All variables in group_structure must be contained in the x.")

  x_names <- setdiff(names(x), group_structure)
  # Consider geometric mean for some variables in the future:
  # across(aggregate_what_into_geometric_means, ~ list(exp(reduce(log(.x), `+`) / length(.x))))
  while(length(group_structure) > 0) {
    group_structure = group_structure[-1]
    x %<>%
      group_by(!!! syms(group_structure), category) %>%
      summarise(
        across(
          intersect(names(x), c("kappa", "nu", "prior", "lapse_rate", "lapse_bias")),
          ~ mean(.x, na.rm = T)),
        across(
          intersect(names(x), c("m", "mu", "S", "Sigma", "Sigma_noise")),
          ~ if (any(unlist(is.null(.x)))) { NULL } else { list(reduce(.x, `+`) / length(.x)) } ))
  }

  x %<>% relocate(!!! syms(x_names))
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
#'
#' @importFrom rlang .data
#' @importFrom tidyselect everything
#' @importFrom dplyr slice_sample
#' @export
sample_MVG_data <- function(
  Ns, mus, Sigmas,
  category.labels = NULL,
  cue.labels = NULL,
  randomize.order = F,
  keep.input_parameters = F
) {
  # Binding variables that RMD Check gets confused about otherwise
  # (since they are in non-standard evaluations)
  category <- n <- mu <- Sigma <- NULL

  assert_that(!is.null(mus), !is.null(Sigmas))
  assert_that(is.null(category.labels) | length(mus) == length(category.labels),
              msg = "Number of category labels mismatch number of mus.")
  assert_that(is.null(cue.labels) | length(mus[[1]]) == length(cue.labels),
              msg = "Number of cue labels mismatches dimensionality of mus.")

  if (is.null(category.labels)) category.labels = 1:length(mus)
  if (is.null(cue.labels)) cue.labels = paste0("cue", 1:length(mus[[1]]))

  x <-
    tibble("category" = category.labels, "n" = Ns, "mu" = mus, "Sigma" = Sigmas) %>%
    mutate(data = pmap(.l = list(.data$n, .data$mu, .data$Sigma), .f = rmvnorm)) %>%
    mutate(data = map(.data$data, ~ .x %>% as.data.frame() %>% rename_all(~ cue.labels))) %>%
    unnest(data)

  if (randomize.order)
    x %<>% slice_sample(prop = 1, replace = FALSE)

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
    return(sample_MVG_data(
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
