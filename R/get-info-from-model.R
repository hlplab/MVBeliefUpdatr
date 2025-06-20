#' Functions that are shared across the various types of likelihood and model types.

check_compatibility_between_input_and_model <- function(x, model) {
  cue_object <- get_cue_representation_from_model(model)

  # mvtnorm::dmvt expects means to be vectors, and x to be either a vector or
  # a matrix. In the latter case, each *row* of the matrix is an input.
  assert_that(is.list(x) | is.vector(x) | is.matrix(x) | is_tibble(x))

  # Do not reorder these conditionals (go from more to less specific)
  if (is_tibble(x)) x %<>% as.matrix() else
    if (is.list(x)) x %<>% reduce(rbind) %>% as.matrix() else
      if (is.vector(x)) x %<>% matrix(nrow = 1)

  stop("NOT YET COMPLETED.")
}


infer_default_noise_treatment <- function(Sigma_noise) {
  noise_treatment <-
    if (
    any(
      is.null(Sigma_noise),
      any(is.null(Sigma_noise)),
      any(map_lgl(Sigma_noise, is.null)))) "no_noise" else "marginalize"

  return(noise_treatment)
}

get_cue_representation_from_model <- function(x) {
  if (is.MVG(x)) {
    x <- x$mu
  } else if (is.NIW_belief(x)) {
    x <- x$m
  } else if (is.exemplars(x)) {
    x <- x$exemplars[[1]]$cues
  } else {
    stop(paste("Model object", deparse(substitute(x)), "not recognized."))
  }
  x <- names(if (is.list(x)) x[[1]] else if (is.numeric(x)) x[1] else stop(paste("Unable to extract cue representation from", deparse(substitute(x)))))

  return(x)
}


#' Get cue dimensionality from likelihood or model
#'
#' Get the number of cues from a likelihood (e.g., MVG or NIW_belief) or model (e.g., an MVG ideal observer
#' or NIW ideal adaptor) object.
#'
#' @param x  A likelihood or model object.
#'
#' @return A numeric.
#'
#' @export
get_cue_dimensionality_from_model <- function(x, indices = NULL) {
  x <- get_cue_representation_from_model(x)
  d <- if (is.null(dim(x))) length(x) else dim(x)[2]
  return(d)
}


#' Get cue labels from likelihood or model
#'
#' Get the names for all cues from a likelihood (e.g., MVG or NIW_belief) or model (e.g., an MVG ideal observer
#' or NIW ideal adaptor) object.
#'
#' @param x  A likelihood or model object.
#' @param indeces A vector of indices that should be turned into the original cue labels corresponding to those
#' indices, or `NULL` if all cue labels should be returned. (default: `NULL`)
#'
#' @return A character vector.
#'
#' @export
get_cue_labels_from_model <- function(x, indices = NULL) {
  x <- get_cue_representation_from_model(x)
  if (is.null(x)) x <- "cue"
  if (is.null(indices)) return(x) else return(x[indices])
  return(x)
}


#' Get category labels from likelihood or model
#'
#' Get the unique labels for all categories from a likelihood (e.g., MVG or NIW_belief) or model (e.g., an MVG ideal observer
#' or NIW ideal adaptor) object.
#'
#' @param x A likelihood or model object.
#'
#' @export
get_category_labels_from_model <- function(x) {
  if (is.MVBU_representation(x) | is.MVBU_model(x)) {
    return(sort(unique(x$category)))
  } else {
    error("Object not recognized.")
  }
}


#' Get number of categories from likelihood or model
#'
#' Get the number of unique category labels from a likelihood (e.g., MVG or NIW_belief) or model (e.g., an MVG ideal observer
#' or NIW ideal adaptor) object.
#'
#' @param x A likelihood or model object.
#'
#' @export
get_nlevels_of_category_labels_from_model <- function(x) {
  if (is.MVBU_representation(x) | is.MVBU_model(x)) {
    return(length(unique(x$category)))
  } else {
    error("Object not recognized.")
  }
}


#' Get priors from model
#'
#' @param model A model object.
#' @param categories A vector of category values.
#'
#' @return A vector of prior values of the same length as \code{categories}.
#'
#' @export
get_priors_from_model <- function(model, categories = model$category) {
  assert_that("prior" %in% names(model),
              msg = "No prior found in model.")

  prior <-
    model %>%
    left_join(tibble(category = categories), by = "category") %>%
    pull(prior)

  return(prior)
}

#' Get lapse rate from model
#'
#' @param model A model object.
#'
#' @export
get_lapse_rate_from_model <- function(model) {
  assert_that("lapse_rate" %in% names(model),
              msg = "No lapse_rate found in model.")

  lapse_rate <- unique(model$lapse_rate)

  assert_that(length(lapse_rate) == 1,
              msg = "More than one lapse_rate found in model.")
  return(lapse_rate)
}

#' Get lapse bias from model
#'
#' @param model A model object.
#' @param categories A vector of category values.
#'
#' @return A vector of lapse bias values of the same length as \code{categories}.
#'
#' @export
get_lapse_biases_from_model <- function(model, categories = model$category) {
  assert_that("lapse_bias" %in% names(model),
              msg = "No lapse_bias found in model.")

  lapse_bias <-
    model %>%
    left_join(tibble(category = categories), by = "category") %>%
    pull(lapse_bias)

  return(lapse_bias)
}


#' Get noise covariance matrix from model
#'
#' @param model A model object.
#'
#' @export
get_perceptual_noise_from_model <- function(model) {
  assert_that("Sigma_noise" %in% names(model),
              msg = "No Sigma_noise found in model.")

  Sigma_noise <- unique(model$Sigma_noise)

  assert_that(length(Sigma_noise) == 1,
              msg = "More than one Sigma_noise found in model.")
  return(Sigma_noise[[1]])
}


#' Nest/unnest the cue information in a model
#'
#' Take the centrality (e.g., mu or m) and scatter parameters (e.g., Sigma or S) of the model and nest them
#' into vector/matrix form or unnest them. The unnested format introduces two cue columns 'cue' and 'cue2',
#' which indicate the cue label for the cue value in the mu/m and Sigma/S columns.
#'
#' @param model A model object.
#'
#' @rdname nest_model
#' @export
nest_cue_information_in_model <- function(model) {
  if (is.MVG(model) | is.MVG_ideal_observer(model)) {
    m <- "mu"
    S <- "Sigma"
  } else if (is.NIW_belief(model) | is.NIW_ideal_adaptor(model) | is.NIW_ideal_adaptor_MCMC(model)) {
    m <- "m"
    S <- "S"
  } else {
    stop("Object not recognized.")
  }

  assert_that(all(c("cue", "cue2") %in% names(model)),
              msg = "cue and cue2 columns not found. There is nothing to nest.")
  model %>%
    group_by(across(-c(cue, cue2, !! sym(m), !! sym(S)))) %>%
    arrange(cue, cue2, .by_group = T) %>%
    summarise(
      !! sym(m) := list(make_named_vector(unique(!! sym(m)), unique(cue))),
      !! sym(S) := list(make_named_square_matrix(!! sym(S), unique(cue)))) %>%
    relocate(starts_with(c("lapse_", "prior")), .after = !! sym(S))
}

#' @rdname nest_model
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr across transmute
#' @export
unnest_cue_information_in_model <- function(model) {
  # Binding variables that RMD Check gets confused about otherwise
  # (since they are in non-standard evaluations)
  cue <- cue2 <- NULL

  if (is.MVG(model) | is.MVG_ideal_observer(model)) {
    m <- "mu"
    S <- "Sigma"
  } else if (is.NIW_belief(model) | is.NIW_ideal_adaptor(model) | is.NIW_ideal_adaptor_MCMC(model)) {
    m <- "m"
    S <- "S"
  } else {
    stop("Object not recognized.")
  }

  assert_that(all(c("cue", "cue2") %nin% names(model)),
              msg = "Cannot create cue and cue2 columns since they already exist in the model.")

  cue.labels <- get_cue_labels_from_model(model)

  model %<>%
    unnest(c(!! sym(m), !! sym(S))) %>%
    group_by(across(-c(!! sym(m), !! sym(S)))) %>%
    mutate(cue = cue.labels)

  for (i in 1:length(cue.labels)) {
      model %<>% mutate(., !! sym(cue.labels[i]) := (!! sym(S))[,i])
  }

  model %>%
    select(-S) %>%
    pivot_longer(cols = all_of(cue.labels), values_to = S, names_to = "cue2") %>%
    ungroup() %>%
    relocate(cue, cue2, .after = nu)
}


format_input_for_likelihood_calculation <- function(x, dim = 1) {
  assert_that(is.vector(x) | is.matrix(x) | is_tibble(x) | is.list(x))
  if (is.list(x)) x %<>% reduce(x, .f = ~ rbind(.x, format_input_for_likelihood_calculation(.y, dim = dim)))
  if (is_tibble(x)) x %<>% as.matrix() else
    if (is.vector(x)) {
      assert_that(length(x) %% dim == 0,
                  msg = paste("x cannot be coerced into matrix of observations with dimensionality", dim))
      x %<>% matrix(ncol = dim)
    }

  return(x)
}


#' Get posterior from model
#'
#' Categorize a single observation based a model. The decision rule can be specified to be either the
#' criterion choice rule, proportional matching (Luce's choice rule), or the sampling-based interpretation of
#' Luce's choice rule.
#'
#' @param x A vector of observations.
#' @param model A model object.
#' @param noise_treatment Determines whether and how multivariate Gaussian noise is added to the input.
#' See \code{\link[=get_MVG_likelihood]{get_MVG_likelihood}}. (default: "sample" if decision_rule is
#' "sample"; "marginalize" otherwise).
#' @param lapse_treatment Determines whether and how lapses will be treated. Can be "no_lapses", "sample" or "marginalize".
#' If "sample", whether a trial is lapsing or not will be sampled for each observations. If a trial is sampled to be
#' a lapsing trial the lapse biases are used as the posterior for that trial. If "marginalize", the posterior probability
#' will be adjusted based on the lapse formula lapse_rate * lapse_bias + (1 - lapse_rate) * posterior probability from
#' perceptual model. (default: "sample" if decision_rule is "sample"; "marginalize" otherwise).
#'
#' @return A tibble of observations with posterior probabilities for each category (in long format).
#'
#' @seealso TBD
#' @keywords TBD
#' @rdname get_posterior_from_model
#' @export
get_posterior_from_model <- function(model, ...) {
  if (is.MVG_ideal_observer(model)) {
    c <- get_posterior_from_MVG_ideal_observer(model = model, ...)
  } else if (is.NIW_ideal_adaptor(model)) {
    c <- get_posterior_from_NIW_ideal_adaptor(model = model, ...)
  } else {
    stop(
      paste(
        "get_categorization_from_* function for model type",
        class(model),
        "does not yet exist."))
  }

  return(c)
}


#' Get categorization from model
#'
#' Categorize a single observation based a model. The decision rule can be specified to be either the
#' criterion choice rule, proportional matching (Luce's choice rule), or the sampling-based interpretation of
#' Luce's choice rule.
#'
#' @param x A vector of observations.
#' @param model A model object.
#' @param decision_rule Must be one of "criterion", "proportional", or "sampling". (default: "sampling")
#' @param noise_treatment Determines whether and how multivariate Gaussian noise is added to the input.
#'   See \code{\link[=get_MVG_likelihood]{get_MVG_likelihood}}. (default: "sample" if decision_rule is
#'   "sample"; "marginalize" otherwise).
#' @param lapse_treatment Determines whether and how lapses will be treated. Can be "no_lapses", "sample" or "marginalize".
#'   If "sample", whether a trial is lapsing or not will be sampled for each observations. If a trial is sampled to be
#'   a lapsing trial the lapse biases are used as the posterior for that trial. If "marginalize", the posterior probability
#'   will be adjusted based on the lapse formula lapse_rate \emph{lapse_bias + (1 - lapse_rate)} posterior probability from
#'   perceptual model. (default: "sample" if decision_rule is "sample"; "marginalize" otherwise).
#' @param simplify Should the output be simplified, and just the label of the selected category be returned? This
#'   option is only available for the criterion and sampling decision rules. (default: `FALSE`)
#'
#' @return Either a tibble of observations with posterior probabilities for each category (in long format), or a
#'   character vector indicating the chosen category in the same order as the observations in x (if simplify = `TRUE`).
#'
#' @seealso TBD
#' @keywords TBD
#' @rdname get_categorization_from_model
#' @export
get_categorization_from_model <- function(model, decision_rule = "sampling", ...) {
  if (is.MVG_ideal_observer(model)) {
    c <- get_categorization_from_MVG_ideal_observer(model = model, decision_rule = decision_rule, ...)
  } else if (is.NIW_ideal_adaptor(model)) {
    c <- get_categorization_from_NIW_ideal_adaptor(model = model, decision_rule = decision_rule, ...)
  } else {
    stop(
      paste(
        "get_categorization_from_* function for model type",
        class(model),
        "does not yet exist."))
  }

  return(c)
}

#' Evaluate the fit of a model
#'
#' Evaluate the fit of a categorization model against a ground truth (e.g., human responses or the category intended
#' by a talker).
#'
#' @param model A model object.
#' @param x A vector of inputs (cue values).
#' @param response_category A vector of category responses corresponding to each input. The model is evaluated against
#'   this ground truth. Must be of the same length as the list of inputs, and each element of the `response_category` must
#'   be one of the `category` levels of `model`.
#' @param method Method for evaluating the model. Can be "accuracy", "likelihood", or "likelihood-up-to-constant". The
#'   latter two return the log-likelihood or the log-likelihood up to a constant that only depends on the data (rather than
#'   the model). This option calculates only sum_i sum_j n_ij log p_ij, where i = 1 ... M and M is the number of unique
#'   inputs x, and j = 1 ... K and K is the number of distinct categories, n_ij is the number of observed responses for
#'   category j at unique input x_i, and p_ij is the predicted posterior probability of response j at unique input x_i.
#'   Calculation of the "likelihood-up-to-constant" is thus particularly fast. Since this quantity is sufficient to compare
#'   models against each other on the same data, it makes the "likelihood-up-to-constant" method particularly useful for
#'   model fitting. It can, however, be important to keep in mind that this log-likelihood will be much smaller than the true
#'   data log-likelihood of the model under the assumption that the order of responses in the data should not be considered
#'   for the evaluation of the model.
#' @inheritParams get_categorization_from_model
#' @param decision_rule For details, see \code{\link{get_categorization_from_model}}. However, `evaluate_model` uses
#'   sensible defaults depending on the value of `method`. If `method` is (only) "accuracy", then defaults to "criterion".
#'   If `method` is "likelihood" or "likelihood-up-to-constant", then defaults to "proportional". Otherwise defaults to
#'   `NULL`, which will result in an error (so users must manually specify the decision rule).
#' @param return_by_x Should results be returned separately for each unique `x`? (default: `FALSE`)
#'
#' @return If `return_by_x`, the accuracy and/or log-likelihood of each unique input `x`. Otherwise, the overall accuracy
#' and/or log-likelihood of the all observations. Note that the overall log-likelihood is *not* simply the sum of the
#' log-likelihoods of all unique inputs `x`.
#'
#' @seealso TBD
#' @keywords TBD
#' @rdname evaluate_model
#' @export
evaluate_model <- function(
    model,
    x,
    response_category,
    method = "likelihood-up-to-constant",
    decision_rule = if (all("accuracy" == method)) "criterion" else if (any(c("likelihood", "likelihood-up-to-constant") %in% method)) "proportional" else NULL,
    ...,
    return_by_x = F
) {
  assert_that(all(method %in% c("likelihood", "likelihood-up-to-constant", "accuracy")))
  # When the input isn't a list, that's ambiguous between the input being a single input or a set of
  # 1D inputs. Use the model's cue dimensionality to disambiguate between the two cases.
  if (!is.list(x)) {
    x <- if (get_cue_dimensionality_from_model(model) == 1) as.list(x) else list(x)
  }
  assert_that(length(x) == length(response_category),
              msg = "Input x and response_category must be lists of the same length.")

  # Get counts of all k possible responses at all *unique* stimulus locations
  # (unique cue combinations)
  d.unique.observations <-
    tibble(
      x = .env$x,
      response_category = .env$response_category) %>%
    group_by(x, response_category) %>%
    tally() %>%
    ungroup()

  # Get predicted posterior probabilities of all k possible responses at all
  # *unique* stimulus locations (unique cue combinations)
  posterior <-
    d.unique.observations %>%
    distinct(x) %>%
    summarise(categorization = list(get_categorization_from_model(x = .data$x, model = .env$model, decision_rule = decision_rule, ...))) %>%
    unnest(categorization) %>%
    rename(posterior = response)

  r <- list()
  if ("accuracy" %in% method) {
    r[["accuracy"]] <-
      d.unique.observations %>%
      left_join(posterior, by = join_by(x == x, response_category == category))

    if (return_by_x) {
      r[["accuracy"]] %<>%
        group_by(x) %>%
        summarise(
          N = sum(.data$n),
          accuracy = sum(.data$posterior * .data$n) / sum(.data$n))
    } else {
      r[["accuracy"]] %<>%
        summarise(accuracy = sum(.data$posterior * .data$n) / sum(.data$n))
    }
  }
  if ("likelihood-up-to-constant" %in% method) {
    # Let n_ij be the number of observed responses for category j = 1 ... K for the
    # stimulus i = 1 ... M. Then the overall log data likelihood of a model that
    # assigns posterior probabilities p_ij to response j at stimulus i is:
    #
    #       log(L(p)) = Sum_i[ Sum_j[n_ij log(p_ij)] ] + constant
    #
    # where the constant only depends on the data (specifically, on both the individual
    # n_ij and on N = Sum_i Sum_j n_ij, but it is invariant for a given data set). To
    # evaluate a model against other models on the same data, it is thus sufficient to
    # compare the non-constant part of the log likelihood Sum_i[ Sum_j[n_ij log(p_ij)] ].
    # That's what is returned here. This makes the computation tractable but also means
    # that the true likelihood is (potentially much) larger (but by a constant amount,
    # as long as the data is held constant).
    r[["likelihood-up-to-constant"]] <-
      # Complete the count of responses to contain also the unobserved responses
      # (n = 0) at each stimulus location. Then join in the predicted posterior
      # probabilities p for each stimulus location.
      d.unique.observations %>%
      left_join(posterior, by = join_by(x == x, response_category == category)) %>%
      group_by(x) %>%
      summarise(
        N = sum(.data$n),
        log_likelihood = sum(.data$n * log(.data$posterior)))

    if (!return_by_x) {
      r[["likelihood-up-to-constant"]] %<>%
        summarise(log_likelihood = sum(.data$log_likelihood))
    }
  }
  if ("likelihood" %in% method) {
    # Derivation of multinomial density: https://statproofbook.github.io/P/mult-pmf.html
    # Proof for multinomial coefficient: https://math.stackexchange.com/questions/548027/prove-multinomial-coefficient-probability-theory
    #
    # Let c_i be the category response to input x_i for 1 ... M observations. For a model
    # that predicts the category response c_i to have a posterior probability p_i, the
    # data likelihood of the observed category responses is the product of all the p_i
    # times the number of possible permutations of all the p_i (since we do not care
    # about the order of responses), --------------------- CONTINUE HERE
    #
    # Let n_ij be the number of observed responses for category j = 1 ... M for the
    # stimulus i = 1 ... N. Then the overall log data likelihood for all observations
    # at stimulus i of a model that assigns posterior probabilities p_ij to response j
    # at stimulus i is:
    #
    #       log(L(p)) = log(N_i!) + Sum(n_ij log(p_ij)) - Sum(log(n_ij!))
    #
    # where N is the total number of responses at location i (i.e., sum_i n_ij), p_ij
    # is the vector of K posterior probabilities (summing to 1) at location i, and n_ij
    # is the number of occurrences of the j-th outcome (out of the 1 ... K possible
    # outcomes) at location i. dmultinom(n, p, log = T) gives us this log likelihood.
    #
    # Computationally, it is most efficient to calculate log-likelihoods for each
    # unique stimulus location (i.e., unique cue combination), and then to aggregate
    # the resulting stimulus-specific log-likelihoods into the overall log-likelihood
    # of the data. The first step of this is well-formed since  all responses to a
    # unique cue combination are predicted to have the same posterior distribution p.
    # So we can use dmultinom(x_[at stimulus location], p_[at stimulus location]).
    #
    # Unfortunately, we cannot simply *sum* the different log-likelihoods of the
    # different stimulus locations since the N! should be based on the aggregate counts
    # *across* stimulus positions, and the n_ij have to be combined in ways I haven't
    # yet figured out.
    #
    warning('method = "likelihood" is not yet working properly. DO NOT TRUST THESE RESULTS!')
    if (return_by_x) {
      # Note that these by-x likelihoods cannot simply be summed up to get the overall
      # likelihood. That would fail to correct for the total number of permutations.
      r[["likelihood"]] <-
        d.unique.observations %>%
        # Complete the count of responses to contain also the unobserved responses
        # (n = 0) at each stimulus location. Then join in the predicted posterior
        # probabilities p for each stimulus location. This is required because we're
        # using dmultinom below.
        complete(x, response_category) %>%
        replace_na(list(n = 0)) %>%
        left_join(posterior, by = join_by(x == x, response_category == category)) %>%
        group_by(x) %>%
        summarise(
          N = sum(n),
          log_likelihood = dmultinom(x = n, prob = posterior, log = T))
    } else {
      r[["likelihood"]] <-
        d.unique.observations %>%
        left_join(posterior, by = join_by(x == x, response_category == category)) %>%
        summarise(log_likelihood = lfactorial(sum(n)) + sum(n * log(.data$posterior)) - sum(lfactorial(n)))
    }
  }

  # Simplify return as much as possible
  if (length(r) == 1) {
    r <- r[[1]]
    if (nrow(r) <= 1) r <- as.numeric(r)
  }
  return(r)
}
