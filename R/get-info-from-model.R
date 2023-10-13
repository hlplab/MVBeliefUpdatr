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

# Just for internal use
get_cue_representation_from_model <- function(x) {
  if (is.MVG(x) | is.MVG_ideal_observer(x)) {
    x <- x$mu
  } else if (is.NIW_belief(x) | is.NIW_ideal_adaptor(x) | is.NIW_ideal_adaptor_MCMC(x)) {
    x <- x$m
  } else if (is.exemplars(x) | is.exemplar_model(x)) {
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
  if (is.MVG(x) | is.NIW_belief(x) | is.exemplars(x) | is.model(x)) {
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
  if (is.MVG(x) | is.NIW_belief(x) | is.exemplars(x) | is.model(x)) {
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
      !! sym(S) := list(make_named_square_matrix(!! sym(S), unique(cue))))
}

#' @rdname nest_model
#' @export
unnest_cue_information_in_model <- function(model) {
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

  model %>%
    unnest(c(!! sym(m), !! sym(S))) %>%
    group_by(across(-c(!! sym(m), !! sym(S)))) %>%
    mutate(cue = cue.labels) %>%
    group_by(across(-c(!! sym(S)))) %>%
    transmute(!! sym(cue.labels[1]) := (!! sym(S))[,1], !! sym(cue.labels[2]) := (!! sym(S))[,2]) %>%
    pivot_longer(cols = all_of(cue.labels), values_to = S, names_to = "cue2") %>%
    ungroup() %>%
    select(cue, cue2, everything())
}


format_input_for_likelihood_calculation <- function(x) {
  assert_that(is.vector(x) | is.matrix(x) | is_tibble(x) | is.list(x))
  if (is.list(x)) x %<>% reduce(x, .f = ~ rbind(.x, format_input_for_likelihood_calculation(.y)))
  if (is_tibble(x)) x %<>% as.matrix() else
    if (is.vector(x)) x %<>% matrix(nrow = 1)

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
#' @examples
#' TBD
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
#' @param decision_rule Must be one of "criterion", "proportional", or "sampling".
#' @param noise_treatment Determines whether and how multivariate Gaussian noise is added to the input.
#' See \code{\link[=get_MVG_likelihood]{get_MVG_likelihood}}. (default: "sample" if decision_rule is
#' "sample"; "marginalize" otherwise).
#' @param lapse_treatment Determines whether and how lapses will be treated. Can be "no_lapses", "sample" or "marginalize".
#' If "sample", whether a trial is lapsing or not will be sampled for each observations. If a trial is sampled to be
#' a lapsing trial the lapse biases are used as the posterior for that trial. If "marginalize", the posterior probability
#' will be adjusted based on the lapse formula lapse_rate * lapse_bias + (1 - lapse_rate) * posterior probability from
#' perceptual model. (default: "sample" if decision_rule is "sample"; "marginalize" otherwise).
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
#' @rdname get_categorization_from_model
#' @export
get_categorization_from_model <- function(model, ...) {
  if (is.MVG_ideal_observer(model)) {
    c <- get_categorization_from_MVG_ideal_observer(model = model, ...)
  } else if (is.NIW_ideal_adaptor(model)) {
    c <- get_categorization_from_NIW_ideal_adaptor(model = model, ...)
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
#' this ground truth. Must be of the same length as the list of inputs, and each element of the `response_category` must
#' be one of the `category` levels of `model`.
#' @param method Method for evaluating the model. Can be "accuracy" or "likelihood". The latter returns the log-likelihood.
#' @inheritParams get_categorization_from_model
#' @param return_by_x Should results be returned separately for each unique `x`? (default: `FALSE`)
#'
#' @return If `return_by_x`, the accuracy and/or log-likelihood of each unique input `x`. Otherwise, the overall accuracy
#' and/or log-likelihood of the all observations. Note that the overall log-likelihood is *not* simply the sum of the
#' log-likelihoods of all unique inputs `x`.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @rdname evaluate_model
#' @export
evaluate_model <- function(model, x, response_category, method = "likelihood", ..., return_by_x = F) {
  assert_that(all(method %in% c("likelihood", "accuracy")))
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
    tally()

  # Get predicted posterior probabilities of all k possible responses at all
  # *unique* stimulus locations (unique cue combinations)
  posterior <-
    d.unique.observations %>%
    ungroup() %>%
    distinct(x) %>%
    summarise(categorization = list(get_categorization_from_model(x = .data$x, model = .env$model, ...))) %>%
    unnest(categorization) %>%
    rename(posterior = response)

  r <- list()
  if ("accuracy" %in% method) {
    r[["accuracy"]] <-
      d.unique.observations %>%
      ungroup() %>%
      left_join(posterior, by = join_by(x == x, response_category == category)) %>%
      { if (return_by_x) group_by(., x) else . } %>%
      summarise(accuracy = sum(.data$posterior * .data$n) / sum(.data$n))
  }
  if ("likelihood" %in% method) {
    # There seem to be two ways to calculate the log-likelihood of the data for a multinomial regression:
    #
    # 1) the sum of all log probabilities for each trial. this can be simplified by calculating log p for
    #    each unique stimulus location, and then summing up those log ps across stimulus locations.
    #
    # 2) however, that would seem to consider one particular order of outcomes. the alternative that is
    #    considering all possible order of outcomes that yield the overall outcome is given by
    #    dmultinom(x, prob, log =T) but this can*not* be simply added up across stimulus locations


    # The multinomial log-likelihood is:
    #
    #       log(L(p)) = log N! + Sum(n_j log(p_j)) - Sum(log(n_j!))
    #
    # where N is the total number of responses, p is the vector of k posterior
    # probabilities (summing to 1), and n_j is the number of occurrences of the j-th
    # outcome (out of the 1 ... k possible outcomes). dmultinom(n, p, log = T) gives
    # us this log likelihood.
    #
    # Computationally, it is most efficient to calculate log-likelihoods for each
    # unique stimulus location (i.e., unique cue combination), and then to aggregate
    # the resulting stimulus-specific log-likelihoods into the overall log-likelihood
    # of the data. The first step of this is well-formed since  all responses to a
    # unique cue combination are predicted to have the same posterior distribution p.
    # So we can use dmultinom(x_[at stimulus location], p_[at stimulus location]).
    #
    # Unfortunately, we cannot simply *sum* the different log-likelihoods of the
    # different stimulus locations since the N! and n_j! should be based on the
    # aggregate counts *across* stimulus positions.
    #
    # However, luckily, the only component of the above equation that is *not* a
    # constant of the data is Sum(n_j log(p_j)), which can be added across locations.
    #
    r[["likelihood"]] <-
      # Complete the count of responses to contain also the unobserved responses
      # (n = 0) at each stimulus location. Then join in the predicted posterior
      # probabilities p for each stimulus location.
      d.unique.observations %>%
      complete(x, response_category) %>%
      replace_na(list(n = 0)) %>%
      left_join(posterior, by = join_by(x == x, response_category == category)) %>%
      group_by(x)

    if (return_by_x) {
      r[["likelihood"]] %<>%
        summarise(
          x = first(x),
          n = sum(n),
          log_likelihood = dmultinom(x = n, prob = posterior, log = T))
    } else {
      r[["likelihood"]] %<>%
        summarise(
          # log-likelihood for x up to constant (so that the components can be correctly summed below)
          log_likelihood = sum(n * log(.data$posterior)),
          N = sum(n),
          n_responses_at_x = list(cbind(correct_category, n))) %>%
        summarise(
          log_likelihood =
            sum(log_likelihood) +
            log(factorial(sum(n))) -
            sum(log(factorial(reduce(n_responses, `+`)))))
      }
  }

  # Simplify return as much as possible
  if (length(r) == 1) {
    r <- r[[1]]
    if (nrow(r) <= 1) r <- as.numeric(r)
  }
  return(r)
}
