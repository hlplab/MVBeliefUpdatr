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
#'
#' @export
get_priors_from_model <- function(model, categories = model$category) {
  assert_that("prior" %in% names(model),
              msg = "No prior found in model.")

  prior <- model %>%
    filter(category %in% categories) %>%
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
#'
#' @export
get_lapse_biases_from_model <- function(model, categories = model$category) {
  assert_that("lapse_bias" %in% names(model),
              msg = "No lapse_bias found in model.")

  lapse_bias <- model %>%
    filter(category %in% categories) %>%
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
    pivot_longer(cols = cue.labels, values_to = S, names_to = "cue2") %>%
    ungroup() %>%
    select(cue, cue2, everything())
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
#' @param x A vector of observations.
#' @param model A model object.
#' @param cues A list of cue values.
#' @param correct_category A list of category labels that is taken to be the ground truth. Must be of the same length as the
#' list of cue values.
#' @param method Method for evaluating the model. Can be "accuracy" or "likelihood".
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
#' @rdname evaluate_model
#' @export
evaluate_model <- function(model, cues, correct_category, method = "likelihood", ...) {
  # Get posterior for all *unique* combinations of cues
  d.unique.observations <-
    tibble(
      x = .env$cues,
      correct_category = .env$correct_category) %>%
    group_by(x, correct_category) %>%
    tally()

  posterior <-
    d.unique.observations %>%
    ungroup() %>%
    summarise(categorization = list(get_categorization_from_model(x = .data$x, model = .env$model, ...))) %>%
    unnest(categorization)

  r <- list()
  if ("accuracy" %in% method) {
    r <-
      append(
        r,
        d.unique.observations %>%
          ungroup() %>%
          left_join(posterior, by = join_by(x == x, correct_category == category)) %>%
          summarise(accuracy = sum(.data$response * .data$n) / sum(.data$n)))
  }
  if ("likelihood" %in% method) {
    r <-
      append(
        r,
        # Get all unique combinations of cues and *possible* responses and fill in 0 as
        # count n for all combinations that aren't observed
        crossing(x = .env$cues, correct_category = .env$model$category) %>%
          left_join(d.unique.observations, by = join_by(x, correct_category)) %>%
          replace_na(list(n = 0)) %>%
          left_join(posterior, by = join_by(x == x, correct_category == category)) %>%
          # Since dmultinom already takes into account the number of observations (size),
          # no need to carry through the number of observations.
          group_by(x) %>%
          summarise(log_likelihood = dmultinom(x = .data$n, prob = .data$response, log = T)) %>%
          summarise(log_likelihood = sum(log_likelihood)))
  }

  if (length(r) == 1) r <- as.numeric(r[[1]])
  return(r)
}
