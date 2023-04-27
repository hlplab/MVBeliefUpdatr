#' Functions that are shared across the various types of likelihood and model types.


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
  if (is.MVG(x) | is.MVG_ideal_observer(x)) {
    x <- x$mu
  } else if (is.NIW_belief(x) | is.NIW_ideal_adaptor(x) | is.NIW_ideal_adaptor_MCMC(x)) {
    x <- x$m
  } else if (is.exemplars(x) | is.exemplar_model(x)) {
    x <- x$exemplars[[1]]$cues
  } else {
    stop("Object not recognized.")
  }

  x <- names(if (is.list(x)) x[[1]] else if (is.numeric(x)) x[1] else error("No suitable information found."))

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


get_categorization_from_model <- function(model, ...) {
  if (is.NIW_ideal_adaptor(model)) {
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
