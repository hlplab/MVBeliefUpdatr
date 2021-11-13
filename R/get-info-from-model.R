#' Functions that are shared across the various types of likelihood and model types.


#' Get cue labels from likelihood or model
#'
#' Get the names for all cues from a likelihood (e.g., MVG or NIW_belief) or model (e.g., an MVG ideal observer
#' or NIW ideal adaptor) object.
#'
#' @param x  A likelihood or model object.
#'
#' @export
get_cue_labels_from_model = function(x) {
  if (is.MVG(x) | is.MVG_ideal_observer(x)) {
    x <- x$mu
  } else if (is.NIW_belief(x) | is.NIW_ideal_adaptor(x)) {
    x <- x$m
  } else {
    error("Object not recognized.")
  }

  x <- names(if (is.list(x)) x[[1]] else if (is.numeric(x)) x[1] else error("No suitable information found."))

  if (is.null(x)) x <- "cue"
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
get_category_labels_from_model = function(x) {
  if (is.MVG(x) | is.MVG_ideal_observer(x) | is.NIW_belief(x) | is.NIW_ideal_adaptor(x)) {
   return(levels(x$category))
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
get_nlevels_of_category_labels_from_model = function(x) {
  if (is.MVG(x) | is.MVG_ideal_observer(x) | is.NIW_belief(x) | is.NIW_ideal_adaptor(x)) {
    return(length(unique(x$category)))
  } else {
    error("Object not recognized.")
  }
}


get_lapse_rate_from_model <- function(model) {
  if ("lapse_rate" %in% names(model)) {
    lapse_rate <- unique(model$lapse_rate)

    asser_that(length(lapse_rate) == 1,
               msg = "More than one lapse_rate found in model.")
    return(lapse_rate)
  } else return(NA)
}

get_lapse_biases_from_model <- function(model) {
  if ("lapse_bias" %in% names(model)) {
    lapse_bias <- model$lapse_bias

    asser_that(sum(lapse_bias) == 1,
               msg = paste("The lapse_biases do not add up to 1:", sum(lapse_bias)))
    return(lapse_bias)
  } else return(NA)
}


get_perceptual_noise_from_model <- function(model) {
  if ("Sigma_noise" %in% names(model)) {
    Sigma_noise <- unique(model$Sigma_noise)

    asser_that(length(Sigma_noise) == 1,
               msg = "More than one Sigma_noise found in model.")
    return(Sigma_noise)
  } else return(NA)
}
