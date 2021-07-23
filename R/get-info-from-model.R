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
