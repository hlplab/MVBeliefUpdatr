#' Get cue labels from NIW belief object
#'
#' Get the names for all cues from an NIW belief object.
#'
#' @param x An NIW belief object.
#'
#' @export
get_cue_labels_from_NIW_belief = function(x) {
  return(names(x$M[[1]]))
}
