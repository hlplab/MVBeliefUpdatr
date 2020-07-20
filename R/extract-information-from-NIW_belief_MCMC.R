NULL

#' Get cue labels from an NIW belief MCMC object
#'
#' Returns the cue labels for the specified indices extracted from the M column of the NIW belief MCMC object.
#' If no indices are provided, then all cue labels are returned.
#'
#' @param x mv_ibbu_stanfit object.
#' @param indeces A vector of indices that should be turned into the original cue labels corresponding to those
#' indices, or `NULL` if all cue labels should be returned. (default: `NULL`)
#'
#' @return A character vector.
#'
#' @rdname get_cue_labels
#' @export
get_cue_labels = function(x, indices = NULL) {
  message("This function can probably be integrated with get_original_levels(), which should be renamed.")
  assert_that(is.NIW_belief_MCMC(x))

  names = names(x$M[[1]])

  if (is.null(indices)) return(names) else return(names(indices))
}
