NULL

#' Get cue labels from an NIW belief MCMC object [DEPRECATED: integrated with get_cue_labels_from_model]
#'
#' Returns the cue labels for the specified indices extracted from the m column of the NIW belief MCMC object.
#' If no indices are provided, then all cue labels are returned.
#'
#' @param x mv_ibbu_stanfit object.
#' @param indices A vector of indices that should be turned into the original cue labels corresponding to those
#' indices, or `NULL` if all cue labels should be returned. (default: `NULL`)
#'
#' @return A character vector.
#'
#' @rdname get_cue_labels_fom_MCMC
#' @export
get_cue_labels_fom_MCMC = function(x, indices = NULL) {
  assert_that(is.NIW_ideal_adaptor_MCMC(x))

  names = names(x$m[[1]])
  if (is.null(names)) rownames(x$S[[1]])

  if (is.null(indices)) return(names) else return(names(indices))
}
