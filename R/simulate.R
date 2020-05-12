#' @import assertthat
#' @importFrom mvtnorm rmvnorm
#' @importFrom purrr pmap
NULL

#' Make multivariate Gaussian exposure data.
#'
#' Returns a tibble of observations drawn from multivariate Gaussians, with one observation per row. Each row
#' provides the category label and cue values. If \code{keep.parameters = T} then the parameters (\code{N, mean, sigma})
#' are also returned.
#'
#' The input is expected to be lists of parameters which the n-th element of each list specifying the n-th Gaussian and the
#' number of observations to be drawn from that Gaussian.
#'
#' @param Ns A list of integers, with each number specifying the number of observations to be drawn from the corresponding
#' Gaussian.
#' @param means List of mean vectors, each specifying the mean of a multivariate Gaussian.
#' @param sigmas List of covariance matrices, each specifying the covariance of a multivariate Gaussian.
#' @param category.labels List of category names, each specifying the category label of a multivariate Gaussian. If \code{NULL}
#' (default) then Gaussians will be numbered from 1:N.
#' @param cue.labels List of cue names. If \code{NULL} (default) then the cues will be numbered cue1, cue2, ...
#' @param keep.parameters Should the parameters handed to this function be included in the output? (default: FALSE)
#'
#' @return A tibble.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
make_MV_exposure_data = function(Ns, means, Sigmas, category.labels = NULL, cue.labels = NULL, keep.parameters = F) {
  assert_that(!is.null(means), !is.null(Sigmas))
  assert_that(is.null(category.labels) | length(means) == length(category.labels),
              msg = "Number of category labels mismatch number of means.")
  assert_that(is.null(cue.labels) | length(means[[1]]) == length(cue.labels),
              msg = "Number of cue labels mismatches dimensionality of means.")

  if (is.null(category.labels)) category.labels = 1:length(means)
  if (is.null(cue.labels)) cue.labels = paste0("cue", 1:length(means[[1]]))

  x = tibble(n = unlist(Ns), mean = means, sigma = Sigmas) %>%
    mutate(data = pmap(.l = list(n, mean, sigma), .f = rmvnorm)) %>%
    mutate(category = unlist(category.labels)) %>%
    mutate(data = map(data, ~ .x %>% as.data.frame() %>% rename_all(~ cue.labels))) %>%
    unnest(data)

  if (keep.parameters) return(x) else return(x %>% select(category, everything(), -c(n, mean, sigma)))
}


