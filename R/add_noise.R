#' Add Gaussian noise to cues in a data set.
#'
#' Adds Gaussian noise to cues in a data set.
#'
#' @param data \code{data.frame} or \code{tibble}.
#' @param cues Column names of cues in \code{data}.
#' @param method If "sample" then a sample of
#' noise is added to the input. If "marginalize" then each observation is transformed into the marginal distribution
#' that results from convolving the input with noise. Currently, only "sample" is implemented (default: "sample")
#' @param Sigma Covariance matrix of noise.
#'
#' @return Same as \code{data}.
#'
#' @export
add_noise = function(
  data,
  cues,
  method = c("sample", "marginalize"),
  Sigma
) {
  assert_that(is_tibble(data) | is.data.frame(data))
  assert_that(nrow(data),
              msg = paste("There must be at least one observation in the data. Found", nrow(data), "observations."))
  assert_cols_in_data(data, cues, scalar = F)
  assert_that(method %in% c("sample", "marginalize"),
              msg = 'method must be one of "sample" or "marginalize"')
  assert_that(!is.null(Sigma),
              msg = "No Sigma provided.")
  asswer_that(is.matrix(Sigma),
              msg = "Sigma must be a matrix.")
  assert_that(all(length(dim(Sigma)) == 2, dim(Sigma)[1] == dim(Sigma)[2], dim(Sigma)[1] == length(cues)),
              msg = "Sigma must be a square matrix of the same dimensionality as the length of cues.")

  data %<>%
    make_vector_column(cues, "cues")

  if (method == "sample") {
    data %<>%
      mutate(cues = map(cues, ~ .x + rmvnorm(n = 1, sigma = Sigma)))
  } else stop("This method is not yet implemented.")

  return(data)
}
