#' Add Gaussian noise to cues in a data set
#'
#' Adds Gaussian noise to cues in a data set.
#'
#' @param data \code{data.frame} or \code{tibble}.
#' @param cues Column names of cues in \code{data}.
#' @param method If \code{"sample"} (default) then a sample of noise is added
#'   to the input. If \code{"marginalize"} then each observation is transformed
#'   into the marginal distribution that results from convolving the input with
#'   noise. Currently, only \code{"sample"} is implemented.
#' @param Sigma Covariance matrix of noise.
#'
#' @return Same as \code{data} with noise added to cues.
#'
#' @importFrom dplyr mutate
#' @export
#' @examples
#' df <- data.frame(vot = c(10, 20, 30), f0 = c(150, 180, 210))
#' noise_cov <- matrix(c(4, 0, 0, 9), nrow = 2)
#' add_noise(df, cues = c("vot", "f0"), Sigma = noise_cov)
add_noise <- function(
  data,
  cues,
  method = c("sample", "marginalize"),
  Sigma
) {
  method <- match.arg(method)
  .assert_data_contains_cols(data, cues)
  .assert_true(
    nrow(data) >= 1,
    msg = paste(
      "There must be at least one observation in the data. Found",
      nrow(data), "observations."
    )
  )
  .assert_matrix(Sigma, msg = "Sigma must be a matrix.")
  .assert_all(
    length(dim(Sigma)) == 2,
    dim(Sigma)[1] == dim(Sigma)[2],
    dim(Sigma)[1] == length(cues),
    msg = paste(
      "Sigma must be a square matrix of the same dimensionality",
      "as the length of cues."
    )
  )

  data %<>%
    .make_vector_column(cues, "cues")

  if (method == "sample") {
    data <- data %>%
      mutate(cues = map(cues, ~ .x + .rmvnorm(n = 1, sigma = Sigma)))
  } else {
    stop("This method is not yet implemented.")
  }

  return(data)
}

