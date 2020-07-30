NULL

#' Get covariance matrix from correlation matrix and standard deviations
#'
#' Get covariance matrix from correlation matrix and vector of standard deviations. Based on code recommended
#' by Ben Bolder for efficiency.
#'
#' @param omega A correlation matrix.
#' @param tau A vector of standard deviations.
#'
#' @return A square matrix.
#'
#' @seealso \code{\link[stats]{cov2cor}}, \code{\link{cov2tau}}
#' @references \url{https://stats.stackexchange.com/questions/62850/obtaining-covariance-matrix-from-correlation-matrix}
#' @keywords TBD
#' @examples
#' TBD
#' @export
cor2cov = function(omega, tau) {
  assert_that(is.matrix(omega))
  assert_that(is.numeric(tau))
  assert_that(nrow(omega) == ncol(omega),
              msg = "omega must be a square matrix.")
  assert_that(length(tau) == nrow(omega),
              msg = "tau must have as many elements as omega has rows (and columns).")

  outer(tau,tau) * omega
}

#' Get vector of standard deviations from covariance matrix
#'
#' Get vector of standard deviations from covariance matrix.
#'
#' @param v A covariance matrix.
#'
#' @return A numeric vector.
#'
#' @seealso \code{\link[stats]{cov2cor}}, \code{\link{cor2cov}}
#' @keywords TBD
#' @examples
#' TBD
#' @export
cov2tau = function(v) {
  assert_that(is.matrix(v))
  assert_that(nrow(v) == ncol(v),
              msg = "omega must be a square matrix.")
  sqrt(diag(v))
}

#' Get sum of uncentered squares
#'
#' Get sum of uncentered squares. This quantity is a sufficient statistic for, for example, multivariate Gaussian
#' belief-updating under an Normal-Inverse-Wishart prior.
#'
#' @param data A `tibble`, `data.frame`, or `matrix`. If data is a `tibble` or `data.frame`, the columns for
#' specified variables are extracted and (together) converted into a matrix with as many colums as there are
#' variables.
#' @param variables Only required if data is not already a `matrix`.
#'
#' @return A matrix.
#'
#' @seealso
#' @keywords TBD
#' @examples
#' TBD
#' @export
get_sum_of_uncentered_squares <- function(data, variables, verbose = F) {
  assert_that(is.tibble(data) | is.data.frame(data) | is.matrix(data))
  if (is.tibble(data) | is.data.frame(data))
    assert_that(variables %in% data,
                msg = paste("Variable column(s)", variables[which(variables %nin% names(data))], "not found in data."))

  data.matrix = if (is_tibble(data) | is.data.frame(data))
    data[[variables]] %>%
      as.matrix() else data

  k = dim(data.matrix)[2]
  m = matrix(ncol = k, nrow = k)

  for (i in 1:k) {
    m[i,i] = sum(data.matrix[,i]**2)
    if (i < k) for (j in (i + 1):k) {
      m[j,i] = sum(data.matrix[,i] * data.matrix[,j])
      m[i,j] = m[j,i]
    }
  }

  return(m)
}
