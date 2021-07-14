#' @importFrom stats rmultinom
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom extraDistr rlst dlst
NULL

# Added here to handle the case of univariate categories
#' @export
colMeans <- function(x, ...) {
  if (!is.array(x) || length(dim(x)) < 2L)
    mean(x) else base::colMeans(x)
}

# Added here to handle the case of univariate categories
#' @export
cov <- function(x, ...) {
  if (!is.array(x) || length(dim(x)) < 2L)
    var(x) else stats::cov(x)
}

#' Overrides stats::rmultinom, reformatting its output. Specifically, the output of stats:rmultinom is transpased, so
#' that the returned value is a matrix in which each row (rather than column) corresponds to one observation and each
#' column (rather than row) corresponds to one of the categorical outcomes. Each cell represents the counts observed
#' for each outcomes on that observation.
#' @export
rmultinom <- function(n, size, prob) {
  return(t(stats::rmultinom(n, size, prob)))
}


# Added here to handle the case of univariate categories
#' @export
rmvnorm <- function(n, mean = rep(0, length(sigma)^.5), sigma = diag(length(mean)), ...) {
  if (length(mean) == 1)
    return(rnorm(n = n, mean = as.vector(mean), sd = as.vector(sigma^.5), ...)) else
      return(mvtnorm::rmvnorm(n, mean, sigma, ...))
}

# Added here to handle the case of univariate categories
#' @export
rmvt <- function (n, delta = rep(0, length(sigma)^.5), sigma = diag(length(mean)), df, ...) {
  if (length(mean) == 1)
    # Using rlst instead of rt since rt is for standardized t distribution (no scale parameter)
    return(rlst(n = n, df = df, mu = as.vector(delta), sigma = as.vector(sigma^.5), ...)) else
      return(mvtnorm::rmvt(n = n, delta = delta, sigma = sigma, df = df, ...))
}

# Added here to handle the case of univariate categories
#' @export
dmvnorm <- function (x, mean = rep(0, length(sigma)^.5), sigma = diag(length(mean)), ...) {
  if (length(mean) == 1)
    return(dnorm(x = x, mean = as.vector(mean), sd = as.vector(sigma^.5))) else
      return(mvtnorm::dmvnorm(x = x, mean, sigma, ...))
}

# Added here to handle the case of univariate categories
#' @export
dmvt <- function (x, delta = rep(0, length(sigma)^.5), sigma = diag(length(mean)), df, ...) {
  if (length(mean) == 1)
    # Using dlst instead of dt since dt is for standardized t distribution (no scale parameter)
    return(dlst(x = x, df = df, mu = as.vector(delta), sigma = as.vector(sigma^.5), ...)) else
      return(mvtnorm::dmvt(x = x, delta = delta, sigma = sigma, df = df, ...))
}

