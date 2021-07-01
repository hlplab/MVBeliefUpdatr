#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom extraDistr rlst dlst
NULL

# Added here to handle the case of univariate categories
#' @export
colMeans <- function(x, ...) {
  if (!is.array(x) || length(dim(x)) < 2L)
    mean(x) else colMeans(x)
}

# Added here to handle the case of univariate categories
#' @export
cov <- function(x, ...) {
  if (!is.array(x) || length(dim(x)) < 2L)
    var(x) else cov(x)
}

# Added here to handle the case of univariate categories
#' @export
rmvnorm <- function (n, mean = rep(0, length(sigma)^.5), sigma = diag(length(mean)), ...) {
  if (length(mean) == 1)
    return(rnorm(n = n, mean = as.vector(mean), sd = as.vector(sigma^.5), ...)) else
      return(rmvnorm(n, mean, sigma, ...))
}

# Added here to handle the case of univariate categories
#' @export
rmvt <- function (n, delta = rep(0, length(sigma)^.5), sigma = diag(length(mean)), df, ...) {
  if (length(mean) == 1)
    # Using rlst instead of rt since rt is for standardized t distribution (no scale parameter)
    return(rlst(n = n, df = df, mu = as.vector(delta), sigma = as.vector(sigma^.5), ...)) else
      return(rmvt(n = n, delta = delta, sigma = sigma, df = df, ...))
}

# Added here to handle the case of univariate categories
#' @export
dmvnorm <- function (x, mean = rep(0, length(sigma)^.5), sigma = diag(length(mean)), ...) {
  if (length(mean) == 1)
    return(dnorm(x = x, mean = as.vector(mean), sd = as.vector(sigma^.5))) else
      return(dmvnorm(x = x, mean, sigma, ...))
}

# Added here to handle the case of univariate categories
#' @export
dmvt <- function (x, delta = rep(0, length(sigma)^.5), sigma = diag(length(mean)), df, ...) {
  if (length(mean) == 1)
    # Using dlst instead of dt since dt is for standardized t distribution (no scale parameter)
    return(dlst(x = x, df = df, mu = as.vector(delta), sigma = as.vector(sigma^.5), ...)) else
      return(dmvt(x = x, delta = delta, sigma = sigma, df = df, ...))
}

