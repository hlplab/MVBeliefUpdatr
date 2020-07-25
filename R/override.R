#' @importFrom mvtnorm rmvnorm
NULL

# Added here to handle the case of univariate categories
rmvnorm = function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), ...) {
  if (length(mean) == 1)
    return(rnorm(n = n, mean = as.vector(mean), sd = as.vector(sigma^.5))) else
      return(mvtnorm::rmvnorm(n, mean, sigma, ...))
}
