#' @import assertthat
NULL

#' Get expected category mean mu or covariance matrix sigma
#'
#' Returns the expected value of posterior marginal distribution over category means mu or
#' category covariance matrix Sigma, marginalized over all MCMC samples.
#'
#' Each MCMC samples' expected value \code{E[mu] = M_n}
#' (i.e, the posterior/updated mean of the mulativariate Normal over category means \code{mu}).
#' Marginalizing across all MCMC samples (representing uncertainty in the true value of
#' \code{M_n}), we get \code{E\[E\[mu\]\] = mean(M_n)}.
#'
#' Each MCMC samples' expected value
#' \code{E[Sigma] = S_n / (nu_n - D - 1)}, where \code{S_n} is the posterior/updated scatter matrix,
#' \code{nu_n} is the posterior/updated pseudocount representing the strength of the posterior/updated
#' beliefs over category covariance matrices sigma (i.e., the inverse-Wishart), and \code{D} is
#' the dimension of the multivariate Normal. Marginalizing across all MCMC samples
#' (representing uncertainty in the true value of \code{S_n}), we get
#' \code{E[E[Sigma]] = mean(S_n / (nu_n - D - 1))}.
#'
#' @seealso TBD
#' @keywords TBD
#' @references Murphy, K. P. (2012). Machine learning: a probabilistic perspective. MIT press.
#' @examples
#' TBD
#' @rdname get_expected_category_statistics
#' @export
get_expected_mu = function(x, category, group) {
  assert_that(is.mv_ibbu_stanfit(x) | is.mv_ibbu_MCMC(x, nested = T, long = T))
  assert_that(is.character(category) | is.numeric(category))
  assert_that(is.character(group) | is.numeric(group))

  if (is.mv_ibbu_stanfit(x))
    x = add_ibbu_draws(x, which = "both", wide = F, nest = T)

  x %<>%
    filter(group == group, category == category) %>%
    summarise(mu.mean = list(M %>% purrr::reduce(`+`) / length(M)))

  return(x$mu.mean[[1]])
}

#' @rdname get_expected_category_statistics
#' @export
get_expected_sigma = function(x, category, group) {
  assert_that(is.mv_ibbu_stanfit(x) | is.mv_ibbu_MCMC(x, nested = T, long = T))
  assert_that(is.character(category) | is.numeric(category))
  assert_that(is.character(group) | is.numeric(group))

  if (is.mv_ibbu_stanfit(x))
    x = add_ibbu_draws(x, which = "both", wide = F, nest = T)

  D = dim(g$S[[1]])[1]
  x %<>%
    filter(group == group, category == category) %>%
    mutate(Sigma = map2(S, nu, ~ .x / (.y - D - 1))) %>%
    summarise(Sigma.mean = list(Sigma %>% purrr::reduce(`+`) / length(Sigma)))

  return(x$Sigma.mean[[1]])
}
