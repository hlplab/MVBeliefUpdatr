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
#' @param x An mv_ibbu_stanfit or mv_ibbu_MCMC object.
#' @param category Character vector with categories (or category) for which category statistics are to be
#' returned.
#' @param group Character vector with groups (or group) for which category statistics are to be
#' returned.
#' @param statistic Either `mu` for category mean or `Sigma` for category covariance matrix.
#'
#' @seealso TBD
#' @keywords TBD
#' @references Murphy, K. P. (2012). Machine learning: a probabilistic perspective. MIT press.
#' @examples
#' TBD
#' @rdname get_expected_category_statistics
#' @export
get_expected_category_statistics = function(x, category, group, statistic = c("mu", "Sigma")) {
  assert_that(statistic %in% c("mu", "Sigma"))
  assert_that(is.mv_ibbu_stanfit(x) | is.mv_ibbu_MCMC(x, nested = T, long = T))
  if (is.mv_ibbu_stanfit(x))
    x = add_ibbu_draws(x, which = "both", wide = F, nest = T)

  if (is.null(category)) category = unique(x$category)
  if (is.null(group)) group = unique(x$group)
  assert_that(is.character(category) | is.numeric(category))
  assert_that(is.character(group) | is.numeric(group))

  x %<>%
    filter(group %in% !! group, category %in% !! category)

  if (statistic == "mu") {
    x %<>%
      group_by(group, category) %>%
      summarise(mu.mean = list(M %>% purrr::reduce(`+`) / length(M)))
  } else if (statistic == "Sigma") {
    D = dim(x$S[[1]])[1]
    x %<>%
      mutate(Sigma = map2(S, nu, ~ .x / (.y - D - 1))) %>%
      group_by(group, category) %>%
      summarise(Sigma.mean = list(Sigma %>% purrr::reduce(`+`) / length(Sigma)))
  } else stop("statistic must be either mu or Sigma.")

  # If just one category and group was requested, just return that object, rather
  # than the tibble
  if (nrow(x) == 1) x = x[,paste0(statistic, ".mean")][[1]][[1]]
  return(x)
}

#' @rdname get_expected_category_statistics
#' @export
get_expected_mu = function(x, category, group) {
  return(get_expected_category_statistics(x, category, group, statistic = "mu"))
}

#' @rdname get_expected_category_statistics
#' @export
get_expected_sigma = function(x, category, group) {
  return(get_expected_category_statistics(x, category, group, statistic = "Sigma"))
}
