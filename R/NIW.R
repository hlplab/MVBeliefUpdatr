#' @import assertthat
#' @importFrom mvtnorm dmvt
NULL


#' Get posterior predictive
#'
#' Get posterior predictive of observation x given the NIW parameters M, S, kappa, and nu. This is
#' a multivariate Student-T distribution (Murphy, 2012, p. 135).
#'
#' @param x Input (observation).
#' @param M The mean of the multivariate Normal distribution of the category mean mu.
#' @param S The scatter matrix of the inverse-Wishart distribution over the category covariance
#' matrix Sigma.
#' @param kappa The strength of the beliefs over the category mean (pseudocounts).
#' @param nu The strength of the beliefs over the category covariance matrix (pseudocounts).
#' @param log Should the log-transformed density be returned (`TRUE`)? (default: `TRUE`)
#'
#' @seealso TBD
#' @keywords TBD
#' @references Murphy, K. P. (2012). Machine learning: a probabilistic perspective. MIT press.
#' @examples
#' TBD
#' @export
get_posterior_predictive = function(x, M, S, kappa, nu, log = T) {
  # mvtnorm::dmvt now expects means to be vectors, and x to be either a vector or a matrix.
  # in the latter case, each *row* of the matrix is an input.
  assert_that(all(is.vector(x) | is.matrix(M), is.vector(M) | is.matrix(M), is.matrix(S)))
  if (is.vector(x)) x = matrix(M, nrow = 1)
  if (is.matrix(M)) M = as.vector(M)

  assert_that(all(is.number(kappa), is.number(nu)))
  assert_that(is.flag(log))
  D = dim(S)[1]
  assert_that(nu >= D,
              msg = "nu must be at least as large as the number of dimensions of the multivariate
              Normal.")
  assert_that(dim(S)[2] == D,
              msg = "S is not a square matrix, and thus not a Scatter matrix")
  assert_that(length(M) == dim(x)[2],
              msg = paste("M and input are not of compatible dimensions. M is of length", length(M), "but input has", dim(x)[2], "columns."))
  assert_that(length(M) == D,
              msg = "S and M are not of compatible dimensions.")

  mvtnorm::dmvt(x,
                delta = M,
                sigma = S * (kappa + 1) / (kappa * (nu - D + 1)),
                df = nu - D + 1,
                log = log)
}

get_posterior_predictive.pmap = function(x, M, S, kappa, nu, ...) {
  get_posterior_predictive(x, M, S, kappa, nu, log = F)
}


#' Get expected category covariance from Scatter matrix and nu
#' @export
get_Sigma_from_S = function(S, nu) {
  return(S / (nu - dim(S)[1] - 1))
}


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
#' returned.  If `NULL` then all categories are included. (default: `NULL`)
#' @param group Character vector with groups (or group) for which category statistics are to be
#' returned. If `NULL` then all groups are included. (default: `NULL`)
#' @param statistic Which category statistic should be returned? `mu` for category mean or `Sigma` for category
#' covariance matrix, or `c("mu", "Sigma")` for both. (default: both)
#'
#' @seealso TBD
#' @keywords TBD
#' @references Murphy, K. P. (2012). Machine learning: a probabilistic perspective. MIT press.
#' @examples
#' TBD
#' @rdname get_expected_category_statistic
#' @export
get_expected_category_statistic = function(x, category = NULL, group = NULL,
                                           statistic = c("mu", "Sigma")) {
  assert_that(all(statistic %in% c("mu", "Sigma")))
  assert_that(is.mv_ibbu_stanfit(x) | is.mv_ibbu_MCMC(x, nested = T, long = T))
  if (is.mv_ibbu_stanfit(x))
    x = add_ibbu_draws(x, which = "both", wide = F, nest = T)

  assert_that(any(is.null(category), is.character(category), is.numeric(category)))
  assert_that(any(is.null(group), is.character(group), is.numeric(group)))
  if (is.null(category)) category = unique(x$category)
  if (is.null(group)) group = unique(x$group)

  x %<>%
    filter(group %in% !! group, category %in% !! category) %>%
    mutate(Sigma = map2(S, nu, get_Sigma_from_S)) %>%
    group_by(group, category) %>%
    summarise(
      mu.mean = list(M %>% purrr::reduce(`+`) / length(M)),
      Sigma.mean = list(Sigma %>% purrr::reduce(`+`) / length(Sigma))) %>%
    select(group, category, !!! rlang::syms(paste0(statistic, ".mean")))

  if (!all(sort(unique(as.character(x$group))) == sort(as.character(group))))
    warning("Not all groups were found in x.")

  # If just one category and group was requested, just return that object, rather
  # than the tibble
  if (nrow(x) == 1) x = x[,paste0(statistic, ".mean")][[1]][[1]]
  return(x)
}

#' @rdname get_expected_category_statistic
#' @export
get_expected_mu = function(x, category, group) {
  return(get_expected_category_statistic(x, category, group, statistic = "mu"))
}

#' @rdname get_expected_category_statistic
#' @export
get_expected_sigma = function(x, category, group) {
  return(get_expected_category_statistic(x, category, group, statistic = "Sigma"))
}
