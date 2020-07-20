#' @importFrom mvtnorm dmvt
#' @importFrom foreach foreach %do%
NULL


#' Get posterior predictive
#'
#' Get posterior predictive of observations x given the NIW parameters M, S, kappa, and nu. This is
#' a multivariate Student-T distribution (Murphy, 2012, p. 135).
#'
#' @param x Observations. Can be a vector with k elements for a single observation or a matrix with k
#' columns and n rows, in which case each row of the matrix is taken to be one observation. If x is a
#' tibble with k columns or a list of vectors of length k, it is reduced into a matrix with k columns.
#' @param M The mean of the multivariate Normal distribution of the category mean mu. Should be a
#' matrix or vector of length k.
#' @param S The scatter matrix of the inverse-Wishart distribution over the category covariance
#' matrix Sigma. Should be a square k x k matrix.
#' @param kappa The strength of the beliefs over the category mean (pseudocounts).
#' @param nu The strength of the beliefs over the category covariance matrix (pseudocounts).
#' @param log Should the log-transformed density be returned (`TRUE`)? (default: `TRUE`)
#'
#' @seealso TBD
#' @keywords TBD
#' @references Murphy, K. P. (2012). Machine learning: a probabilistic perspective. MIT press.
#' @examples
#' TBD
#' @rdname get_posterior_predictive
#' @export
get_posterior_predictive = function(x, M, S, kappa, nu, log = T) {
  # mvtnorm::dmvt now expects means to be vectors, and x to be either a vector or a matrix.
  # in the latter case, each *row* of the matrix is an input.
  assert_that(all(is.vector(x) | is.matrix(x) | is_tibble(x) | is.list(x), is.vector(M) | is.matrix(M), is.matrix(S)))
  # do not reorder these conditionals (go from more to less specific)
  if (is_tibble(x)) x %<>% as.matrix() else
    if (is.list(x)) x %<>% reduce(rbind) else
      if (is.vector(x)) x %<>% matrix(nrow = 1)

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

  dmvt(x,
       delta = M,
       sigma = S * (kappa + 1) / (kappa * (nu - D + 1)),
       df = nu - D + 1,
       log = log)
}


#' @rdname get_posterior_predictive
#' @export
get_posterior_predictive.pmap = function(x, M, S, kappa, nu, ...) {
  get_posterior_predictive(x, M, S, kappa, nu, log = F)
}

#' @export
get_posterior_predictive_from_NIW_belief = function(
  x,
  belief,
  log = T,
  category = "category",
  category.label = NULL,
  wide = FALSE
) {
  assert_that(is.NIW_belief(belief))
  assert_that(any(is.null(category.label) | is.character(category.label)))

  if (is.null(category.label)) {
    belief %<>%
      droplevels()
    category.label = belief %>%
      pull(!! sym(category)) %>%
      unique()
  }

  pp = foreach(c = category.label) %do% {
    b = belief %>%
      filter(!! sym(category) == c)

    get_posterior_predictive(
      x = x,
      M = b$M[[1]], S = b$S[[1]], kappa = b$kappa[[1]], nu = b$nu[[1]], log = log) %>%
      as_tibble() %>%
      rename_all(~ if (log) "lpp" else "pp") %>%
      mutate(!! sym(category) := c)
  }

  pp = reduce(pp, rbind)
  if (wide)
    pp %<>%
    pivot_wider(
      values_from = if (log) "lpp" else "pp",
      names_from = !! sym(category),
      names_prefix = if (log) "lpp." else "pp.") %>%
    unnest()

  return(pp)
}


#' Get expected category covariance from Scatter matrix S and pseudocount nu
#' @export
get_Sigma_from_S = function(S, nu) {
  return(S / (nu - dim(S)[1] - 1))
}

#' Get Scatter matrix S from expected category covariance Sigma and pseudocount nu
#' @export
get_S_from_Sigma = function(Sigma, nu) {
  return(Sigma * (nu - dim(Sigma)[1] - 1))
}

