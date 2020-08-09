#' @importFrom mvtnorm dmvt
#' @importFrom foreach foreach %do%
NULL


#' Get posterior predictive
#'
#' Get posterior predictive of observations x given the NIW parameters m, S, kappa, and nu. This is
#' a multivariate Student-T distribution \insertCite{@see @murphy2012 p. 134}{MVBeliefUpdatr}.
#'
#' @param x Observations. Can be a vector with k elements for a single observation or a matrix with k
#' columns and n rows, in which case each row of the matrix is taken to be one observation. If x is a
#' tibble with k columns or a list of vectors of length k, it is reduced into a matrix with k columns.
#' @param m The mean of the multivariate Normal distribution of the category mean mu. Should be a
#' matrix or vector of length k.
#' @param S The scatter matrix of the inverse-Wishart distribution over the category covariance
#' matrix Sigma. Should be a square k x k matrix.
#' @param kappa The strength of the beliefs over the category mean (pseudocounts).
#' @param nu The strength of the beliefs over the category covariance matrix (pseudocounts).
#' @param log Should the log-transformed density be returned (`TRUE`)? (default: `TRUE`)
#'
#' @seealso TBD
#' @keywords TBD
#' @references \insertRef{murphy2012}{MVBeliefUpdatr}
#' @examples
#' TBD
#' @rdname get_posterior_predictive
#' @export
get_posterior_predictive = function(x, m, S, kappa, nu, log = T) {
  # mvtnorm::dmvt now expects means to be vectors, and x to be either a vector or a matrix.
  # in the latter case, each *row* of the matrix is an input.
  assert_that(all(is.vector(x) | is.matrix(x) | is_tibble(x) | is.list(x), is.vector(m) | is.matrix(m), is.matrix(S)))
  # do not reorder these conditionals (go from more to less specific)
  if (is_tibble(x)) x %<>% as.matrix() else
    if (is.list(x)) x %<>% reduce(rbind) else
      if (is.vector(x)) x %<>% matrix(nrow = 1)

  if (is.matrix(m)) m = as.vector(m)

  assert_that(all(is.number(kappa), is.number(nu)))
  assert_that(is.flag(log))
  D = dim(S)[1]
  assert_that(nu >= D,
              msg = "nu must be at least as large as the number of dimensions of the multivariate
              Normal.")
  assert_that(dim(S)[2] == D,
              msg = "S is not a square matrix, and thus not a Scatter matrix")
  assert_that(length(m) == dim(x)[2],
              msg = paste("m and input are not of compatible dimensions. m is of length", length(m), "but input has", dim(x)[2], "columns."))
  assert_that(length(m) == D,
              msg = "S and m are not of compatible dimensions.")

  dmvt(x,
       delta = m,
       sigma = S * (kappa + 1) / (kappa * (nu - D + 1)),
       df = nu - D + 1,
       log = log)
}


#' Get categorization function
#'
#' Returns a categorization function for the first category, based on a set of parameters for the Normal-inverse-wishart (NIW)
#' distribtuion. ms, Ss, kappas, nus, and priors are assumed to be of the same length and sorted the same way, so that the first
#' element of ms is corresponding to the same category as the first element of Ss, kappas, nus, and priors, etc.
#'
#' @param ms List of IBBU-inferred means describing the multivariate normal distribution over category means.
#' @param Ss List of IBBU-inferred scatter matrices describing the inverse Wishart distribution over category
#' covariance matrices.
#' @param kappas List of IBBU-inferred kappas describing the strength of the beliefs into the distribution over catgory means.
#' @param nus List of IBBU-inferred nus describing the strength of the beliefs into the distribution over catgory covariance matrices.
#' @param lapse_rate An IBBU-inferred lapse rate for the categorization responses.
#' @param priors Vector of categories' prior probabilities. (default: uniform prior over categories)
#' @param n.cat Number of categories. Is inferred from the input, but can be set manually.
#' @param logit Should the function that is returned return log-odds (TRUE) or probabilities (FALSE)? (default: TRUE)
#'
#' @return A function that takes as input cue values and returns posterior probabilities of the first category,
#' based on the posterior predictive of the cues given the (IBBU-derived parameters for the) categories' m, S,
#' kappa, nu, and prior, as well as the lapse rate.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @rdname get_categorization_function
#' @export
get_categorization_function = function(
  ms, Ss, kappas, nus,
  lapse_rate = NULL, bias = NULL,
  priors = rep(1 / n.cat, n.cat),
  n.cat = length(ms),
  logit = FALSE
) {
  assert_that(are_equal(length(ms), length(Ss)),
              are_equal(length(ms), length(priors)),
              are_equal(length(ms), length(kappas)),
              are_equal(length(ms), length(nus)),
              msg = "The number of ms, Ss, kappas, nus, and priors must be identical.")
  if (!is.null(lapse_rate)) assert_that(all(between(lapse_rate, 0, 1))) else lapse_rate = rep(0, n.cat)
  if (!is.null(bias)) assert_that(all(between(bias, 0, 1), sum(bias) == lapse_rate[1]),
                                  msg = "biases must sum up to lapse rate.") else bias = lapse_rate / n.cat

  # Get dimensions of multivariate category
  D = length(ms[[1]])
  assert_that(nus[[1]] >= D,
              msg = "Nu must be at least K (number of dimensions of the multivariate Gaussian category).")

  f <- function(x, target_category = 1) {
    log_p = matrix(
      nrow = length(x),
      ncol = n.cat
    )
    for (cat in 1:n.cat) {
      log_p[, cat] = get_posterior_predictive(x, ms[[cat]], Ss[[cat]], kappas[[cat]], nus[[cat]], log = T)
    }

    p_target =
      exp(
        log_p[,target_category] + log(priors[target_category]) -
          log(rowSums(exp(log_p) * priors))) *
      (1 - lapse_rate[target_category]) + bias[target_category]

    if (logit)
      return(qlogis(p_target))
    else
      return(p_target)
  }

  return(f)
}

#' @rdname get_posterior_predictive
#' @export
get_posterior_predictive.pmap = function(x, m, S, kappa, nu, ...) {
  get_posterior_predictive(x, m, S, kappa, nu, log = F)
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

