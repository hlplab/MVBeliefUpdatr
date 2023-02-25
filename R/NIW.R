#' @importFrom mvtnorm dmvt
#' @importFrom foreach foreach %do%
NULL


#' Get dimensionality of mean or covariance matrix
#'
#' @param x A mean or covariance matrix, or lists of means or covariance matrices. If x is a list, then the
#' dimensionality of the first list element is determined. This function does not currently check whether all
#' list elements have the same dimensionality.
#'
#' @export
get_D <- function(x) {
  if (is.list(x)) {
    return(get_D(x[[1]]))
  } else if (is.matrix(x)) {
    # could be mean or cov, so take max of dims
    return(max(dim(x)))
  } else if (is.vector(x)) {
    return(length(x))
  } else if (is.numeric(x) & length(x) == 1) {
    return(1)
  } else stop("Cannot recognize dimensionality of x.")
}

#' Get expected category mu from mean of means m
#'
#' See Murphy (2012, p. 134).
#'
#' @param mu expected category mean mu.
#' @param m Mean of means m.
#'
#' @rdname get_expected_mu_from_m
#' @export
get_expected_mu_from_m = function(m) {
  if (!is.list(m)) m <- list(m) # in case the input is not a list
  if (all(unlist(map(m, ~ length(.x) == 1)))) mu <- unlist(m) else mu <- m
  if (length(mu) == 1) mu <- mu[[1]]
  return(mu)
}

#' Get mean of means from expected category mean mu
#'
#' @rdname get_expected_mu_from_m
#' @export
get_m_from_expected_mu = function(mu) {
  if (!is.list(mu)) mu <- list(mu) # in case the input is not a list
  if (all(unlist(map(mu, ~ length(.x) == 1)))) m <- unlist(mu) else m <- mu
  if (length(m) == 1) m <- m[[1]]
  return(m)
}

#' Get expected category covariance from Scatter matrix S and pseudocount nu
#'
#' See Murphy (2012, p. 134).
#'
#' @param Sigma expected category covariance matrix.
#' @param S Scatter matrix S.
#' @param nu Strength of belief (pseudocount) about Sigma.
#'
#' @rdname get_expected_Sigma_from_S
#' @export
get_expected_Sigma_from_S = function(S, nu) {
  if (!is.list(S)) {
    # in case the input is not a list
    S <- list(S)
    nu <- list(nu)
  }

  Sigma = map2(S, nu, .f = function(S, nu) {
    D = get_D(S)
    return(S / (nu - D - 1))
  })

  if (all(unlist(map(Sigma, ~ length(.x) == 1)))) Sigma <- unlist(Sigma)
  if (length(Sigma) == 1) Sigma <- Sigma[[1]]
  return(Sigma)
}

#' Get Scatter matrix S from expected category covariance Sigma and pseudocount nu
#'
#' @rdname get_expected_Sigma_from_S
#' @export
get_S_from_expected_Sigma = function(Sigma, nu) {
  if (!is.list(Sigma)) {
    # in case the input is not a list
    Sigma <- list(Sigma)
    nu <- list(nu)
  }

  S = map2(Sigma, nu, .f = function(Sigma, nu) {
    D = get_D(Sigma)
    return(Sigma * (nu - D - 1))
  })

  if (all(unlist(map(S, ~ length(.x) == 1)))) S <- unlist(S)
  if (length(S) == 1) S <- S[[1]]
  return(S)
}


#' Get posterior predictive
#'
#' Get posterior predictive of observations x given the NIW parameters m, S, kappa, and nu. This is
#' the density of a multivariate Student-T distribution \insertCite{@see @murphy2012 p. 134}{MVBeliefUpdatr}.
#'
#' @param x Observation(s). Can be a vector with k elements for a single observation or a matrix with k
#' columns and n rows, in which case each row of the matrix is taken to be one observation. If x is a
#' tibble with k columns or a list of vectors of length k, it is reduced into a matrix with k columns.
#' @param m The mean of the multivariate Normal distribution of the category mean mu. Should be a
#' matrix or vector of length k.
#' @param S The scatter matrix of the inverse-Wishart distribution over the category covariance
#' matrix Sigma. Should be a square k x k matrix.
#' @param kappa The strength of the beliefs over the category mean (pseudocounts).
#' @param nu The strength of the beliefs over the category covariance matrix (pseudocounts).
#' @param Sigma_noise Optionally, a covariance matrix describing the perceptual noise to be applied while
#' calculating the posterior predictive. (default: `NULL`)
#' @param noise_treatment Determines whether and how multivariate Gaussian noise is added to the input. Can be "no_noise", "sample"
#' or "marginalize". If "no_noise", no noise will be applied. This describes the idealized categorization function if the percept
#' was known. If "sample" or "marginalize", `Sigma_noise` must be a covariance
#' matrix of appropriate dimensions. If "sample", observations are adjusted by samples drawn from the noise distribution before applying
#' categorization. If "marginalize" then each observation is transformed into the marginal distribution
#' that results from convolving the input with noise. This latter option might be helpful, for example, if one is
#' interested in estimating the consequences of noise across individuals. (default: "no_noise").
#' @param log Should the log-transformed density be returned (`TRUE`)? (default: `TRUE`)
#'
#' @seealso TBD
#' @keywords TBD
#' @references \insertRef{murphy2012}{MVBeliefUpdatr}
#' @examples
#' TBD
#' @rdname get_NIW_posterior_predictive
#' @export
get_NIW_posterior_predictive = function(x, m, S, kappa, nu, Sigma_noise = NULL, noise_treatment = "no_noise", log = T) {
  # mvtnorm::dmvt expects means to be vectors, and x to be either a vector or a matrix.
  # in the latter case, each *row* of the matrix is an input.
  assert_that(is.vector(x) | is.matrix(x) | is_tibble(x) | is.list(x))
  assert_that(is.vector(m) | is.matrix(m) | is_scalar_double(m))
  assert_that(is.matrix(S) | is_scalar_double(S))
  # do not reorder these conditionals (go from more to less specific)
  if (is_tibble(x)) x %<>% as.matrix() else
    if (is.list(x)) x %<>% reduce(rbind) %<>% as.matrix(nrow = 1) else
      if (is.vector(x)) x %<>% matrix(nrow = 1)
  if (is.matrix(m)) m = as.vector(m)

  assert_that(all(is.number(kappa), is.number(nu)))
  assert_that(is.flag(log))
  assert_that(any(noise_treatment %in% c("no_noise", "sample", "marginalize")),
              msg = "noise_treatment must be one of 'no_noise', 'sample' or 'marginalize'.")
  if (noise_treatment != "no_noise") {
    assert_that(is.Sigma(Sigma_noise))
    assert_that(all(dim(S) == dim(Sigma_noise)),
                msg = 'If noise_treatment is not "no_noise", Sigma_noise must be a covariance matrix of appropriate dimensions, matching those of the scatter matrices S.')
  }

  D = get_D(S)
  assert_that(nu >= D,
              msg = "nu must be at least as large as the number of dimensions of the multivariate
              Normal.")

  if (D == 1) {
    assert_that(is_scalar_double(m), msg = "S and m are not of compatible dimensions.")
  } else {
    assert_that(dim(S)[2] == D,
                msg = "S is not a square matrix, and thus not a Scatter matrix")
    assert_that(length(m) == dim(x)[2],
                msg = paste("m and input are not of compatible dimensions. m is of length", length(m), "but input has", dim(x)[2], "columns."))
    assert_that(length(m) == D,
                msg = "S and m are not of compatible dimensions.")
  }

  # How should noise be treated?
  if (noise_treatment == "sample") {
    assert_that(
      is_weakly_greater_than(length(x), 1),
      msg = "For noise sampling, x must be of length 1 or longer.")

    x <- map(x, ~ rmvnorm(n = 1, mean = .x, sigma = Sigma_noise))
  } else if (noise_treatment == "marginalize") {
    warning("noise_treatment == 'marginalize' is experimental. The math has not yet been verified. Use with caution.")
    S = get_S_from_expected_Sigma(get_expected_Sigma_from_S(S, nu) + Sigma_noise, nu)
  }

  dmvt(x,
       delta = m,
       sigma = S * ((kappa + 1) / (kappa * (nu - D + 1))),
       df = nu - D + 1,
       log = log)
}


#' @rdname get_NIW_posterior_predictive
#' @export
get_NIW_posterior_predictive.pmap = function(x, m, S, kappa, nu, ...) {
  get_NIW_posterior_predictive(x = x, m = m, S = S, kappa = kappa, nu = nu, ...)
}


#' Get NIW categorization function
#'
#' Returns a categorization function for the first category, based on a set of parameters for the Normal-Inverse-Wishart (NIW)
#' distribution. ms, Ss, kappas, nus, and priors are assumed to be of the same length and sorted the same way, so that the first
#' element of ms is corresponding to the same category as the first element of Ss, kappas, nus, and priors, etc.
#'
#' @param ms Means of the multivariate normal distributions over category means.
#' @param Ss Scatter matrices of the inverse Wishart distribution over category covariance matrices.
#' @param kappas Strength of the beliefs into the distribution over category means.
#' @param nus Strength of the beliefs into the distribution over category covariance matrices.
#' @param priors Vector of categories' prior probabilities. (default: uniform prior over categories)
#' @param lapse_rate A lapse rate for the categorization responses.
#' @param lapse_biases A lapse bias for the categorization responses. (default: uniform bias over categories)
#' @param Sigma_noise A noise matrix. (default: a 0-matrix)
#' @param noise_treatment How should the noise specified in \code{Sigma_noise} be considered in the categorization function?
#' For details, see \code{\link{get_NIW_posterior_predictive}} (default: "no_noise")
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
#' @rdname get_NIW_categorization_function
#' @export
get_NIW_categorization_function = function(
  ms, Ss, kappas, nus,
  priors = rep(1 / length(ms), length(ms)),
  lapse_rate = NULL,
  lapse_biases = rep(1 / length(ms), length(ms)),
  Sigma_noise = matrix(
    0,
    nrow = if (is.null(dim(Ss[[1]]))) 1 else max(dim(Ss[[1]])),
    ncol = if (is.null(dim(Ss[[1]]))) 1 else max(dim(Ss[[1]]))),
  noise_treatment = "no_noise",
  logit = FALSE
) {
  tolerance = 1e-5
  assert_that(are_equal(length(ms), length(Ss)),
              are_equal(length(ms), length(priors)),
              are_equal(length(ms), length(kappas)),
              are_equal(length(ms), length(nus)),
              msg = "The number of ms, Ss, kappas, nus, and priors must be identical.")
  n.cat = length(ms)

  assert_that(all(between(priors, 0, 1), between(sum(priors), 1 - tolerance, 1 + tolerance)),
              msg = "priors must sum to 1.")
  if (!is.null(lapse_rate)) assert_that(all(between(lapse_rate, 0, 1))) else lapse_rate = rep(0, n.cat)
  if (!is.null(lapse_biases)) {
    assert_that(all(between(lapse_biases, 0, 1), between(sum(lapse_biases), 1 - tolerance, 1 + tolerance)),
                msg = "biases must sum to 1.")
  } else lapse_biases <- 1 / n.cat

  # Get dimensions of multivariate category
  D = get_D(ms)
  assert_that(
    nus[[1]] >= D,
    msg = "Nu must be at least K (number of dimensions of the multivariate Gaussian category).")

  f <- function(x, target_category = 1) {
    log_p = matrix(
      nrow = length(x),
      ncol = n.cat
    )
    for (cat in 1:n.cat) {
      log_p[, cat] = get_NIW_posterior_predictive(
        x,
        ms[[cat]], Ss[[cat]], kappas[[cat]], nus[[cat]],
        Sigma_noise = Sigma_noise[[cat]], noise_treatment = noise_treatment,
        log = T)
    }

    p_target <-
      (1 - lapse_rate[target_category]) * exp(log_p[,target_category] + log(priors[target_category]) - log(rowSums(exp(log_p) * priors))) +
      lapse_rate[target_category] * lapse_biases[target_category]

    if (logit)
      return(qlogis(p_target))
    else
      return(p_target)
  }

  return(f)
}


