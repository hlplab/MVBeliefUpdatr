#' @importFrom mvtnorm dmvt
#' @importFrom foreach foreach %do%
#' @importFrom rlang is_scalar_double
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
#' @importFrom purrr map_lgl
#' @rdname get_expected_mu_from_m
#' @export
get_expected_mu_from_m = function(m) {
  if (!is.list(m)) m <- list(m) # in case the input is not a list
  if (all(map_lgl(m, ~ length(.x) == 1))) mu <- unlist(m) else mu <- m
  if (length(mu) == 1) mu <- mu[[1]]
  return(mu)
}

#' Get mean of means from expected category mean mu
#'
#' @rdname get_expected_mu_from_m
#' @export
get_m_from_expected_mu = function(mu) {
  if (!is.list(mu)) mu <- list(mu) # in case the input is not a list
  if (all(map_lgl(mu, ~ length(.x) == 1))) m <- unlist(mu) else m <- mu
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
#' @importFrom purrr map_lgl
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

  if (all(map_lgl(Sigma, ~ length(.x) == 1))) Sigma <- unlist(Sigma)
  if (length(Sigma) == 1) Sigma <- Sigma[[1]]
  return(Sigma)
}

#' Get Scatter matrix S from expected category covariance Sigma and pseudocount nu
#'
#' @rdname get_expected_Sigma_from_S
#' @export
get_S_from_expected_Sigma <- function(Sigma, nu) {
  if (!is.list(Sigma)) {
    # in case the input is not a list
    Sigma <- list(Sigma)
    nu <- list(nu)
  }

  S <-  map2(Sigma, nu, .f = function(Sigma, nu) {
    D = get_D(Sigma)
    return(Sigma * (nu - D - 1))
  })

  # if (all(map_lgl(S, ~ length(.x) == 1))) {
  #   S <- map(
  #     S,
  #     function(x) {
  #       x <- as.vector(x)
  #       names(x) <- dimnames(Sigma[[1]])[1][[1]]
  #       return(x)
  #     })
  # }
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
#' @param noise_treatment Determines whether perceptual noise is considered during categorization, and how.
#' Can be "no_noise", "sample", or "marginalize". If "no_noise", no noise will be applied to the input,
#' and no noise will be assumed during categorization. If "marginalize", average noise (i.e., no noise)
#' will be added to the stimulus, and `Sigma_noise` is added to Sigma when calculating the likelihood.
#' This simulates the expected consequences for perceptual noise on categorization *in the limit*, i.e,
#' if the input was categorized infinitely many times. If "sample", then noise is sampled and applied to
#' the input, and `Sigma_noise` is added to Sigma when calculating the likelihood. This simulates the
#' consequence of perceptual noise *on a particular observation*. If "sample" or "marginalize" are chosen,
#' `Sigma_noise` must be a covariance matrix of appropriate dimensions. (default: "no_noise" if Sigma_noise
#' is NULL, "marginalize" otherwise).
#' @param log Should the log-transformed density be returned (`TRUE`)? (default: `TRUE`)
#'
#' @seealso TBD
#' @keywords TBD
#' @references \insertRef{murphy2012}{MVBeliefUpdatr}
#'
#' @rdname get_NIW_posterior_predictive
#' @export
get_NIW_posterior_predictive = function(
    x, m, S, kappa, nu, Sigma_noise = NULL,
    noise_treatment = if (any(is.null(Sigma_noise), all(is.null(Sigma_noise)), all(map_lgl(Sigma_noise, is.null)))) "no_noise" else "marginalize",
    log = T
) {
  # mvtnorm::dmvt expects means to be vectors, and x to be either a vector or a matrix.
  # in the latter case, each *row* of the matrix is an input.
  assert_that(is.vector(m) | is.matrix(m) | is_scalar_double(m))
  assert_that(is.matrix(S) | is_scalar_double(S))
  # do not reorder these conditionals (go from more to less specific)
  if (is.matrix(m)) m <- as.vector(m)

  x %<>% format_input_for_likelihood_calculation(dim = length(m))
  assert_that(dim(x)[2] == length(m),
              msg = "Input x and m are not of compatible dimensions.")

  assert_that(all(is.number(kappa), is.number(nu)))
  assert_that(is.flag(log))
  assert_that(any(noise_treatment %in% c("no_noise", "sample", "marginalize")),
              msg = "noise_treatment must be one of 'no_noise', 'sample' or 'marginalize'.")
  if (noise_treatment != "no_noise") {
    assert_that(is.Sigma(Sigma_noise))
    assert_that(all(dim(S) == dim(Sigma_noise)),
                msg = 'If noise_treatment is not "no_noise", Sigma_noise must be a covariance matrix of appropriate dimensions, matching those of the scatter matrices S.')
  }

  D <- get_D(S)
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
  }


  if (noise_treatment %in% c("sample", "marginalize")) {
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


#' @rdname get_NIW_posterior_predictive
#' @export
get_posterior_predictive_from_NIW_belief = function(
  x,
  model,
  noise_treatment = if (is.NIW_ideal_adaptor(model)) { if (!is.null(first(model$Sigma_noise))) "marginalize" else "no_noise" } else "no_noise",
  log = T,
  category = "category",
  category.label = NULL,
  wide = FALSE
) {
  assert_that(is.NIW_belief(model))
  assert_that(any(is.null(category.label) | is.character(category.label)))
  assert_that(any(noise_treatment == "no_noise", is.NIW_ideal_adaptor(model)),
              msg = 'No noise matrix Sigma_noise found. If noise_treatment is not "no_noise", then model must be an NIW_ideal_adaptor.')

  if (is.null(category.label)) {
    model %<>%
      droplevels()
    category.label = model %>%
      pull(!! sym(category)) %>%
      unique()
  }

  posterior_predictive <- foreach(c = category.label) %do% {
    m <-
      model %>%
      filter(!! sym(category) == c)

    get_NIW_posterior_predictive(
      x = x,
      m = m$m[[1]],
      S = m$S[[1]],
      kappa = m$kappa[[1]],
      nu = m$nu[[1]],
      log = log,
      noise_treatment = noise_treatment,
      Sigma_noise = if (noise_treatment == "no_noise") NULL else m$Sigma_noise[[1]]) %>%
      as_tibble(.name_repair = "minimal") %>%
      rename_with(~ if (log) { "log_posterior_predictive" } else { "posterior_predictive" }) %>%
      mutate(!! sym(category) := c)
  }

  posterior_predictive %<>% reduce(rbind)
  if (wide)
    posterior_predictive %<>%
    pivot_wider(
      values_from = if (log) "log_posterior_predictive" else "posterior_predictive",
      names_from = !! sym(category),
      names_prefix = if (log) "log_posterior_predictive." else "posterior_predictive.") %>%
    unnest()

  return(posterior_predictive)
}

# If there's a grouping variable extract the posterior predictive for each level of that grouping variable
get_posterior_predictives_from_NIW_beliefs = function(
  x,
  model,
  noise_treatment = if (is.NIW_ideal_adaptor(model)) { if (!is.null(first(model$Sigma_noise))) "marginalize" else "no_noise" } else "no_noise",
  log = T,
  category = "category",
  category.label = NULL,
  grouping.var,
  wide = FALSE
) {
  if (is.null(grouping.var)) {
    return(get_posterior_predictive_from_NIW_belief(
      x,
      model,
      log = log,
      category = category,
      category.label = category.label,
      wide = wide))
  } else {
    assert_that(grouping.var %in% names(x),
                msg = "Grouping variable not found in the NIW belief object.")

    foreach (i = unique(x[[grouping.var]])) %do% {
      posterior_predictive <-
        get_posterior_predictive_from_NIW_belief(
          x,
          model %>% filter(!! sym(grouping.var) == i),
          log = log,
          category = category,
          category.label = category.label,
          wide = wide
        ) %>%
        mutate(!! sym(grouping.var) := i)
    }

    posterior_predictive %<>% reduce(rbind)
    return(posterior_predictive)
  }
}



