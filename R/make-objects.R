#' Make multivariate Gaussian(s) from data.
#'
#' Constructs an \code{\link[=is.MVG]{MVG}} or \code{\link[=is.MVG_ideal_observer]{MVG_ideal_observer}} object with category
#' information for all categories found in the data.
#'
#'
#' @param data The tibble or data.frame from which to construct the prior.
#' @param groups Optionally, a group variable can be specified. If group is not NULL, one set of multivariate
#' Gaussian categories will be derived for each level of group. (default: NULL)
#' @param category Name of variable in \code{data} that contains the category information. (default: "category")
#' @param cues Name(s) of variables in \code{data} that contain the cue information.
#' @param prior Optionally specify a prior probability for each category (in each group). (default: a uniform
#' prior over all categories).
#' @param lapse_rate Optionally specify a lapse rate. (default: \code{NA})
#' @param bias Optionally specify a response bias. (default: \code{NA})
#' @param Sigma_noise Optionally specify a (multivariate Gaussian) noise covariance matrix. This argument will be
#' @param keep.category_parameters Should categories' mu and Sigma be included in the output (in addition to m
#' and S of the prior)? (default: FALSE)
#' @param verbose If true provides more information. (default: FALSE)
#'
#' @return A tibble that is an MVG object.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
#'
make_MVG_from_data = function(
  data,
  groups = NULL,
  category = "category",
  cues,
  keep.category_parameters = F,
  verbose = F
) {
  assert_that(is.data.frame(data) | is_tibble(data))
  assert_that(all(is.null(groups) | all(is.character(groups) | is_symbol(groups), length(groups) == 1)))
  assert_that(all(is.character(category) | is_symbol(category), length(category) == 1))
  assert_that(all(is.character(cues) | is_symbol(cues), length(cues) > 0))

  if (is.character(groups)) groups = syms(groups)
  if (is.character(category)) category = sym(category)
  if (is.character(cues)) cues = syms(cues)

  assert_that(as_name(category) %in% names(data),
              msg = paste0("Category variable (", as_name(category), ") not found in data."))
  assert_that(all(as_name(cues) %in% names(data)),
              msg = paste0("Some cues not found in data: ", setdiff(as_name(cues), intersect(as_name(cues), names(data)))))

  data %<>%
    select(!! category, !!! cues, !!! groups) %>%
    mutate(cues = pmap(list(!!! cues),
                       .f = function(...) {
                         x = c(...)
                         names(x) = as.character(cues)
                         return(x) })) %>%
    { if (is.null(groups)) group_by(., !! category) else group_by(., !!! groups, !! category) } %>%
    summarise(
      mu = list(reduce(cues, `+`) / length(cues)),
      Sigma = list(cov(cbind(!!! cues))))

  while(length(groups) > 1) {
    groups = groups[2:length(groups)]
    data %<>%
      group_by(., !!! groups, !! category) %>%
      summarise_at(vars(starts_with(c("mu", "Sigma"))),
                   ~ list(reduce(.x, `+`) / length(.x)))
  }

  data %<>%
    mutate(
      !! category := factor(!! category)) %>%
    ungroup()

  if (!keep.category_parameters) data %<>% select(-c(mu, Sigma))
  if (!is.MVG(data, category = as_name(category), verbose = verbose)) warning("Something went wrong. The returned object is not an MVG. Try again with verbose = T?")

  return(data)
}

#' @export
#' @rdname make_MVG_from_data
make_MVG_ideal_observer_from_data = function(
  data,
  groups = NULL,
  category = "category",
  cues,
  prior = NA_real_,
  lapse_rate = NA_real_,
  bias = NA_real_,
  Sigma_noise = NULL,
  keep.category_parameters = F,
  verbose = F
) {
  data <- make_MVG_from_data(data, groups = groups, category = category, cues = cues, keep.category_parameters = keep.category_parameters, verbose = verbose)

  if (is.character(groups)) groups = syms(groups)
  if (is.character(category)) category = sym(category)
  if (is.character(cues)) cues = syms(cues)
  assert_that(all(is.numeric(lapse_rate), is.numeric(bias), is.numeric(prior)),
              msg = "The category prior, lapse rate, and bias must be numeric.")

  n.cat = length(unique(data[[!! category]]))
  if (is.na(prior) | is.null(prior)) {
    message(paste0("No prior specified. Defaulting to uniform prior over the ", n.cat, " categories found in the data."))
    prior = rep(1 / n.cat, n.cat)
  }

  assert_that(all(
    is.na(lapse_rate) | between(lapse_rate, 0, 1),
    is.na(bias) | between(bias, 0, 1),
    between(prior, 0, 1)),
              msg = "If not NA, the category prior, lapse rate and bias must have values between 0 and 1.")
  assert_that(sum(prior) == 1,
              msg = paste0("Priors must add up to 1. (instead: ", sum(prior), ")."))

  assert_that(is.na(Sigma_noise) | is.matrix(Sigma_noise),
              msg = "If not NULL, Sigma_noise must be a matrix.")
  if (!is.na(Sigma_noise)) {
    assert_that(all(dim(Sigma_noise) == dim(first(data$Sigma))),
                msg = "If not NULL, Sigma_noise must be a matrix of the same dimensionality as Sigma.")
    assert_that(all(dimnames(Sigma_noise) == dimnames(first(data$Sigma))),
                msg = "If Sigma_noise is not NULL, the dimnames of Sigma_noise and Sigma must match.")
  }

  data %<>%
    mutate(
      prior = prior,
      lapse_rate = lapse_rate,
      bias = bias,
      Sigma_noise = list(Sigma_noise))

  return(data)
}



#' Make NIW belief from data.
#'
#' Constructs an \code{\link[=NIW_belief]{NIW_belief}} object, representing Normal-Inverse Wishart (NIW) parameters
#' for all categories found in the data. This object can be used as a prior for functions like
#' \code{\link{update_NIW_beliefs_incrementally}}.
#'
#' Currently, \code{make_NIW_prior_from_data()} does not infer kappa or nu, nor does it fit hierarchical data. Rather
#' the function simply estimates the category mean and covariance matrix from the sample (\code{data}), assumes them
#' to be the expected category mean (mu) and covariance (Sigma), and derives the m and S parameters
#' of the NIW from mu and Sigma based on the user-provided kappa and nu. That means m = mu and S = Sigma * (nu - D -1),
#' where D is the dimensionality of the data.
#'
#' @param data The tibble or data.frame from which to construct the prior.
#' @param groups Optionally, a group variable can be specified. If group is not NULL, one prior will be derived for
#' each level of group. If groups is a vector, then the data will be aggregated in the order of the specified
#' grouping variables: the first group variable will be used to obtain a mean and covariance matrix; the next group
#' variable will be used to obtain the average of those means and covariance matrices; and so on. For this process
#' to have a reasonable outcome, it is important that the group variables are hierarchically organized: the second
#' group variable should be a superset of the first group variable; the third group variable should be a superset
#' of the second group variable; etc. (default: NULL)
#' @param category Name of variable in \code{data} that contains the category information. (default: "category")
#' @param cues Name(s) of variables in \code{data} that contain the cue information.
#' @param kappa The strength of the beliefs over the category mean (pseudocounts). (default: same as nu)
#' @param nu The strength of the beliefs over the category covariance matrix (pseudocounts). (default: number of
#' cues + 2)
#' @param lapse_rate Optionally specify a lapse rate. (default: \code{NA})
#' @param Sigma_noise Optionally specify a (multivariate Gaussian) noise covariance matrix. This argument will be
#' ignored if `NULL`. (default: NULL)
#' @param keep.category_parameters Should categories' mu and Sigma be included in the output (in addition to m
#' and S of the prior)? (default: FALSE)
#'
#' @return A tibble that is an NIW_belief object.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
#'
make_NIW_belief_from_data = function(
  data,
  groups = NULL,
  category = "category",
  cues,
  kappa = nu,
  nu = length(cues) + 2,
  lapse_rate = NA_real_,
  Sigma_noise = NULL,
  keep.category_parameters = F
) {
  assert_that(is.data.frame(data) | is_tibble(data))
  assert_that(all(is.null(groups) | all(is.character(groups) | is_symbol(groups))))
  assert_that(all(is.character(category) | is_symbol(category), length(category) == 1))
  assert_that(all(is.character(cues) | is_symbol(cues), length(cues) > 0))
  assert_that(all(is.numeric(kappa), is.numeric(nu), !is.na(kappa), !is.na(nu)))
  assert_that(nu > length(cues) + 1,
              msg = paste0("nu must be larger than dimensionality of cues + 1 (>", length(cues) + 1, ")."))

  if (is.character(groups)) groups = syms(groups)
  if (is.character(category)) category = sym(category)
  if (is.character(cues)) cues = syms(cues)

  data %<>%
    select(!! category, !!! cues, !!! groups) %>%
    mutate(cues = pmap(list(!!! cues),
                       .f = function(...) {
                         x = c(...)
                         names(x) = as.character(cues)
                         return(x) })) %>%
    { if (is.null(groups)) group_by(., !! category) else group_by(., !!! groups, !! category) } %>%
    summarise(
      mu = list(reduce(cues, `+`) / length(cues)),
      Sigma = list(cov(cbind(!!! cues))))

  while(length(groups) > 1) {
    groups = groups[2:length(groups)]
    data %<>%
      group_by(., !!! groups, !! category) %>%
      summarise_at(vars(starts_with(c("mu", "Sigma"))),
                   ~ list(reduce(.x, `+`) / length(.x)))
  }

  message("S is set so that the expected category covariance matrix Sigma matches the category covariance in the sample (given nu). It might be safer to fit an Inverse-Wishart distribution to the entire set of covariance matrices.")
  data %<>%
    mutate(
      !! category := factor(!! category),
      kappa = kappa,
      nu = nu,
      m = mu,
      S = map2(Sigma, nu, get_S_from_Sigma),
      lapse_rate = lapse_rate) %>%
    { if (!is.null(Sigma_noise)) mutate(., Sigma_noise = list(Sigma_noise)) else . } %>%
    ungroup()

  if (!keep.category_parameters) data %<>% select(-c(mu, Sigma))
  if (!is.NIW_belief(data, category = as_name(category))) warning("Something went wrong. The returned object is not an NIW belief.")

  return(data)
}

#' @rdname make_NIW_belief_from_data
#' @export
make_NIW_prior_from_data <- make_NIW_belief_from_data


