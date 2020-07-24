#' @importFrom dplyr sample_frac
#' @importFrom magrittr is_weakly_greater_than
#' @importFrom rlang is_scalar_integerish
#' @importFrom mvtnorm rmvnorm
NULL

#' Example NIW priors.
#'
#' @export
example_NIW_prior = function(example = 1) {
  if (example == 1) {
    message("An example belief for two categories in a 2D cue continuum that differ in means and correlatation, but not standard deviations.")
    tibble(
      category = c("A", "B"),
      kappa = 10,
      nu = 30,
      M = list(c("cue1" = -2, "cue2" = -2), c("cue1" = 2, "cue2" = 2)),
      S = list(matrix(c(1, .3, .3, 1), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2"))),
               matrix(c(1, -.3, -.3, 1), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2")))),
      lapse_rate = .05
    ) %>%
    mutate(category = factor(category))
  } else if (example == 2) {
    message("An example belief for two categories in a 2D cue continuum that differ in means and correlatation, but not standard deviations.
            Same as Example 1, but with independent perceptual noise along both cue dimensions.")
    tibble(
      category = c("A", "B"),
      kappa = 10,
      nu = 30,
      M = list(c("cue1" = -2, "cue2" = -2), c("cue1" = 2, "cue2" = 2)),
      S = list(matrix(c(1, .3, .3, 1), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2"))),
               matrix(c(1, -.3, -.3, 1), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2")))),
      lapse_rate = .05,
      Sigma_noise = list(matrix(c(1, 0, 0, .25), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2"))))
    ) %>%
      mutate(category = factor(category))
  }
}


#' Make NIW prior from data.
#'
#' Constructs an NIW_belief object, representing Normal-Inverse Wishart (NIW) parameters for all categories found in
#' the data. This object can be used as a prior for functions like \code{\link{update_NIW_beliefs}}.
#'
#' Currently, \code{make_NIW_prior_from_data()} does not infer kappa or nu, nor does it fit hierarchical data. Rather
#' the function simply estimates the category mean and covariance matrix from the sample (\code{data}), assumes them
#' to be the expected category mean (mu) and covariance (Sigma), and derives the M and S parameters
#' of the NIW from mu and Sigma based on the user-provided kappa and nu. That means M = mu and S = Sigma * (nu - D -1),
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
#' @param kappa The strength of the beliefs over the category mean (pseudocounts).
#' @param nu The strength of the beliefs over the category covariance matrix (pseudocounts).
#' @param lapse_rate Optionally specify a lapse rate. (default: \code{NA})
#' @param Sigma_noise Optionally specify a (multivariate Gaussian) noise covariance matrix. This argument will be
#' ignored if `NULL`. (default: NULL)
#' @param keep.category_parameters Should categories' mu and Sigma be included in the output (in addition to M
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
make_NIW_prior_from_data = function(
  data,
  groups = NULL,
  category = "category",
  cues,
  kappa = NA,
  nu = NA,
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
              msg = "nu must be larger than D (dimensionality of cues) + 1.")

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

  data %<>%
    mutate(
      !! category := factor(!! category),
      kappa = kappa,
      nu = nu,
      M = mu,
      S = map2(Sigma, nu, get_S_from_Sigma),
      lapse_rate = lapse_rate) %>%
    { if (!is.null(Sigma_noise)) mutate(., Sigma_noise = list(Sigma_noise)) else . } %>%
    ungroup()

  if (!keep.category_parameters) data %<>% select(-c(mu, Sigma))
  if (!is.NIW_belief(data, category = as_name(category))) warning("Something went wrong. The returned object is not an NIW belief.")

  return(data)
}

#' Make multivariate Gaussian exposure data.
#'
#' Returns a tibble of observations drawn from multivariate Gaussians, with one observation per row. Each row
#' provides the category label and cue values. If \code{keep.parameters = T} then the parameters (\code{N, mean, sigma})
#' are also returned.
#'
#' The input is expected to be lists/vectors of parameters with the n-th element of each list/vector specifying the
#' category label, number of observations, \code{mu}, and \code{Sigma} of the n-th Gaussian.
#'
#' @param Ns Integer vector, with each number specifying the number of observations to be drawn from the corresponding
#' Gaussian.
#' @param mus List of mean vectors, each specifying the mean of a multivariate Gaussian.
#' @param sigmas List of covariance matrices, each specifying the covariance of a multivariate Gaussian.
#' @param category.labels Character vector of category names, each specifying the category label of a multivariate Gaussian. If \code{NULL}
#' (default) then Gaussians will be numbered from 1:N.
#' @param cue.labels Character vector of cue names. If \code{NULL} (default) then the cues will be numbered cue1, cue2, ...
#' @param randomize.order Should the order of the data be randomized? (default: FALSE) This won't affect the final outcome of
#' NIW belief updating, but it will change the incremental updates (and thus, for example, visualizations of the update process).
#' @param keep.input_parameters Should the parameters handed to this function be included in the output? (default: FALSE)
#'
#' @return A tibble.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
#'
make_MV_exposure_data = function(
  Ns, mus, Sigmas,
  category.labels = NULL,
  cue.labels = NULL,
  randomize.order = F,
  keep.input_parameters = F
) {
  assert_that(!is.null(mus), !is.null(Sigmas))
  assert_that(is.null(category.labels) | length(mus) == length(category.labels),
              msg = "Number of category labels mismatch number of mus.")
  assert_that(is.null(cue.labels) | length(mus[[1]]) == length(cue.labels),
              msg = "Number of cue labels mismatches dimensionality of mus.")

  if (is.null(category.labels)) category.labels = 1:length(mus)
  if (is.null(cue.labels)) cue.labels = paste0("cue", 1:length(mus[[1]]))

  x = tibble(category = category.labels, n = Ns, mu = mus, Sigma = Sigmas) %>%
    mutate(data = pmap(.l = list(n, mu, Sigma), .f = rmvnorm)) %>%
    mutate(data = map(data, ~ .x %>% as.data.frame() %>% rename_all(~ cue.labels))) %>%
    unnest(data)

  if (randomize.order)
    x = sample_frac(x, 1)

  if (keep.input_parameters) return(x) else return(x %>% select(category, everything(), -c(n, mu, Sigma)))
}


#' Update parameters of NIW prior beliefs about multivariate Gaussian category.
#'
#' Returns updated/posterior M, S, kappa, or nu based on \insertCite{@see @murphy2012 p. 134;textual}{MVBeliefUpdatr}.
#'
#' @rdname update_NIW_parameters
#' @export
update_NIW_belief_kappa = function(kappa_0, x_N) { kappa_0 + x_N }
update_NIW_belief_nu = function(nu_0, x_N) { nu_0 + x_N }
update_NIW_belief_M = function(kappa_0, M_0, x_N, x_mean) { (kappa_0 / (kappa_0 + x_N)) * M_0 + x_N / (kappa_0 + x_N) * x_mean }
update_NIW_belief_S = function(kappa_0, M_0, S_0, x_N, x_mean, x_S) { S_0 + x_S + (kappa_0 * x_N) / (kappa_0 + x_N) * (x_mean - M_0) %*% t(x_mean - M_0) }



#' Update NIW prior beliefs about multivariate Gaussian category based on sufficient statistics of observations.
#'
#' Returns updated/posterior beliefs about the Gaussian categories based on conjugate NIW prior.
#'
#' The priors for the categories are specified through the \code{priors} argument, an
#' \code{\link[=is.NIW_belief]{NIW_belief}} object. Which category is updated is determined by x_category.
#' Updating proceeds as in \insertCite{@see @murphy2012 p. 134;textual}{MVBeliefUpdatr}. The prior kappa
#' and nu will be incremented by the number of observations (x_N). The prior M and S will be updated based on the
#' prior kappa, prior nu, x_N and, of course, the sample mean (x_mean) and sum of squares (x_S) of the observations.
#'
#' A number of different updating schemes are supported, including supervised updating based on labeled data and
#' unsupervised updating based on unlabeled data.
#' \itemize{
#'   \item "no-updating" doesn't update the prior. Combined with keep_history = T, this allows the creation of baseline
#'   against which to compare the updated beliefs. This option is likely most useful when used as part of a call to
#'   \code{\link{update_NIW_beliefs}}.
#'   \item "label-certain" assumes that the label is provided and known to the observer with 100% certainty, resulting
#'   in fully Bayesian supervised belief-updating.
#'   \item "nolabel-criterion" implements a winner-takes-all update based on the prior beliefs. The input is attributed
#'   to the category with the highest posterior probability (calculated based on the prior beliefs), and this category
#'   is updated using the "label-certain" method.
#' }
#' This functionality could be extended with additional proposals. Please feel free to suggest additional features.
#'
#' @param prior An \code{\link[=is.NIW_belief]{NIW_belief}} object, specifying the prior beliefs.
#' @param x_category Character. The label of the category that is to be updated.
#' @param x A single observation.
#' @param x_mean Mean of the observations.
#' @param x_S Centered sum of squares matrix of the observations.
#' @param x_N Number of observations that went into the mean and sum of squares matrix.
#' @param add_noise Determines whether multivariate Gaussian noise is added to the input.
#' If `NULL`, no noise is added during the updating. If "sample" then a sample of
#' noise is added to the input. If "marginalize" then each observation is transformed into the marginal distribution
#' that results from convolving the input with noise. This latter option might be helpful, for example, if one is
#' interested in estimating the consequences of noise across individuals. If add_noise is not `NULL` a Sigma_noise
#' column must be present in the NIW_belief object specified as the priors argument. (default: `NULL`)
#' @param method Which updating method should be used? See details. (default: "supervised-certain")
#'
#' @return A tibble.
#'
#' @seealso TBD
#' @keywords TBD
#' @references \insertRef{murphy2012}{MVBeliefUpdatr}
#' @examples
#' TBD
#' @rdname update_NIW_belief
#' @export
update_NIW_belief_by_sufficient_statistics = function(
  prior, x_category, x_mean, x_S, x_N,
  add_noise = NULL,
  method = "label-certain"
) {
  # TO DO: check match between dimensionality of belief and of input, check that input category is part of belief, etc.
  assert_NIW_belief(prior)
  assert_that(all(is_scalar_double(x_N), x_N >= 0), msg = "x_N must be >= 0.")
  assert_that(method %in% c("no-updating", "label-certain", "nolabel-criterion"),
              msg = paste0(method, "is not an acceptable updating method. See details section of help page."))
  if (method %nin% c("no-updating", "label-certain"))
    assert_that(x_N > 1,
                msg = "For this updating method, only incremental updating (one observations at a time) is implemented.")
  assert_that(any(is.null(add_noise), add_noise %in% c("sample", "marginalize")),
              msg = 'add_noise must be one of "sample" or "marginalize"')
  if (!is.null(add_noise))
    assert_that(any(add_noise %in% "marginalize",
                    all(add_noise %in% c("sample"), is_scalar_integerish(x_N), is_weakly_greater_than(x_N, 1))),
              msg = "For noise sampling, x_N must be a positive integer")
  assert_that(any(is.null(add_noise), "Sigma_noise" %in% names(prior)),
              msg = "Can't add noise: argument priors does not have column Sigma_noise.")

  x_Ns = as.list(rep(0, length(prior$category)))
  # Determine how observation should be distributed across categories
  if (method == "no-updating") return(prior) else
    if (method == "label-certain") x_Ns[[which(prior$category == x_category)]] = x_N else
      if (method == "nolabel-uniform") x_Ns = as.list(1 / length(prior$category)) * x_N else {
        decision_rule = case_when(
          method == "nolabel-criterion" ~ "criterion",
          method == "nolabel-posterior" ~ "proportional",
          T ~ NA_character_
        )
        message("get_categorization_from_NIW_belief still needs to be written. simplify = F is meant to return a vector of of
                posterior probabilities. the function should also have an option 'sampling' which allows to sample based on luce's
                choice rule. if simplify = T only the label of the category that is chosen will be displayed. incompatible with
                decision_rule = 'proportional'")
        x_Ns = as.list(get_categorization_from_NIW_belief(x = x_mean, belief = prior, decision_rule = decision_rule,
                                                  simplify = F) * x_N)
      }

  # Handle noise
  if (is.null(add_noise)) add_noise = ""
  if (add_noise == "sample") {
    x = rmvnorm(n = x_N,
                sigma = prior$Sigma_noise[[1]])
    x_mean = x_mean + colMeans(x)
    if (x_N > 1) x_S = x_S + cov(x)
  } else if (add_noise == "marginalize") x_S = x_S + prior$Sigma_noise[[1]]

  x_mean = replicate(length(prior$category), x_mean)
  x_S = replicate(length(prior$category), x_S)

  prior %<>%
    mutate(
      M = pmap(.l = list(kappa, M, x_Ns, x_mean), update_NIW_belief_M),
      kappa = map2(kappa, x_Ns, update_NIW_belief_kappa),
      nu = map2(nu, x_Ns, update_NIW_belief_nu),
      S = pmap(.l = list(kappa, M, S, x_Ns, x_mean, x_S), update_NIW_belief_M))

  # if (method %in% c("label-certain", "nolabel-criterion")) {
  #   M_0 = prior[prior$category == x_category,]$M[[1]]
  #   kappa_0 = prior[prior$category == x_category,]$kappa
  #   nu_0 = prior[prior$category == x_category,]$nu
  #   S_0 = prior[prior$category == x_category,]$S[[1]]
  #
  #   prior %<>%
  #     mutate(
  #       M = list((kappa_0 / (kappa_0 + x_N)) * M_0 + x_N / (kappa_0 + x_N) * x_mean),
  #       kappa = kappa_0 + x_N,
  #       nu = nu_0 + x_N,
  #       S = list(S_0 + x_S + (kappa_0 * x_N) / (kappa_0 + x_N) * (x_mean - M_0) %*% t(x_mean - M_0)))
  # } else if (method == "nolabel-criterion") {
  #
  #   prior %<>%
  #     mutate(
  #       M = list((kappa_0 / (kappa_0 + x_N)) * M_0 + x_N / (kappa_0 + x_N) * x_mean),
  #       kappa = kappa_0 + x_N,
  #       nu = nu_0 + x_N,
  #       S = list(S_0 + x_S + (kappa_0 * x_N) / (kappa_0 + x_N) * (x_mean - M_0) %*% t(x_mean - M_0)))
  # }

  return(prior)
}


#' @rdname update_NIW_belief
#' @export
update_NIW_belief_by_one_observation = function(
  prior, x_category, x,
  add_noise = NULL,
  method = "label-certain"
) {
  update_NIW_belief_by_sufficient_statistics(prior, x_category = x_category, x_mean = x, x_S = 0L, x_N = 1L,
                                             add_noise = add_noise, method = method)
}


#' Update NIW prior beliefs about multivariate Gaussian category based on exposure data.
#'
#' Returns updated/posterior beliefs about the Gaussian categories based on conjugate NIW prior.
#'
#' The priors for the categories are specified through the \code{priors} argument. This is expected to be a tibble
#' of the same format as the posterior draws stored in an MV IBBU stanfit object. Each row of the tibble specifies
#' the prior for one category (specified in the \code{category} column). The four parameters of the NIW are the
#' pseudocounts indicating the strength of the prior beliefs into the mean (\code{kappa}) and covariance (\code{
#' nu}), as well as the prior mean of means (\code{M}, same as \code{M_0} in Murphy, 2012) and the prior scatter
#' matrix (\code{S}, same as \code{S_0} in Murphy, 2012).
#'
#' @param prior A \code{\link[=is.NIW_belief]{NIW_belief}} object, specifying the prior beliefs.
#' @param exposure \code{data.frame} or \code{tibble} with exposure data. Each row is assumed to contain one observation.
#' @param category Name of variable in \code{data} that contains the category information. (default: "category")
#' @param cues Name(s) of variables in \code{data} that contain the cue information. By default these cue names are
#' extracted from the prior object.
#' @param exposure.order Name of variable in \code{data} that contains the order of the exposure data. If `NULL` the
#' exposure data is assumed to be in the order in which it should be presented.
#' @param add_noise Determines whether multivariate Gaussian noise is added to the input. See \code{\link{update_NIW_belief}}.
#' (default: `NULL`)
#' @param method Which updating method should be used? See \code{\link{update_NIW_belief}}. (default: "supervised-certain")
#' @param keep.update_history Should the history of the belief-updating be stored and returned? If so, the output is
#' tibble with the one set of NIW beliefs for each exposure observation. This is useful, for example, if one wants to
#' visualize the changes in the category parameters, posterior predictive, categorization function, or alike across time.
#' (default: `TRUE`)
#' @param keep.exposure_data Should the input data be included in the output? If `FALSE` then only the category and cue
#' columns will be kept. If `TRUE` then all columns will be kept. (default: `FALSE`)
#'
#' @return A tibble.
#'
#' @seealso TBD
#' @keywords TBD
#' @references \insertRef{murphy2012}{MVBeliefUpdatr}
#' @examples
#' TBD
#' @export
#'
update_NIW_beliefs <- function(
  prior,
  exposure,
  category = "category",
  cues = names(prior$M[[1]]),
  exposure.order = NULL,
  add_noise = NULL,
  method = "label-certain",
  keep.update_history = TRUE,
  keep.exposure_data = FALSE
){
  assert_NIW_belief(prior)
  assert_that(any(is_tibble(exposure), is.data.frame(exposure)))
  assert_that(all(is.flag(keep.update_history), is.flag(keep.exposure_data)))
  assert_that(any(is.null(exposure.order), exposure.order %in% names(exposure)),
              msg = paste0("exposure.order variable not found: ", exposure.order, " must be a column in the exposure data."))
  assert_that(any(is.null(exposure.order), if (!is.null(exposure.order)) is.numeric(exposure[[exposure.order]]) else T),
              msg = "exposure.order variable must be numeric.")

  # Number of dimensions/cues
  D = length(cues)
  if (any(prior$nu < D + 2))
    message(paste0("Prior for at least one category had nu smaller than allowed (", D, ").\n"))

  # Prepare exposure data
  exposure %<>%
    { if (!is.null(exposure.order)) arrange(., !! sym(exposure.order)) else . } %>%
    mutate(cues = pmap(.l = list(!!! syms(cues)), .f = ~ c(...)))

  if (keep.update_history)
    prior %<>%
    mutate(observation.n = 0)

  for (i in 1:nrow(exposure)) {
    posterior = if (keep.update_history) prior %>% filter(observation.n == i - 1) else prior
    posterior =
      update_NIW_belief_by_one_observation(
        prior = posterior,
        x = unlist(exposure[i, "cues"]),
        x_category = exposure[i,]$category,
        add_noise = add_noise,
        method = method)

    if (keep.update_history) {
      posterior %<>%
        mutate(observation.n = i)
      prior = rbind(prior, posterior)
    } else prior = posterior
  }

  if (keep.exposure_data) {
    exposure %<>%
      { if (!is.null(exposure.order))
        rename(., observation.n = !! sym(exposure.order)) else
          mutate(., observation.n = 1:nrow(exposure)) } %>%
      rename_with(~ paste0("observation.", .x), !starts_with("observation.n"))

    prior %<>%
      left_join(exposure)
  }

  return(prior)
}
