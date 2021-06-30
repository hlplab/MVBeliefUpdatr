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
      m = list(c("cue1" = -2, "cue2" = -2), c("cue1" = 2, "cue2" = 2)),
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
      m = list(c("cue1" = -2, "cue2" = -2), c("cue1" = 2, "cue2" = 2)),
      S = list(matrix(c(1, .3, .3, 1), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2"))),
               matrix(c(1, -.3, -.3, 1), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2")))),
      lapse_rate = .05,
      Sigma_noise = list(matrix(c(1, 0, 0, .25), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2"))))
    ) %>%
      mutate(category = factor(category))
  }
}



#' Update parameters of NIW prior beliefs about multivariate Gaussian category.
#'
#' Returns updated/posterior m, S, kappa, or nu based on \insertCite{@see @murphy2012 p. 134;textual}{MVBeliefUpdatr}.
#' These functions are called, for example, by \code{\link{update_NIW_belief_by_sufficient_statistics_of_one_category}},
#' \code{\link{update_NIW_belief_by_one_observation}}, and \code{\link{update_NIW_beliefs_incrementally}}.
#'
#' @return Depending on the updated parameter, a numeric scalar (kappa, nu), vector (m) or square matrix (S).
#'
#' @rdname update_NIW_parameters
#' @export
update_NIW_belief_kappa = function(kappa_0, x_N) { kappa_0 + x_N }
#' @rdname update_NIW_parameters
#' @export
update_NIW_belief_nu = function(nu_0, x_N) { nu_0 + x_N }
#' @rdname update_NIW_parameters
#' @export
update_NIW_belief_m = function(kappa_0, m_0, x_N, x_mean) { (kappa_0 / (kappa_0 + x_N)) * m_0 + x_N / (kappa_0 + x_N) * x_mean }
#' @rdname update_NIW_parameters
#' @export
update_NIW_belief_S = function(kappa_0, m_0, S_0, x_N, x_mean, x_S) { S_0 + x_S + (kappa_0 * x_N) / (kappa_0 + x_N) * (x_mean - m_0) %*% t(x_mean - m_0) }



#' Update NIW prior beliefs about multivariate Gaussian category based on sufficient statistics of observations.
#'
#' Returns updated/posterior beliefs about the Gaussian categories based on conjugate NIW prior.
#'
#' The priors for the categories are specified through the \code{priors} argument, an
#' \code{\link[=is.NIW_belief]{NIW_belief}} object. Which category is updated is determined by x_category.
#' Updating proceeds as in \insertCite{@see @murphy2012 p. 134;textual}{MVBeliefUpdatr}. The prior kappa
#' and nu will be incremented by the number of observations (x_N). The prior m and S will be updated based on the
#' prior kappa, prior nu, x_N and, of course, the sample mean (x_mean) and sum of squares (x_S) of the observations.
#'
#' A number of different updating schemes are supported, including supervised updating based on labeled data and
#' unsupervised updating based on unlabeled data. Except for "nolabel-sampling", the unsupervised updating rules
#' were originally presented in \insertCite{yan:jaeger2018;textual}{MVBeliefUpdatr}.
#' \itemize{
#'   \item "no-updating" doesn't update the prior. Combined with keep_history = T, this allows the creation of baseline
#'   against which to compare the updated beliefs. This option is likely most useful when used as part of a call to
#'   \code{\link{update_NIW_beliefs_incrementally}}.
#'   \item "label-certain" assumes that the label is provided and known to the observer with 100% certainty, resulting
#'   in fully Bayesian supervised belief-updating. This update method is \emph{not order sensitive}. Under this update
#'   method, the order of a batch of observations does not affect the final posterior belief after all observations
#'   have been made (it does, of course, effect the beliefs held in the interim). This reflects the "exchangeability"
#'   assumption of Bayesian belief-updating under the assumption that the observations are drawn from a stationary
#'   random process.
#'   \item "nolabel-criterion" implements a winner-takes-all update based on the prior beliefs. The
#'   posterior probability of all categories under the prior is calculated, and only the category with the
#'   highest posterior probability is updated (see decision_rule = "criterion" in \code{\link{get_categorization_from_NIW_ideal_adaptor}}).
#'   The update for this category proceed the same ways as under the "label-certain" method. method "nolabel-criterion" is
#'   thus \emph{order sensitive}, but deterministic (i.e., it yields the same result on each repeated run), provided
#'   there is no noise (specifically, no noise = "sample").
#'   \item "nolabel-sampling" implements another winner-takes-all update based on the prior beliefs. The
#'   posterior probability of all categories under the prior is calculated. Then a single category is \emph{sampled}
#'   based on this posterior (see decision_rule = "sampling" in \code{\link{get_categorization_from_NIW_ideal_adaptor}}), and only this
#'   sampled category is  updated with the observation.
#'   Like "nolabel-criterion", this updating method is \emph{order sensitive}. Unlike "nolabel-criterion", this method had a
#'   random element and thus does \emph{not replicate the same result on each run}.
#'   \item "nolabel-proportional" captures the uncertainty about the category label (and is this sense fully Bayesian). The
#'   posterior probability of all categories under the prior is calculated (see decision_rule = "proportional" in
#'   \code{\link{get_categorization_from_NIW_ideal_adaptor}}). The input is then distributed across all categories
#'   based on their posterior probability under the prior beliefs. That is, \emph{all categories are updated}, but to different
#'   degrees. Like "nolabel-criterion" and "nolabel-sampling", this updating method is \emph{order sensitive}. Like "nolabel-criterion",
#'   it is \emph{deterministic}.
#'   \item "nolabel-uniform" implements unsupervised updating under maximal \emph{uninformed} uncertainty. The input is attributed
#'   to equal parts to all categories. Like all the "nolabel-*" methods, this updating method is \emph{order sensitive}. Like
#'   "nolabel-criterion" and "nolabel-proportional, it is \emph{deterministic}. Unlike any of these methods, it completely ignores
#'   the knowledge listeners have (i.e., the prior).
#' }
#' \strong{Please consider the following important limitations of these functions:}
#' \itemize{
#'   \item the category priors---p(category)---are currently assumed to be uniform. This affects the categorization of the input
#' based on the NIW prior and thus the way that the input observations are distributed across the categories during
#' unsupervised updating.
#'   \item the updating models introduced here treat the cue distributions of the different categories as \emph{independent} of
#'   each other. That is, the updating does not consider knowledge about, for example, the covariance of the category means, although
#'   there is evidence that listeners have knowledge of this covariance and draw on it during categorization and adaptation.
#' }
#' Please feel free to suggest additional features.
#'
#' @param prior An \code{\link[=is.NIW_belief]{NIW_belief}} object, specifying the prior beliefs.
#' @param x_category Character. The label of the category that is to be updated.
#' @param x A single observation.
#' @param x_mean mean of the observations.
#' @param x_S Centered sum of squares matrix of the observations.
#' @param x_N Number of observations that went into the mean and sum of squares matrix.
#' @param add_noise Determines whether multivariate Gaussian noise is added to the input.
#' If `NULL`, no noise is added during the updating. If "sample" then a sample of
#' noise is added to the input. If "marginalize" then each observation is transformed into the marginal distribution
#' that results from convolving the input with noise. This latter option might be helpful, for example, if one is
#' interested in estimating the consequences of noise across individuals. If add_noise is not `NULL` a Sigma_noise
#' column must be present in the NIW_belief object specified as the priors argument. (default: `NULL`)
#' @param method Which updating method should be used? See details. (default: "supervised-certain")
#' @param verbose Should more informative output be provided?
#'
#' @return An NIW_belief object.
#'
#' @seealso \code{\link{update_NIW_belief_kappa}}, \code{\link{update_NIW_belief_nu}}, \code{\link{update_NIW_belief_m}},
#' \code{\link{update_NIW_belief_S}}, all of which are called by \code{update_NIW_belief_by_sufficient_statistics_of_one_category}.
#' @keywords belief-updating NIW Normal-Inverse Wishart
#' @references \insertRef{murphy2012}{MVBeliefUpdatr}
#' @examples
#' TBD
#' @rdname update_NIW_belief
#' @export
update_NIW_belief_by_sufficient_statistics_of_one_category = function(
  prior, x_category, x_mean, x_S, x_N,
  add_noise = NULL,
  method = "label-certain",
  verbose = FALSE
) {
  # TO DO: check match between dimensionality of belief and of input, check that input category is part of belief, etc.
  assert_NIW_belief(prior)
  assert_that(all(is.scalar(x_N), is.numeric(x_N)), msg = "x_N must be a scalar numeric.")
  assert_that(x_N >= 0, msg = paste("x_N is", x_N, "but must be >= 0."))
  # if there's nothing to update, return prior.
  if (any(is.na(x_mean), x_N == 0)) {
    if (verbose) message("Missing observation (x_N == 0 or x_mean is NA). Returning prior as posterior.")
    return(prior)
  }

  assert_that(method %in% c("no-updating",
                            "label-certain",
                            "nolabel-criterion", "nolabel-sampling", "nolabel-proportional",
                            "nolabel-uniform"),
              msg = paste(method, "is not an acceptable updating method. See details section of help page."))
  if (method %nin% c("no-updating", "label-certain"))
    assert_that(x_N <= 1,
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
    if (method == "label-certain") x_Ns[[which(prior$category == x_category)]] <- x_N else
      if (method == "nolabel-uniform") x_Ns <- as.list(1 / length(prior$category) * x_N) else
        if (method %in% c("nolabel-criterion", "nolabel-sampling", "nolabel-proportional")) {
        x_Ns = as.list(get_categorization_from_NIW_ideal_adaptor(
          x = x_mean, belief = prior,
          decision_rule = gsub("nolabel-", "", method),
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

  x_mean = list(x_mean)
  x_S = list(x_S)
  prior %<>%
    mutate(
      m = pmap(.l = list(kappa, m, x_Ns, x_mean), update_NIW_belief_m),
      kappa = unlist(map2(kappa, x_Ns, update_NIW_belief_kappa)),
      nu = unlist(map2(nu, x_Ns, update_NIW_belief_nu)),
      S = pmap(.l = list(kappa, m, S, x_Ns, x_mean, x_S), update_NIW_belief_S))

  return(prior)
}


#' @rdname update_NIW_belief
#' @export
update_NIW_belief_by_one_observation = function(
  prior, x_category, x,
  add_noise = NULL,
  method = "label-certain",
  verbose = FALSE
) {
  update_NIW_belief_by_sufficient_statistics_of_one_category(prior, x_category = x_category, x_mean = x, x_S = 0L, x_N = 1L,
                                             add_noise = add_noise, method = method, verbose = verbose)
}


#' Update NIW prior beliefs about multivariate Gaussian category based on exposure data.
#'
#' Returns updated/posterior beliefs about the Gaussian categories based on conjugate NIW prior.
#'
#' The priors for the categories are specified through the \code{priors} argument. This is expected to be a tibble
#' of the same format as the posterior draws stored in an MV IBBU stanfit object. Each row of the tibble specifies
#' the prior for one category (specified in the \code{category} column). The four parameters of the NIW are the
#' pseudocounts indicating the strength of the prior beliefs into the mean (\code{kappa}) and covariance (\code{
#' nu}), as well as the prior mean of means (\code{m}, same as \code{m_0} in Murphy, 2012) and the prior scatter
#' matrix (\code{S}, same as \code{S_0} in Murphy, 2012).
#'
#' @param prior An \code{\link[=is.NIW_belief]{NIW_belief}} object, specifying the prior beliefs.
#' @param exposure \code{data.frame} or \code{tibble} with exposure data. Each row is assumed to contain one observation.
#' @param exposure.category Name of variable in \code{data} that contains the category information. (default: "category")
#' @param exposure.cues Name(s) of variables in \code{data} that contain the cue information. By default these cue names are
#' extracted from the prior object.
#' @param exposure.order Name of variable in \code{data} that contains the order of the exposure data. If `NULL` the
#' exposure data is assumed to be in the order in which it should be presented.
#' @param add_noise Determines whether multivariate Gaussian noise is added to the input. See \code{
#' \link{update_NIW_belief_by_sufficient_statistics_of_one_category}}. (default: `NULL`)
#' @param add_lapse Determines the proportion of trials on which the listener lapses, not updating beliefs. Must be a number
#' between 0 and 1, or NULL if lapses are to be ignored. (default: `NULL`)
#' @param method Which updating method should be used? See \code{\link{update_NIW_belief_by_sufficient_statistics_of_one_category}}.
#' The length of this argument should either be 1 (in which case it is recycled for each observation) or the same as
#' the number of rows in \code{expsure}. (default: "label-certain").
#' @param keep.update_history Should the history of the belief-updating be stored and returned? If so, the output is
#' tibble with the one set of NIW beliefs for each exposure observation. This is useful, for example, if one wants to
#' visualize the changes in the category parameters, posterior predictive, categorization function, or alike across time.
#' (default: `TRUE`)
#' @param keep.exposure_data Should the input data be included in the output? If `FALSE` then only the category and cue
#' columns will be kept. If `TRUE` then all columns will be kept. (default: `FALSE`)
#' @param verbose Should more informative output be provided?
#'
#' @return An NIW_belief object.
#'
#' @seealso \code{\link{update_NIW_belief_by_one_observation}}, which is called by \code{update_NIW_beliefs_incrementally}
#' @keywords belief-updating NIW Normal-Inverse Wishart
#' @references \insertRef{murphy2012}{MVBeliefUpdatr}
#' @examples
#' TBD
#' @export
update_NIW_beliefs_incrementally <- function(
  prior,
  exposure,
  exposure.category = "category",
  exposure.cues = get_cue_labels_from_model(prior),
  exposure.order = NULL,
  add_noise = NULL,
  add_lapse = NULL,
  method = "label-certain",
  keep.update_history = TRUE,
  keep.exposure_data = FALSE,
  verbose = FALSE
){
  if (verbose) message("Assuming that category variable in NIW belief/ideal adaptor is called category.")

  assert_NIW_belief(prior)
  assert_that(any(is.null(add_noise), is_scalar_character(add_noise)),
              msg = "add_noise must be NULL or a scalar character.")
  assert_that(any(is.null(add_lapse), is_scalar_double(add_lapse)),
              msg = "add_lapse must be NULL or a scalar double (between 0 or 1).")
  assert_that(between(add_lapse, 0, 1),
              msg = "If not NULL, add_noise must be between 0 to 1.")
  assert_that(all(is.flag(keep.update_history), is.flag(keep.exposure_data)))
  assert_that(any(is_tibble(exposure), is.data.frame(exposure)))
  assert_that(exposure.category %in% names(exposure),
              msg = paste0("exposure.category variable not found: ", exposure.category, " must be a column in the exposure data."))
  assert_that(any(is.null(exposure.order), exposure.order %in% names(exposure)),
              msg = paste0("exposure.order variable not found: ", exposure.order, " must be a column in the exposure data."))
  assert_that(any(is.null(exposure.order), if (!is.null(exposure.order)) is.numeric(exposure[[exposure.order]]) else T),
              msg = "exposure.order variable must be numeric.")
  assert_that(length(method) == 1 | length(method) == nrow(exposure),
              msg = paste0("Length of method argument must be either 1 or the number of rows in exposure (", nrow(exposure),")."))
  if (length(method) == 1) method = rep(method, nrow(exposure))

  # Number of dimensions/cues
  D = length(exposure.cues)
  if (any(prior$nu <= D + 1))
    message(paste0("Prior for at least one category had nu smaller than allowed (is ", min(prior$nu), "; should be >", D+1, ").\n"))

  # Prepare exposure data
  exposure %<>%
    { if (!is.null(exposure.order)) arrange(., !! sym(exposure.order)) else . } %>%
    make_vector_column(exposure.cues, "cues") %>%
    # Apply lapses
    { if (!is.null(add_lapse)) mutate(., cues = ifelse(rbinom(nrow(.), 1, add_lapse), NA, cues)) else . }

  if (keep.update_history)
    prior %<>%
    mutate(observation.n = 0)

  for (i in 1:nrow(exposure)) {
    posterior = if (keep.update_history) prior %>% filter(observation.n == i - 1) else prior
    posterior =
      update_NIW_belief_by_one_observation(
        prior = posterior,
        x = unlist(exposure[i, "cues"]),
        x_category = exposure[i,][[exposure.category]],
        add_noise = add_noise,
        method = method[i],
        verbose = verbose)

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
