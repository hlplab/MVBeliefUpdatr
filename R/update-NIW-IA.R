#' @importFrom dplyr sample_frac
#' @importFrom magrittr is_weakly_greater_than
#' @importFrom rlang is_scalar_integerish
#' @importFrom mvtnorm rmvnorm
NULL

#' Update parameters of NIW prior beliefs about multivariate Gaussian category.
#'
#' Returns updated/posterior m, S, kappa, or nu based on \insertCite{@see @murphy2012 p. 134;textual}{MVBeliefUpdatr}.
#' These functions are called, for example, by \code{\link{update_NIW_belief_by_sufficient_statistics_of_one_category}},
#' \code{\link{update_NIW_belief_by_one_observation}}, and \code{\link{update_NIW_beliefs_incrementally}}.
#'
#' @param kappa_0,nu_0,m_0,S_0 Prior parameter values of NIW.
#' @param x_n,x_mean,x_SS Number of observations, mean of those observations, and *centered* sum-of-square matrix of those
#' observations based on which the parameters of the NIW should be updated.
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
update_NIW_belief_m = function(kappa_0, m_0, x_N, x_mean) { (kappa_0 / (kappa_0 + x_N)) * m_0 + (x_N / (kappa_0 + x_N)) * x_mean }
#' @rdname update_NIW_parameters
#' @export
update_NIW_belief_S = function(kappa_0, m_0, S_0, x_N, x_mean, x_SS) { S_0 + x_SS + ((kappa_0 * x_N) / (kappa_0 + x_N)) * (x_mean - m_0) %*% t(x_mean - m_0) }


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
#' @param x The cues of single observation.
#' @param x_category The category label(s) of one or more observations.
#' @param x_mean The cue mean of observations.
#' @param x_SS *Centered* The sum of squares matrix of observations.
#' @param x_N Number of observations that went into the mean and sum of squares matrix.
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
#' @param lapse_treatment Determines whether attentional lapses can occur during which no updating occurs.
#' Can be "no_lapses", "sample", or "marginalize". If "no_lapses", no lapses occur (even if the ideal adaptor specifies
#' a non-zero `lapse_rate`), and all observations lead to updating. If "sample" or "marginalize", the lapse rate from '
#' the ideal adaptor will be used. For "sample", sampling determines for each observation whether it was a lapse or not.
#' If an observation was a lapse no updating occurs. For "marginalize", 1 - lapse_rate is the proportion of observations
#' that are assumed to be lapsing trials (default: "no_lapses" for NIW_beliefs; "sample" for NIW_ideal_adaptors)
#' @param method Which updating method should be used? See details. (default: "supervised-certain")
#' @param verbose Should more informative output be provided?
#'
#' @return An NIW_belief object.
#'
#' @seealso \code{\link{update_NIW_belief_kappa}}, \code{\link{update_NIW_belief_nu}}, \code{\link{update_NIW_belief_m}},
#' \code{\link{update_NIW_belief_S}}, all of which are called by \code{update_NIW_belief_by_sufficient_statistics_of_one_category}.
#' @keywords belief-updating NIW Normal-Inverse Wishart
#' @references \insertRef{murphy2012}{MVBeliefUpdatr}
#'
#' @importFrom purrr map2_dbl
#' @rdname update_NIW_belief
#' @export
update_NIW_belief_by_sufficient_statistics_of_one_category <- function(
  prior, x_category, x_mean, x_SS, x_N,
  noise_treatment = if (is.NIW_ideal_adaptor(prior)) { if (!is.null(first(prior$Sigma_noise))) "marginalize" else "no_noise" } else "no_noise",
  lapse_treatment = if (is.NIW_ideal_adaptor(prior)) "sample" else "no_lapses",
  method = "label-certain",
  verbose = FALSE
) {
  # Binding variables that RMD Check gets confused about otherwise
  # (since they are in non-standard evaluations)
  response <- NULL

  # TO DO: check match between dimensionality of belief and of input, check that input category is part of belief, etc.
  assert_that(all(is_scalar_character(noise_treatment)), is_scalar_character(lapse_treatment))
  if (any(noise_treatment != "no_noise", lapse_treatment != "no_lapses"))
    assert_NIW_ideal_adaptor(prior, verbose = verbose) else assert_NIW_belief(prior, verbose = verbose)

  assert_that(all(is.scalar(x_N), is.numeric(x_N)), msg = "x_N must be a scalar numeric.")
  assert_that(x_N >= 0, msg = paste("x_N is", x_N, "but must be >= 0."))

  # Handle lapses
  assert_that(lapse_treatment %in% c("no_lapses", "sample", "marginalize"),
              msg = paste(lapse_treatment, "is not an acceptable lapse_treatment. See details section of help page."))
  if (lapse_treatment == "sample") {
    x_N <- rbinom(1, x_N, 1 - get_lapse_rate_from_model(prior))
  } else if (lapse_treatment == "marginalize") {
    warning("Using lapse_treatment == 'marginalize' can result in updating by *fractions* of observations, which might not be wellformed.", call. = FALSE)
    x_N <- (1 - get_lapse_rate_from_model(prior)) * x_N
  }

  # If there's nothing to update (x_N == 0), return prior.
  if (any(is.na(x_mean), x_N == 0)) {
    if (verbose) message("No observations to update on (x_N == 0 or x_mean is NA). This can happen, for example, if observations are missing or because model was lapsing during all observations. Returning prior as posterior.")
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

  x_Ns <- as.list(rep(0, length(prior$category)))
  # Determine how observations should be distributed across categories
  if (method == "no-updating") return(prior) else
    if (method == "label-certain") x_Ns[[which(prior$category == x_category)]] <- x_N else
      if (method == "nolabel-uniform") x_Ns <- as.list(1 / length(prior$category) * x_N) else
        if (method %in% c("nolabel-criterion", "nolabel-proportional")) {
          x_Ns <-
            get_categorization_from_NIW_ideal_adaptor(
              x = x_mean, model = prior,
              noise_treatment = noise_treatment,
              lapse_treatment = lapse_treatment,
              decision_rule = gsub("nolabel-", "", method),
              simplify = F) %>%
            pull(response) * x_N %>%
            as.list()
        } else if (method %in% c("nolabel-sampling")) {
          x_Ns <-
            get_categorization_from_NIW_ideal_adaptor(
              x = x_mean, model = prior,
              noise_treatment = noise_treatment,
              lapse_treatment = lapse_treatment,
              decision_rule = "proportional", # perhaps decision_rule = "sampling" could be used here with simplify = T, pre-empting the remaining rows?
              simplify = F) %>%
            pull(response) %>%
            rmultinom(1, x_N, .) %>%
            as.list()
        }

  # Handle noise
  assert_that(noise_treatment %in% c("no_noise", "sample", "marginalize"))
  Sigma_noise <- get_perceptual_noise_from_model(prior)
  if (noise_treatment == "sample") {
    assert_that(all(is_scalar_integerish(x_N), is_weakly_greater_than(x_N, 1)),
                msg = "If noise_treatment is 'sample', x_N must be a positive integer.")
    if (verbose) message("Sampling perceptual noise and adding it to each observation")
    warning("Updating while including noise_treatment = sample has not yet been thoroughly tested. If noise is included in perception but not in the prior beliefs, it should be discounted during the updating. This implementation has not been tested. You might want to construct the model while adding the perceptual noise to the category beliefs and use categorization that does not add the noise again (noise_treatment = 'no_noise').")
    x <- rmvnorm(n = x_N, sigma = Sigma_noise)
    x_mean <- x_mean + colMeans(x)
    # Add sampled stimulus-level noise and subtract expected noise
    if (x_N > 1) x_SS <- x_SS + (ss(x, center = T) - cov2css(Sigma_noise, n = x_N))
    # Handle (hopefuly very rare) case where subtraction results in negative diagonal values
    # (this is now handled below)
    # if (any(diag(x_SS) < 0)) diag(x_SS) <- pmax(diag(x_SS), rep(0, length(diag(x_SS))))
  } else if (noise_treatment == "marginalize") {
    # As long the actual perceptual noise experienced during the input is on average the same as the noise
    # assumed by the input, then the two noise sources cancel each other out. I.e.,
    # x_SS <- x_SS +
    #     cov2css(Sigma_noise, n = x_N) - # average stimulus noise experience
    #     cov2css(Sigma_noise, n = x_N)   # explaining away expected noise
  }

  x_mean <- list(x_mean)
  x_SS <- list(x_SS)
  prior %<>%
    mutate(
      # Order of application matters here since all of these update functions assume inputs (kappa, nu, m, S) that are not yet updated.
      S = pmap(.l = list(.data$kappa, .data$m, .data$S, x_Ns, x_mean, x_SS), update_NIW_belief_S),
      m = pmap(.l = list(.data$kappa, .data$m, x_Ns, x_mean), update_NIW_belief_m),
      kappa = map2_dbl(.data$kappa, x_Ns, update_NIW_belief_kappa),
      nu = map2_dbl(.data$nu, x_Ns, update_NIW_belief_nu))

  if (noise_treatment == "sample") {
    # To correct of expected consequence of perceptual noise on the uncertainty about the mean (the additional
    # component in the updating of S, cf. update_NIW_belief_S), we need to correct S by subtracting
    #
    #   ((kappa * x_N) / (kappa + x_N)) * Sigma_noise / x_N =
    #   ((kappa * x_N) / (kappa + x_N)) * Sigma_noise / x_N =
    #   ((kappa) / (kappa + x_N)) * Sigma_noise
    prior %<>%
      mutate(
        # Since kappa has already been updated  ((kappa) / (kappa + x_N)) --> (kappa - x_N) / kappa
        S = pmap(.l = list(.data$kappa, x_Ns, .data$S), ~ ..3 - (..1 - ..2) / (..1) * .env$Sigma_noise),
        # Handle (hopefuly very rare) case where subtraction results in negative diagonal values
        S = .data$S - diag(diag(.data$S)) + diag(pmax(diag(.data$S), rep(0, length(diag(.data$S)))))
      )
  }

  return(prior)
}


#' @rdname update_NIW_belief
#' @export
update_NIW_belief_by_one_observation = function(
  prior, x_category, x,
  noise_treatment = if (is.NIW_ideal_adaptor(prior)) { if (!is.null(first(prior$Sigma_noise))) "marginalize" else "no_noise" } else "no_noise",
  lapse_treatment = if (is.NIW_ideal_adaptor(prior)) "sample" else "no_lapses",
  method = "label-certain",
  verbose = FALSE
) {
  update_NIW_belief_by_sufficient_statistics_of_one_category(
    prior, x_category = x_category, x_mean = x, x_SS = 0L, x_N = 1,
    noise_treatment = noise_treatment, lapse_treatment = lapse_treatment, method = method, verbose = verbose)
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
#' @param noise_treatment Determines whether and how multivariate Gaussian noise is considered during categorization.
#' See \code{\link{update_NIW_belief_by_sufficient_statistics_of_one_category}}.
#' @param lapse_treatment Determines whether attentional lapses can occur during which no updating occurs.
#' See \code{\link{update_NIW_belief_by_sufficient_statistics_of_one_category}}.
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
#' @export
update_NIW_ideal_adaptor_incrementally <- function(
  prior,
  exposure,
  exposure.category = "category",
  exposure.cues = get_cue_labels_from_model(prior),
  exposure.order = NULL,
  noise_treatment = if (is.NIW_ideal_adaptor(prior)) { if (!is.null(first(prior$Sigma_noise))) "marginalize" else "no_noise" } else "no_noise",
  lapse_treatment = if (is.NIW_ideal_adaptor(prior)) "sample" else "no_lapses",
  method = "label-certain",
  keep.update_history = TRUE,
  keep.exposure_data = FALSE,
  verbose = FALSE
){
  if (verbose) message("Assuming that category variable in NIW belief/ideal adaptor is called category.")
  if (lapse_treatment == "marginalize")
    warning("Using lapse_treatment == 'marginalize' can result in updating by *fractions* of observations, which might not be wellformed.", call. = FALSE)

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

  # Prepare exposure data
  exposure %<>%
    { if (!is.null(exposure.order)) arrange(., !! sym(exposure.order)) else . } %>%
    make_vector_column(exposure.cues, "cues")

  if (keep.update_history)
    prior %<>%
    mutate(observation.n = 0)

  for (i in 1:nrow(exposure)) {
    posterior <- if (keep.update_history) prior %>% filter(.data$observation.n == i - 1) else prior
    posterior <-
      suppressWarnings(
        update_NIW_belief_by_one_observation(
          prior = posterior,
          x = unlist(exposure[i, "cues"]),
          x_category = exposure[i,][[exposure.category]],
          noise_treatment = noise_treatment,
          lapse_treatment = lapse_treatment,
          method = method[i],
          verbose = verbose))

    if (keep.update_history) {
      posterior %<>%
        mutate(observation.n = i)
      prior <- rbind(prior, posterior)
    } else prior <- posterior
  }

  if (keep.exposure_data) {
    exposure %<>%
      { if (!is.null(exposure.order))
        rename(., observation.n = !! sym(exposure.order)) else
          mutate(., observation.n = 1:nrow(exposure)) } %>%
      rename_with(~ paste0("observation.", .x), !starts_with("observation.n"))

    prior %<>%
      left_join(exposure, by = "observation.n")
  }

  return(prior)
}


#' @rdname update_NIW_ideal_adaptor_incrementally
#' @export
update_NIW_ideal_adaptor_batch <- function(
  prior,
  exposure,
  exposure.category = "category",
  exposure.cues = get_cue_labels_from_model(prior),
  # Could add lapse treatment here, though it would only make sense to interpret it in the limit
  # as the fraction of trials that will be missed, changing only N, without affecting the other
  # sufficient statistics.
  noise_treatment = if (is.NIW_ideal_adaptor(prior)) { if (!is.null(first(prior$Sigma_noise))) "marginalize" else "no_noise" } else "no_noise",
  verbose = FALSE
){
  if (verbose) message("Assuming that category variable in NIW belief/ideal adaptor is called category.")

  assert_that(any(is_tibble(exposure), is.data.frame(exposure)))
  assert_that(exposure.category %in% names(exposure),
              msg = paste0("exposure.category variable not found: ", exposure.category, " must be a column in the exposure data."))
  assert_that(noise_treatment %in% c("no_noise", "marginalize"))

  # Prepare exposure data
  exposure %<>%
    make_vector_column(exposure.cues, "cues") %>%
    group_by(!! sym(exposure.category)) %>%
    summarise(
      x_N = length(!! sym(exposure.category)),
      x_mean = list(colMeans(cbind(!!! syms(exposure.cues)))),
      x_SS = list(get_sum_of_centered_squares_from_df(cbind(!!! syms(exposure.cues)), verbose = verbose)))

  categories = unique(exposure[[exposure.category]])
  for (c in categories) {
    posterior <-
      suppressWarnings(
        update_NIW_belief_by_sufficient_statistics_of_one_category(
          prior = prior,
          x_category = c,
          x_mean = exposure[exposure$category == c,]$x_mean[[1]],
          x_SS = exposure[exposure$category == c,]$x_SS[[1]],
          x_N = exposure[exposure$category == c,]$x_N[[1]],
          noise_treatment = noise_treatment,
          lapse_treatment = "no_lapses",
          method = "label-certain",
          verbose = verbose))

    prior <- posterior
  }

  return(prior)
}


#' @rdname update_NIW_ideal_adaptor_incrementally
#' @export
update_NIW_beliefs_incrementally <- function(...){
  dots <- list(...)
  assert_that(all(is.null(dots[["noise_treatment"]]), is.null(dots[["lapse_treatment"]])),
              msg = "NIW beliefs do not have noise or lapse rates. Perhaps you meant to use update_NIW_ideal_adaptor_incrementally()?")

  update_NIW_ideal_adaptor_incrementally(..., noise_treatment = "no_noise", lapse_treatment = "no_lapses")
}


