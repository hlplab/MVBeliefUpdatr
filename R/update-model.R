# Functions for incremental bias change model
probability2logit <- function(p, refcat = 1)
  if (length(p) == 2)
    qlogis(p) else
      message("not yet defined")

logit2probability <- function(l, refcat = 1)
  if (length(l) == 2)
    plogis(l) else
      ifelse(1:length(l) == refcat, 1, exp(l)) / (1 + sum(exp(l[-refcat])))


#' Update model's decision biases based on labelled input
#'
#' Returns model with updated decision biases. This updating process is order-sensitive. For each input observation x,
#' a posterior expectation for each category is obtained based on the specified noise and lapse treatment. This expectation
#' is used to calculate the surprisal (in bits) of the category label x_category. The surprisal is weighted by the learning
#' rate beta to increase the (log-odds of the) decision bias for that category.
#'
#' @param model A model object with decision biases.
#' @param beta The learning rate with which decision biases change. Might be integrated into model objects in the future.
#' @param x_category Character. The label of the category that is to be updated.
#' @param x The cues of single observation.
#' @param x_category The category label(s) of one or more observations.
#' @param noise_treatment Determines whether multivariate Gaussian noise is added to the input.
#' If `NULL`, no noise is added during the updating. If "sample" then a sample of
#' noise is added to the input. If "marginalize" then each observation is transformed into the marginal distribution
#' that results from convolving the input with noise. This latter option might be helpful, for example, if one is
#' interested in estimating the consequences of noise across individuals. If add_noise is not `NULL` a Sigma_noise
#' column must be present in the NIW_belief object specified as the priors argument.
#' (default: "no_noise")
#' @param lapse_treatment Determines whether attentional lapses can occur during which no updating occurs.
#' Can be "no_lapses", "sample", or "marginalize". If "no_lapses", no lapses occur (even if the model specifies
#' a non-zero `lapse_rate`), and all observations lead to updating. If "sample" or "marginalize", the lapse rate from '
#' the model will be used. For "sample", sampling determines for each observation whether it was a lapse or not.
#' If an observation was a lapse no updating occurs. For "marginalize", 1 - lapse_rate is the proportion of observations
#' that are assumed to be lapsing trials (default: "no_lapses")
#' @param update_prior Should the prior probability of each category be updated along with the decision bias?
#' @param verbose Should more informative output be provided?
#'
#' @return A model object.
#'
#' @keywords updating, decision bias, response bias
#' @examples
#' TBD
#' @rdname update_model_decision_bias
#' @export
update_model_decision_bias_by_one_observation <- function(
    model,
    beta,
    x,
    x_category,
    noise_treatment = "no_noise",
    lapse_treatment = "no_lapses",
    update_prior = T,
    verbose = F
) {
  assert_that(all(is_scalar_character(noise_treatment)), is_scalar_character(lapse_treatment))
  if (any(noise_treatment != "no_noise", lapse_treatment != "no_lapses")) {
    # implement check that this is a model
  }

  model %>%
    left_join(
      # The response variable that is returned by the get_categorization function tells us *how often* each category response is
      # expected to be observed. The x_category argument tells what the (supervising) category label is.
      get_categorization_from_model(x = x, model = model, decision_rule = "proportional", noise_treatment = noise_treatment, lapse_treatment = lapse_treatment) %>%
        # Calculate the amount of change (in log-odds) that the observed category label (x_category) causes.
        # Here we are assuming that each observation leads to change proportional to its surprisal (the
        # surprisal experienced when seeing the category label, x_category):
        mutate(delta_logodds = beta * sum(ifelse(category == x_category, -log2(response) * response, 0))),
      # If one wants to use an error signal that is either 0 or 1 (rather than the prediction error) then
      # delta_logodds would simply be the sum of all error signals multiplied by beta:
      # mutate(delta_logodds = beta * sum(ifelse(category == x_category, 0 * response, 1 * response))),
      by = "category") %>%
    # For the category that was actually observed (x_category) increase the response bias by delta_logodds.
    # Subtract that same *total* amount from the response bias of all other categories (which were not observed),
    # by decreasing the response bias of each of those categories by delta_logodds divided by the number of those
    # categories. Prior to rounding errors, this keeps the sum of the response bias across all categories at 1.
    mutate(
      lapse_bias = logit2probability(probability2logit(lapse_bias) + ifelse(category == x_category, +delta_logodds, -delta_logodds / (length(category) - 1))),
      # correct for rounding errors by re-normalizing
      lapse_bias = lapse_bias / sum(lapse_bias)) %>%
    { if (update_prior) {
      mutate(
        .,
        prior = logit2probability(probability2logit(prior) + ifelse(category == x_category, +delta_logodds, -delta_logodds / (length(category) - 1))),
        prior = prior / sum(prior))
    } else . } %>%
    select(-c(observationID, x, response, delta_logodds)) %>%
    ungroup()
}

#' @rdname update_model_decision_bias
#' @export
update_model_decision_bias_incrementally <- function(
    model,
    beta,
    exposure,
    exposure.category = "category",
    exposure.cues = get_cue_labels_from_model(model),
    exposure.order = NULL,
    noise_treatment = "no_noise",
    lapse_treatment = "no_lapses",
    keep.update_history = TRUE,
    keep.exposure_data = FALSE,
    verbose = FALSE
){
  if (verbose) message("Assuming that category variable in model is called category.")
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

  # Prepare exposure data
  exposure %<>%
    { if (!is.null(exposure.order)) arrange(., !! sym(exposure.order)) else . } %>%
    make_vector_column(exposure.cues, "cues")

  if (keep.update_history)
    prior %<>%
    mutate(observation.n = 0)

  for (i in 1:nrow(exposure)) {
    posterior <- if (keep.update_history) model %>% filter(observation.n == i - 1) else model

    posterior <-
      suppressWarnings(
        update_model_decision_bias_by_one_observation(
          model = posterior,
          beta = beta,
          x = matrix(unlist(exposure[i,][["cues"]]), nrow = 1),
          x_category = exposure[i,][[exposure.category]],
          noise_treatment = noise_treatment,
          lapse_treatment = lapse_treatment,
          verbose = verbose))

    if (keep.update_history) {
      posterior %<>%
        mutate(observation.n = i)
      model <- rbind(model, posterior)
    } else model <- posterior
  }

  if (keep.exposure_data) {
    exposure %<>%
      { if (!is.null(exposure.order))
        rename(., observation.n = !! sym(exposure.order)) else
          mutate(., observation.n = 1:nrow(exposure)) } %>%
      rename_with(~ paste0("observation.", .x), !starts_with("observation.n"))

    model %<>%
      left_join(exposure)
  }

  return(model)
}
