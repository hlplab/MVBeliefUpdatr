#' @importFrom dplyr arrange first rename_with

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
#' @rdname update_model_decision_bias
#' @importFrom rlang .data .env
#' @export
update_model_decision_bias_by_one_observation <- function(
    model,
    beta,
    x,
    x_category,
    noise_treatment = if (!is.null(first(model$Sigma_noise))) "marginalize" else "no_noise",
    lapse_treatment = "no_lapses",
    update_prior = T,
    verbose = F
) {
  # Binding variables that RMD Check gets confused about otherwise
  # (since they are in non-standard evaluations)
  observationID <- response <- delta_logodds <- NULL

  assert_that(all(is_scalar_character(noise_treatment)), is_scalar_character(lapse_treatment))
  if (any(noise_treatment != "no_noise", lapse_treatment != "no_lapses")) {
    # implement check that this is a model
  }

  model %>%
    left_join(
      # The response variable that is returned by the get_categorization function tells us *how often* each category response is
      # expected to be observed. The x_category argument tells what the (supervising) category label is.
      get_categorization_from_model(
        x = x, model = model,
        decision_rule = "proportional", noise_treatment = noise_treatment, lapse_treatment = lapse_treatment) %>%
        # Calculate the amount of change (in log-odds) that the observed category label (x_category) causes.
        # Here we are assuming that each observation leads to change proportional to its surprisal (the
        # surprisal experienced when seeing the category label, x_category):
        mutate(delta_logodds = .env$beta * sum(ifelse(.data$category == .env$x_category, -log2(.data$response) * .data$response, 0))),
      # If one wants to use an error signal that is either 0 or 1 (rather than the prediction error) then
      # delta_logodds would simply be the sum of all error signals multiplied by beta:
      # mutate(delta_logodds = beta * sum(ifelse(category == x_category, 0 * response, 1 * response))),
      by = "category") %>%
    # For the category that was actually observed (x_category) increase the response bias by delta_logodds.
    # Subtract that same *total* amount from the response bias of all other categories (which were not observed),
    # by decreasing the response bias of each of those categories by delta_logodds divided by the number of those
    # categories. Prior to rounding errors, this keeps the sum of the response bias across all categories at 1.
    mutate(
      lapse_bias =
        logit2probability(
          probability2logit(.data$lapse_bias) +
            ifelse(.data$category == .env$x_category, +.data$delta_logodds, -.data$delta_logodds / (length(.data$category) - 1))),
      # correct for rounding errors by re-normalizing
      lapse_bias = .data$lapse_bias / sum(.data$lapse_bias)) %>%
    { if (update_prior) {
      mutate(
        .,
        prior =
          logit2probability(
            probability2logit(.data$prior) +
              ifelse(.data$category == .env$x_category, +.data$delta_logodds, -.data$delta_logodds / (length(.data$category) - 1))),
        prior = .data$prior / sum(.data$prior))
    } else . } %>%
    select(-c(observationID, x, response, delta_logodds)) %>%
    ungroup()
}

#' Update model's decision biases based on exposure data.
#'
#' Returns the model with updated decision biases.
#'
#' @param model A \code{\link[=is.model]{model}} object with decision biases.
#' @param exposure \code{data.frame} or \code{tibble} with exposure data. Each row is assumed to contain one observation.
#' @param exposure.category Name of variable in \code{data} that contains the category information. (default: "category")
#' @param exposure.cues Name(s) of variables in \code{data} that contain the cue information. By default these cue names are
#' extracted from the prior object.
#' @param exposure.order Name of variable in \code{data} that contains the order of the exposure data. If `NULL` the
#' exposure data is assumed to be in the order in which it should be presented.
#' @param noise_treatment Determines whether and how multivariate Gaussian noise is considered during categorization.
#' See \code{\link{update_model_decision_bias_by_one_observation}}.
#' @param lapse_treatment Determines whether attentional lapses can occur during which no updating occurs.
#' See \code{\link{update_model_decision_bias_by_one_observation}}.
#' @param keep.update_history Should the history of the updating be stored and returned? If so, the output is
#' tibble with the one model for each exposure observation. This is useful, for example, if one wants to
#' visualize the changes in the category parameters, posterior predictive, categorization function, or alike across time.
#' (default: `TRUE`)
#' @param keep.exposure_data Should the input data be included in the output? If `FALSE` then only the category and cue
#' columns will be kept. If `TRUE` then all columns will be kept. (default: `FALSE`)
#' @param verbose Should more informative output be provided?
#'
#' @return An model object.
#'
#' @seealso \code{\link{update_model_decision_bias_by_one_observation}}, which is called by \code{update_model_decision_bias_incrementally}
#' @export
update_model_decision_bias_incrementally <- function(
    model,
    beta,
    exposure,
    exposure.category = "category",
    exposure.cues = get_cue_labels_from_model(model),
    exposure.order = NULL,
    noise_treatment = if (!is.null(first(model$Sigma_noise))) "marginalize" else "no_noise",
    lapse_treatment = "no_lapses",
    keep.update_history = TRUE,
    keep.exposure_data = FALSE,
    verbose = FALSE
){
  if (verbose) message("Assuming that category variable in model is called category.")
  if (lapse_treatment == "marginalize")
    warning("Using lapse_treatment == 'marginalize' can result in updating by *fractions* of observations, which might not be wellformed.\n", call. = FALSE)

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
    model %<>%
    mutate(observation.n = 0)

  posterior <- model
  for (i in 1:nrow(exposure)) {
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

    if (keep.update_history)
      model <- rbind(model, posterior %>% mutate(observation.n = .env$i))
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
