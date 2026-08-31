#' @include S7-update-model.R S7-generics.R
NULL

# Deprecated NIW update compatibility wrappers now require S7 representations or models.

#' @name update_NIW_parameters
#' @title Deprecated NIW parameter helpers
#' @param kappa_0 Prior kappa value.
#' @param nu_0 Prior nu value.
#' @param m_0 Prior mean vector.
#' @param S_0 Prior scatter matrix.
#' @param x_N Number of observations.
#' @param x_mean Observation mean.
#' @param x_SS Centered observation sum-of-squares matrix.
#' @keywords internal
NULL

#' Deprecated: update_NIW_belief_kappa
#' @rdname update_NIW_parameters
#' @export
update_NIW_belief_kappa <- function(kappa_0, x_N) {
  lifecycle::deprecate_warn("0.0.9", "update_NIW_belief_kappa()", "update_NIW_category_representation_kappa()")
  .update_NIW_category_representation_kappa(kappa_0, x_N)
}

#' Deprecated: update_NIW_belief_nu
#' @rdname update_NIW_parameters
#' @export
update_NIW_belief_nu <- function(nu_0, x_N) {
  lifecycle::deprecate_warn("0.0.9", "update_NIW_belief_nu()", "update_NIW_category_representation_nu()")
  .update_NIW_category_representation_nu(nu_0, x_N)
}

#' Deprecated: update_NIW_belief_m
#' @rdname update_NIW_parameters
#' @export
update_NIW_belief_m <- function(kappa_0, m_0, x_N, x_mean) {
  lifecycle::deprecate_warn("0.0.9", "update_NIW_belief_m()", "update_NIW_category_representation_m()")
  .update_NIW_category_representation_m(kappa_0, m_0, x_N, x_mean)
}

#' Deprecated: update_NIW_belief_S
#' @rdname update_NIW_parameters
#' @export
update_NIW_belief_S <- function(kappa_0, m_0, S_0, x_N, x_mean, x_SS) {
  lifecycle::deprecate_warn("0.0.9", "update_NIW_belief_S()", "update_NIW_category_representation_S()")
  .update_NIW_category_representation_S(kappa_0, m_0, S_0, x_N, x_mean, x_SS)
}

#' Deprecated: update_NIW_belief_by_sufficient_statistics_of_one_category
#' @param prior_model An NIW S7 category representation.
#' @param x_category Deprecated and ignored for S7 representations.
#' @param x_mean Observation mean.
#' @param x_SS Centered observation sum-of-squares matrix.
#' @param x_N Number of observations.
#' @param ... Deprecated compatibility arguments.
#' @return An updated NIW category representation.
#' @export
update_NIW_belief_by_sufficient_statistics_of_one_category <- function(prior_model, x_category = NULL, x_mean, x_SS, x_N, ...) {
  lifecycle::deprecate_warn("0.0.9", "update_NIW_belief_by_sufficient_statistics_of_one_category()", "update_category_representation()")
  .assert_true(S7::S7_inherits(prior_model, NIW_CategoryRepresentation), msg = "prior_model must be an NIW_CategoryRepresentation; legacy tibble inputs are no longer supported.")
  .update_NIW_category_representation_by_sufficient_statistics(prior_model, x_mean, x_SS, x_N)
}

#' Deprecated: update_NIW_belief_by_one_observation
#' @param prior_model An NIW category representation or ideal-adaptor model.
#' @param x_category Category label for the observation.
#' @param x Numeric observation vector.
#' @param noise_treatment Noise treatment.
#' @param lapse_treatment Lapse treatment.
#' @param method Updating method.
#' @param verbose Whether to print additional output.
#' @return An updated S7 representation or model.
#' @export
update_NIW_belief_by_one_observation <- function(prior_model, x_category, x, noise_treatment = "no_noise", lapse_treatment = "no_lapses", method = "label-certain", verbose = FALSE) {
  lifecycle::deprecate_warn("0.0.9", "update_NIW_belief_by_one_observation()", "update_category_representation()")
  if (S7::S7_inherits(prior_model, NIW_CategoryRepresentation)) {
    observation <- as.numeric(x)
    return(update_category_representation(
      prior_model,
      x_N = 1,
      x_mean = observation,
      x_SS = matrix(0, nrow = length(observation), ncol = length(observation))
    ))
  }
  .assert_true(S7::S7_inherits(prior_model, NIW_IdealAdaptor), msg = "prior_model must be an NIW S7 representation or ideal-adaptor model; legacy tibble inputs are no longer supported.")
  observations <- as.data.frame(as.list(as.numeric(x)))
  names(observations) <- get_cue_labels(prior_model)
  observations$category <- x_category
  observations <- observations[, c("category", get_cue_labels(prior_model)), drop = FALSE]
  update_template(
    prior_model,
    observations,
    updating = "incremental",
    keep_history = FALSE,
    lapse_treatment = lapse_treatment,
    noise_treatment = noise_treatment,
    update_method = method
  )
}

#' Deprecated: update_NIW_ideal_adaptor_incrementally
#' @param prior_model An NIW S7 ideal-adaptor model.
#' @param exposure Observation data.
#' @param exposure.category Name of the category column.
#' @param exposure.cues Cue columns.
#' @param exposure.order Deprecated ordering column.
#' @param noise_treatment Noise treatment.
#' @param lapse_treatment Lapse treatment.
#' @param method Updating method.
#' @param keep.update_history Whether to return intermediate models.
#' @param keep.exposure_data Deprecated and ignored.
#' @param verbose Whether to print additional output.
#' @return An updated S7 model or model history.
#' @export
update_NIW_ideal_adaptor_incrementally <- function(prior_model, exposure, exposure.category = "category", exposure.cues = get_cue_labels(prior_model), exposure.order = NULL, noise_treatment = "no_noise", lapse_treatment = "no_lapses", method = "label-certain", keep.update_history = TRUE, keep.exposure_data = FALSE, verbose = FALSE) {
  lifecycle::deprecate_warn("0.0.9", "update_NIW_ideal_adaptor_incrementally()", "update_template()")
  .assert_true(S7::S7_inherits(prior_model, NIW_IdealAdaptor), msg = "prior_model must be an NIW_IdealAdaptor; legacy tibble inputs are no longer supported.")
  observations <- exposure
  if (exposure.category != "category") names(observations)[names(observations) == exposure.category] <- "category"
  update_template(prior_model, observations, updating = "incremental", keep_history = keep.update_history, lapse_treatment = lapse_treatment, noise_treatment = noise_treatment, update_method = method)
}

#' Deprecated: update_NIW_ideal_adaptor_batch
#' @rdname update_NIW_ideal_adaptor_incrementally
#' @export
update_NIW_ideal_adaptor_batch <- function(prior_model, exposure, exposure.category = "category", exposure.cues = get_cue_labels(prior_model), noise_treatment = "no_noise", verbose = FALSE) {
  lifecycle::deprecate_warn("0.0.9", "update_NIW_ideal_adaptor_batch()", "update_template()")
  .assert_true(S7::S7_inherits(prior_model, NIW_IdealAdaptor), msg = "prior_model must be an NIW_IdealAdaptor; legacy tibble inputs are no longer supported.")
  observations <- exposure
  if (exposure.category != "category") names(observations)[names(observations) == exposure.category] <- "category"
  update_template(prior_model, observations, updating = "batch", lapse_treatment = "no_lapses", noise_treatment = noise_treatment, update_method = "label-certain")
}

#' Deprecated: update_NIW_beliefs_incrementally
#' @rdname update_NIW_ideal_adaptor_incrementally
#' @export
update_NIW_beliefs_incrementally <- function(prior_model, exposure, ...) {
  lifecycle::deprecate_warn("0.0.9", "update_NIW_beliefs_incrementally()", "update_template()")
  update_template(prior_model, exposure, updating = "incremental", ...)
}
