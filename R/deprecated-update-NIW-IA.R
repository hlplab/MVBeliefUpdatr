#' @include S7-update-model.R S7-generics.R
NULL

# Deprecated NIW update compatibility wrappers now require S7 representations or models.

#' Deprecated: update_NIW_belief_kappa
#'
#' @description `r lifecycle::badge("deprecated")`
#' `update_NIW_belief_kappa()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [update_category_representation()] instead.
#'
#' @param kappa_0 Prior kappa value.
#' @param x_N Number of observations.
#' @seealso [update_category_representation()]
#' @keywords internal
#' @export
update_NIW_belief_kappa <- function(kappa_0, x_N) {
  lifecycle::deprecate_warn("0.1.0", "update_NIW_belief_kappa()", with = "update_category_representation()")
  .update_NIW_category_representation_kappa(kappa_0, x_N)
}

#' Deprecated: update_NIW_belief_nu
#'
#' @description `r lifecycle::badge("deprecated")`
#' `update_NIW_belief_nu()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [update_category_representation()] instead.
#'
#' @param nu_0 Prior nu value.
#' @param x_N Number of observations.
#' @seealso [update_category_representation()]
#' @keywords internal
#' @export
update_NIW_belief_nu <- function(nu_0, x_N) {
  lifecycle::deprecate_warn("0.1.0", "update_NIW_belief_nu()", with = "update_category_representation()")
  .update_NIW_category_representation_nu(nu_0, x_N)
}

#' Deprecated: update_NIW_belief_m
#'
#' @description `r lifecycle::badge("deprecated")`
#' `update_NIW_belief_m()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [update_category_representation()] instead.
#'
#' @param kappa_0 Prior kappa value.
#' @param m_0 Prior mean vector.
#' @param x_N Number of observations.
#' @param x_mean Observation mean.
#' @seealso [update_category_representation()]
#' @keywords internal
#' @export
update_NIW_belief_m <- function(kappa_0, m_0, x_N, x_mean) {
  lifecycle::deprecate_warn("0.1.0", "update_NIW_belief_m()", with = "update_category_representation()")
  .update_NIW_category_representation_m(kappa_0, m_0, x_N, x_mean)
}

#' Deprecated: update_NIW_belief_S
#'
#' @description `r lifecycle::badge("deprecated")`
#' `update_NIW_belief_S()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [update_category_representation()] instead.
#'
#' @param kappa_0 Prior kappa value.
#' @param m_0 Prior mean vector.
#' @param S_0 Prior scatter matrix.
#' @param x_N Number of observations.
#' @param x_mean Observation mean.
#' @param x_SS Centered observation sum-of-squares matrix.
#' @seealso [update_category_representation()]
#' @keywords internal
#' @export
update_NIW_belief_S <- function(kappa_0, m_0, S_0, x_N, x_mean, x_SS) {
  lifecycle::deprecate_warn("0.1.0", "update_NIW_belief_S()", with = "update_category_representation()")
  .update_NIW_category_representation_S(kappa_0, m_0, S_0, x_N, x_mean, x_SS)
}


#' Deprecated: update_NIW_belief_by_sufficient_statistics_of_one_category
#'
#' @description `r lifecycle::badge("deprecated")`
#' `update_NIW_belief_by_sufficient_statistics_of_one_category()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [update_category_representation()] instead.
#'
#' @param prior_model An NIW S7 category representation.
#' @param x_category Deprecated and ignored for S7 representations.
#' @param x_mean Observation mean.
#' @param x_SS Centered observation sum-of-squares matrix.
#' @param x_N Number of observations.
#' @param ... Deprecated compatibility arguments.
#' @return An updated NIW category representation.
#' @seealso [update_category_representation()]
#' @keywords internal
#' @export
update_NIW_belief_by_sufficient_statistics_of_one_category <- function(prior_model, x_category = NULL, x_mean, x_SS, x_N, ...) {
  lifecycle::deprecate_warn("0.1.0", "update_NIW_belief_by_sufficient_statistics_of_one_category()", with = "update_category_representation()")
  .assert_true(S7::S7_inherits(prior_model, NIW_CategoryRepresentation), msg = "prior_model must be an NIW_CategoryRepresentation; legacy tibble inputs are no longer supported.")
  .update_NIW_category_representation_by_sufficient_statistics(prior_model, x_mean, x_SS, x_N)
}

#' Deprecated: update_NIW_belief_by_one_observation
#'
#' @description `r lifecycle::badge("deprecated")`
#' `update_NIW_belief_by_one_observation()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [update_category_representation()] or [update_template()] instead.
#'
#' @param prior_model An NIW category representation or ideal-adaptor model.
#' @param x_category Category label for the observation.
#' @param x Numeric observation vector.
#' @param noise_treatment Noise treatment.
#' @param lapse_treatment Lapse treatment.
#' @param method Updating method.
#' @param verbose Whether to print additional output.
#' @return An updated S7 representation or model.
#' @seealso [update_category_representation()], [update_template()]
#' @keywords internal
#' @export
update_NIW_belief_by_one_observation <- function(prior_model, x_category, x, noise_treatment = "no_noise", lapse_treatment = "no_lapses", method = "label-certain", verbose = FALSE) {
  lifecycle::deprecate_warn("0.1.0", "update_NIW_belief_by_one_observation()", with = "update_category_representation()")
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
#'
#' @description `r lifecycle::badge("deprecated")`
#' `update_NIW_ideal_adaptor_incrementally()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [update_template()] instead.
#'
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
#' @seealso [update_template()]
#' @keywords internal
#' @export
update_NIW_ideal_adaptor_incrementally <- function(prior_model, exposure, exposure.category = "category", exposure.cues = get_cue_labels(prior_model), exposure.order = NULL, noise_treatment = "no_noise", lapse_treatment = "no_lapses", method = "label-certain", keep.update_history = TRUE, keep.exposure_data = FALSE, verbose = FALSE) {
  lifecycle::deprecate_warn("0.1.0", "update_NIW_ideal_adaptor_incrementally()", with = "update_template()")
  .assert_true(S7::S7_inherits(prior_model, NIW_IdealAdaptor), msg = "prior_model must be an NIW_IdealAdaptor; legacy tibble inputs are no longer supported.")
  observations <- exposure
  if (exposure.category != "category") names(observations)[names(observations) == exposure.category] <- "category"
  update_template(prior_model, observations, updating = "incremental", keep_history = keep.update_history, lapse_treatment = lapse_treatment, noise_treatment = noise_treatment, update_method = method)
}

#' Deprecated: update_NIW_ideal_adaptor_batch
#'
#' @description `r lifecycle::badge("deprecated")`
#' `update_NIW_ideal_adaptor_batch()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [update_template()] instead.
#'
#' @seealso [update_template()]
#' @keywords internal
#' @rdname update_NIW_ideal_adaptor_incrementally
#' @export
update_NIW_ideal_adaptor_batch <- function(prior_model, exposure, exposure.category = "category", exposure.cues = get_cue_labels(prior_model), noise_treatment = "no_noise", verbose = FALSE) {
  lifecycle::deprecate_warn("0.1.0", "update_NIW_ideal_adaptor_batch()", with = "update_template()")
  .assert_true(S7::S7_inherits(prior_model, NIW_IdealAdaptor), msg = "prior_model must be an NIW_IdealAdaptor; legacy tibble inputs are no longer supported.")
  observations <- exposure
  if (exposure.category != "category") names(observations)[names(observations) == exposure.category] <- "category"
  update_template(prior_model, observations, updating = "batch", lapse_treatment = "no_lapses", noise_treatment = noise_treatment, update_method = "label-certain")
}

#' Deprecated: update_NIW_beliefs_incrementally
#'
#' @description `r lifecycle::badge("deprecated")`
#' `update_NIW_beliefs_incrementally()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [update_template()] instead.
#'
#' @seealso [update_template()]
#' @keywords internal
#' @rdname update_NIW_ideal_adaptor_incrementally
#' @export
update_NIW_beliefs_incrementally <- function(prior_model, exposure, ...) {
  lifecycle::deprecate_warn("0.1.0", "update_NIW_beliefs_incrementally()", with = "update_template()")
  update_template(prior_model, exposure, updating = "incremental", ...)
}

