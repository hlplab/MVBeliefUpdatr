#' @include S7-family-coercion.R
#' @include S7-core-uvg-classes.R
#' @include S7-core-nix-classes.R
#' @include S7-core-muvg-classes.R
#' @include S7-core-mnix-classes.R
#' @include S7-core-mvg-classes.R
#' @include S7-core-niw-classes.R
#' @include S7-core-exemplar-classes.R
NULL

# -----------------------------------------------------------------------------
# Dynamic Object Construction from Data
# -----------------------------------------------------------------------------

#' Construct S7 Category Representations, Templates, and Models from Data
#'
#' Dynamic constructors that inspect the input `type` (or `family`) argument
#' and dispatch to the appropriate family-specific constructor from data.
#'
#' @name new-objects-from-data
#' @rdname new-objects-from-data
#' @param data A data frame containing observation cues and category labels.
#'   For individual category representations, data must contain observations
#'   from exactly one category.
#' @param type Family type string: `"UVG"`, `"NIX"`, `"MUVG"`, `"MNIX"`,
#'   `"MVG"`, `"NIW"`, or `"EXEMPLAR"`.
#' @param category Character string giving the category column name in `data`.
#'   Defaults to `"category"`.
#' @param cues Character vector of cue column names in `data`.
#' @param category_prior Optional numeric vector of prior category
#'   probabilities summing to 1.
#' @param decision_rule Categorization decision rule: `"sampling"` or
#'   `"argmax"`. Defaults to `"sampling"`.
#' @param lapse_rate Numeric scalar lapse probability in `[0, 1]`. Defaults
#'   to 0.
#' @param lapse_bias Optional numeric vector of lapse category probabilities.
#' @param Sigma_noise Optional perceptual noise covariance matrix.
#' @param noise_treatment Treatment of noise: `"no_noise"`, `"sample"`, or
#'   `"marginalize"`.
#' @param lapse_treatment Treatment of lapses: `"no_lapses"`, `"sample"`, or
#'   `"marginalize"`.
#' @param ... Additional family-specific parameters (such as `kappa`, `nu`,
#'   `m_0`, `S_0`, `c`, or `bandwidth`).
#'
#' @return An S7 representation, template, or cognitive model object.
#' @seealso [family-uvg], [family-nix], [family-muvg], [family-mnix],
#'   [family-mvg], [family-niw], [family-exemplar]
#' @export
new_category_representation_from_data <- function(
  data,
  type,
  category = "category",
  cues,
  ...
) {
  type <- toupper(type)
  constructor <- switch(
    type,
    UVG = new_uvg_category_representation_from_data,
    NIX = new_nix_category_representation_from_data,
    MUVG = new_muvg_category_representation_from_data,
    MNIX = new_mnix_category_representation_from_data,
    MVG = new_mvg_category_representation_from_data,
    NIW = new_niw_category_representation_from_data,
    EXEMPLAR = new_exemplar_category_representation_from_data,
    .stop("type must be one of UVG, NIX, MUVG, MNIX, MVG, NIW, or EXEMPLAR.")
  )
  constructor(data = data, category = category, cues = cues, ...)
}

#' @rdname new-objects-from-data
#' @export
new_category_representation_template_from_data <- function(
  data,
  type,
  category = "category",
  cues,
  ...
) {
  type <- toupper(type)
  constructor <- switch(
    type,
    UVG = new_uvg_category_representation_template_from_data,
    NIX = new_nix_category_representation_template_from_data,
    MUVG = new_muvg_category_representation_template_from_data,
    MNIX = new_mnix_category_representation_template_from_data,
    MVG = new_mvg_category_representation_template_from_data,
    NIW = new_niw_category_representation_template_from_data,
    EXEMPLAR = new_exemplar_category_representation_template_from_data,
    .stop("type must be one of UVG, NIX, MUVG, MNIX, MVG, NIW, or EXEMPLAR.")
  )
  constructor(data = data, category = category, cues = cues, ...)
}

#' @rdname new-objects-from-data
#' @export
new_model_from_data <- function(
  data,
  type,
  category = "category",
  cues,
  ...
) {
  type <- toupper(type)
  constructor <- switch(
    type,
    UVG = new_uvg_ideal_observer_from_data,
    NIX = new_nix_ideal_adaptor_from_data,
    MUVG = new_muvg_ideal_observer_from_data,
    MNIX = new_mnix_ideal_adaptor_from_data,
    MVG = new_mvg_ideal_observer_from_data,
    NIW = new_niw_ideal_adaptor_from_data,
    EXEMPLAR = new_exemplar_model_from_data,
    .stop("type must be one of UVG, NIX, MUVG, MNIX, MVG, NIW, or EXEMPLAR.")
  )
  constructor(data = data, category = category, cues = cues, ...)
}
