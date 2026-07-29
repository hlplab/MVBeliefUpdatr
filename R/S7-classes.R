# S7 foundation: base classes and core generic scaffolding.
# This file is intentionally minimal and non-breaking.

#' @importFrom S7 new_class new_generic method S7_inherits
NULL

# -------------------------
# Base class hierarchy
# -------------------------

#' MVBU Core S7 Class Architecture
#'
#' The MVBeliefUpdatr S7 organizes model objects around an explicit
#' compositional structure:
#'
#' - [MVBU_CategoryRepresentation] defines the structure for one category-level
#'   representational object. This could be a uni- or multivariate Gaussian, 
#'   a mixture of Gaussians, or an exemplar-based representation.
#' - [MVBU_CategoryRepresentationTemplate] collects one or more
#'   category-representation objects into a validated set used by a model. 
#'   This corresponds to the notion of "templates" in e.g., \cite{nearey-assmann07}. 
#' - [MVBU_CognitiveModel] combines category-template structure with model-level
#'   decision behavior (`decision_rule`, `category_prior`, `lapse_rate`, `lapse_bias`).
#' - [MVBU_ModelDistribution] represents distributions over models, which are 
#'   needed, for example, to represent researchers' uncertainty about a model as 
#'   resulting, for example, from fitting models to behavioral data. 
#'
#' Family-level semantics are layered on top of these base classes. Ideal observer/adaptor pairings follow:
#'
#' - [UVG_IdealObserver] `UVG` <-> [NIX_IdealAdaptor] `NIX`
#' - [MUVG_IdealObserver] `MUVG` <-> [MNIX_IdealAdaptor] `MNIX`
#' - [MVG_IdealObserver] `MVG` <-> [NIW_IdealAdaptor] `NIW`
#'
#' while [Exemplar_Model] `EXEMPLAR` is treated as a standalone family.
#'
#' @name MVBU-core-classes
#' @keywords internal
#' @references
#' \insertAllCited{}
NULL

#' @rdname MVBU-core-classes
#' @section MVBU_Object:
#' Root S7 base class for MVBeliefUpdatr objects.
MVBU_Object <- S7::new_class("MVBU_Object", package = NULL)

#' @rdname MVBU-core-classes
#' @section MVBU_CategoryRepresentation:
#' Abstract category-level representational class.
#'
#' Expected properties:
#' - `category_labels`: one or more category labels associated with the object.
#' - `cue_labels`: one or more cue-dimension labels.
#' - `category_likelihood_function`: family-specific likelihood function placeholder.
#' - `metadata`: optional auxiliary metadata list.
MVBU_CategoryRepresentation <- S7::new_class(
  "MVBU_CategoryRepresentation",
  package = NULL,
  parent = MVBU_Object,
  properties = list(
    category_likelihood_function = S7::class_function,
    metadata = S7::class_list
  ),
  validator = function(self) {
    if (!is.list(self@metadata)) {
      return("metadata must be a list.")
    }

    label_information <- .mvbu_label_information(self@metadata)
    if (is.null(label_information$category)) {
      label_information$category <- character(0)
    }
    if (is.null(label_information$cue)) {
      label_information$cue <- character(0)
    }

    if (length(label_information$category) < 1) {
      return("metadata$label_information$category must contain at least one label.")
    }
    if (length(label_information$cue) < 1) {
      return("metadata$label_information$cue must contain at least one label.")
    }
    NULL
  }
)

#' @rdname MVBU-core-classes
#' @section MVBU_CategoryRepresentationTemplate:
#' Container class for per-category representation objects.
#'
#' Expected properties:
#' - `representations`: list of objects inheriting from
#'   `MVBU_CategoryRepresentation`.
#' - `metadata`: optional template-level metadata.
#'
#' Validators enforce at least one representation and (if named) unique,
#' non-empty representation names.
MVBU_CategoryRepresentationTemplate <- S7::new_class(
  "MVBU_CategoryRepresentationTemplate",
  package = NULL,
  parent = MVBU_Object,
  properties = list(
    representations = S7::class_list,
    metadata = S7::class_list
  ),
  validator = function(self) {
    if (length(self@representations) < 1) {
      return("representations must contain at least one category representation object.")
    }

    if (!all(vapply(self@representations, function(r) S7::S7_inherits(r, MVBU_CategoryRepresentation), logical(1)))) {
      return("all representations entries must inherit from MVBU_CategoryRepresentation.")
    }

    repr_names <- names(self@representations)
    if (!is.null(repr_names)) {
      if (length(repr_names) != length(self@representations) || any(repr_names == "") || anyDuplicated(repr_names) > 0) {
        return("if representations are named, names must be non-empty and unique.")
      }
    }

    NULL
  }
)

#' @rdname MVBU-core-classes
#' @section MVBU_CognitiveModel:
#' Abstract model class binding a category-representation template to
#' model-level decision and uncertainty parameters.
#'
#' Expected properties:
#' - `category_template`: `MVBU_CategoryRepresentationTemplate`.
#' - `category_posterior_functions`: named list of treatment-specific posterior functions.
#' - `decision_rule`: scalar character rule label.
#' - `category_prior`: numeric probability vector over represented categories.
#' - `lapse_behavior`: list with `lapse_rate`, `lapse_bias`, and `lapse_treatment`.
#' - `noise_behavior`: list with `Sigma_noise` and `noise_treatment`.
#'   `Sigma_noise` is stored as a square matrix and may be supplied as a vector
#'   (which is converted to a diagonal matrix) or as a matrix.
#' - `metadata`: optional model metadata.
MVBU_CognitiveModel <- S7::new_class(
  "MVBU_CognitiveModel",
  package = NULL,
  parent = MVBU_Object,
  properties = list(
    category_template = MVBU_CategoryRepresentationTemplate,
    category_posterior_functions = S7::class_list,
    decision_rule = S7::class_character,
    category_prior = S7::class_numeric,
    lapse_behavior = S7::class_list,
    noise_behavior = S7::class_list,
    metadata = S7::class_list
  ),
  validator = function(self) {
    n_repr <- length(self@category_template@representations)

    if (length(self@decision_rule) != 1) {
      return("decision_rule must be a scalar character value.")
    }

    lapse_behavior <- self@lapse_behavior
    if (!is.list(lapse_behavior) || !all(c("lapse_rate", "lapse_bias", "lapse_treatment") %in% names(lapse_behavior))) {
      return("lapse_behavior must be a list containing lapse_rate, lapse_bias, and lapse_treatment.")
    }

    lapse_rate <- as.numeric(lapse_behavior$lapse_rate)
    lapse_bias <- as.numeric(lapse_behavior$lapse_bias)
    lapse_treatment <- as.character(lapse_behavior$lapse_treatment)
    if (length(lapse_rate) != 1 || lapse_rate < 0 || lapse_rate > 1) {
      return("lapse_rate must be a scalar numeric in [0, 1].")
    }

    if (length(lapse_bias) != n_repr) {
      return("lapse_bias length must match the number of category representations.")
    }

    if (length(lapse_bias) > 0) {
      if (any(lapse_bias < 0) || any(lapse_bias > 1)) {
        return("lapse_bias entries must be in [0, 1].")
      }
      if (abs(sum(lapse_bias) - 1) > MVBU_PROB_TOL) {
        return("lapse_bias entries must sum to 1.")
      }
    }

    if (!lapse_treatment %in% c("no_lapses", "sample", "marginalize")) {
      return("lapse_treatment must be one of 'no_lapses', 'sample', or 'marginalize'.")
    }

    noise_behavior <- self@noise_behavior
    if (!is.list(noise_behavior) || !all(c("Sigma_noise", "noise_treatment") %in% names(noise_behavior))) {
      return("noise_behavior must be a list containing Sigma_noise and noise_treatment.")
    }

    Sigma_noise <- noise_behavior$Sigma_noise
    noise_treatment <- as.character(noise_behavior$noise_treatment)
    if (!is.null(Sigma_noise)) {
      if (!is.matrix(Sigma_noise)) {
        return("Sigma_noise must be NULL or a matrix.")
      }
      if (!is.numeric(Sigma_noise)) {
        return("Sigma_noise must be NULL or a numeric matrix.")
      }
      if (any(!is.finite(Sigma_noise))) {
        return("Sigma_noise entries must be finite.")
      }
      if (any(Sigma_noise < 0)) {
        return("Sigma_noise entries must be non-negative.")
      }
      cue_labels <- get_cue_labels(self@category_template)
      if (nrow(Sigma_noise) != length(cue_labels) || ncol(Sigma_noise) != length(cue_labels)) {
        return("Sigma_noise dimensions must match the number of cue labels.")
      }
    }

    if (!noise_treatment %in% c("no_noise", "sample", "marginalize")) {
      return("noise_treatment must be one of 'no_noise', 'sample', or 'marginalize'.")
    }

    if (length(self@category_prior) != n_repr) {
      return("category_prior length must match the number of category representations.")
    }

    if (any(self@category_prior < 0) || any(self@category_prior > 1)) {
      return("category_prior entries must be in [0, 1].")
    }

    if (abs(sum(self@category_prior) - 1) > MVBU_PROB_TOL) {
      return("category_prior entries must sum to 1.")
    }

    # Category association can be by explicit names or by order.
    repr_names <- names(self@category_template@representations)
    if (!is.null(repr_names)) {
      if (length(repr_names) != n_repr || any(repr_names == "") || anyDuplicated(repr_names) > 0) {
        return("if representations are named, names must be non-empty and unique.")
      }

      prior_names <- names(self@category_prior)
      if (!is.null(prior_names)) {
        if (length(prior_names) != n_repr || any(prior_names == "") || anyDuplicated(prior_names) > 0) {
          return("if category_prior is named, names must be non-empty and unique.")
        }
        if (!setequal(prior_names, repr_names)) {
          return("if category_prior is named, names must match category_template names.")
        }
      }

      lapse_names <- names(lapse_bias)
      if (!is.null(lapse_names)) {
        if (length(lapse_names) != n_repr || any(lapse_names == "") || anyDuplicated(lapse_names) > 0) {
          return("if lapse_bias is named, names must be non-empty and unique.")
        }
        if (!setequal(lapse_names, repr_names)) {
          return("if lapse_bias is named, names must match category_template names.")
        }
      }
    }

    NULL
  }
)

#' @rdname MVBU-core-classes
#' @section MVBU_ModelDistribution:
#' Abstract inferred-model/distribution class for family-specific fitted or
#' posterior-oriented objects.
#'
#' Expected properties:
#' - `model_family`: scalar family identifier.
#' - `cache`: list for optional computed intermediates.
#' - `metadata`: list for provenance/schema tags.
#' - `group_label`: scalar grouping label for grouped workflows.
MVBU_ModelDistribution <- S7::new_class(
  "MVBU_ModelDistribution",
  package = NULL,
  parent = MVBU_Object,
  properties = list(
    model_family = S7::class_character,
    cache = S7::class_list,
    metadata = S7::class_list,
    group_label = S7::class_character
  ),
  validator = function(self) {
    if (length(self@model_family) != 1) {
      return("model_family must be a scalar character value.")
    }
    if (length(self@group_label) != 1) {
      return("group_label must be a scalar character value.")
    }
    NULL
  }
)

# -------------------------
# Family extension hooks
# -------------------------

#' MVBU S7 Family Glossary
#'
#' 1) `UVG` = Univariate Gaussian ([UVG_IdealObserver])
#' 2) `NIX` = Normal-Inverse-chi^2 over UVG parameters ([NIX_IdealAdaptor])
#' 3) `MUVG` = Cue integration over multiple Univariate Gaussian ([MUVG_IdealObserver])
#' 4) `MNIX` = Cue integration over multiple Normal-Inverse-chi^2 ([MNIX_IdealAdaptor])
#' 5) `MVG` = Multivariate Gaussian ([MVG_IdealObserver])
#' 6) `NIW` = Normal-Inverse-Wishart over MVG parameters ([NIW_IdealAdaptor])
#' 7) `EXEMPLAR` = exemplar-based model ([Exemplar_Model])
#'
#' "Ideal Observer" and "Ideal Adaptor" are model names used in this
#' package's terminology.
#' @keywords internal

.mvbu_family_registry <- new.env(parent = emptyenv())
.mvbu_family_registry$families <- list(
  UVG = list(
    category_representation = "UVG_CategoryRepresentation",
    cognitive_model = "UVG_IdealObserver",
    model_distribution = "UVG_IdealObserverDistribution"
  ),
  NIX = list(
    category_representation = "NIX_CategoryRepresentation",
    cognitive_model = "NIX_IdealAdaptor",
    model_distribution = "NIX_IdealAdaptorDistribution"
  ),
  MUVG = list(
    category_representation = "MUVG_CategoryRepresentation",
    cognitive_model = "MUVG_IdealObserver",
    model_distribution = "MUVG_IdealObserverDistribution"
  ),
  MNIX = list(
    category_representation = "MNIX_CategoryRepresentation",
    cognitive_model = "MNIX_IdealAdaptor",
    model_distribution = "MNIX_IdealAdaptorDistribution"
  ),
  MVG = list(
    category_representation = "MVG_CategoryRepresentation",
    cognitive_model = "MVG_IdealObserver",
    model_distribution = "MVG_IdealObserverDistribution"
  ),
  NIW = list(
    category_representation = "NIW_CategoryRepresentation",
    cognitive_model = "NIW_IdealAdaptor",
    model_distribution = "NIW_IdealAdaptorDistribution"
  ),
  EXEMPLAR = list(
    category_representation = "Exemplar_CategoryRepresentation",
    cognitive_model = "Exemplar_Model",
    model_distribution = "Exemplar_ModelDistribution"
  )
)

.normalize_family_name <- function(family) {
  if (!is.character(family) || length(family) != 1 || nchar(family) == 0) {
    stop("family must be a non-empty scalar character value.", call. = FALSE)
  }
  toupper(family)
}

#' Register an MVBU model family
#' @keywords internal
register_model_family <- function(family, category_representation_class, cognitive_model_class, model_distribution_class) {
  family <- .normalize_family_name(family)

  class_fields <- list(
    category_representation = category_representation_class,
    cognitive_model = cognitive_model_class,
    model_distribution = model_distribution_class
  )

  if (!all(vapply(class_fields, function(x) is.character(x) && length(x) == 1 && nchar(x) > 0, logical(1)))) {
    stop("category_representation_class, cognitive_model_class, and model_distribution_class must be non-empty scalar character values.", call. = FALSE)
  }

  .mvbu_family_registry$families[[family]] <- class_fields
  invisible(TRUE)
}

#' Get registered MVBU model families
#' @keywords internal
get_registered_model_families <- function() {
  sort(names(.mvbu_family_registry$families))
}

#' Get class registration for an MVBU model family
#' @keywords internal
get_model_family_registration <- function(family) {
  family <- .normalize_family_name(family)
  registration <- .mvbu_family_registry$families[[family]]
  if (is.null(registration)) {
    stop(paste0("No model family registered for '", family, "'."), call. = FALSE)
  }
  registration
}

#' Register Stan-family extension hooks
#' @keywords internal
register_stan_family_hooks <- function(
    family,
    stanfit_class = NULL,
    staninput_class = NULL,
    bridge_methods = character(0),
    dependency_rationale = character(0)
) {
  family <- .normalize_family_name(family)

  if (!is.null(stanfit_class) && (!is.character(stanfit_class) || length(stanfit_class) != 1)) {
    stop("stanfit_class must be NULL or a scalar character value.", call. = FALSE)
  }
  if (!is.null(staninput_class) && (!is.character(staninput_class) || length(staninput_class) != 1)) {
    stop("staninput_class must be NULL or a scalar character value.", call. = FALSE)
  }
  if (!is.character(bridge_methods)) {
    stop("bridge_methods must be a character vector.", call. = FALSE)
  }
  if (!is.character(dependency_rationale)) {
    stop("dependency_rationale must be a character vector.", call. = FALSE)
  }

  if (is.null(.mvbu_family_registry$stan_hooks)) {
    .mvbu_family_registry$stan_hooks <- list()
  }

  .mvbu_family_registry$stan_hooks[[family]] <- list(
    stanfit_class = stanfit_class,
    staninput_class = staninput_class,
    bridge_methods = bridge_methods,
    dependency_rationale = dependency_rationale
  )

  invisible(TRUE)
}

#' Get Stan-family extension hooks
#' @keywords internal
get_stan_family_hooks <- function(family = NULL) {
  hooks <- .mvbu_family_registry$stan_hooks
  if (is.null(hooks)) {
    return(list())
  }

  if (is.null(family)) {
    return(hooks)
  }

  family <- .normalize_family_name(family)
  hooks[[family]]
}

# Baseline bridge-watchlist hooks (Phase 1 contract).
register_stan_family_hooks(
  family = "NIW",
  bridge_methods = c("get_stanfit", "as_stanfit", "get_draws", "summary", "print", "loo", "posterior::as_draws_df"),
  dependency_rationale = "Align with rstan/tidybayes workflows without adding dependencies beyond demonstrated usage."
)

# -------------------------
# Family-specific representation classes
# -------------------------

#' UVG Category Representation
#'
#' Category-level representation for a univariate Gaussian (`UVG`) family.
#'
#' Assumed parameterization:
#' - `mu`: scalar mean.
#' - `sigma2`: scalar variance, constrained to be positive.
#'
#' Structural assumptions:
#' - exactly one cue dimension (`length(cue_labels) == 1`).
#'
#' Category likelihood:
#' for observation \eqn{x},
#'
#' \deqn{p(x\mid c) = \mathcal{N}(x;\mu_c,\sigma_c^2).}
#'
#' With `log = TRUE`, the implementation returns
#' \eqn{\log p(x\mid c)}.
#'
#' @seealso [UVG_IdealObserver], [NIX_CategoryRepresentation]
#' @name UVG-CategoryRepresentation-class
#' @keywords internal
UVG_CategoryRepresentation <- S7::new_class(
  "UVG_CategoryRepresentation",
  parent = MVBU_CategoryRepresentation,
  properties = list(
    mu = S7::class_numeric,
    sigma2 = S7::class_numeric
  ),
  validator = function(self) {
    label_information <- .mvbu_label_information(self@metadata)
    if (length(label_information$cue) != 1) {
      return("UVG representation must describe a single cue dimension.")
    }
    if (length(self@mu) != 1) {
      return("UVG mu must be a scalar.")
    }
    if (length(self@sigma2) != 1 || self@sigma2 <= 0) {
      return("UVG sigma2 must be a positive scalar.")
    }

    NULL
  }
)

#' NIX Category Representation
#'
#' Category-level representation for a Normal-Inverse-chi^2 (`NIX`) family,
#' used as an adaptor-style uncertainty parameterization over univariate
#' Gaussian category structure.
#'
#' Assumed parameterization:
#' - `m`: scalar location hyperparameter.
#' - `kappa`: positive scalar precision scaling.
#' - `nu`: positive scalar degrees-of-freedom-like quantity.
#' - `sigma2`: positive scalar scale/variance parameter.
#'
#' Structural assumptions:
#' - exactly one cue dimension (`length(cue_labels) == 1`).
#'
#' Category likelihood:
#' the posterior predictive induced by NIX hyperparameters is Student-\eqn{t}.
#' With
#' \eqn{s_c = \sqrt{\sigma_c^2(\kappa_c+1)/\kappa_c}},
#'
#' \deqn{p(x\mid c) = \frac{1}{s_c} \, t_{\nu_c}\!\left(\frac{x-m_c}{s_c}\right).}
#'
#' With `log = TRUE`, the implementation returns
#' \eqn{\log p(x\mid c)}.
#'
#' @seealso [NIX_IdealAdaptor], [UVG_CategoryRepresentation]
#' @name NIX-CategoryRepresentation-class
#' @keywords internal
NIX_CategoryRepresentation <- S7::new_class(
  "NIX_CategoryRepresentation",
  parent = MVBU_CategoryRepresentation,
  properties = list(
    m = S7::class_numeric,
    kappa = S7::class_numeric,
    nu = S7::class_numeric,
    sigma2 = S7::class_numeric
  ),
  validator = function(self) {
    label_information <- .mvbu_label_information(self@metadata)
    if (length(label_information$cue) != 1) {
      return("NIX representation must describe a single cue dimension.")
    }
    if (length(self@m) != 1) {
      return("NIX m must be a scalar.")
    }
    if (length(self@kappa) != 1 || self@kappa <= 0) {
      return("NIX kappa must be a positive scalar.")
    }
    if (length(self@nu) != 1 || self@nu <= 0) {
      return("NIX nu must be a positive scalar.")
    }
    if (length(self@sigma2) != 1 || self@sigma2 <= 0) {
      return("NIX sigma2 must be a positive scalar.")
    }

    NULL
  }
)

#' MUVG Category Representation
#'
#' Category-level representation for a multi-cue univariate Gaussian integration
#' (`MUVG`) family.
#'
#' Assumed parameterization:
#' - `component_mu`: per-cue means.
#' - `component_sigma2`: per-cue positive variances.
#' - `component_weights`: per-cue integration weights in `[0,1]` summing to 1.
#'
#' Structural assumptions:
#' - all per-cue vectors have equal length.
#' - cue dimensionality equals number of cue-level parameters.
#'
#' If `component_weights` are omitted at construction, ideal independent-cue
#' precision weights are used:
#'
#' \deqn{\lambda_i = \frac{1}{\sigma_i^2},\qquad w_i = \frac{\lambda_i}{\sum_j \lambda_j}.}
#'
#' Under conditional cue independence and integrated cue
#' \eqn{z = \sum_i w_i x_i}, the category likelihood is
#'
#' \deqn{p(z\mid c) = \mathcal{N}\!\left(z;\; \sum_i w_i\mu_i,\; \sum_i w_i^2\sigma_i^2\right).}
#'
#' With `log = TRUE`, the implementation returns
#' \eqn{\log p(z\mid c)}.
#'
#' @seealso [MUVG_IdealObserver], [MNIX_CategoryRepresentation]
#' @name MUVG-CategoryRepresentation-class
#' @keywords internal
MUVG_CategoryRepresentation <- S7::new_class(
  "MUVG_CategoryRepresentation",
  parent = MVBU_CategoryRepresentation,
  properties = list(
    component_mu = S7::class_numeric,
    component_sigma2 = S7::class_numeric,
    component_weights = S7::class_numeric
  ),
  validator = function(self) {
    n_comp <- length(self@component_mu)

    if (n_comp < 1) {
      return("MUVG representation must contain at least one cue component.")
    }
    if (length(self@component_sigma2) != n_comp || length(self@component_weights) != n_comp) {
      return("MUVG component parameter vectors must all have equal length.")
    }
    label_information <- .mvbu_label_information(self@metadata)
    if (length(label_information$cue) != n_comp) {
      return("MUVG cue_labels length must match the number of cue components.")
    }
    if (any(is.na(self@component_mu))) {
      return("MUVG component_mu entries must be numeric values.")
    }
    if (any(is.na(self@component_sigma2))) {
      return("MUVG component_sigma2 entries must be numeric values.")
    }
    if (any(self@component_sigma2 <= 0)) {
      return("MUVG component_sigma2 entries must be positive.")
    }
    if (any(self@component_weights < 0) || any(self@component_weights > 1)) {
      return("MUVG component_weights entries must be in [0, 1].")
    }
    if (abs(sum(self@component_weights) - 1) > MVBU_PROB_TOL) {
      return("MUVG component_weights entries must sum to 1.")
    }

    NULL
  }
)

#' MNIX Category Representation
#'
#' Category-level representation for a Mixture of Normal-Inverse-chi^2
#' (`MNIX`) family.
#'
#' Assumed parameterization:
#' - `component_m`: per-component location hyperparameters.
#' - `component_kappa`: per-component positive precision scalings.
#' - `component_nu`: per-component positive degrees-of-freedom-like values.
#' - `component_sigma2`: per-component positive scale/variance parameters.
#' - `component_weights`: per-component mixture weights in `[0,1]` summing to 1.
#'
#' Structural assumptions:
#' - exactly one cue dimension (`length(cue_labels) == 1`).
#' - all component parameter vectors have equal length.
#'
#' Category likelihood:
#' a weighted finite mixture of Student-\eqn{t} predictive components.
#' For component \eqn{i}, let
#' \eqn{s_i = \sqrt{\sigma_i^2(\kappa_i+1)/\kappa_i}}.
#' Then
#'
#' \deqn{p(x\mid c) = \sum_{i=1}^{K} w_i \frac{1}{s_i} t_{\nu_i}\!\left(\frac{x-m_i}{s_i}\right).}
#'
#' With `log = TRUE`, the implementation returns
#' \eqn{\log p(x\mid c)} computed in log-space via log-sum-exp.
#'
#' @seealso [MNIX_IdealAdaptor], [MUVG_CategoryRepresentation]
#' @name MNIX-CategoryRepresentation-class
#' @keywords internal
MNIX_CategoryRepresentation <- S7::new_class(
  "MNIX_CategoryRepresentation",
  parent = MVBU_CategoryRepresentation,
  properties = list(
    component_m = S7::class_numeric,
    component_kappa = S7::class_numeric,
    component_nu = S7::class_numeric,
    component_sigma2 = S7::class_numeric,
    component_weights = S7::class_numeric
  ),
  validator = function(self) {
    n_comp <- length(self@component_m)

    if (n_comp < 1) {
      return("MNIX representation must contain at least one mixture component.")
    }
    if (length(self@component_kappa) != n_comp ||
        length(self@component_nu) != n_comp ||
        length(self@component_sigma2) != n_comp ||
        length(self@component_weights) != n_comp) {
      return("MNIX component parameter vectors must all have equal length.")
    }
    label_information <- .mvbu_label_information(self@metadata)
    if (length(label_information$cue) != 1) {
      return("MNIX representation must describe a single cue dimension.")
    }
    if (any(self@component_kappa <= 0)) {
      return("MNIX component_kappa entries must be > 0.")
    }
    if (any(self@component_nu <= 0)) {
      return("MNIX component_nu entries must be > 0.")
    }
    if (any(self@component_sigma2 <= 0)) {
      return("MNIX component_sigma2 entries must be > 0.")
    }
    if (any(self@component_weights < 0) || any(self@component_weights > 1)) {
      return("MNIX component_weights entries must be in [0, 1].")
    }
    if (abs(sum(self@component_weights) - 1) > MVBU_PROB_TOL) {
      return("MNIX component_weights entries must sum to 1.")
    }

    NULL
  }
)

#' MVG Category Representation
#'
#' Category-level representation for a Multivariate Gaussian (`MVG`) family.
#'
#' Assumed parameterization:
#' - `mu`: mean vector.
#' - `Sigma`: covariance matrix.
#'
#' Structural assumptions:
#' - `length(mu)` equals cue dimensionality.
#' - `Sigma` is a numeric, symmetric, square matrix with dimensions matching
#'   cue dimensionality.
#'
#' Category likelihood:
#' for cue vector \eqn{\mathbf{x}},
#'
#' \deqn{p(\mathbf{x}\mid c) = \mathcal{N}(\mathbf{x};\boldsymbol{\mu}_c,\mathbf{\Sigma}_c).}
#'
#' With `log = TRUE`, the implementation returns
#' \eqn{\log p(\mathbf{x}\mid c)}.
#'
#' @seealso [MVG_IdealObserver], [NIW_CategoryRepresentation]
#' @name MVG-CategoryRepresentation-class
#' @keywords internal
MVG_CategoryRepresentation <- S7::new_class(
  "MVG_CategoryRepresentation",
  parent = MVBU_CategoryRepresentation,
  properties = list(
    mu = S7::class_numeric,
    Sigma = S7::class_any
  ),
  validator = function(self) {
    label_information <- .mvbu_label_information(self@metadata)
    d <- length(label_information$cue)

    if (length(self@mu) != d) {
      return("MVG mu must have length equal to cue dimensionality.")
    }
    if (!is.matrix(self@Sigma) || !is.numeric(self@Sigma) || nrow(self@Sigma) != d || ncol(self@Sigma) != d) {
      return("MVG Sigma must be a numeric square matrix with cue dimensionality.")
    }
    if (!isTRUE(all.equal(self@Sigma, t(self@Sigma), tolerance = MVBU_PROB_TOL))) {
      return("MVG Sigma must be symmetric.")
    }

    NULL
  }
)

#' NIW Category Representation
#'
#' Category-level representation for a Normal-Inverse-Wishart (`NIW`) family,
#' used as an adaptor-style uncertainty parameterization over multivariate
#' Gaussian category structure.
#'
#' Assumed parameterization:
#' - `m`: location vector.
#' - `kappa`: positive scalar precision scaling.
#' - `nu`: scalar greater than cue dimensionality minus one.
#' - `S`: symmetric numeric square scale matrix.
#'
#' Structural assumptions:
#' - `length(m)` equals cue dimensionality.
#' - `S` dimensions match cue dimensionality.
#'
#' Category likelihood:
#' NIW induces a multivariate Student-\eqn{t} posterior predictive.
#' For cue dimensionality \eqn{D},
#' \eqn{\nu_t = \nu - D + 1} and
#' \eqn{\mathbf{\Sigma}_t = \frac{\kappa+1}{\kappa\nu_t}\mathbf{S}}.
#' Then
#'
#' \deqn{p(\mathbf{x}\mid c) = t_{\nu_t}(\mathbf{x};\mathbf{m},\mathbf{\Sigma}_t).}
#'
#' With `log = TRUE`, the implementation returns
#' \eqn{\log p(\mathbf{x}\mid c)}.
#'
#' @seealso [NIW_IdealAdaptor], [MVG_CategoryRepresentation]
#' @name NIW-CategoryRepresentation-class
#' @keywords internal
NIW_CategoryRepresentation <- S7::new_class(
  "NIW_CategoryRepresentation",
  parent = MVBU_CategoryRepresentation,
  properties = list(
    m = S7::class_numeric,
    kappa = S7::class_numeric,
    nu = S7::class_numeric,
    S = S7::class_any
  ),
  validator = function(self) {
    label_information <- .mvbu_label_information(self@metadata)
    d <- length(label_information$cue)

    if (length(self@m) != d) {
      return("NIW m must have length equal to cue dimensionality.")
    }
    if (length(self@kappa) != 1 || self@kappa <= 0) {
      return("NIW kappa must be a positive scalar.")
    }
    if (length(self@nu) != 1 || self@nu <= (d - 1)) {
      return("NIW nu must be a scalar greater than cue dimensionality minus one.")
    }
    if (!is.matrix(self@S) || !is.numeric(self@S) || nrow(self@S) != d || ncol(self@S) != d) {
      return("NIW S must be a numeric square matrix with cue dimensionality.")
    }
    if (!isTRUE(all.equal(self@S, t(self@S), tolerance = MVBU_PROB_TOL))) {
      return("NIW S must be symmetric.")
    }

    NULL
  }
)

#' Exemplar Category Representation
#'
#' Category-level exemplar-storage representation.
#'
#' Assumed parameterization:
#' - `exemplars`: numeric matrix of stored exemplars (rows are exemplars;
#'   columns correspond to cue dimensions).
#' - `exemplar_weights`: exemplar weighting vector in `[0,1]` summing to 1.
#'
#' Structural assumptions:
#' - at least one exemplar row.
#' - exemplar column count equals cue dimensionality.
#'
#' Category likelihood:
#' implemented as a kernel mixture over stored exemplars,
#'
#' \deqn{p(\mathbf{x}\mid c) = \sum_{i=1}^{N} w_i\,\mathcal{N}(\mathbf{x};\mathbf{e}_i,\mathbf{\Sigma}_{\text{kern}}).}
#'
#' Here \eqn{\mathbf{e}_i} are exemplar rows and
#' \eqn{\mathbf{\Sigma}_{\text{kern}}} is estimated from exemplar covariance,
#' regularized by adding \eqn{\varepsilon\mathbf{I}} (with identity fallback
#' when covariance is unavailable).
#'
#' With `log = TRUE`, the implementation returns
#' \eqn{\log p(\mathbf{x}\mid c)} computed in log-space via log-sum-exp.
#'
#' @seealso [Exemplar_Model]
#' @name Exemplar-CategoryRepresentation-class
#' @keywords internal
Exemplar_CategoryRepresentation <- S7::new_class(
  "Exemplar_CategoryRepresentation",
  parent = MVBU_CategoryRepresentation,
  properties = list(
    exemplars = S7::class_any,
    exemplar_weights = S7::class_numeric
  ),
  validator = function(self) {
    if (!is.matrix(self@exemplars) || !is.numeric(self@exemplars)) {
      return("Exemplar exemplars must be a numeric matrix.")
    }

    n_ex <- nrow(self@exemplars)
    d <- ncol(self@exemplars)
    if (n_ex < 1) {
      return("Exemplar exemplars must contain at least one row.")
    }
    label_information <- .mvbu_label_information(self@metadata)
    if (d != length(label_information$cue)) {
      return("Exemplar exemplar column count must match cue dimensionality.")
    }
    if (length(self@exemplar_weights) != n_ex) {
      return("Exemplar exemplar_weights length must match number of exemplars.")
    }
    if (any(self@exemplar_weights < 0) || any(self@exemplar_weights > 1)) {
      return("Exemplar exemplar_weights entries must be in [0, 1].")
    }
    if (abs(sum(self@exemplar_weights) - 1) > MVBU_PROB_TOL) {
      return("Exemplar exemplar_weights entries must sum to 1.")
    }

    NULL
  }
)

# -------------------------
# Family-specific cognitive model classes
# -------------------------

#' UVG Ideal Observer Model Class
#'
#' Cognitive model class for univariate Gaussian observation. 
#' Inherits all base cognitive-model machinery and constrains
#' category representations to [UVG_CategoryRepresentation] at constructor level.
#'
#' @seealso [NIX_IdealAdaptor], [UVG_IdealObserverDistribution]
#' @name UVG-IdealObserver-class
#' @keywords internal
UVG_IdealObserver <- S7::new_class("UVG_IdealObserver", parent = MVBU_CognitiveModel)

#' NIX Ideal Adaptor Model Class
#'
#' Cognitive model class for NIX-based adaptor inference over univariate
#' Gaussian category structure. Inherits base cognitive-model behavior and
#' expects [NIX_CategoryRepresentation] entries in templates.
#'
#' @seealso [UVG_IdealObserver], [NIX_IdealAdaptorDistribution]
#' @name NIX-IdealAdaptor-class
#' @keywords internal
NIX_IdealAdaptor <- S7::new_class("NIX_IdealAdaptor", parent = MVBU_CognitiveModel)

#' MUVG Ideal Observer Model Class
#'
#' Ideal observer model class for independent-cue Gaussian integration.
#' Category-specific cue integration is represented inside each
#' [MUVG_CategoryRepresentation], so no additional model-level component
#' weights are required.
#'
#' @seealso [MNIX_IdealAdaptor], [MUVG_IdealObserverDistribution]
#' @name MUVG-IdealObserver-class
#' @keywords internal
MUVG_IdealObserver <- S7::new_class("MUVG_IdealObserver", parent = MVBU_CognitiveModel)

#' MNIX Ideal Adaptor Model Class
#'
#' Cognitive model class for mixture adaptor inference. Cue integration is
#' represented inside each [MNIX_CategoryRepresentation], so no additional
#' model-level component weights are required.
#'
#' @seealso [MUVG_IdealObserver], [MNIX_IdealAdaptorDistribution]
#' @name MNIX-IdealAdaptor-class
#' @keywords internal
MNIX_IdealAdaptor <- S7::new_class("MNIX_IdealAdaptor", parent = MVBU_CognitiveModel)

#' MVG Ideal Observer Model Class
#'
#' Cognitive model class for multivariate Gaussian observers. 
#' Inherits all base cognitive-model machinery and is paired with
#' [MVG_CategoryRepresentation] templates.
#'
#' @seealso [NIW_IdealAdaptor], [MVG_IdealObserverDistribution]
#' @name MVG-IdealObserver-class
#' @keywords internal
MVG_IdealObserver <- S7::new_class("MVG_IdealObserver", parent = MVBU_CognitiveModel)

#' NIW Ideal Adaptor Model Class
#'
#' Cognitive model class for NIW-based adaptor inference over multivariate
#' Gaussian category structure. Inherits base cognitive-model behavior and
#' expects [NIW_CategoryRepresentation] entries in templates.
#'
#' @seealso [MVG_IdealObserver], [NIW_IdealAdaptorDistribution]
#' @name NIW-IdealAdaptor-class
#' @keywords internal
NIW_IdealAdaptor <- S7::new_class("NIW_IdealAdaptor", parent = MVBU_CognitiveModel)

#' Exemplar Model Class
#'
#' Cognitive model class for exemplar-based category models. This family is
#' currently treated as standalone rather than an observer/adaptor conjugate
#' pair.
#'
#' @seealso [Exemplar_CategoryRepresentation], [Exemplar_ModelDistribution]
#' @name Exemplar-Model-class
#' @keywords internal
Exemplar_Model <- S7::new_class("Exemplar_Model", parent = MVBU_CognitiveModel)

# -------------------------
# Family-specific model distribution classes
# -------------------------

#' UVG Ideal Observer Distribution Class
#'
#' Family-typed model-distribution class requiring `model_family = "UVG"`.
#'
#' @seealso [UVG_IdealObserver]
#' @name UVG-IdealObserverDistribution-class
#' @keywords internal
UVG_IdealObserverDistribution <- S7::new_class(
  "UVG_IdealObserverDistribution",
  parent = MVBU_ModelDistribution,
  validator = function(self) {
    if (self@model_family != "UVG") {
      return("UVG_IdealObserverDistribution requires model_family = 'UVG'.")
    }
    NULL
  }
)

#' NIX Ideal Adaptor Distribution Class
#'
#' Family-typed model-distribution class requiring `model_family = "NIX"`.
#'
#' @seealso [NIX_IdealAdaptor]
#' @name NIX-IdealAdaptorDistribution-class
#' @keywords internal
NIX_IdealAdaptorDistribution <- S7::new_class(
  "NIX_IdealAdaptorDistribution",
  parent = MVBU_ModelDistribution,
  validator = function(self) {
    if (self@model_family != "NIX") {
      return("NIX_IdealAdaptorDistribution requires model_family = 'NIX'.")
    }
    NULL
  }
)

#' MUVG Ideal Observer Distribution Class
#'
#' Family-typed model-distribution class requiring `model_family = "MUVG"`.
#'
#' @seealso [MUVG_IdealObserver]
#' @name MUVG-IdealObserverDistribution-class
#' @keywords internal
MUVG_IdealObserverDistribution <- S7::new_class(
  "MUVG_IdealObserverDistribution",
  parent = MVBU_ModelDistribution,
  validator = function(self) {
    if (self@model_family != "MUVG") {
      return("MUVG_IdealObserverDistribution requires model_family = 'MUVG'.")
    }
    NULL
  }
)

#' MNIX Ideal Adaptor Distribution Class
#'
#' Family-typed model-distribution class requiring `model_family = "MNIX"`.
#'
#' @seealso [MNIX_IdealAdaptor]
#' @name MNIX-IdealAdaptorDistribution-class
#' @keywords internal
MNIX_IdealAdaptorDistribution <- S7::new_class(
  "MNIX_IdealAdaptorDistribution",
  parent = MVBU_ModelDistribution,
  validator = function(self) {
    if (self@model_family != "MNIX") {
      return("MNIX_IdealAdaptorDistribution requires model_family = 'MNIX'.")
    }
    NULL
  }
)

#' MVG Ideal Observer Distribution Class
#'
#' Family-typed model-distribution class requiring `model_family = "MVG"`.
#'
#' @seealso [MVG_IdealObserver]
#' @name MVG-IdealObserverDistribution-class
#' @keywords internal
MVG_IdealObserverDistribution <- S7::new_class(
  "MVG_IdealObserverDistribution",
  parent = MVBU_ModelDistribution,
  validator = function(self) {
    if (self@model_family != "MVG") {
      return("MVG_IdealObserverDistribution requires model_family = 'MVG'.")
    }
    NULL
  }
)

#' NIW Ideal Adaptor Distribution Class
#'
#' Family-typed model-distribution class requiring `model_family = "NIW"`.
#'
#' @seealso [NIW_IdealAdaptor]
#' @name NIW-IdealAdaptorDistribution-class
#' @keywords internal
NIW_IdealAdaptorDistribution <- S7::new_class(
  "NIW_IdealAdaptorDistribution",
  parent = MVBU_ModelDistribution,
  validator = function(self) {
    if (self@model_family != "NIW") {
      return("NIW_IdealAdaptorDistribution requires model_family = 'NIW'.")
    }
    NULL
  }
)

#' Exemplar Model Distribution Class
#'
#' Family-typed model-distribution class requiring `model_family =
#' "EXEMPLAR"`.
#'
#' @seealso [Exemplar_Model]
#' @name Exemplar-ModelDistribution-class
#' @keywords internal
Exemplar_ModelDistribution <- S7::new_class(
  "Exemplar_ModelDistribution",
  parent = MVBU_ModelDistribution,
  validator = function(self) {
    if (self@model_family != "EXEMPLAR") {
      return("Exemplar_ModelDistribution requires model_family = 'EXEMPLAR'.")
    }
    NULL
  }
)

# -------------------------
# Constructors and helpers
# -------------------------

#' Construct a base MVBU object
#' @keywords internal
new_mvbu_object <- function() {
  MVBU_Object()
}

#' Construct a base inferred-model object
#' @keywords internal
new_model_distribution <- function(model_family, cache = list(), metadata = list(), group_label = "") {
  MVBU_ModelDistribution(
    model_family = .normalize_family_name(model_family),
    cache = cache,
    metadata = metadata,
    group_label = as.character(group_label)
  )
}

#' Construct a UVG model distribution object
#' @keywords internal
new_uvg_model_distribution <- function(cache = list(), metadata = list(), group_label = "") {
  UVG_IdealObserverDistribution(
    model_family = "UVG",
    cache = cache,
    metadata = metadata,
    group_label = as.character(group_label)
  )
}

#' Construct a NIX model distribution object
#' @keywords internal
new_nix_model_distribution <- function(cache = list(), metadata = list(), group_label = "") {
  NIX_IdealAdaptorDistribution(
    model_family = "NIX",
    cache = cache,
    metadata = metadata,
    group_label = as.character(group_label)
  )
}

#' Construct a MUVG model distribution object
#' @keywords internal
new_muvg_model_distribution <- function(cache = list(), metadata = list(), group_label = "") {
  MUVG_IdealObserverDistribution(
    model_family = "MUVG",
    cache = cache,
    metadata = metadata,
    group_label = as.character(group_label)
  )
}

#' Construct a MNIX model distribution object
#' @keywords internal
new_mnix_model_distribution <- function(cache = list(), metadata = list(), group_label = "") {
  MNIX_IdealAdaptorDistribution(
    model_family = "MNIX",
    cache = cache,
    metadata = metadata,
    group_label = as.character(group_label)
  )
}

#' Construct a MVG model distribution object
#' @keywords internal
new_mvg_model_distribution <- function(cache = list(), metadata = list(), group_label = "") {
  MVG_IdealObserverDistribution(
    model_family = "MVG",
    cache = cache,
    metadata = metadata,
    group_label = as.character(group_label)
  )
}

#' Construct a NIW model distribution object
#' @keywords internal
new_niw_model_distribution <- function(cache = list(), metadata = list(), group_label = "") {
  NIW_IdealAdaptorDistribution(
    model_family = "NIW",
    cache = cache,
    metadata = metadata,
    group_label = as.character(group_label)
  )
}

#' Construct an exemplar model distribution object
#' @keywords internal
new_exemplar_model_distribution <- function(cache = list(), metadata = list(), group_label = "") {
  Exemplar_ModelDistribution(
    model_family = "EXEMPLAR",
    cache = cache,
    metadata = metadata,
    group_label = as.character(group_label)
  )
}

#' Construct a base cognitive model object
#'
#' By default, lapse_bias equals category_prior.
#' Category association is by order unless explicit names are supplied.
#' @keywords internal
.mvbu_normalize_sigma_noise <- function(Sigma_noise, cue_labels) {
  if (is.null(Sigma_noise)) {
    return(NULL)
  }

  if (is.matrix(Sigma_noise)) {
    Sigma_noise <- as.matrix(Sigma_noise)
  } else if (is.numeric(Sigma_noise) && length(Sigma_noise) > 0 && is.null(dim(Sigma_noise))) {
    if (length(Sigma_noise) != length(cue_labels)) {
      stop("Sigma_noise length must match the number of cue labels.", call. = FALSE)
    }
    Sigma_noise <- diag(Sigma_noise, nrow = length(Sigma_noise), ncol = length(Sigma_noise))
  } else {
    stop("Sigma_noise must be NULL, a numeric vector, or a matrix.", call. = FALSE)
  }

  if (!is.numeric(Sigma_noise)) {
    stop("Sigma_noise must be numeric.", call. = FALSE)
  }
  if (any(!is.finite(Sigma_noise))) {
    stop("Sigma_noise entries must be finite.", call. = FALSE)
  }
  if (any(Sigma_noise < 0)) {
    stop("Sigma_noise entries must be non-negative.", call. = FALSE)
  }
  if (nrow(Sigma_noise) != length(cue_labels) || ncol(Sigma_noise) != length(cue_labels)) {
    stop("Sigma_noise dimensions must match the number of cue labels.", call. = FALSE)
  }

  Sigma_noise
}

new_cognitive_model <- function(
    category_template = NULL,
    decision_rule = "sampling",
    category_prior = NULL,
    lapse_rate = 0,
    lapse_bias = NULL,
    Sigma_noise = NULL,
    noise_treatment = "no_noise",
    lapse_treatment = "no_lapses",
    metadata = list()
) {
  if (is.null(category_template)) {
    stop("category_template must be supplied.", call. = FALSE)
  }

  if (!S7::S7_inherits(category_template, MVBU_CategoryRepresentationTemplate)) {
    stop("category_template must be an MVBU_CategoryRepresentationTemplate.", call. = FALSE)
  }

  n_repr <- length(category_template@representations)
  repr_names <- names(category_template@representations)

  if (!noise_treatment %in% c("no_noise", "sample", "marginalize")) {
    stop("noise_treatment must be one of 'no_noise', 'sample', or 'marginalize'.", call. = FALSE)
  }
  if (!lapse_treatment %in% c("no_lapses", "sample", "marginalize")) {
    stop("lapse_treatment must be one of 'no_lapses', 'sample', or 'marginalize'.", call. = FALSE)
  }

  if (is.null(category_prior)) {
    category_prior <- rep(1 / n_repr, n_repr)
  }

  if (is.null(lapse_bias)) {
    lapse_bias <- category_prior
  }

  .mvbu_align_probability_vector <- function(values, target_names, arg_name) {
    value_names <- names(values)
    values <- as.numeric(values)

    if (is.null(value_names)) {
      if (!is.null(target_names) && length(target_names) == length(values)) {
        if (any(target_names == "") || anyDuplicated(target_names) > 0) {
          stop(paste0(arg_name, " names cannot be validated because category_template names are missing or invalid."), call. = FALSE)
        }
        names(values) <- target_names
      }
      return(values)
    }

    if (length(value_names) != length(values) || any(value_names == "") || anyDuplicated(value_names) > 0) {
      stop(paste0(arg_name, " names must be non-empty and unique when provided."), call. = FALSE)
    }

    if (is.null(target_names) || length(target_names) != length(values) || any(target_names == "") || anyDuplicated(target_names) > 0) {
      stop(paste0(arg_name, " names cannot be validated because category_template names are missing or invalid."), call. = FALSE)
    }

    if (!setequal(value_names, target_names)) {
      stop(paste0(arg_name, " names must match category_template names."), call. = FALSE)
    }

    values <- values[match(target_names, value_names)]
    names(values) <- target_names
    values
  }

  category_prior <- .mvbu_align_probability_vector(category_prior, repr_names, "category_prior")
  lapse_bias <- .mvbu_align_probability_vector(lapse_bias, repr_names, "lapse_bias")

  if (length(category_prior) != n_repr) {
    stop("category_prior length must match the number of category representations.", call. = FALSE)
  }
  if (length(lapse_bias) != n_repr) {
    stop("lapse_bias length must match the number of category representations.", call. = FALSE)
  }

  cue_labels <- get_cue_labels(category_template)
  Sigma_noise <- .mvbu_normalize_sigma_noise(Sigma_noise, cue_labels)
  noise_behavior <- list(
    Sigma_noise = Sigma_noise,
    noise_treatment = as.character(noise_treatment)
  )
  lapse_behavior <- list(
    lapse_rate = as.numeric(lapse_rate),
    lapse_bias = lapse_bias,
    lapse_treatment = as.character(lapse_treatment)
  )

  model <- MVBU_CognitiveModel(
    category_template = category_template,
    category_posterior_functions = list(),
    decision_rule = as.character(decision_rule),
    category_prior = category_prior,
    lapse_behavior = lapse_behavior,
    noise_behavior = noise_behavior,
    metadata = metadata
  )

  noise_treatments <- if (!is.null(model@noise_behavior$Sigma_noise)) c("no_noise", "sample", "marginalize") else "no_noise"
  lapse_treatments <- c("no_lapses", "sample", "marginalize")
  posterior_functions <- list()
  for (noise_treatment_i in noise_treatments) {
    for (lapse_treatment_i in lapse_treatments) {
      key <- paste(noise_treatment_i, lapse_treatment_i, sep = "__")
      posterior_functions[[key]] <- function(new_data, categories = NULL, .noise_treatment = noise_treatment_i, .lapse_treatment = lapse_treatment_i) {
        .mvbu_posterior_matrix(model, new_data, categories = categories, noise_treatment = .noise_treatment, lapse_treatment = .lapse_treatment)
      }
    }
  }
  model@category_posterior_functions <- posterior_functions

  model
}

#' Construct a category representation set
#' @keywords internal
.mvbu_label_information <- function(metadata = list()) {
  if (!is.list(metadata)) {
    metadata <- as.list(metadata)
  }

  if (!is.null(metadata$label_information) && is.list(metadata$label_information)) {
    label_information <- metadata$label_information
  } else {
    label_information <- list()
  }

  if (is.null(label_information$category) && !is.null(metadata$category)) {
    label_information$category <- metadata$category
  }
  if (is.null(label_information$cue) && !is.null(metadata$cue)) {
    label_information$cue <- metadata$cue
  }

  label_information
}

.mvbu_label_metadata <- function(category_labels = character(), cue_labels = character(), metadata = list()) {
  if (!is.list(metadata)) {
    metadata <- as.list(metadata)
  }

  label_information <- .mvbu_label_information(metadata)
  label_information$category <- as.character(category_labels)
  label_information$cue <- as.character(cue_labels)

  metadata$label_information <- label_information
  metadata$category <- NULL
  metadata$cue <- NULL
  metadata
}

.mvbu_extract_label_metadata <- function(x) {
  if (is.null(x)) {
    return(list(category = character(0), cue = character(0), group = character(0)))
  }

  if (S7::S7_inherits(x, MVBU_CategoryRepresentation)) {
    metadata <- x@metadata
  } else if (S7::S7_inherits(x, MVBU_CategoryRepresentationTemplate)) {
    metadata <- x@metadata
  } else if (S7::S7_inherits(x, MVBU_CognitiveModel)) {
    likelihood <- S7::method(get_category_likelihood, MVBU_CognitiveModel)(x)
    return(.mvbu_extract_label_metadata(likelihood))
  } else {
    return(list(category = character(0), cue = character(0), group = character(0)))
  }

  if (!is.list(metadata)) {
    metadata <- list()
  }

  label_information <- .mvbu_label_information(metadata)

  list(
    category = if (!is.null(label_information$category)) as.character(label_information$category) else character(0),
    cue = if (!is.null(label_information$cue)) as.character(label_information$cue) else character(0),
    group = if (!is.null(label_information$group)) as.character(label_information$group) else character(0)
  )
}

.mvbu_validate_cue_consistency <- function(representations) {
  if (length(representations) < 2) {
    return(invisible(NULL))
  }

  reference_cues <- .mvbu_extract_label_metadata(representations[[1]])$cue
  for (i in seq_along(representations)[-1]) {
    rep_cues <- .mvbu_extract_label_metadata(representations[[i]])$cue
    if (!identical(as.character(rep_cues), as.character(reference_cues))) {
      stop("cue labels must be consistent across all representations in a template.", call. = FALSE)
    }
  }

  invisible(NULL)
}

.mvbu_template_metadata <- function(representations, metadata = list()) {
  if (!is.list(metadata)) {
    metadata <- as.list(metadata)
  }

  category_labels <- unlist(lapply(representations, function(rep) .mvbu_extract_label_metadata(rep)$category), use.names = FALSE)
  cue_labels <- .mvbu_extract_label_metadata(representations[[1]])$cue

  metadata$label_information <- list(
    category = category_labels,
    cue = cue_labels
  )
  metadata$category <- NULL
  metadata$cue <- NULL
  metadata
}

add_category_representation <- function(template, representation, name = NULL) {
  if (!S7::S7_inherits(template, MVBU_CategoryRepresentationTemplate)) {
    stop("template must be an MVBU_CategoryRepresentationTemplate.", call. = FALSE)
  }
  if (!S7::S7_inherits(representation, MVBU_CategoryRepresentation)) {
    stop("representation must be an MVBU_CategoryRepresentation.", call. = FALSE)
  }
  if (!is.null(name) && (length(name) != 1 || !nzchar(name))) {
    stop("name must be a non-empty scalar character value.", call. = FALSE)
  }

  representations <- template@representations
  representations[[length(representations) + 1]] <- representation

  if (!is.null(name)) {
    names(representations)[length(representations)] <- name
  }

  .mvbu_validate_cue_consistency(representations)

  metadata <- as.list(template@metadata)
  rep_labels <- .mvbu_extract_label_metadata(representation)
  template_labels <- .mvbu_extract_label_metadata(template)
  if (length(template_labels$cue) > 0 && length(rep_labels$cue) > 0 && !identical(as.character(template_labels$cue), as.character(rep_labels$cue))) {
    stop("cue labels must be consistent across all representations in a template.", call. = FALSE)
  }
  if (length(template_labels$cue) == 0 && length(rep_labels$cue) > 0) {
    metadata$label_information$cue <- rep_labels$cue
  } else if (length(template_labels$cue) > 0) {
    metadata$label_information$cue <- template_labels$cue
  }

  metadata$label_information$category <- c(template_labels$category, rep_labels$category)
  metadata$category <- NULL
  metadata$cue <- NULL
  MVBU_CategoryRepresentationTemplate(
    representations = representations,
    metadata = metadata
  )
}

new_category_representation_template <- function(representations, metadata = list()) {
  if (!is.list(representations) || length(representations) < 1) {
    stop("representations must contain at least one category representation object.", call. = FALSE)
  }
  if (!all(vapply(representations, function(r) S7::S7_inherits(r, MVBU_CategoryRepresentation), logical(1)))) {
    stop("all representations entries must inherit from MVBU_CategoryRepresentation.", call. = FALSE)
  }

  rep_names <- names(representations)
  if (!is.null(rep_names) && (length(rep_names) != length(representations) || any(rep_names == "") || anyDuplicated(rep_names) > 0)) {
    stop("if representations are named, names must be non-empty and unique.", call. = FALSE)
  }

  template <- NULL
  rep_names <- names(representations)
  for (i in seq_along(representations)) {
    rep_name <- if (!is.null(rep_names)) rep_names[i] else NULL
    if (is.null(template)) {
      template_representations <- list(representations[[i]])
      if (!is.null(rep_name)) {
        names(template_representations) <- rep_name
      }
      template <- MVBU_CategoryRepresentationTemplate(
        representations = template_representations,
        metadata = .mvbu_template_metadata(template_representations, metadata)
      )
    } else {
      template <- add_category_representation(template, representations[[i]], name = rep_name)
    }
  }

  template
}

#' Construct a base representation object
#' @keywords internal
new_category_representation <- function(category_labels, cue_labels, category_likelihood_function = NULL, metadata = list()) {
  if (is.null(category_likelihood_function)) {
    category_likelihood_function <- function(...) stop("category_likelihood not implemented.", call. = FALSE)
  }

  MVBU_CategoryRepresentation(
    category_likelihood_function = category_likelihood_function,
    metadata = .mvbu_label_metadata(as.character(category_labels), as.character(cue_labels), metadata)
  )
}

.as_observation_matrix <- function(x, d, arg_name = "x") {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (is.matrix(x)) {
    if (!is.numeric(x)) {
      stop(arg_name, " must be numeric.", call. = FALSE)
    }
    if (ncol(x) != d) {
      stop(arg_name, " must have ", d, " column(s).", call. = FALSE)
    }
    return(x)
  }

  if (is.atomic(x) && is.numeric(x)) {
    if (d == 1) {
      return(matrix(as.numeric(x), ncol = 1))
    }
    if (length(x) %% d != 0) {
      stop(arg_name, " length must be a multiple of ", d, ".", call. = FALSE)
    }
    return(matrix(as.numeric(x), ncol = d, byrow = TRUE))
  }

  stop(arg_name, " must be a numeric vector, matrix, or data frame.", call. = FALSE)
}

.logsumexp_rows <- function(log_mat) {
  row_max <- apply(log_mat, 1, max)
  row_max + log(rowSums(exp(log_mat - row_max)))
}

.dmvnorm_density <- function(x, mean, Sigma, log = FALSE) {
  x <- .as_observation_matrix(x, d = length(mean), arg_name = "x")
  centered <- sweep(x, 2, mean, "-")
  inv_sigma <- solve(Sigma)
  log_det <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  qf <- rowSums((centered %*% inv_sigma) * centered)
  log_d <- -0.5 * (length(mean) * log(2 * pi) + log_det + qf)
  if (isTRUE(log)) {
    log_d
  } else {
    exp(log_d)
  }
}

.dmvt_density <- function(x, mean, Sigma, df, log = FALSE) {
  x <- .as_observation_matrix(x, d = length(mean), arg_name = "x")
  centered <- sweep(x, 2, mean, "-")
  inv_sigma <- solve(Sigma)
  log_det <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  qf <- rowSums((centered %*% inv_sigma) * centered)
  d <- length(mean)
  log_coef <- lgamma((df + d) / 2) -
    lgamma(df / 2) -
    0.5 * (d * log(df * pi) + log_det)
  log_d <- log_coef + (-(df + d) / 2) * log1p(qf / df)
  if (isTRUE(log)) {
    log_d
  } else {
    exp(log_d)
  }
}

#' Construct a UVG representation object
#' @keywords internal
new_uvg_category_representation <- function(
    category_labels,
    cue_labels,
    mu,
    sigma2,
    metadata = list()
) {
  UVG_CategoryRepresentation(
    category_likelihood_function = {
      mu0 <- as.numeric(mu)
      sigma20 <- as.numeric(sigma2)
      function(x, log = FALSE, noise_treatment = "no_noise", Sigma_noise = NULL) {
        x <- .as_observation_matrix(x, d = 1, arg_name = "x")
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(n = nrow(x), mean = rep(0, ncol(x)), sigma = Sigma_noise)
        }
        noise_variance <- if (!is.null(Sigma_noise) && (identical(noise_treatment, "sample") || identical(noise_treatment, "marginalize"))) as.numeric(Sigma_noise[1, 1]) else 0
        stats::dnorm(x[, 1], mean = mu0, sd = sqrt(sigma20 + noise_variance), log = log)
      }
    },
    metadata = .mvbu_label_metadata(as.character(category_labels), as.character(cue_labels), metadata),
    mu = as.numeric(mu),
    sigma2 = as.numeric(sigma2)
  )
}

#' Construct a NIX representation object
#' @keywords internal
new_nix_category_representation <- function(
    category_labels,
    cue_labels,
    m,
    kappa,
    nu,
    sigma2,
    metadata = list()
) {
  NIX_CategoryRepresentation(
    category_likelihood_function = {
      m0 <- as.numeric(m)
      kappa0 <- as.numeric(kappa)
      nu0 <- as.numeric(nu)
      sigma20 <- as.numeric(sigma2)
      function(x, log = FALSE, noise_treatment = "no_noise", Sigma_noise = NULL) {
        x <- .as_observation_matrix(x, d = 1, arg_name = "x")
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(n = nrow(x), mean = rep(0, ncol(x)), sigma = Sigma_noise)
        }
        noise_variance <- if (!is.null(Sigma_noise) && (identical(noise_treatment, "sample") || identical(noise_treatment, "marginalize"))) as.numeric(Sigma_noise[1, 1]) else 0
        scale_eff <- sqrt((sigma20 * (kappa0 + 1) / kappa0) + noise_variance)
        z <- (x[, 1] - m0) / scale_eff
        if (isTRUE(log)) {
          stats::dt(z, df = nu0, log = TRUE) - log(scale_eff)
        } else {
          stats::dt(z, df = nu0) / scale_eff
        }
      }
    },
    metadata = .mvbu_label_metadata(as.character(category_labels), as.character(cue_labels), metadata),
    m = as.numeric(m),
    kappa = as.numeric(kappa),
    nu = as.numeric(nu),
    sigma2 = as.numeric(sigma2)
  )
}

#' Construct a MUVG representation object
#' @keywords internal
new_muvg_category_representation <- function(
    category_labels,
    cue_labels,
    component_mu,
    component_sigma2,
    component_weights = NULL,
    metadata = list()
) {
  n_comp <- length(component_mu)
  if (n_comp < 1) {
    stop("component_mu must contain at least one element.", call. = FALSE)
  }
  if (!is.numeric(component_mu) || !is.numeric(component_sigma2)) {
    stop("component_mu and component_sigma2 must be numeric.", call. = FALSE)
  }
  if (is.null(component_weights)) {
    precision <- 1 / as.numeric(component_sigma2)
    component_weights <- precision / sum(precision)
  }
  if (!is.numeric(component_weights)) {
    stop("component_weights must be numeric.", call. = FALSE)
  }

  MUVG_CategoryRepresentation(
    category_likelihood_function = {
      mu0 <- as.numeric(component_mu)
      sigma20 <- as.numeric(component_sigma2)
      w0 <- as.numeric(component_weights)
      z_mean <- sum(w0 * mu0)
      z_var <- sum((w0^2) * sigma20)
      function(x, log = FALSE, noise_treatment = "no_noise", Sigma_noise = NULL) {
        x <- .as_observation_matrix(x, d = length(w0), arg_name = "x")
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(n = nrow(x), mean = rep(0, ncol(x)), sigma = Sigma_noise)
        }
        noise_variance <- if (!is.null(Sigma_noise) && (identical(noise_treatment, "sample") || identical(noise_treatment, "marginalize"))) sum((w0^2) * diag(Sigma_noise)) else 0
        z <- as.numeric(x %*% w0)
        stats::dnorm(z, mean = z_mean, sd = sqrt(z_var + noise_variance), log = log)
      }
    },
    metadata = .mvbu_label_metadata(as.character(category_labels), as.character(cue_labels), metadata),
    component_mu = as.numeric(component_mu),
    component_sigma2 = as.numeric(component_sigma2),
    component_weights = as.numeric(component_weights)
  )
}

#' Construct a MNIX representation object
#' @keywords internal
new_mnix_category_representation <- function(
    category_labels,
    cue_labels,
    component_m,
    component_kappa,
    component_nu,
    component_sigma2,
    component_weights = NULL,
    metadata = list()
) {
  n_comp <- length(component_m)
  if (n_comp < 1) {
    stop("component_m must contain at least one element.", call. = FALSE)
  }
  if (!is.numeric(component_m) || !is.numeric(component_kappa) || !is.numeric(component_nu) || !is.numeric(component_sigma2)) {
    stop("component_m, component_kappa, component_nu, and component_sigma2 must be numeric.", call. = FALSE)
  }
  if (is.null(component_weights)) {
    component_weights <- rep(1 / n_comp, n_comp)
  }
  if (!is.numeric(component_weights)) {
    stop("component_weights must be numeric.", call. = FALSE)
  }

  MNIX_CategoryRepresentation(
    category_likelihood_function = {
      m0 <- as.numeric(component_m)
      kappa0 <- as.numeric(component_kappa)
      sigma20 <- as.numeric(component_sigma2)
      nu0 <- as.numeric(component_nu)
      pred_var <- sigma20 * (kappa0 + 1) / kappa0
      w0 <- as.numeric(component_weights)
      function(x, log = FALSE, noise_treatment = "no_noise", Sigma_noise = NULL) {
        x <- .as_observation_matrix(x, d = 1, arg_name = "x")
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(n = nrow(x), mean = rep(0, ncol(x)), sigma = Sigma_noise)
        }
        noise_variance <- if (!is.null(Sigma_noise) && (identical(noise_treatment, "sample") || identical(noise_treatment, "marginalize"))) as.numeric(Sigma_noise[1, 1]) else 0
        component_logdens <- sapply(seq_along(w0), function(i) {
          scale_i <- sqrt(pred_var[i] + noise_variance)
          stats::dt((x[, 1] - m0[i]) / scale_i, df = nu0[i], log = TRUE) - log(scale_i)
        })
        if (is.vector(component_logdens)) {
          component_logdens <- matrix(component_logdens, ncol = length(w0))
        }
        weighted_logdens <- sweep(component_logdens, 2, log(w0), "+")
        log_mix <- .logsumexp_rows(weighted_logdens)
        if (isTRUE(log)) {
          log_mix
        } else {
          exp(log_mix)
        }
      }
    },
    metadata = .mvbu_label_metadata(as.character(category_labels), as.character(cue_labels), metadata),
    component_m = as.numeric(component_m),
    component_kappa = as.numeric(component_kappa),
    component_nu = as.numeric(component_nu),
    component_sigma2 = as.numeric(component_sigma2),
    component_weights = as.numeric(component_weights)
  )
}

#' Construct a MVG representation object
#' @keywords internal
new_mvg_category_representation <- function(
    category_labels,
    cue_labels,
    mu,
    Sigma,
    metadata = list()
) {
  MVG_CategoryRepresentation(
    category_likelihood_function = {
      mu0 <- as.numeric(mu)
      Sigma0 <- as.matrix(Sigma)
      function(x, log = FALSE, noise_treatment = "no_noise", Sigma_noise = NULL) {
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(n = nrow(x), mean = rep(0, ncol(x)), sigma = Sigma_noise)
        }
        if (!is.null(Sigma_noise) && (identical(noise_treatment, "sample") || identical(noise_treatment, "marginalize"))) {
          Sigma_eff <- Sigma0 + Sigma_noise
        } else {
          Sigma_eff <- Sigma0
        }
        .dmvnorm_density(x, mean = mu0, Sigma = Sigma_eff, log = log)
      }
    },
    metadata = .mvbu_label_metadata(as.character(category_labels), as.character(cue_labels), metadata),
    mu = as.numeric(mu),
    Sigma = as.matrix(Sigma)
  )
}

#' Construct a NIW representation object
#' @keywords internal
new_niw_category_representation <- function(
    category_labels,
    cue_labels,
    m,
    kappa,
    nu,
    S,
    metadata = list()
) {
  NIW_CategoryRepresentation(
    category_likelihood_function = {
      m0 <- as.numeric(m)
      kappa0 <- as.numeric(kappa)
      nu0 <- as.numeric(nu)
      S0 <- as.matrix(S)
      d0 <- length(m0)
      df0 <- nu0 - d0 + 1
      Sigma0 <- ((kappa0 + 1) / (kappa0 * df0)) * S0
      function(x, log = FALSE, noise_treatment = "no_noise", Sigma_noise = NULL) {
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(n = nrow(x), mean = rep(0, ncol(x)), sigma = Sigma_noise)
        }
        if (!is.null(Sigma_noise) && (identical(noise_treatment, "sample") || identical(noise_treatment, "marginalize"))) {
          Sigma_eff <- Sigma0 + Sigma_noise
        } else {
          Sigma_eff <- Sigma0
        }
        .dmvt_density(x, mean = m0, Sigma = Sigma_eff, df = df0, log = log)
      }
    },
    metadata = .mvbu_label_metadata(as.character(category_labels), as.character(cue_labels), metadata),
    m = as.numeric(m),
    kappa = as.numeric(kappa),
    nu = as.numeric(nu),
    S = as.matrix(S)
  )
}

#' Construct an Exemplar representation object
#' @keywords internal
new_exemplar_category_representation <- function(
    category_labels,
    cue_labels,
    exemplars,
    exemplar_weights = NULL,
    metadata = list()
) {
  exemplars <- as.matrix(exemplars)
  if (!is.numeric(exemplars)) {
    stop("exemplars must be numeric.", call. = FALSE)
  }
  n_ex <- nrow(exemplars)
  if (is.null(exemplar_weights)) {
    exemplar_weights <- rep(1 / n_ex, n_ex)
  }
  if (!is.numeric(exemplar_weights)) {
    stop("exemplar_weights must be numeric.", call. = FALSE)
  }

  Exemplar_CategoryRepresentation(
    category_likelihood_function = {
      ex0 <- exemplars
      w0 <- as.numeric(exemplar_weights)
      d0 <- ncol(ex0)
      n0 <- nrow(ex0)
      if (n0 > 1) {
        Sigma0 <- stats::cov(ex0)
      } else {
        Sigma0 <- diag(1, d0)
      }
      if (!is.matrix(Sigma0) || any(!is.finite(Sigma0))) {
        Sigma0 <- diag(1, d0)
      }
      Sigma0 <- Sigma0 + diag(MVBU_PROB_TOL, d0)

      function(x, log = FALSE, noise_treatment = "no_noise", Sigma_noise = NULL) {
        x <- .as_observation_matrix(x, d = d0, arg_name = "x")
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(n = nrow(x), mean = rep(0, ncol(x)), sigma = Sigma_noise)
        }
        if (!is.null(Sigma_noise) && (identical(noise_treatment, "sample") || identical(noise_treatment, "marginalize"))) {
          Sigma_eff <- Sigma0 + Sigma_noise
        } else {
          Sigma_eff <- Sigma0
        }
        log_dens_by_exemplar <- sapply(seq_len(n0), function(i) {
          .dmvnorm_density(x, mean = ex0[i, ], Sigma = Sigma_eff, log = TRUE)
        })
        if (is.vector(log_dens_by_exemplar)) {
          log_dens_by_exemplar <- matrix(log_dens_by_exemplar, ncol = n0)
        }
        weighted_logdens <- sweep(log_dens_by_exemplar, 2, log(w0), "+")
        log_mix <- .logsumexp_rows(weighted_logdens)
        if (isTRUE(log)) {
          log_mix
        } else {
          exp(log_mix)
        }
      }
    },
    metadata = .mvbu_label_metadata(as.character(category_labels), as.character(cue_labels), metadata),
    exemplars = exemplars,
    exemplar_weights = as.numeric(exemplar_weights)
  )
}

.new_family_cognitive_model <- function(
    category_template = NULL,
    category_representation_class,
    model_class,
    family_label,
    decision_rule = "sampling",
    category_prior = NULL,
    lapse_rate = 0,
    lapse_bias = NULL,
    Sigma_noise = NULL,
    noise_treatment = "no_noise",
    lapse_treatment = "no_lapses",
    metadata = list()
) {
  if (is.null(category_template)) {
    stop("category_template must be supplied.", call. = FALSE)
  }

  if (!S7::S7_inherits(category_template, MVBU_CategoryRepresentationTemplate)) {
    stop("category_template must be an MVBU_CategoryRepresentationTemplate.", call. = FALSE)
  }

  family_representations <- category_template@representations

  if (!all(vapply(family_representations, function(r) S7::S7_inherits(r, category_representation_class), logical(1)))) {
    stop(
      paste0("All category_likelihood_template entries must inherit from ", family_label, " category representation class."),
      call. = FALSE
    )
  }

  base_model <- new_cognitive_model(
    category_template = category_template,
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment,
    metadata = metadata
  )

  model <- model_class(
    category_template = base_model@category_template,
    decision_rule = base_model@decision_rule,
    category_prior = base_model@category_prior,
    lapse_behavior = base_model@lapse_behavior,
    noise_behavior = base_model@noise_behavior,
    metadata = base_model@metadata
  )

  model@category_posterior_functions <- base_model@category_posterior_functions

  model
}

#' Construct a UVG ideal observer object
#' @keywords internal
new_uvg_ideal_observer <- function(
    category_template = NULL,
    decision_rule = "sampling",
    category_prior = NULL,
    lapse_rate = 0,
    lapse_bias = NULL,
    Sigma_noise = NULL,
    noise_treatment = "no_noise",
    lapse_treatment = "no_lapses",
    metadata = list()
) {
  .new_family_cognitive_model(
    category_template = category_template,
    category_representation_class = UVG_CategoryRepresentation,
    model_class = UVG_IdealObserver,
    family_label = "UVG",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    metadata = metadata
  )
}

#' Construct a NIX ideal adaptor object
#' @keywords internal
new_nix_ideal_adaptor <- function(
    category_template = NULL,
    decision_rule = "sampling",
    category_prior = NULL,
    lapse_rate = 0,
    lapse_bias = NULL,
    Sigma_noise = NULL,
    noise_treatment = "no_noise",
    lapse_treatment = "no_lapses",
    metadata = list()
) {
  .new_family_cognitive_model(
    category_template = category_template,
    category_representation_class = NIX_CategoryRepresentation,
    model_class = NIX_IdealAdaptor,
    family_label = "NIX",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    metadata = metadata
  )
}

#' Construct a MUVG ideal observer object
#' @keywords internal
new_muvg_ideal_observer <- function(
    category_template = NULL,
    decision_rule = "sampling",
    category_prior = NULL,
    lapse_rate = 0,
    lapse_bias = NULL,
    Sigma_noise = NULL,
    noise_treatment = "no_noise",
    lapse_treatment = "no_lapses",
    metadata = list()
) {
  .new_family_cognitive_model(
    category_template = category_template,
    category_representation_class = MUVG_CategoryRepresentation,
    model_class = MUVG_IdealObserver,
    family_label = "MUVG",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    metadata = metadata
  )
}

#' Construct a MNIX ideal adaptor object
#' @keywords internal
new_mnix_ideal_adaptor <- function(
    category_template = NULL,
    decision_rule = "sampling",
    category_prior = NULL,
    lapse_rate = 0,
    lapse_bias = NULL,
    Sigma_noise = NULL,
    noise_treatment = "no_noise",
    lapse_treatment = "no_lapses",
    metadata = list()
) {
  .new_family_cognitive_model(
    category_template = category_template,
    category_representation_class = MNIX_CategoryRepresentation,
    model_class = MNIX_IdealAdaptor,
    family_label = "MNIX",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    metadata = metadata
  )
}

#' Construct a MVG ideal observer object
#' @keywords internal
new_mvg_ideal_observer <- function(
    category_template = NULL,
    decision_rule = "sampling",
    category_prior = NULL,
    lapse_rate = 0,
    lapse_bias = NULL,
    Sigma_noise = NULL,
    noise_treatment = "no_noise",
    lapse_treatment = "no_lapses",
    metadata = list()
) {
  .new_family_cognitive_model(
    category_template = category_template,
    category_representation_class = MVG_CategoryRepresentation,
    model_class = MVG_IdealObserver,
    family_label = "MVG",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    metadata = metadata
  )
}

#' Construct a NIW ideal adaptor object
#' @keywords internal
new_niw_ideal_adaptor <- function(
    category_template = NULL,
    decision_rule = "sampling",
    category_prior = NULL,
    lapse_rate = 0,
    lapse_bias = NULL,
    Sigma_noise = NULL,
    noise_treatment = "no_noise",
    lapse_treatment = "no_lapses",
    metadata = list()
) {
  .new_family_cognitive_model(
    category_template = category_template,
    category_representation_class = NIW_CategoryRepresentation,
    model_class = NIW_IdealAdaptor,
    family_label = "NIW",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    metadata = metadata
  )
}

#' Construct an Exemplar model object
#' @keywords internal
new_exemplar_model <- function(
    category_template = NULL,
    decision_rule = "sampling",
    category_prior = NULL,
    lapse_rate = 0,
    lapse_bias = NULL,
    Sigma_noise = NULL,
    noise_treatment = "no_noise",
    lapse_treatment = "no_lapses",
    metadata = list()
) {
  .new_family_cognitive_model(
    category_template = category_template,
    category_representation_class = Exemplar_CategoryRepresentation,
    model_class = Exemplar_Model,
    family_label = "EXEMPLAR",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    metadata = metadata
  )
}

#' Family-specific model constructors are canonical in v1
#' @keywords internal
new_model <- function(...) {
  stop(
    "Use family-specific constructors in v1 (e.g., new_niw_ideal_adaptor_model_from_data()).",
    call. = FALSE
  )
}

#' Validate a model object
#' @keywords internal
validate_object <- function(x) {
  if (!tryCatch(S7::S7_inherits(x, MVBU_Object), error = function(e) FALSE)) {
    stop("x must be an S7 object.", call. = FALSE)
  }
  TRUE
}

#' Safe validation check
#' @keywords internal
is_valid <- function(x) {
  tryCatch({
    validate_object(x)
    TRUE
  }, error = function(e) FALSE)
}

.mvbu_not_implemented <- function(generic_name, class_name) {
  stop(
    paste0("Method for ", generic_name, "() not yet implemented for class ", class_name, "."),
    call. = FALSE
  )
}

