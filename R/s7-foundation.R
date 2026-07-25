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
#' Root S7 base class for MVBeliefUpdatr scaffold objects.
MVBU_Object <- S7::new_class("MVBU_Object")

#' @rdname MVBU-core-classes
#' @section MVBU_CategoryRepresentation:
#' Abstract category-level representational class.
#'
#' Expected properties:
#' - `category_labels`: one or more category labels associated with the object.
#' - `cue_labels`: one or more cue-dimension labels.
#' - `category_likelihood`: family-specific likelihood function placeholder.
#' - `metadata`: optional auxiliary metadata list.
MVBU_CategoryRepresentation <- S7::new_class(
  "MVBU_CategoryRepresentation",
  parent = MVBU_Object,
  properties = list(
    category_labels = S7::class_character,
    cue_labels = S7::class_character,
    category_likelihood = S7::class_function,
    metadata = S7::class_list
  ),
  validator = function(self) {
    if (length(self@category_labels) < 1) {
      return("category_labels must contain at least one label.")
    }
    if (length(self@cue_labels) < 1) {
      return("cue_labels must contain at least one label.")
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
#' - `category_likelihood`: `MVBU_CategoryRepresentationTemplate`.
#' - `decision_rule`: scalar character rule label.
#' - `category_prior`: numeric probability vector over represented categories.
#' - `lapse_rate`: scalar probability in `[0,1]`.
#' - `lapse_bias`: numeric probability vector over represented categories.
#' - `metadata`: optional model metadata.
MVBU_CognitiveModel <- S7::new_class(
  "MVBU_CognitiveModel",
  parent = MVBU_Object,
  properties = list(
    category_likelihood = MVBU_CategoryRepresentationTemplate,
    decision_rule = S7::class_character,
    category_prior = S7::class_numeric,
    lapse_rate = S7::class_numeric,
    lapse_bias = S7::class_numeric,
    metadata = S7::class_list
  ),
  validator = function(self) {
    n_repr <- length(self@category_likelihood@representations)

    if (length(self@decision_rule) != 1) {
      return("decision_rule must be a scalar character value.")
    }

    if (length(self@lapse_rate) != 1 || self@lapse_rate < 0 || self@lapse_rate > 1) {
      return("lapse_rate must be a scalar numeric in [0, 1].")
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

    if (length(self@lapse_bias) != n_repr) {
      return("lapse_bias length must match the number of category representations.")
    }

    if (length(self@lapse_bias) > 0) {
      if (any(self@lapse_bias < 0) || any(self@lapse_bias > 1)) {
        return("lapse_bias entries must be in [0, 1].")
      }
      if (abs(sum(self@lapse_bias) - 1) > MVBU_PROB_TOL) {
        return("lapse_bias entries must sum to 1.")
      }
    }

    # Category association can be by explicit names or by order.
    repr_names <- names(self@category_likelihood@representations)
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
          return("if category_prior is named, names must match category_likelihood names.")
        }
      }

      lapse_names <- names(self@lapse_bias)
      if (!is.null(lapse_names)) {
        if (length(lapse_names) != n_repr || any(lapse_names == "") || anyDuplicated(lapse_names) > 0) {
          return("if lapse_bias is named, names must be non-empty and unique.")
        }
        if (!setequal(lapse_names, repr_names)) {
          return("if lapse_bias is named, names must match category_likelihood names.")
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
# Family-typed subclasses
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
    if (length(self@cue_labels) != 1) {
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
    if (length(self@cue_labels) != 1) {
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
    if (length(self@cue_labels) != n_comp) {
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
    if (length(self@cue_labels) != 1) {
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
    d <- length(self@cue_labels)

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
    d <- length(self@cue_labels)

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
    if (d != length(self@cue_labels)) {
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
new_cognitive_model <- function(
    category_likelihood,
    decision_rule = "sampling",
    category_prior = NULL,
    lapse_rate = 0,
    lapse_bias = NULL,
  metadata = list()
) {
  if (!S7::S7_inherits(category_likelihood, MVBU_CategoryRepresentationTemplate)) {
    stop("category_likelihood must be an MVBU_CategoryRepresentationTemplate.", call. = FALSE)
  }

  n_repr <- length(category_likelihood@representations)
  if (n_repr < 1) {
    stop("category_likelihood must contain at least one element.", call. = FALSE)
  }

  if (is.null(category_prior)) {
    category_prior <- rep(1 / n_repr, n_repr)
  }

  if (is.null(lapse_bias)) {
    lapse_bias <- category_prior
  }

  # If category_likelihood is named, use those names for category_prior/lapse_bias unless
  # already explicitly named.
  repr_names <- names(category_likelihood@representations)
  if (!is.null(repr_names) && length(repr_names) == n_repr && all(repr_names != "")) {
    if (is.null(names(category_prior)) && length(category_prior) == n_repr) names(category_prior) <- repr_names
    if (is.null(names(lapse_bias)) && length(lapse_bias) == n_repr) names(lapse_bias) <- repr_names
  }

  prior_names <- names(category_prior)
  lapse_names <- names(lapse_bias)

  category_prior <- as.numeric(category_prior)
  lapse_bias <- as.numeric(lapse_bias)

  if (!is.null(prior_names)) names(category_prior) <- prior_names
  if (!is.null(lapse_names)) names(lapse_bias) <- lapse_names

  MVBU_CognitiveModel(
    category_likelihood = category_likelihood,
    decision_rule = as.character(decision_rule),
    lapse_rate = as.numeric(lapse_rate),
    category_prior = category_prior,
    lapse_bias = lapse_bias,
    metadata = metadata
  )
}

#' Construct a category representation set
#' @keywords internal
new_category_representation_template <- function(representations, metadata = list()) {
  MVBU_CategoryRepresentationTemplate(
    representations = representations,
    metadata = metadata
  )
}

#' Construct a base representation object
#' @keywords internal
new_category_representation <- function(category_labels, cue_labels, category_likelihood = function(...) stop("category_likelihood not implemented.", call. = FALSE), metadata = list()) {
  MVBU_CategoryRepresentation(
    category_labels = as.character(category_labels),
    cue_labels = as.character(cue_labels),
    category_likelihood = category_likelihood,
    metadata = metadata
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
    category_labels = as.character(category_labels),
    cue_labels = as.character(cue_labels),
    category_likelihood = {
      mu0 <- as.numeric(mu)
      sigma20 <- as.numeric(sigma2)
      function(x, log = FALSE) {
        x <- .as_observation_matrix(x, d = 1, arg_name = "x")
        stats::dnorm(x[, 1], mean = mu0, sd = sqrt(sigma20), log = log)
      }
    },
    metadata = metadata,
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
    category_labels = as.character(category_labels),
    cue_labels = as.character(cue_labels),
    category_likelihood = {
      m0 <- as.numeric(m)
      kappa0 <- as.numeric(kappa)
      nu0 <- as.numeric(nu)
      sigma20 <- as.numeric(sigma2)
      scale0 <- sqrt(sigma20 * (kappa0 + 1) / kappa0)
      function(x, log = FALSE) {
        x <- .as_observation_matrix(x, d = 1, arg_name = "x")
        z <- (x[, 1] - m0) / scale0
        if (isTRUE(log)) {
          stats::dt(z, df = nu0, log = TRUE) - log(scale0)
        } else {
          stats::dt(z, df = nu0) / scale0
        }
      }
    },
    metadata = metadata,
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
    category_labels = as.character(category_labels),
    cue_labels = as.character(cue_labels),
    category_likelihood = {
      mu0 <- as.numeric(component_mu)
      sigma20 <- as.numeric(component_sigma2)
      w0 <- as.numeric(component_weights)
      z_mean <- sum(w0 * mu0)
      z_var <- sum((w0^2) * sigma20)
      function(x, log = FALSE) {
        x <- .as_observation_matrix(x, d = length(w0), arg_name = "x")
        z <- as.numeric(x %*% w0)
        stats::dnorm(z, mean = z_mean, sd = sqrt(z_var), log = log)
      }
    },
    metadata = metadata,
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
    category_labels = as.character(category_labels),
    cue_labels = as.character(cue_labels),
    category_likelihood = {
      m0 <- as.numeric(component_m)
      kappa0 <- as.numeric(component_kappa)
      sigma20 <- as.numeric(component_sigma2)
      nu0 <- as.numeric(component_nu)
      pred_var <- sigma20 * (kappa0 + 1) / kappa0
      w0 <- as.numeric(component_weights)
      function(x, log = FALSE) {
        x <- .as_observation_matrix(x, d = 1, arg_name = "x")
        component_logdens <- sapply(seq_along(w0), function(i) {
          scale_i <- sqrt(pred_var[i])
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
    metadata = metadata,
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
    category_labels = as.character(category_labels),
    cue_labels = as.character(cue_labels),
    category_likelihood = {
      mu0 <- as.numeric(mu)
      Sigma0 <- as.matrix(Sigma)
      function(x, log = FALSE) {
        .dmvnorm_density(x, mean = mu0, Sigma = Sigma0, log = log)
      }
    },
    metadata = metadata,
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
    category_labels = as.character(category_labels),
    cue_labels = as.character(cue_labels),
    category_likelihood = {
      m0 <- as.numeric(m)
      kappa0 <- as.numeric(kappa)
      nu0 <- as.numeric(nu)
      S0 <- as.matrix(S)
      d0 <- length(m0)
      df0 <- nu0 - d0 + 1
      Sigma0 <- ((kappa0 + 1) / (kappa0 * df0)) * S0
      function(x, log = FALSE) {
        .dmvt_density(x, mean = m0, Sigma = Sigma0, df = df0, log = log)
      }
    },
    metadata = metadata,
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
    category_labels = as.character(category_labels),
    cue_labels = as.character(cue_labels),
    category_likelihood = {
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

      function(x, log = FALSE) {
        x <- .as_observation_matrix(x, d = d0, arg_name = "x")
        log_dens_by_exemplar <- sapply(seq_len(n0), function(i) {
          .dmvnorm_density(x, mean = ex0[i, ], Sigma = Sigma0, log = TRUE)
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
    metadata = metadata,
    exemplars = exemplars,
    exemplar_weights = as.numeric(exemplar_weights)
  )
}

.new_family_cognitive_model <- function(
    category_likelihood,
    category_representation_class,
    model_class,
    family_label,
    decision_rule = "sampling",
    category_prior = NULL,
    lapse_rate = 0,
    lapse_bias = NULL,
    metadata = list()
) {
  if (!S7::S7_inherits(category_likelihood, MVBU_CategoryRepresentationTemplate)) {
    stop("category_likelihood must be an MVBU_CategoryRepresentationTemplate.", call. = FALSE)
  }

  family_representations <- category_likelihood@representations

  if (!all(vapply(family_representations, function(r) S7::S7_inherits(r, category_representation_class), logical(1)))) {
    stop(
      paste0("All category_likelihood entries must inherit from ", family_label, " category representation class."),
      call. = FALSE
    )
  }

  base_model <- new_cognitive_model(
    category_likelihood = category_likelihood,
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    metadata = metadata
  )

  model_class(
    category_likelihood = base_model@category_likelihood,
    decision_rule = base_model@decision_rule,
    category_prior = base_model@category_prior,
    lapse_rate = base_model@lapse_rate,
    lapse_bias = base_model@lapse_bias,
    metadata = base_model@metadata
  )
}

#' Construct a UVG ideal observer object
#' @keywords internal
new_uvg_ideal_observer <- function(
    category_likelihood,
    decision_rule = "sampling",
    category_prior = NULL,
    lapse_rate = 0,
    lapse_bias = NULL,
    metadata = list()
) {
  .new_family_cognitive_model(
    category_likelihood = category_likelihood,
    category_representation_class = UVG_CategoryRepresentation,
    model_class = UVG_IdealObserver,
    family_label = "UVG",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    metadata = metadata
  )
}

#' Construct a NIX ideal adaptor object
#' @keywords internal
new_nix_ideal_adaptor <- function(
    category_likelihood,
    decision_rule = "sampling",
    category_prior = NULL,
    lapse_rate = 0,
    lapse_bias = NULL,
    metadata = list()
) {
  .new_family_cognitive_model(
    category_likelihood = category_likelihood,
    category_representation_class = NIX_CategoryRepresentation,
    model_class = NIX_IdealAdaptor,
    family_label = "NIX",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    metadata = metadata
  )
}

#' Construct a MUVG ideal observer object
#' @keywords internal
new_muvg_ideal_observer <- function(
    category_likelihood,
    decision_rule = "sampling",
    category_prior = NULL,
    lapse_rate = 0,
    lapse_bias = NULL,
    metadata = list()
) {
  .new_family_cognitive_model(
    category_likelihood = category_likelihood,
    category_representation_class = MUVG_CategoryRepresentation,
    model_class = MUVG_IdealObserver,
    family_label = "MUVG",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    metadata = metadata
  )
}

#' Construct a MNIX ideal adaptor object
#' @keywords internal
new_mnix_ideal_adaptor <- function(
    category_likelihood,
    decision_rule = "sampling",
    category_prior = NULL,
    lapse_rate = 0,
    lapse_bias = NULL,
    metadata = list()
) {
  .new_family_cognitive_model(
    category_likelihood = category_likelihood,
    category_representation_class = MNIX_CategoryRepresentation,
    model_class = MNIX_IdealAdaptor,
    family_label = "MNIX",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    metadata = metadata
  )
}

#' Construct a MVG ideal observer object
#' @keywords internal
new_mvg_ideal_observer <- function(
    category_likelihood,
    decision_rule = "sampling",
    category_prior = NULL,
    lapse_rate = 0,
    lapse_bias = NULL,
    metadata = list()
) {
  .new_family_cognitive_model(
    category_likelihood = category_likelihood,
    category_representation_class = MVG_CategoryRepresentation,
    model_class = MVG_IdealObserver,
    family_label = "MVG",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    metadata = metadata
  )
}

#' Construct a NIW ideal adaptor object
#' @keywords internal
new_niw_ideal_adaptor <- function(
    category_likelihood,
    decision_rule = "sampling",
    category_prior = NULL,
    lapse_rate = 0,
    lapse_bias = NULL,
    metadata = list()
) {
  .new_family_cognitive_model(
    category_likelihood = category_likelihood,
    category_representation_class = NIW_CategoryRepresentation,
    model_class = NIW_IdealAdaptor,
    family_label = "NIW",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    metadata = metadata
  )
}

#' Construct an Exemplar model object
#' @keywords internal
new_exemplar_model <- function(
    category_likelihood,
    decision_rule = "sampling",
    category_prior = NULL,
    lapse_rate = 0,
    lapse_bias = NULL,
    metadata = list()
) {
  .new_family_cognitive_model(
    category_likelihood = category_likelihood,
    category_representation_class = Exemplar_CategoryRepresentation,
    model_class = Exemplar_Model,
    family_label = "EXEMPLAR",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
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

# -------------------------
# Core generics
# -------------------------

construct_mvbu <- S7::new_generic("construct_mvbu", "x")
validate_mvbu <- S7::new_generic("validate_mvbu", "x")
summarize_mvbu <- S7::new_generic("summarize_mvbu", "x")
print_mvbu <- S7::new_generic("print_mvbu", "x")
categorize_mvbu <- S7::new_generic("categorize_mvbu", c("x", "new_data"))
predict_mvbu <- S7::new_generic("predict_mvbu", c("x", "new_data"))
posterior_mvbu <- S7::new_generic("posterior_mvbu", "x")
plot_prep_mvbu <- S7::new_generic("plot_prep_mvbu", "x")

get_model_family <- S7::new_generic("get_model_family", "x")
get_category_likelihood <- S7::new_generic("get_category_likelihood", "x")
get_parameters <- S7::new_generic("get_parameters", "x")
get_category_prior <- S7::new_generic("get_category_prior", "x")
get_cue_labels <- S7::new_generic("get_cue_labels", "x")
get_category_labels <- S7::new_generic("get_category_labels", "x")
get_group_labels <- S7::new_generic("get_group_labels", "x")

get_categorization <- S7::new_generic("get_categorization", c("x", "new_data"))
get_category_prediction <- S7::new_generic("get_category_prediction", c("x", "new_data"))
get_posterior_prediction <- S7::new_generic("get_posterior_prediction", c("x", "new_data"))

update_model <- S7::new_generic("update_model", c("x", "data"))
update_category_likelihood <- S7::new_generic("update_category_likelihood", c("x", "data"))

get_posterior <- S7::new_generic("get_posterior", "x")
get_expected_category <- S7::new_generic("get_expected_category", "x")

plot_categories <- S7::new_generic("plot_categories", "x")
plot_parameters <- S7::new_generic("plot_parameters", "x")
plot_diagnostics <- S7::new_generic("plot_diagnostics", "x")

# -------------------------
# Minimal base methods
# -------------------------

S7::method(get_model_family, MVBU_Object) <- function(x) {
  class(x)[1]
}

S7::method(construct_mvbu, MVBU_Object) <- function(x) {
  x
}

S7::method(validate_mvbu, MVBU_Object) <- function(x) {
  validate_object(x)
}

S7::method(summarize_mvbu, MVBU_Object) <- function(x) {
  list(
    class = class(x)[1],
    model_family = get_model_family(x)
  )
}

S7::method(print_mvbu, MVBU_Object) <- function(x) {
  cat("<", class(x)[1], ">\n", sep = "")
  invisible(x)
}

S7::method(categorize_mvbu, list(MVBU_CognitiveModel, S7::class_any)) <- function(x, new_data) {
  get_categorization(x, new_data)
}

S7::method(predict_mvbu, list(MVBU_CognitiveModel, S7::class_any)) <- function(x, new_data) {
  get_category_prediction(x, new_data)
}

S7::method(posterior_mvbu, MVBU_Object) <- function(x) {
  get_posterior(x)
}

S7::method(plot_prep_mvbu, MVBU_Object) <- function(x) {
  .mvbu_not_implemented("plot_prep_mvbu", class(x)[1])
}

S7::method(get_model_family, MVBU_ModelDistribution) <- function(x) {
  x@model_family
}

S7::method(get_category_likelihood, MVBU_CognitiveModel) <- function(x) {
  x@category_likelihood
}

S7::method(get_parameters, MVBU_Object) <- function(x) {
  .mvbu_not_implemented("get_parameters", class(x)[1])
}

S7::method(get_parameters, UVG_CategoryRepresentation) <- function(x) {
  list(
    mu = x@mu,
    sigma2 = x@sigma2
  )
}

S7::method(get_parameters, NIX_CategoryRepresentation) <- function(x) {
  list(
    m = x@m,
    kappa = x@kappa,
    nu = x@nu,
    sigma2 = x@sigma2
  )
}

S7::method(get_parameters, MUVG_CategoryRepresentation) <- function(x) {
  list(
    component_mu = x@component_mu,
    component_sigma2 = x@component_sigma2,
    component_weights = x@component_weights
  )
}

S7::method(get_parameters, MNIX_CategoryRepresentation) <- function(x) {
  list(
    component_m = x@component_m,
    component_kappa = x@component_kappa,
    component_nu = x@component_nu,
    component_sigma2 = x@component_sigma2,
    component_weights = x@component_weights
  )
}

S7::method(get_parameters, MVG_CategoryRepresentation) <- function(x) {
  list(
    mu = x@mu,
    Sigma = x@Sigma
  )
}

S7::method(get_parameters, NIW_CategoryRepresentation) <- function(x) {
  list(
    m = x@m,
    kappa = x@kappa,
    nu = x@nu,
    S = x@S
  )
}

S7::method(get_parameters, Exemplar_CategoryRepresentation) <- function(x) {
  list(
    exemplars = x@exemplars,
    exemplar_weights = x@exemplar_weights
  )
}

S7::method(get_category_prior, MVBU_Object) <- function(x) {
  .mvbu_not_implemented("get_category_prior", class(x)[1])
}

S7::method(get_category_prior, MVBU_CognitiveModel) <- function(x) {
  x@category_prior
}

S7::method(get_cue_labels, MVBU_CategoryRepresentation) <- function(x) {
  x@cue_labels
}

S7::method(get_category_labels, MVBU_CategoryRepresentation) <- function(x) {
  x@category_labels
}

S7::method(get_group_labels, MVBU_Object) <- function(x) {
  character(0)
}

# NOTE: group labels currently identify model instances in combinations of
# models. Revisit later for richer grouped-model containers.
S7::method(get_group_labels, MVBU_ModelDistribution) <- function(x) {
  x@group_label
}
