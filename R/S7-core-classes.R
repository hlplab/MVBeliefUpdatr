#' @include asserts.R
#' @importFrom S7 new_class new_generic method S7_inherits
NULL

# S7 foundation: base classes and core generic scaffolding.
# This file is intentionally minimal and non-breaking.

# -------------------------
# Base class hierarchy
# -------------------------

#' MVBU Core S7 Class Architecture
#'
#' The MVBeliefUpdatr S7 organizes model objects around an explicit
#' compositional structure. These classes define the canonical S7 API; any
#' older list/data-frame-style inputs are supported only through transitional
#' compatibility methods and should be phased out over time.
#'
#' - [MVBU_CategoryRepresentation] defines the structure for one category-level
#'   representational object. This could be a uni- or multivariate Gaussian,
#'   a mixture of Gaussians, or an exemplar-based representation.
#' - [MVBU_CategoryRepresentationTemplate] collects one or more
#'   category-representation objects into a validated set used by a model.
#'   This corresponds to the notion of "templates" in e.g.,
#'   \cite{nearey-assmann07}.
#' - [MVBU_CognitiveModel] combines category-template structure with model-level
#'   decision behavior (`decision_rule`, `category_prior`,
#'   `lapse_rate`, `lapse_bias`).
#'
#' Family-level semantics are layered on top of these base classes.
#' Ideal observer/adaptor pairings follow:
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
#' - `category_likelihood_function`: family-specific likelihood function
#'   placeholder.
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
      return(
        "metadata$label_information$category must contain at least one label."
      )
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
      return(
        paste0(
          "representations must contain at least one ",
          "category representation object."
        )
      )
    }

    if (!all(vapply(
      self@representations,
      function(r) S7::S7_inherits(r, MVBU_CategoryRepresentation),
      logical(1)
    ))) {
      return(
        paste0(
          "all representations entries must inherit from ",
          "MVBU_CategoryRepresentation."
        )
      )
    }

    repr_names <- names(self@representations)
    if (!is.null(repr_names)) {
      if (length(repr_names) != length(self@representations) ||
          any(repr_names == "") || anyDuplicated(repr_names) > 0) {
        return(
          "if representations are named, names must be non-empty and unique."
        )
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
#' - `category_posterior_functions`: named list of treatment-specific
#'   posterior functions.
#' - `decision_rule`: scalar character rule label.
#' - `category_prior`: numeric probability vector over represented categories.
#' - `lapse_behavior`: list with `lapse_rate`, `lapse_bias`, and
#'   `lapse_treatment`.
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
    if (!is.list(lapse_behavior) ||
        !all(c("lapse_rate", "lapse_bias", "lapse_treatment") %in%
             names(lapse_behavior))) {
      return(
        paste0(
          "lapse_behavior must be a list containing lapse_rate, ",
          "lapse_bias, and lapse_treatment."
        )
      )
    }

    lapse_rate <- as.numeric(lapse_behavior$lapse_rate)
    lapse_bias <- as.numeric(lapse_behavior$lapse_bias)
    lapse_treatment <- as.character(lapse_behavior$lapse_treatment)
    if (length(lapse_rate) != 1 || lapse_rate < 0 || lapse_rate > 1) {
      return("lapse_rate must be a scalar numeric in [0, 1].")
    }

    if (length(lapse_bias) != n_repr) {
      return(
        "lapse_bias length must match the number of category representations."
      )
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
      return(
        paste0(
          "lapse_treatment must be one of 'no_lapses', 'sample', ",
          "or 'marginalize'."
        )
      )
    }

    noise_behavior <- self@noise_behavior
    if (!is.list(noise_behavior) ||
        !all(c("Sigma_noise", "noise_treatment") %in%
             names(noise_behavior))) {
      return(
        paste0(
          "noise_behavior must be a list containing Sigma_noise ",
          "and noise_treatment."
        )
      )
    }

    Sigma_noise <- noise_behavior$Sigma_noise
    noise_treatment <- as.character(noise_behavior$noise_treatment)
    if (!is.null(Sigma_noise)) {
      if (!is.matrix(Sigma_noise) || !is.numeric(Sigma_noise)) {
        return("Sigma_noise must be a numeric matrix when supplied.")
      }
      if (nrow(Sigma_noise) != ncol(Sigma_noise)) {
        return("Sigma_noise must be a square matrix.")
      }
      if (
        !isTRUE(
          all.equal(Sigma_noise, t(Sigma_noise), tolerance = MVBU_PROB_TOL)
        )
      ) {
        return("Sigma_noise must be symmetric.")
      }
      if (any(diag(Sigma_noise) < 0)) {
        return("Sigma_noise diagonal entries must be non-negative.")
      }
      if (any(Sigma_noise < 0)) {
        return("Sigma_noise entries must be non-negative.")
      }
      cue_labels <- get_cue_labels(self@category_template)
      if (nrow(Sigma_noise) != length(cue_labels) ||
          ncol(Sigma_noise) != length(cue_labels)) {
        return("Sigma_noise dimensions must match the number of cue labels.")
      }
    }

    if (!noise_treatment %in% c("no_noise", "sample", "marginalize")) {
      return(
        "noise_treatment must be one of 'no_noise', 'sample', or 'marginalize'."
      )
    }

    if (length(self@category_prior) != n_repr) {
      return(
        paste0(
          "category_prior length must match the number of ",
          "category representations."
        )
      )
    }

    if (any(self@category_prior < 0) || any(self@category_prior > 1)) {
      return("category_prior entries must be in [0, 1].")
    }
    if (abs(sum(self@category_prior) - 1) > MVBU_PROB_TOL) {
      return("category_prior entries must sum to 1.")
    }

    repr_names <- names(self@category_template@representations)
    if (!is.null(repr_names)) {
      if (length(repr_names) != n_repr || any(repr_names == "") ||
          anyDuplicated(repr_names) > 0) {
        return(
          "if representations are named, names must be non-empty and unique."
        )
      }

      prior_names <- names(self@category_prior)
      if (!is.null(prior_names)) {
        if (length(prior_names) != n_repr || any(prior_names == "") ||
            anyDuplicated(prior_names) > 0) {
          return(
            "if category_prior is named, names must be non-empty and unique."
          )
        }
        if (!setequal(prior_names, repr_names)) {
          return(
            paste0(
              "if category_prior is named, names must match ",
              "category_template names."
            )
          )
        }
      }

      lapse_names <- names(lapse_bias)
      if (!is.null(lapse_names)) {
        if (length(lapse_names) != n_repr || any(lapse_names == "") ||
            anyDuplicated(lapse_names) > 0) {
          return("if lapse_bias is named, names must be non-empty and unique.")
        }
        if (!setequal(lapse_names, repr_names)) {
          return(
            "if lapse_bias is named, names must match category_template names."
          )
        }
      }
    }

    NULL
  }
)

# -------------------------
# Family registration framework
# -------------------------

.mvbu_family_registry <- new.env(parent = emptyenv())
.mvbu_family_registry$families <- list(
  UVG = list(
    category_representation = "UVG_CategoryRepresentation",
    cognitive_model = "UVG_IdealObserver"
  ),
  NIX = list(
    category_representation = "NIX_CategoryRepresentation",
    cognitive_model = "NIX_IdealAdaptor"
  ),
  MUVG = list(
    category_representation = "MUVG_CategoryRepresentation",
    cognitive_model = "MUVG_IdealObserver"
  ),
  MNIX = list(
    category_representation = "MNIX_CategoryRepresentation",
    cognitive_model = "MNIX_IdealAdaptor"
  ),
  MVG = list(
    category_representation = "MVG_CategoryRepresentation",
    cognitive_model = "MVG_IdealObserver"
  ),
  NIW = list(
    category_representation = "NIW_CategoryRepresentation",
    cognitive_model = "NIW_IdealAdaptor"
  ),
  EXEMPLAR = list(
    category_representation = "Exemplar_CategoryRepresentation",
    cognitive_model = "Exemplar_Model"
  )
)

#' Normalize a model-family name to a canonical uppercase form.
#'
#' @param family Model family name.
#' @return Uppercase character string representing the family name.
#' @keywords internal
.normalize_family_name <- function(family) {
  .assert_non_NA_scalar_character(family)
  toupper(trimws(family))
}

#' Register a model family in the MVBU family registry.
#'
#' @param family Model family name.
#' @param category_representation_class Name of the category-representation
#'   class.
#' @param cognitive_model_class Name of the cognitive-model class.
#' @return Invisibly TRUE.
#' @keywords internal
.register_model_family <- function(
  family,
  category_representation_class,
  cognitive_model_class
) {
  family <- .normalize_family_name(family)

  .assert_non_NA_scalar_character(category_representation_class)
  .assert_non_NA_scalar_character(cognitive_model_class)

  .mvbu_family_registry$families[[family]] <- list(
    category_representation = category_representation_class,
    cognitive_model = cognitive_model_class
  )
  invisible(TRUE)
}

.get_registered_model_families <- function() {
  sort(names(.mvbu_family_registry$families))
}

#' Get the registration information for a model family.
#'
#' @param family Model family name.
#' @return A list with the registered class names.
#' @keywords internal
.get_model_family_registration <- function(family) {
  family <- .normalize_family_name(family)
  registration <- .mvbu_family_registry$families[[family]]
  if (is.null(registration)) {
    .stop(paste0("No model family registered for '", family, "'."))
  }
  registration
}

#' Register Stan-family extension hooks.
#'
#' @param family Model family name.
#' @param stanfit_class Optional Stanfit class name.
#' @param staninput_class Optional Staninput class name.
#' @param bridge_methods Character vector of bridge methods.
#' @param dependency_rationale Character vector describing the dependency
#'   rationale.
#' @return Invisibly TRUE.
#' @keywords internal
.register_stan_family_hooks <- function(
  family,
  stanfit_class = NULL,
  staninput_class = NULL,
  bridge_methods = character(),
  dependency_rationale = character()
) {
  family <- .normalize_family_name(family)

  .assert_true(
    is.null(stanfit_class) || .is_non_NA_scalar_character(stanfit_class),
    msg = "stanfit_class must be NULL or a non-NA scalar character value."
  )
  .assert_true(
    is.null(staninput_class) || .is_non_NA_scalar_character(staninput_class),
    msg = "staninput_class must be NULL or a non-NA scalar character value."
  )
  .assert_non_NA_character(bridge_methods)
  .assert_non_NA_character(dependency_rationale)

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

.get_stan_family_hooks <- function(family = NULL) {
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
.register_stan_family_hooks(
  family = "NIW",
  bridge_methods = c(
    "get_stanfit", "as_stanfit", "get_draws", "summary", "print", "loo",
    "posterior::as_draws_df"
  ),
  dependency_rationale = paste0(
    "Align with rstan/tidybayes workflows without adding dependencies ",
    "beyond demonstrated usage."
  )
)

#' Construct a base MVBU object
#'
#' @return A new instance of `MVBU_Object`.
#' @keywords internal
.new_mvbu_object <- function() {
  MVBU_Object()
}

#' Construct a base MVBU object
#'
#' @return A new instance of `MVBU_Object`.
#' @export
new_mvbu_object <- function() {
  .new_mvbu_object()
}

#' Register a model family in the MVBU family registry.
#'
#' @param family Model-family name.
#' @param category_representation_class Name of the category-representation
#'   class.
#' @param cognitive_model_class Name of the cognitive-model class.
#' @return Invisibly TRUE.
#' @export
register_model_family <- function(
  family,
  category_representation_class,
  cognitive_model_class
) {
  .register_model_family(
    family = family,
    category_representation_class = category_representation_class,
    cognitive_model_class = cognitive_model_class
  )
}

#' Get registered model families in MVBU.
#'
#' @return Character vector of registered family names.
#' @export
get_registered_model_families <- function() {
  .get_registered_model_families()
}

#' List registered model families in MVBU.
#'
#' @return Character vector of registered family names.
#' @export
list_model_families <- function() {
  .get_registered_model_families()
}

#' Get the registration information for a model family.
#'
#' @param family Model-family name.
#' @return A list with the registered class names.
#' @export
get_model_family_registration <- function(family) {
  .get_model_family_registration(family)
}

#' Register Stan-family extension hooks.
#'
#' @param family Model-family name.
#' @param stanfit_class Optional Stanfit class name.
#' @param staninput_class Optional Staninput class name.
#' @param bridge_methods Character vector of bridge methods.
#' @param dependency_rationale Character vector describing the dependency
#'   rationale.
#' @return Invisibly TRUE.
#' @export
register_stan_family_hooks <- function(
  family,
  stanfit_class = NULL,
  staninput_class = NULL,
  bridge_methods = character(),
  dependency_rationale = character()
) {
  .register_stan_family_hooks(
    family = family,
    stanfit_class = stanfit_class,
    staninput_class = staninput_class,
    bridge_methods = bridge_methods,
    dependency_rationale = dependency_rationale
  )
}

#' Get Stan-family extension hooks.
#'
#' @param family Optional model-family name. If NULL, returns all hooks.
#' @return A list with hook definitions.
#' @export
get_stan_family_hooks <- function(family = NULL) {
  .get_stan_family_hooks(family = family)
}

# -------------------------
# Internal Shared Model and Representation Construction Helpers
# -------------------------

#' Normalize perceptual noise covariance matrix.
#' @keywords internal
.mvbu_normalize_sigma_noise <- function(Sigma_noise, cue_labels) {
  if (is.null(Sigma_noise)) {
    return(NULL)
  }

  if (is.matrix(Sigma_noise)) {
    Sigma_noise <- as.matrix(Sigma_noise)
  } else if (is.numeric(Sigma_noise) && length(Sigma_noise) > 0 &&
             is.null(dim(Sigma_noise))) {
    .assert_true(
      length(Sigma_noise) == length(cue_labels),
      msg = "Sigma_noise length must match the number of cue labels."
    )
    Sigma_noise <- diag(
      Sigma_noise,
      nrow = length(Sigma_noise),
      ncol = length(Sigma_noise)
    )
  } else {
    .stop("Sigma_noise must be NULL, a numeric vector, or a matrix.")
  }

  .assert_true(is.numeric(Sigma_noise), msg = "Sigma_noise must be numeric.")
  .assert_true(
    all(is.finite(Sigma_noise)),
    msg = "Sigma_noise entries must be finite."
  )
  .assert_true(
    all(Sigma_noise >= 0),
    msg = "Sigma_noise entries must be non-negative."
  )
  .assert_true(
    nrow(Sigma_noise) == length(cue_labels) &&
      ncol(Sigma_noise) == length(cue_labels),
    msg = "Sigma_noise dimensions must match the number of cue labels."
  )

  Sigma_noise
}

.mvbu_cognitive_model_constructor <- function(
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
  .assert_true(
    !is.null(category_template),
    msg = "category_template must be supplied."
  )
  .assert_true(
    S7::S7_inherits(category_template, MVBU_CategoryRepresentationTemplate),
    msg = "category_template must be an MVBU_CategoryRepresentationTemplate."
  )

  n_repr <- length(category_template@representations)
  repr_names <- names(category_template@representations)

  .assert_non_NA_scalar_character(decision_rule)
  .assert_true(
    noise_treatment %in% c("no_noise", "sample", "marginalize"),
    msg = paste0(
      "noise_treatment must be one of 'no_noise', 'sample', ",
      "or 'marginalize'."
    )
  )
  .assert_true(
    lapse_treatment %in% c("no_lapses", "sample", "marginalize"),
    msg = paste0(
      "lapse_treatment must be one of 'no_lapses', 'sample', ",
      "or 'marginalize'."
    )
  )

  if (is.null(category_prior)) {
    category_prior <- rep(1 / n_repr, n_repr)
    if (!is.null(repr_names)) {
      names(category_prior) <- repr_names
    }
  }
  .assert_non_NA_numeric(category_prior)

  .assert_non_NA_scalar_numeric(lapse_rate)
  .assert_true(
    lapse_rate >= 0 && lapse_rate <= 1,
    msg = "lapse_rate must be a scalar in [0, 1]."
  )

  if (is.null(lapse_bias)) {
    lapse_bias <- rep(1 / n_repr, n_repr)
    if (!is.null(repr_names)) {
      names(lapse_bias) <- repr_names
    }
  }
  .assert_non_NA_numeric(lapse_bias)

  .mvbu_align_probability_vector <- function(values, target_names, arg_name) {
    if (is.null(values)) {
      return(values)
    }

    value_names <- names(values)
    if (is.null(value_names)) {
      if (!is.null(target_names) && length(target_names) == length(values)) {
        if (any(target_names == "") || anyDuplicated(target_names) > 0) {
          .stop(
            paste0(
              arg_name,
              " names cannot be validated because category_template ",
              "names are missing or invalid."
            )
          )
        }
        names(values) <- target_names
      }
      return(values)
    }

    if (length(value_names) != length(values) ||
        any(value_names == "") || anyDuplicated(value_names) > 0) {
      .stop(
        paste0(
          arg_name,
          " names must be non-empty and unique when provided."
        )
      )
    }

    if (is.null(target_names) || length(target_names) != length(values) ||
        any(target_names == "") || anyDuplicated(target_names) > 0) {
      .stop(
        paste0(
          arg_name,
          " names cannot be validated because category_template ",
          "names are missing or invalid."
        )
      )
    }

    if (!setequal(value_names, target_names)) {
      .stop(paste0(arg_name, " names must match category_template names."))
    }

    values <- values[target_names]
    values
  }

  category_prior <- .mvbu_align_probability_vector(
    category_prior, repr_names, "category_prior"
  )
  lapse_bias <- .mvbu_align_probability_vector(
    lapse_bias, repr_names, "lapse_bias"
  )

  if (length(category_prior) != n_repr) {
    .stop(
      paste0(
        "category_prior length must match the number of ",
        "category representations."
      )
    )
  }
  if (length(lapse_bias) != n_repr) {
    .stop(
      "lapse_bias length must match the number of category representations."
    )
  }

  cue_labels <- get_cue_labels(category_template)
  Sigma_noise <- .mvbu_normalize_sigma_noise(Sigma_noise, cue_labels)

  lapse_behavior <- list(
    lapse_rate = as.numeric(lapse_rate),
    lapse_bias = lapse_bias,
    lapse_treatment = as.character(lapse_treatment)
  )

  noise_behavior <- list(
    Sigma_noise = Sigma_noise,
    noise_treatment = as.character(noise_treatment)
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

  noise_treatments <- if (!is.null(model@noise_behavior$Sigma_noise)) {
    c("no_noise", "sample", "marginalize")
  } else {
    "no_noise"
  }
  lapse_treatments <- c("no_lapses", "sample", "marginalize")
  posterior_functions <- list()
  for (noise_treatment_i in noise_treatments) {
    for (lapse_treatment_i in lapse_treatments) {
      key <- paste(noise_treatment_i, lapse_treatment_i, sep = "__")
      posterior_functions[[key]] <- function(
        new_data,
        categories = NULL,
        .noise_treatment = noise_treatment_i,
        .lapse_treatment = lapse_treatment_i
      ) {
        .mvbu_posterior_matrix(
          model,
          new_data,
          categories = categories,
          noise_treatment = .noise_treatment,
          lapse_treatment = .lapse_treatment
        )
      }
    }
  }
  model@category_posterior_functions <- posterior_functions

  model
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
  .mvbu_cognitive_model_constructor(
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
}

.mvbu_label_information <- function(metadata = list()) {
  if (is.null(metadata)) {
    metadata <- list()
  }
  if (!is.list(metadata)) {
    metadata <- as.list(metadata)
  }

  if (!is.null(metadata$label_information) &&
      is.list(metadata$label_information)) {
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
  if (is.null(label_information$group) && !is.null(metadata$group)) {
    label_information$group <- metadata$group
  }

  label_information
}

.mvbu_label_metadata <- function(
  category_labels = character(),
  cue_labels = character(),
  metadata = list()
) {
  if (!is.list(metadata)) {
    metadata <- as.list(metadata)
  }

  metadata$label_information <- list(
    category = as.character(category_labels),
    cue = as.character(cue_labels),
    group = if (!is.null(metadata$group)) {
      as.character(metadata$group)
    } else {
      character(0)
    }
  )
  metadata$category <- NULL
  metadata$cue <- NULL
  metadata
}

.mvbu_extract_label_metadata <- function(x) {
  if (is.null(x)) {
    return(list(
      category = character(0),
      cue = character(0),
      group = character(0)
    ))
  }

  if (S7::S7_inherits(x, MVBU_CognitiveModel)) {
    label_info <- .mvbu_label_information(x@category_template@metadata)
  } else if (S7::S7_inherits(x, MVBU_Object) &&
             "metadata" %in% names(S7::props(x))) {
    label_info <- .mvbu_label_information(x@metadata)
  } else if (S7::S7_inherits(x, MVBU_Object) &&
             "labels" %in% names(S7::props(x))) {
    label_info <- x@labels
  } else if (is.list(x) && !is.null(x$metadata)) {
    label_info <- .mvbu_label_information(x$metadata)
  } else if (is.list(x) && !is.null(x$labels)) {
    label_info <- x$labels
  } else {
    return(list(
      category = character(0),
      cue = character(0),
      group = character(0)
    ))
  }

  if (!is.list(label_info)) {
    label_info <- list()
  }

  list(
    category = if (!is.null(label_info$category)) {
      as.character(label_info$category)
    } else {
      character(0)
    },
    cue = if (!is.null(label_info$cue)) {
      as.character(label_info$cue)
    } else {
      character(0)
    },
    group = if (!is.null(label_info$group)) {
      as.character(label_info$group)
    } else {
      character(0)
    }
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
      .stop(
        paste0(
          "cue labels must be consistent across all representations ",
          "in a template."
        )
      )
    }
  }

  invisible(NULL)
}

.mvbu_template_metadata <- function(representations, metadata = list()) {
  if (!is.list(metadata)) {
    metadata <- as.list(metadata)
  }

  category_labels <- unlist(
    lapply(representations, function(rep) {
      .mvbu_extract_label_metadata(rep)$category
    }),
    use.names = FALSE
  )
  cue_labels <- .mvbu_extract_label_metadata(representations[[1]])$cue

  metadata$label_information <- list(
    category = category_labels,
    cue = cue_labels,
    group = if (!is.null(metadata$group)) {
      as.character(metadata$group)
    } else {
      character(0)
    }
  )
  metadata$category <- NULL
  metadata$cue <- NULL
  metadata
}

add_category_representation <- function(
  template,
  representation,
  name = NULL
) {
  if (!S7::S7_inherits(template, MVBU_CategoryRepresentationTemplate)) {
    .stop("template must be an MVBU_CategoryRepresentationTemplate.")
  }
  if (!S7::S7_inherits(representation, MVBU_CategoryRepresentation)) {
    .stop("representation must be an MVBU_CategoryRepresentation.")
  }
  if (!is.null(name) && (length(name) != 1 || !nzchar(name))) {
    .stop("name must be a non-empty scalar character value.")
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
  if (length(template_labels$cue) > 0 &&
      length(rep_labels$cue) > 0 &&
      !identical(as.character(template_labels$cue),
                 as.character(rep_labels$cue))) {
    .stop(
      paste0(
        "cue labels must be consistent across all representations ",
        "in a template."
      )
    )
  }
  if (length(template_labels$cue) == 0 && length(rep_labels$cue) > 0) {
    metadata$label_information$cue <- rep_labels$cue
  } else if (length(template_labels$cue) > 0) {
    metadata$label_information$cue <- template_labels$cue
  }

  metadata$label_information$category <- c(
    template_labels$category,
    rep_labels$category
  )
  metadata$category <- NULL
  metadata$cue <- NULL
  MVBU_CategoryRepresentationTemplate(
    representations = representations,
    metadata = metadata
  )
}

new_category_representation_template <- function(
  representations,
  metadata = list()
) {
  if (!is.list(representations) || length(representations) < 1) {
    .stop(
      paste0(
        "representations must contain at least one ",
        "category representation object."
      )
    )
  }
  if (!all(vapply(
    representations,
    function(r) S7::S7_inherits(r, MVBU_CategoryRepresentation),
    logical(1)
  ))) {
    .stop(
      paste0(
        "all representations entries must inherit from ",
        "MVBU_CategoryRepresentation."
      )
    )
  }

  rep_names <- names(representations)
  if (!is.null(rep_names) &&
      (length(rep_names) != length(representations) ||
       any(rep_names == "") || anyDuplicated(rep_names) > 0)) {
    .stop("if representations are named, names must be non-empty and unique.")
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
        metadata = .mvbu_template_metadata(
          template_representations,
          metadata
        )
      )
    } else {
      template <- add_category_representation(
        template,
        representations[[i]],
        name = rep_name
      )
    }
  }

  template
}

#' Construct a base representation object
new_category_representation <- function(
  category_labels,
  cue_labels,
  category_likelihood_function = NULL,
  metadata = list()
) {
  if (is.null(category_likelihood_function)) {
    category_likelihood_function <- function(...) {
      .stop("category_likelihood not implemented.")
    }
  }

  MVBU_CategoryRepresentation(
    category_likelihood_function = category_likelihood_function,
    metadata = .mvbu_label_metadata(
      as.character(category_labels),
      as.character(cue_labels),
      metadata
    )
  )
}

.logsumexp_rows <- function(matrix_log) {
  max_val <- apply(matrix_log, 1, max)
  max_val + log(rowSums(exp(matrix_log - max_val)))
}

.as_observation_matrix <- function(x, d = NULL, arg_name = "x") {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (is.list(x) && !is.data.frame(x)) {
    if (length(x) == 0) {
      return(matrix(numeric(0), nrow = 0, ncol = if (!is.null(d)) d else 0))
    }
    valid_list <- vapply(x, function(el) {
      if (is.matrix(el) || is.data.frame(el)) {
        m <- as.matrix(el)
        return(is.numeric(m) && (is.null(d) || ncol(m) == d))
      }
      if (is.numeric(el) && is.vector(el)) {
        return(is.null(d) || length(el) == d)
      }
      FALSE
    }, logical(1))

    if (all(valid_list)) {
      mats <- lapply(x, function(el) {
        if (is.matrix(el)) return(el)
        if (is.data.frame(el)) return(as.matrix(el))
        ncol_val <- if (!is.null(d)) d else length(el)
        return(matrix(as.numeric(el), ncol = ncol_val, byrow = TRUE))
      })
      return(do.call(rbind, mats))
    }
  }

  if (is.matrix(x)) {
    if (!is.numeric(x)) {
      .stop(arg_name, " must be numeric.")
    }
    if (!is.null(d) && ncol(x) != d) {
      .stop(arg_name, " must have ", d, " column(s).")
    }
    return(x)
  }

  if (is.atomic(x) && is.numeric(x)) {
    if (is.null(d) || d == 1) {
      return(matrix(as.numeric(x), ncol = 1))
    }
    if (length(x) %% d != 0) {
      .stop(arg_name, " length must be a multiple of ", d, ".")
    }
    return(matrix(as.numeric(x), ncol = d, byrow = TRUE))
  }

  .stop(
    arg_name,
    paste0(
      " must be a numeric vector, matrix, data frame, ",
      "or list of numeric vectors."
    )
  )
}

.dmvnorm_density <- function(x, mean, Sigma, log = FALSE) {
  x <- .as_observation_matrix(x, d = length(mean), arg_name = "x")
  d <- ncol(x)
  diff <- sweep(x, 2, mean, "-")
  chol_sigma <- tryCatch(chol(Sigma), error = function(e) NULL)

  if (is.null(chol_sigma)) {
    ridge <- diag(MVBU_PROB_TOL, d)
    chol_sigma <- chol(Sigma + ridge)
  }

  log_det <- 2 * sum(log(diag(chol_sigma)))
  half_sol <- backsolve(chol_sigma, t(diff), transpose = TRUE)
  quad <- colSums(half_sol^2)
  log_d <- -0.5 * (d * log(2 * pi) + log_det + quad)

  if (isTRUE(log)) {
    log_d
  } else {
    exp(log_d)
  }
}

.dmvt_density <- function(x, mean, Sigma, df, log = FALSE) {
  x <- .as_observation_matrix(x, d = length(mean), arg_name = "x")
  d <- ncol(x)
  diff <- sweep(x, 2, mean, "-")
  chol_sigma <- tryCatch(chol(Sigma), error = function(e) NULL)

  if (is.null(chol_sigma)) {
    ridge <- diag(MVBU_PROB_TOL, d)
    chol_sigma <- chol(Sigma + ridge)
  }

  log_det <- 2 * sum(log(diag(chol_sigma)))
  half_sol <- backsolve(chol_sigma, t(diff), transpose = TRUE)
  quad <- colSums(half_sol^2)

  log_const <- lgamma((df + d) / 2) - lgamma(df / 2) -
    0.5 * (d * log(df * pi) + log_det)
  log_d <- log_const - 0.5 * (df + d) * log1p(quad / df)

  if (isTRUE(log)) {
    log_d
  } else {
    exp(log_d)
  }
}

.mvbu_posterior_matrix <- function(
  model,
  new_data,
  categories = NULL,
  noise_treatment = "no_noise",
  lapse_treatment = "no_lapses"
) {
  representations <- model@category_template@representations
  rep_names <- names(representations)
  n_repr <- length(representations)

  new_data_mat <- .as_observation_matrix(new_data, arg_name = "new_data")
  n_obs <- nrow(new_data_mat)

  log_lik <- matrix(NA_real_, nrow = n_obs, ncol = n_repr)
  for (i in seq_len(n_repr)) {
    lik_fun <- representations[[i]]@category_likelihood_function
    if (is.null(lik_fun)) {
      .stop(
        "category representation ", i,
        " is missing category_likelihood_function."
      )
    }
    log_lik[, i] <- as.numeric(
      lik_fun(
        new_data_mat,
        log = TRUE,
        noise_treatment = noise_treatment,
        Sigma_noise = model@noise_behavior$Sigma_noise
      )
    )
  }

  prior <- model@category_prior
  log_prior <- log(prior)
  log_joint <- sweep(log_lik, 2, log_prior, "+")

  decision_rule <- model@decision_rule
  post_probs <- matrix(0, nrow = n_obs, ncol = n_repr)

  if (identical(decision_rule, "sampling")) {
    log_norm <- .logsumexp_rows(log_joint)
    post_probs <- exp(log_joint - log_norm)
  } else if (identical(decision_rule, "argmax")) {
    max_idx <- max.col(log_joint, ties.method = "random")
    for (row_i in seq_len(n_obs)) {
      post_probs[row_i, max_idx[row_i]] <- 1
    }
  } else {
    .stop("Unsupported decision_rule: ", decision_rule)
  }

  lapse_rate <- model@lapse_behavior$lapse_rate
  lapse_bias <- model@lapse_behavior$lapse_bias

  if (identical(lapse_treatment, "marginalize")) {
    if (lapse_rate > 0) {
      post_probs <- (1 - lapse_rate) * post_probs +
        matrix(
          lapse_rate * lapse_bias,
          nrow = n_obs,
          ncol = n_repr,
          byrow = TRUE
        )
    }
  } else if (identical(lapse_treatment, "sample")) {
    if (lapse_rate > 0) {
      is_lapse <- stats::rbinom(n = n_obs, size = 1, prob = lapse_rate) == 1
      if (any(is_lapse)) {
        post_probs[is_lapse, ] <- matrix(
          lapse_bias,
          nrow = sum(is_lapse),
          ncol = n_repr,
          byrow = TRUE
        )
      }
    }
  }

  colnames(post_probs) <- if (!is.null(repr_names)) {
    repr_names
  } else {
    paste0("cat_", seq_len(n_repr))
  }

  if (!is.null(categories)) {
    categories <- as.character(categories)
    missing_cats <- setdiff(categories, colnames(post_probs))
    if (length(missing_cats) > 0) {
      .stop(
        "Requested categories not found in model representations: ",
        paste(missing_cats, collapse = ", ")
      )
    }
    post_probs <- post_probs[, categories, drop = FALSE]
  }

  post_probs
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
    .stop("category_template must be supplied.")
  }

  if (!S7::S7_inherits(
    category_template,
    MVBU_CategoryRepresentationTemplate
  )) {
    .stop("category_template must be an MVBU_CategoryRepresentationTemplate.")
  }

  family_representations <- category_template@representations

  if (!all(vapply(
    family_representations,
    function(r) S7::S7_inherits(r, category_representation_class),
    logical(1)
  ))) {
    .stop(
      paste0(
        "All category_likelihood_template entries must inherit from ",
        family_label, " category representation class."
      )
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

#' Validate a model object
#' @keywords internal
validate_object <- function(x) {
  if (!S7::S7_inherits(x, MVBU_Object)) {
    .stop("x must be an S7 object.")
  }
  TRUE
}

# Safe validation check.
.is_valid <- function(x) {
  tryCatch(
    {
      validate_object(x)
      TRUE
    },
    error = function(e) FALSE
  )
}

.mvbu_not_implemented <- function(generic_name, class_name) {
  .stop(
    paste0(
      "Method for ", generic_name, "() not yet implemented for class ",
      class_name, "."
    )
  )
}
