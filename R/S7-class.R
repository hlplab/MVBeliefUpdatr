#' @include internal-asserts.R
#' @include S7-family.R
#' @importFrom S7 new_class new_generic method S7_inherits
NULL

# S7 foundation: base classes and core generic scaffolding.
# This file is intentionally minimal and non-breaking.

# Retrieve cached computations stored on an object.
#
# Arguments:
#   x: An S7 or S3 object.
#
# Returns:
#   The cached list or an empty list if not set.
.get_cache <- function(x) {
  res <- tryCatch(x@cache, error = function(e) attr(x, "cache"))
  if (is.null(res)) list() else res
}

# Store cached computations on an object.
#
# Arguments:
#   x: An S7 or S3 object.
#   val: A list or value to store in the cache.
#
# Returns:
#   The modified object with cache set.
.set_cache <- function(x, val) {
  tryCatch(
    {
      x@cache <- val
      x
    },
    error = function(e) {
      attr(x, "cache") <- val
      x
    }
  )
}

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
#'   \insertCite{nearey-assmann2007}{MVBeliefUpdatr}.
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
#' @references
#' \insertAllCited{}
#'
#' @name MVBU-core-classes
#' @docType class
#' @keywords internal
NULL

#' @section MVBU_Object:
#' Root S7 base class for MVBeliefUpdatr objects.
#'
#' @rdname MVBU-core-classes
#' @usage NULL
#' @export
MVBU_Object <- S7::new_class("MVBU_Object", package = NULL)

#' Construct a base MVBU object
#'
#' @return A new instance of `MVBU_Object`.
#' @export
#' @examples
#' obj <- new_mvbu_object()
#' inherits(obj, "MVBU_Object")
new_mvbu_object <- function() {
  MVBU_Object()
}

#' @section MVBU_CategoryRepresentation:
#' Abstract category-level representational class.
#'
#' Expected properties:
#' - `category_labels`: one or more category labels associated with the object.
#' - `cue_labels`: one or more cue-dimension labels.
#' - `category_likelihood_function`: family-specific likelihood function
# Core Model Family Class Hierarchy Design Rationale:
# We retain explicit typed classes (e.g., UVG_CategoryRepresentation, NIW_IdealAdaptor)
# for the built-in core families for two primary reasons:
# 1. Interactive tab completion in IDEs (RStudio, VSCode, Positron) for named properties.
# 2. Strict S7 slot typing and early validation at construction time, preventing invalid
#    dimensional or numerical structures prior to downstream computation or Stan fitting.
# Despite retaining typed subclasses for convenience, both built-in core classes and
# newly defined/custom classes follow the exact same underlying family contract:
# their category likelihood functions and computational moments are obtained directly
# from their registered MVBU_Family objects.

#' @rdname MVBU-core-classes
#' @usage NULL
#' @export
MVBU_CategoryRepresentation <- S7::new_class(
  "MVBU_CategoryRepresentation",
  package = NULL,
  parent = MVBU_Object,
  properties = list(
    category_likelihood_function = S7::new_property(
      class = S7::new_union(NULL, S7::class_function),
      default = quote(NULL)
    ),
    metadata = S7::new_property(
      class = S7::class_list,
      default = quote(list(parameters = list())),
      setter = function(self, value) {
        val <- as.list(value)
        if (is.null(val$parameters)) {
          val$parameters <- list()
        }
        self@metadata <- val
        self
      }
    )
  ),
  validator = function(self) {
    if (!is.list(self@metadata)) {
      return("metadata must be a list.")
    }

    if (is.null(self@metadata$label_information) ||
      !is.list(self@metadata$label_information)) {
      return("metadata$label_information must be a list.")
    }

    label_information <- self@metadata$label_information
    if (is.null(label_information$category) ||
      length(label_information$category) < 1) {
      return(
        "metadata$label_information$category must contain at least one label."
      )
    }
    if (is.null(label_information$cue) ||
      length(label_information$cue) < 1) {
      return("metadata$label_information$cue must contain at least one label.")
    }
    NULL
  }
)

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
#'
#' @rdname MVBU-core-classes
#' @usage NULL
#' @export
MVBU_CategoryRepresentationTemplate <- S7::new_class(
  "MVBU_CategoryRepresentationTemplate",
  package = NULL,
  parent = MVBU_Object,
  properties = list(
    representations = S7::class_list,
    metadata = S7::new_property(
      class = S7::class_list,
      default = quote(list(parameters = list())),
      setter = function(self, value) {
        val <- as.list(value)
        if (is.null(val$parameters)) {
          reps <- tryCatch(self@representations, error = function(e) list())
          if (length(reps) > 0) {
            val$parameters <- lapply(reps, function(r) {
              tryCatch(get_parameters(r), error = function(e) list())
            })
            if (!is.null(names(reps))) {
              names(val$parameters) <- names(reps)
            }
          } else {
            val$parameters <- list()
          }
        }
        self@metadata <- val
        self
      }
    )
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
#' - `lapse_behavior`: list bundling attentional lapse parameters:
#'   `lapse_rate` (scalar probability in `[0, 1]`), `lapse_bias` (vector of
#'   lapse response probabilities across categories), and `lapse_treatment`
#'   (`"no_lapses"`, `"sample"`, or `"marginalize"`).
#' - `noise_behavior`: list bundling perceptual noise parameters:
#'   `Sigma_noise` (perceptual noise covariance matrix) and `noise_treatment`
#'   (`"no_noise"`, `"sample"`, or `"marginalize"`). `Sigma_noise` is stored
#'   as a square matrix and may be supplied as a vector (converted to a
#'   diagonal matrix) or as a matrix.
#' - `metadata`: optional model metadata.
#'
#' @references
#' \insertRef{xie2023}{MVBeliefUpdatr}
#'
#' Xie, X., Jaeger, T. F., & Kurumada, C. (2021). What we do (not) know about
#' the mechanisms underlying adaptive speech perception: A computational review.
#'
#' @rdname MVBU-core-classes
#' @usage NULL
#' @export
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
    metadata = S7::new_property(
      class = S7::class_list,
      default = quote(list(parameters = list())),
      setter = function(self, value) {
        val <- as.list(value)
        if (is.null(val$parameters)) {
          tmpl <- tryCatch(self@category_template, error = function(e) NULL)
          if (!is.null(tmpl)) {
            val$parameters <- tryCatch(
              get_parameters(tmpl),
              error = function(e) list()
            )
          } else {
            val$parameters <- list()
          }
        }
        self@metadata <- val
        self
      }
    ),
    cache = S7::class_list
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
# Internal Shared Model and Representation Construction Helpers
# -------------------------

#' Normalize perceptual noise covariance matrix
#' @keywords internal
#' @noRd
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
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
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

  if (is.null(noise_treatment)) {
    noise_treatment <- .infer_noise_treatment(Sigma_noise)
  }
  if (is.null(lapse_treatment)) {
    lapse_treatment <- .infer_lapse_treatment(lapse_rate)
  }

  n_repr <- length(category_template@representations)
  repr_names <- names(category_template@representations)

  .assert_character_scalar(decision_rule)
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

  .assert_numeric_scalar(lapse_rate)
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

  metadata <- set_labels(
    metadata,
    cue = cue_labels,
    category = get_category_labels(category_template),
    response_category = get_category_labels(category_template),
    group = get_group_labels(category_template)
  )
  if (!is.null(category_template)) {
    metadata$parameters <- get_parameters(category_template)
  }

  model <- MVBU_CognitiveModel(
    category_template = category_template,
    category_posterior_functions = list(),
    decision_rule = as.character(decision_rule),
    category_prior = category_prior,
    lapse_behavior = lapse_behavior,
    noise_behavior = noise_behavior,
    metadata = metadata
  )

  noise_treatments <- if (!is.null(get_noise(model))) {
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

#' Create a cognitive model
#'
#' Constructs a cognitive model binding a category representation template to
#' model-level decision and uncertainty parameters.
#'
#' @param category_template An [MVBU_CategoryRepresentationTemplate] object
#'   containing category representations.
#' @param decision_rule Categorization decision rule: `"sampling"`,
#'   `"criterion"`, or `"proportional"`. Defaults to `"sampling"`.
#' @param category_prior Optional numeric vector of prior category
#'   probabilities summing to 1. Defaults to equal priors.
#' @param lapse_rate Numeric scalar lapse probability in `[0, 1]`. Defaults
#'   to 0.
#' @param lapse_bias Optional numeric vector of lapse category probabilities
#'   summing to 1. Defaults to uniform bias across categories.
#' @param Sigma_noise Optional perceptual/measurement noise covariance matrix
#'   or scalar variance. Defaults to `NULL`.
#' @param noise_treatment Treatment of perceptual noise during evaluation and
#'   updating: `"no_noise"`, `"sample"`, or `"marginalize"`.
#'   \itemize{
#'     \item \code{"no_noise"}: No perceptual noise is applied to input cues,
#'       and no noise is assumed during categorization.
#'     \item \code{"sample"}: Simulates sensory noise on a specific observation.
#'       Noise drawn from \eqn{\mathcal{N}(0, \Sigma_{\text{noise}})} is added
#'       to input cues (and projected into composite cue space in early
#'       integration), and the observer accounts for their own noise
#'       uncertainty during categorization by adding \eqn{\Sigma_{\text{noise}}}
#'       to their likelihood scale.
#'     \item \code{"marginalize"}: Evaluates categorization in the limit across
#'       infinitely many noisy encounters (expected consequences). No noise is
#'       sampled onto the input stimulus, but \eqn{\Sigma_{\text{noise}}} is
#'       added to the category covariance scale.
#'   }
#'   Defaults to `.infer_noise_treatment(Sigma_noise)`: `"marginalize"` if
#'   noise is present, `"no_noise"` otherwise.
#' @param lapse_treatment Treatment of attentional lapses during evaluation and
#'   updating: `"no_lapses"`, `"sample"`, or `"marginalize"`.
#'   \itemize{
#'     \item \code{"no_lapses"}: Lapses are ignored; categorization and
#'       updating reflect pure normative evidence evaluation (even if
#'       \code{lapse_rate > 0}).
#'     \item \code{"sample"}: Simulates stochastic attentional lapses on
#'       individual observations. With probability \code{lapse_rate}, a lapse
#'       occurs: during categorization, the posterior is replaced by
#'       \code{lapse_bias}; during updating, the lapsing observation is treated
#'       as uninformative.
#'     \item \code{"marginalize"}: Evaluates categorization in the limit across
#'       infinitely many encounters (expected response probabilities). The
#'       category posterior is the mixture
#'       \eqn{(1 - \lambda) p(c \mid x) + \lambda b_c}, and evidence during
#'       updating is downweighted by \eqn{1 - \lambda}.
#'   }
#'   Defaults to `.infer_lapse_treatment(lapse_rate)`: `"marginalize"` if
#'   `lapse_rate > 0`, `"no_lapses"` otherwise.
#' @param metadata Optional list of metadata. Defaults to `list()`.
#'
#' @return An [MVBU_CognitiveModel] object.
#' @rdname MVBU-core-classes
#' @export
new_cognitive_model <- function(
  category_template = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
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

.mvbu_validate_cue_consistency <- function(representations) {
  if (length(representations) < 2) {
    return(invisible(NULL))
  }

  reference_cues <- get_cue_labels(representations[[1]])
  for (i in seq_along(representations)[-1]) {
    rep_cues <- get_cue_labels(representations[[i]])
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

#' Construct an MVBU Category Representation Template
#'
#' @param representations A named or unnamed list of
#'   \code{MVBU_CategoryRepresentation} objects.
#' @param metadata Optional metadata list for template-level properties.
#'   Defaults to \code{list()}.
#'
#' @return An object of class \code{MVBU_CategoryRepresentationTemplate}.
#' @export
#' @examples
#' rep_b <- new_uvg_category_representation(
#'   category_labels = "/b/", cue_labels = "VOT",
#'   mu = 10, sigma2 = 25
#' )
#' rep_d <- new_uvg_category_representation(
#'   category_labels = "/d/", cue_labels = "VOT",
#'   mu = 50, sigma2 = 25
#' )
#' tmpl <- new_category_representation_template(list(rep_b, rep_d))
#' tmpl
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
  if (is.null(rep_names)) {
    inferred <- vapply(representations, function(r) {
      lbls <- get_category_labels(r)
      if (length(lbls) == 1L && is.character(lbls) && nzchar(lbls)) lbls else ""
    }, character(1))
    if (all(nzchar(inferred)) && anyDuplicated(inferred) == 0) {
      rep_names <- inferred
      names(representations) <- rep_names
    }
  }
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
      first_rep <- template_representations[[1]]
      metadata <- set_labels(
        metadata,
        cue = get_cue_labels(first_rep),
        category = get_category_labels(first_rep),
        response_category = get_category_labels(first_rep),
        group = if (!is.null(metadata$group)) as.character(metadata$group) else character(0)
      )
      template_params <- lapply(template_representations, get_parameters)
      if (!is.null(rep_name)) {
        names(template_params) <- rep_name
      }
      metadata$parameters <- template_params
      template <- MVBU_CategoryRepresentationTemplate(
        representations = template_representations,
        metadata = metadata
      )
    } else {
      template <- add_category_representation(
        template,
        representations[[i]],
        name = rep_name
      )
    }
  }

  template_params <- lapply(template@representations, get_parameters)
  if (!is.null(names(template@representations))) {
    names(template_params) <- names(template@representations)
  }
  template@metadata$parameters <- template_params

  template
}

#' Construct a base category representation object
#'
#' @param category_labels Character vector of category label(s).
#' @param cue_labels Character vector of cue label(s).
#' @param category_likelihood_function Optional category likelihood function.
#'   Default: \code{NULL}.
#' @param parameters Optional named list of representation parameters.
#'   Default: \code{list()}.
#' @param metadata Optional list of metadata. Default: \code{list()}.
#' @return An [MVBU_CategoryRepresentation] object.
#' @keywords internal
#' @export
new_category_representation <- function(
  category_labels,
  cue_labels,
  category_likelihood_function = NULL,
  parameters = list(),
  metadata = list()
) {
  metadata <- as.list(metadata)
  if (is.null(metadata$parameters)) {
    metadata$parameters <- as.list(parameters)
  }
  out <- MVBU_CategoryRepresentation(
    category_likelihood_function = category_likelihood_function,
    metadata = set_labels(
      metadata,
      cue = as.character(cue_labels),
      category = as.character(category_labels),
      response_category = as.character(category_labels),
      group = if (!is.null(metadata$group)) {
        as.character(metadata$group)
      } else {
        character(0)
      }
    )
  )
  if (is.null(out@category_likelihood_function)) {
    fam <- tryCatch(get_family(out), error = function(e) NULL)
    if (!is.null(fam)) {
      out@category_likelihood_function <- .make_family_likelihood_function(
        out,
        fam
      )
    }
  }
  out
}


# Internal constructor helper for family cognitive models (NIX, MNIX, NIW).
#
# Arguments:
#   category_template: Category representation template object.
#   category_representation_class: Expected S7 class of representations.
#   model_class: Target S7 cognitive model class.
#   family_label: Character label of the distribution family (e.g., "NIX").
#   decision_rule: Decision rule strategy ("sampling", "criterion", or
#     "proportional").
#   category_prior: Prior category probability vector.
#   lapse_rate: Probability of a lapse response.
#   lapse_bias: Optional response category bias vector during lapses.
#   Sigma_noise: Perceptual noise covariance matrix.
#   noise_treatment: Treatment mode for perceptual noise.
#   lapse_treatment: Treatment mode for lapses.
#   metadata: List of metadata attributes.
#   ...: Additional arguments forwarded to the model constructor.
#
# Returns:
#   A new instance of `model_class`.
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
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
  metadata = list(),
  ...
) {
  if (is.null(noise_treatment)) {
    noise_treatment <- .infer_noise_treatment(Sigma_noise)
  }
  if (is.null(lapse_treatment)) {
    lapse_treatment <- .infer_lapse_treatment(lapse_rate)
  }
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

  extra_args <- list(...)
  model_args <- c(
    list(
      category_template = base_model@category_template,
      decision_rule = base_model@decision_rule,
      category_prior = base_model@category_prior,
      lapse_behavior = base_model@lapse_behavior,
      noise_behavior = base_model@noise_behavior,
      metadata = base_model@metadata
    ),
    extra_args
  )
  if (!is.null(extra_args$cues)) {
    model_args$metadata$label_information$cue <- as.character(extra_args$cues)
  } else if (!is.null(extra_args$input_cues)) {
    model_args$metadata$label_information$cue <- as.character(extra_args$input_cues)
  }
  model <- do.call(model_class, model_args)

  model@category_posterior_functions <- base_model@category_posterior_functions

  model
}

#' Validate a model object
#' @keywords internal
#' @noRd
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

# Signal that a generic method is not implemented for a given class.
#
# Arguments:
#   generic_name: Name of the generic function.
#   class_name: Name of the class lacking the implementation.
.mvbu_not_implemented <- function(generic_name, class_name) {
  .stop(
    paste0(
      "Method for ", generic_name, "() not yet implemented for class ",
      class_name, "."
    )
  )
}
