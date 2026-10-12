#' @include internal-asserts.R
#' @importFrom S7 new_class new_generic method S7_inherits
NULL

# S7 foundation: base classes and core generic scaffolding.
# This file is intentionally minimal and non-breaking.

.get_cache <- function(x) {
  res <- tryCatch(x@cache, error = function(e) attr(x, "cache"))
  if (is.null(res)) list() else res
}

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
#' - `lapse_behavior`: list with `lapse_rate`, `lapse_bias`, and
#'   `lapse_treatment`.
#' - `noise_behavior`: list with `Sigma_noise` and `noise_treatment`.
#'   `Sigma_noise` is stored as a square matrix and may be supplied as a vector
#'   (which is converted to a diagonal matrix) or as a matrix.
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
# Family registration framework
# -------------------------

#' Model Family Class
#'
#' Encapsulates the mathematical primitives, moment extraction, and naming
#' conventions for a distributional model family in `MVBeliefUpdatr`.
#'
#' @param family_names A list containing `name_short`, `name_long`, and
#'   `name_details`.
#' @param density_function Function `f(x, rep, log = FALSE)` computing density.
#' @param sampling_function Function `r(rep, n = 1L, with_replacement = TRUE)`
#'   drawing sample observations from the category representation.
#' @param expected_moments Function `f(rep)` returning a list with `mu`, `Sigma`,
#'   and optional component `weights`.
#' @param marginal_moments Optional function `f(rep)` returning marginal
#'   predictive moments for adaptors.
#' @param sufficient_statistics_function Optional function extracting
#'   sufficient data statistics.
#'
#' @name MVBU_Family
#' @export
MVBU_Family <- S7::new_class(
  "MVBU_Family",
  properties = list(
    family_names = S7::class_list,
    density_function = S7::class_function,
    sampling_function = S7::class_function,
    expected_moments = S7::class_function,
    marginal_moments = S7::class_any,
    sufficient_statistics_function = S7::class_any
  ),
  validator = function(self) {
    if (is.null(self@family_names$name_short) ||
        !nzchar(self@family_names$name_short)) {
      return("family_names must contain a non-empty name_short character string.")
    }
    NULL
  }
)

#' Construct a new MVBU_Family object
#'
#' @param name_short Short uppercase identifier for the family (e.g., "UVG").
#' @param name_long Long descriptive name (defaults to `name_short`).
#' @param name_details Optional list specifying tier class/type names:
#'   `CategoryRepresentation_name`, `CategoryTemplate_name`,
#'   `CognitiveModel_name`, `CognitiveModelList_name`.
#' @param density_function Density function `function(x, rep, log = FALSE)`.
#' @param sampling_function Sampling function
#'   `function(rep, n = 1L, with_replacement = TRUE, ...)`.
#' @param expected_moments Function `function(rep)` returning
#'   `list(mu, Sigma, weights)`.
#' @param marginal_moments Optional function returning marginal predictive
#'   moments.
#' @param sufficient_statistics_function Optional function returning
#'   sufficient statistics.
#' @return An `MVBU_Family` object.
#' @export
new_mvbu_family <- function(
  name_short,
  name_long = NULL,
  name_details = NULL,
  density_function,
  sampling_function,
  expected_moments,
  marginal_moments = NULL,
  sufficient_statistics_function = NULL
) {
  .assert_character_scalar(name_short)
  name_short <- toupper(trimws(name_short))
  name_long <- if (!is.null(name_long)) {
    .assert_character_scalar(name_long)
    trimws(name_long)
  } else {
    name_short
  }

  defaults_details <- list(
    CategoryRepresentation_name = paste0(name_short, "_CategoryRepresentation"),
    CategoryTemplate_name = paste0(name_short, "_CategoryRepresentationTemplate"),
    CognitiveModel_name = paste0(name_short, "_CognitiveModel"),
    CognitiveModelList_name = paste0(name_short, "_CognitiveModelList")
  )
  if (is.null(name_details)) {
    name_details <- defaults_details
  } else {
    for (k in names(defaults_details)) {
      if (is.null(name_details[[k]])) {
        name_details[[k]] <- defaults_details[[k]]
      }
    }
  }

  MVBU_Family(
    family_names = list(
      name_short = name_short,
      name_long = name_long,
      name_details = name_details
    ),
    density_function = density_function,
    sampling_function = sampling_function,
    expected_moments = expected_moments,
    marginal_moments = marginal_moments,
    sufficient_statistics_function = sufficient_statistics_function
  )
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
        if (is.matrix(el)) {
          return(el)
        }
        if (is.data.frame(el)) {
          return(as.matrix(el))
        }
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

.logsumexp_rows <- function(matrix_log) {
  max_val <- apply(matrix_log, 1, max)
  max_val + log(rowSums(exp(matrix_log - max_val)))
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

.get_family_from_class_name <- function(class_name) {
  class_name <- sub("^.*::", "", class_name)
  families <- tryCatch(list_model_families(), error = function(e) character(0))
  for (fam in families) {
    reg <- tryCatch(.get_model_family_registration(fam), error = function(e) NULL)
    if (!is.null(reg) && class_name %in% c(reg$category_representation, reg$cognitive_model)) {
      return(fam)
    }
  }
  if (grepl("MNIX", class_name)) return("MNIX")
  if (grepl("NIW", class_name)) return("NIW")
  if (grepl("NIX", class_name)) return("NIX")
  if (grepl("MUVG", class_name)) return("MUVG")
  if (grepl("UVG", class_name)) return("UVG")
  if (grepl("MVG", class_name)) return("MVG")
  clean_name <- gsub(
    paste0(
      "_(IdealObserver|IdealAdaptor|IdealAdaptorStanfit|IdealAdaptorStaninput|",
      "CategoryRepresentation|CategoryRepresentationTemplate|Model|CognitiveModelList)$"
    ),
    "",
    class_name
  )
  if (toupper(clean_name) == "EXEMPLAR") {
    return("EXEMPLAR")
  }
  clean_name
}

# Unified likelihood closure generator from family object
.make_family_likelihood_function <- function(rep, fam = NULL) {
  if (is.null(fam)) {
    fam <- if (!is.null(rep@metadata$family) &&
               S7::S7_inherits(rep@metadata$family, MVBU_Family)) {
      rep@metadata$family
    } else {
      fam_name <- .get_family_from_class_name(class(rep)[1])
      .get_model_family(fam_name)
    }
  }
  fn <- fam@density_function
  fn_formals <- names(formals(fn))
  has_noise <- "noise_treatment" %in% fn_formals || "..." %in% fn_formals

  function(x,
           log = FALSE,
           noise_treatment = "no_noise",
           Sigma_noise = NULL) {
    if (has_noise) {
      fn(
        x = x,
        rep = rep,
        log = log,
        noise_treatment = noise_treatment,
        Sigma_noise = Sigma_noise
      )
    } else {
      if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
        x <- as.matrix(x) + mvtnorm::rmvnorm(
          n = nrow(as.matrix(x)),
          mean = rep(0, ncol(as.matrix(x))),
          sigma = Sigma_noise
        )
      }
      fn(x = x, rep = rep, log = log)
    }
  }
}

.mnix_cue_log_densities <- function(x, rep, noise_treatment = "no_noise", Sigma_noise = NULL) {
  k <- length(rep@m)
  x <- .as_observation_matrix(x, d = k, arg_name = "x")
  if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
    x <- x + mvtnorm::rmvnorm(
      n = nrow(x),
      mean = rep(0, ncol(x)),
      sigma = Sigma_noise
    )
  }
  noise_diag <- if (!is.null(Sigma_noise) &&
    (identical(noise_treatment, "sample") ||
      identical(noise_treatment, "marginalize"))) {
    if (is.matrix(Sigma_noise)) diag(Sigma_noise) else as.numeric(Sigma_noise)
  } else {
    rep(0, k)
  }
  kappa_vec <- if (length(rep@kappa) == 1L) rep(rep@kappa, k) else rep@kappa
  nu_vec <- if (length(rep@nu) == 1L) rep(rep@nu, k) else rep@nu
  d_mat <- vapply(seq_len(k), function(i) {
    scale <- sqrt((rep@sigma2[i] + noise_diag[i]) * (kappa_vec[i] + 1) / kappa_vec[i])
    stats::dt(
      (x[, i] - rep@m[i]) / scale,
      df = nu_vec[i],
      log = TRUE
    ) - log(scale)
  }, numeric(nrow(x)))
  matrix(d_mat, nrow = nrow(x), ncol = k)
}

.build_builtin_families <- function() {
  list(
    UVG = new_mvbu_family(
      name_short = "UVG",
      name_long = "Univariate Gaussian Ideal Observer",
      name_details = list(
        CategoryRepresentation_name = "UVG_CategoryRepresentation",
        CategoryTemplate_name = "UVG_CategoryRepresentationTemplate",
        CognitiveModel_name = "UVG_IdealObserver",
        CognitiveModelList_name = "UVG_CognitiveModelList"
      ),
      density_function = function(x, rep, log = FALSE, noise_treatment = "no_noise", Sigma_noise = NULL) {
        x <- .as_observation_matrix(x, d = 1, arg_name = "x")
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(
            n = nrow(x),
            mean = rep(0, ncol(x)),
            sigma = Sigma_noise
          )
        }
        noise_variance <- if (!is.null(Sigma_noise) &&
          (identical(noise_treatment, "sample") ||
            identical(noise_treatment, "marginalize"))) {
          as.numeric(Sigma_noise[1, 1])
        } else {
          0
        }
        stats::dnorm(
          x[, 1],
          mean = rep@mu,
          sd = sqrt(rep@sigma2 + noise_variance),
          log = log
        )
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        matrix(stats::rnorm(n, mean = rep@mu, sd = sqrt(rep@sigma2)), ncol = 1)
      },
      expected_moments = function(rep) {
        list(mu = rep@mu, Sigma = matrix(rep@sigma2, 1, 1), weights = 1)
      }
    ),
    NIX = new_mvbu_family(
      name_short = "NIX",
      name_long = "Normal-Inverse-Chi-Squared Ideal Adaptor",
      name_details = list(
        CategoryRepresentation_name = "NIX_CategoryRepresentation",
        CategoryTemplate_name = "NIX_CategoryRepresentationTemplate",
        CognitiveModel_name = "NIX_IdealAdaptor",
        CognitiveModelList_name = "NIX_CognitiveModelList"
      ),
      density_function = function(x, rep, log = FALSE, noise_treatment = "no_noise", Sigma_noise = NULL) {
        x <- .as_observation_matrix(x, d = 1, arg_name = "x")
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(
            n = nrow(x),
            mean = rep(0, ncol(x)),
            sigma = Sigma_noise
          )
        }
        noise_variance <- if (!is.null(Sigma_noise) &&
          (identical(noise_treatment, "sample") ||
            identical(noise_treatment, "marginalize"))) {
          as.numeric(Sigma_noise[1, 1])
        } else {
          0
        }
        scale <- sqrt((rep@sigma2 + noise_variance) * (rep@kappa + 1) / rep@kappa)
        d <- stats::dt((x[, 1] - rep@m) / scale, df = rep@nu, log = log)
        if (log) d - log(scale) else d / scale
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        scale <- sqrt(rep@sigma2 * (rep@kappa + 1) / rep@kappa)
        matrix(rep@m + scale * stats::rt(n, df = rep@nu), ncol = 1)
      },
      expected_moments = function(rep) {
        list(mu = rep@m, Sigma = matrix(get_expected_sigma(rep), 1, 1),
             weights = 1)
      },
      marginal_moments = function(rep) {
        list(mu = rep@m, Sigma = matrix(get_marginal_sigma(rep), 1, 1),
             weights = 1)
      }
    ),
    MUVG = new_mvbu_family(
      name_short = "MUVG",
      name_long = "Multivariate Independent Gaussian Ideal Observer",
      name_details = list(
        CategoryRepresentation_name = "MUVG_CategoryRepresentation",
        CategoryTemplate_name = "MUVG_CategoryRepresentationTemplate",
        CognitiveModel_name = "MUVG_IdealObserver",
        CognitiveModelList_name = "MUVG_CognitiveModelList"
      ),
      density_function = function(x, rep, log = FALSE, noise_treatment = "no_noise", Sigma_noise = NULL) {
        k <- length(rep@mu)
        x <- .as_observation_matrix(x, d = k, arg_name = "x")
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(
            n = nrow(x),
            mean = rep(0, ncol(x)),
            sigma = Sigma_noise
          )
        }
        noise_diag <- if (!is.null(Sigma_noise) &&
          (identical(noise_treatment, "sample") ||
            identical(noise_treatment, "marginalize"))) {
          diag(Sigma_noise)
        } else {
          rep(0, k)
        }
        d_mat <- vapply(seq_len(k), function(i) {
          stats::dnorm(
            x[, i],
            mean = rep@mu[i],
            sd = sqrt(rep@sigma2[i] + noise_diag[i]),
            log = TRUE
          )
        }, numeric(nrow(x)))
        tot <- rowSums(matrix(d_mat, nrow = nrow(x), ncol = k))
        if (log) tot else exp(tot)
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        k <- length(rep@mu)
        out <- vapply(seq_len(k), function(i) {
          stats::rnorm(n, mean = rep@mu[i], sd = sqrt(rep@sigma2[i]))
        }, numeric(n))
        matrix(out, nrow = n, ncol = k)
      },
      expected_moments = function(rep) {
        list(mu = rep@mu, Sigma = diag(rep@sigma2, nrow = length(rep@sigma2)),
             weights = rep@weights)
      }
    ),
    MNIX = new_mvbu_family(
      name_short = "MNIX",
      name_long = "Multivariate Independent Normal-Inverse-Chi-Squared Ideal Adaptor",
      name_details = list(
        CategoryRepresentation_name = "MNIX_CategoryRepresentation",
        CategoryTemplate_name = "MNIX_CategoryRepresentationTemplate",
        CognitiveModel_name = "MNIX_IdealAdaptor",
        CognitiveModelList_name = "MNIX_CognitiveModelList"
      ),
      density_function = function(x, rep, log = FALSE, noise_treatment = "no_noise", Sigma_noise = NULL) {
        if (length(rep@kappa) == 1L) {
          k <- length(rep@m)
          x <- .as_observation_matrix(x, d = k, arg_name = "x")
          if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
            x <- x + mvtnorm::rmvnorm(
              n = nrow(x),
              mean = rep(0, ncol(x)),
              sigma = Sigma_noise
            )
          }
          w <- rep@weights
          z <- as.numeric(x %*% w)
          m_z <- sum(w * rep@m)
          s2_z <- sum(w^2 * rep@sigma2)
          noise_var <- if (!is.null(Sigma_noise) &&
            (identical(noise_treatment, "sample") ||
              identical(noise_treatment, "marginalize"))) {
            if (is.matrix(Sigma_noise)) as.numeric(t(w) %*% Sigma_noise %*% w) else sum(w^2 * as.numeric(Sigma_noise))
          } else {
            0
          }
          scale <- sqrt((s2_z + noise_var) * (rep@kappa + 1) / rep@kappa)
          d <- stats::dt((z - m_z) / scale, df = rep@nu, log = TRUE) - log(scale)
          if (log) d else exp(d)
        } else {
          d_mat <- .mnix_cue_log_densities(x, rep, noise_treatment = noise_treatment, Sigma_noise = Sigma_noise)
          tot <- rowSums(d_mat)
          if (log) tot else exp(tot)
        }
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        k <- length(rep@m)
        kappa_vec <- if (length(rep@kappa) == 1L) rep(rep@kappa, k) else rep@kappa
        nu_vec <- if (length(rep@nu) == 1L) rep(rep@nu, k) else rep@nu
        out <- vapply(seq_len(k), function(i) {
          scale <- sqrt(rep@sigma2[i] * (kappa_vec[i] + 1) / kappa_vec[i])
          rep@m[i] + scale * stats::rt(n, df = nu_vec[i])
        }, numeric(n))
        matrix(out, nrow = n, ncol = k)
      },
      expected_moments = function(rep) {
        list(mu = rep@m, Sigma = get_expected_sigma(rep), weights = rep@weights)
      },
      marginal_moments = function(rep) {
        list(mu = rep@m, Sigma = get_marginal_sigma(rep), weights = rep@weights)
      }
    ),
    MVG = new_mvbu_family(
      name_short = "MVG",
      name_long = "Multivariate Gaussian Ideal Observer",
      name_details = list(
        CategoryRepresentation_name = "MVG_CategoryRepresentation",
        CategoryTemplate_name = "MVG_CategoryRepresentationTemplate",
        CognitiveModel_name = "MVG_IdealObserver",
        CognitiveModelList_name = "MVG_CognitiveModelList"
      ),
      density_function = function(x, rep, log = FALSE, noise_treatment = "no_noise", Sigma_noise = NULL) {
        d <- length(rep@mu)
        x <- .as_observation_matrix(x, d = d, arg_name = "x")
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(
            n = nrow(x),
            mean = rep(0, ncol(x)),
            sigma = Sigma_noise
          )
        }
        Sigma_eff <- if (!is.null(Sigma_noise) &&
            (identical(noise_treatment, "sample") ||
              identical(noise_treatment, "marginalize"))) {
          rep@Sigma + Sigma_noise
        } else {
          rep@Sigma
        }
        mvtnorm::dmvnorm(x, mean = rep@mu, sigma = Sigma_eff, log = log)
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        mvtnorm::rmvnorm(n, mean = rep@mu, sigma = rep@Sigma)
      },
      expected_moments = function(rep) {
        list(mu = rep@mu, Sigma = rep@Sigma, weights = NULL)
      }
    ),
    NIW = new_mvbu_family(
      name_short = "NIW",
      name_long = "Normal-Inverse-Wishart Ideal Adaptor",
      name_details = list(
        CategoryRepresentation_name = "NIW_CategoryRepresentation",
        CategoryTemplate_name = "NIW_CategoryRepresentationTemplate",
        CognitiveModel_name = "NIW_IdealAdaptor",
        CognitiveModelList_name = "NIW_CognitiveModelList"
      ),
      density_function = function(x, rep, log = FALSE, noise_treatment = "no_noise", Sigma_noise = NULL) {
        d <- length(rep@m)
        x <- .as_observation_matrix(x, d = d, arg_name = "x")
        df <- rep@nu - d + 1
        Sigma_pred <- ((rep@kappa + 1) / (rep@kappa * df)) * rep@S
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(
            n = nrow(x),
            mean = rep(0, ncol(x)),
            sigma = Sigma_noise
          )
        }
        Sigma_eff <- if (!is.null(Sigma_noise) &&
            (identical(noise_treatment, "sample") ||
              identical(noise_treatment, "marginalize"))) {
          Sigma_pred + Sigma_noise
        } else {
          Sigma_pred
        }
        .dmvt_density(x, mean = rep@m, Sigma = Sigma_eff, df = df, log = log)
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        d <- length(rep@m)
        df <- rep@nu - d + 1
        Sigma_pred <- ((rep@kappa + 1) / (rep@kappa * df)) * rep@S
        .rmvt(n = n, delta = rep@m, sigma = Sigma_pred, df = df)
      },
      expected_moments = function(rep) {
        list(mu = get_expected_mu_from_m(rep@m),
             Sigma = get_expected_Sigma_from_S(rep@S, rep@nu),
             weights = NULL)
      },
      marginal_moments = function(rep) {
        list(mu = get_marginal_mu(rep), Sigma = get_marginal_sigma(rep),
             weights = NULL)
      }
    ),
    EXEMPLAR = new_mvbu_family(
      name_short = "EXEMPLAR",
      name_long = "Exemplar Model",
      name_details = list(
        CategoryRepresentation_name = "Exemplar_CategoryRepresentation",
        CategoryTemplate_name = "Exemplar_CategoryRepresentationTemplate",
        CognitiveModel_name = "Exemplar_Model",
        CognitiveModelList_name = "Exemplar_CognitiveModelList"
      ),
      density_function = function(x, rep, log = FALSE, noise_treatment = "no_noise", Sigma_noise = NULL) {
        ex0 <- rep@exemplars
        d0 <- ncol(ex0)
        n0 <- nrow(ex0)
        w0 <- as.numeric(rep@exemplar_weights)
        bandwidth_scale <- if (n0 > 1) {
          (4 / ((d0 + 2) * n0))^(2 / (d0 + 4))
        } else {
          1.0
        }
        Sigma0 <- if (n0 > 1) {
          stats::cov(ex0) * bandwidth_scale
        } else {
          diag(1, d0) * bandwidth_scale
        }
        if (!is.matrix(Sigma0) || any(!is.finite(Sigma0))) {
          Sigma0 <- diag(1, d0) * bandwidth_scale
        }
        Sigma0 <- Sigma0 + diag(MVBU_PROB_TOL, d0)

        x <- .as_observation_matrix(x, d = d0, arg_name = "x")
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(
            n = nrow(x),
            mean = rep(0, ncol(x)),
            sigma = Sigma_noise
          )
        }
        Sigma_eff <- if (!is.null(Sigma_noise) &&
            (identical(noise_treatment, "sample") ||
              identical(noise_treatment, "marginalize"))) {
          Sigma0 + Sigma_noise
        } else {
          Sigma0
        }

        chol_sigma <- tryCatch(chol(Sigma_eff), error = function(e) NULL)
        if (is.null(chol_sigma)) {
          chol_sigma <- chol(Sigma_eff + diag(MVBU_PROB_TOL, d0))
        }
        log_det <- 2 * sum(log(diag(chol_sigma)))

        x_std <- t(backsolve(chol_sigma, t(x), transpose = TRUE))
        ex_std <- t(backsolve(chol_sigma, t(ex0), transpose = TRUE))

        sq_x <- rowSums(x_std^2)
        sq_ex <- rowSums(ex_std^2)
        dist_mat <- outer(sq_x, sq_ex, "+") - 2 * tcrossprod(x_std, ex_std)
        dist_mat[dist_mat < 0] <- 0

        log_dens_by_exemplar <- -0.5 * (d0 * log(2 * pi) + log_det + dist_mat)
        weighted_logdens <- sweep(log_dens_by_exemplar, 2, log(w0), "+")
        log_mix <- .logsumexp_rows(weighted_logdens)
        if (isTRUE(log)) log_mix else exp(log_mix)
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        n_avail <- nrow(rep@exemplars)
        if (isTRUE(with_replacement)) {
          idx <- sample.int(n_avail, size = n, replace = TRUE,
                            prob = rep@exemplar_weights)
        } else {
          if (n > n_avail) {
            cat_lbl <- get_category_labels(rep)
            .stop(sprintf(
              paste0("Cannot sample %d observations without replacement ",
                     "from category '%s' which contains only %d exemplars."),
              n, cat_lbl, n_avail
            ))
          }
          idx <- sample.int(n_avail, size = n, replace = FALSE,
                            prob = rep@exemplar_weights)
        }
        rep@exemplars[idx, , drop = FALSE]
      },
      expected_moments = function(rep) {
        list(
          mu = as.numeric(rep@exemplar_weights %*% rep@exemplars),
          Sigma = stats::cov.wt(rep@exemplars, wt = rep@exemplar_weights)$cov,
          weights = NULL
        )
      }
    )
  )
}

.mvbu_family_registry <- new.env(parent = emptyenv())
.mvbu_family_registry$families <- .build_builtin_families()

#' Normalize a model-family name to a canonical uppercase form
#'
#' @param family Model family name.
#' @return Uppercase character string representing the family name.
#' @keywords internal
#' @noRd
.normalize_family_name <- function(family) {
  .assert_character_scalar(family)
  toupper(trimws(family))
}

#' Register a model family in the MVBU family registry
#'
#' Registers a model family in the package registry, allowing either an
#' [MVBU_Family] object or class names.
#'
#' @param family Model family name or [MVBU_Family] object.
#' @param category_representation_class Name of category-representation class.
#' @param cognitive_model_class Name of cognitive-model class.
#' @param family_obj Optional [MVBU_Family] object.
#' @return Invisibly TRUE.
#' @export
register_model_family <- function(
  family,
  category_representation_class = NULL,
  cognitive_model_class = NULL,
  family_obj = NULL
) {
  if (S7::S7_inherits(family, MVBU_Family)) {
    family_obj <- family
    family <- family_obj@family_names$name_short
  } else if (!is.null(family_obj) && S7::S7_inherits(family_obj, MVBU_Family)) {
    family <- if (!missing(family)) .normalize_family_name(family) else family_obj@family_names$name_short
  } else {
    family <- .normalize_family_name(family)
    .assert_character_scalar(category_representation_class)
    .assert_character_scalar(cognitive_model_class)
    family_obj <- new_mvbu_family(
      name_short = family,
      name_details = list(
        CategoryRepresentation_name = category_representation_class,
        CognitiveModel_name = cognitive_model_class
      ),
      density_function = function(x, rep, log = FALSE) numeric(0),
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        matrix(numeric(0), n, 0)
      },
      expected_moments = function(rep) {
        list(mu = numeric(0), Sigma = matrix(numeric(0), 0, 0), weights = NULL)
      }
    )
  }
  .register_model_family(family, family_obj)
}

.register_model_family <- function(family, family_obj) {
  family <- .normalize_family_name(family)
  .mvbu_family_registry$families[[family]] <- family_obj
  invisible(TRUE)
}

.get_registered_model_families <- function() {
  sort(names(.mvbu_family_registry$families))
}

#' List all registered model families
#'
#' @return Character vector of registered family short names.
#' @export
list_model_families <- function() {
  .get_registered_model_families()
}

#' Get an MVBU_Family object from the registry
#'
#' @param family Model family short name.
#' @return The [MVBU_Family] object.
#' @export
get_model_family <- function(family) {
  .get_model_family(family)
}

.get_model_family <- function(family) {
  family <- .normalize_family_name(family)
  fam <- .mvbu_family_registry$families[[family]]
  if (is.null(fam)) {
    .stop(paste0("No model family registered for '", family, "'."))
  }
  fam
}

#' Get the registration information for a model family
#'
#' @param family Model family name.
#' @return A list with registered class names and the family object.
#' @keywords internal
#' @noRd
.get_model_family_registration <- function(family) {
  fam <- .get_model_family(family)
  list(
    category_representation = fam@family_names$name_details$CategoryRepresentation_name,
    cognitive_model = fam@family_names$name_details$CognitiveModel_name,
    names = fam@family_names,
    family = fam
  )
}

#' Register Stan-family extension hooks
#'
#' @param family Model family name.
#' @param stanfit_class Optional Stanfit class name.
#' @param staninput_class Optional Staninput class name.
#' @param bridge_methods Character vector of bridge methods.
#' @param dependency_rationale Character vector describing the dependency
#'   rationale.
#' @return Invisibly TRUE.
#' @keywords internal
#' @noRd
.register_stan_family_hooks <- function(
  family,
  stanfit_class = NULL,
  staninput_class = NULL,
  bridge_methods = character(),
  dependency_rationale = character()
) {
  family <- .normalize_family_name(family)

  .assert_true(
    is.null(stanfit_class) || .is_scalar_character(stanfit_class),
    msg = "stanfit_class must be NULL or a non-NA scalar character value."
  )
  .assert_true(
    is.null(staninput_class) || .is_scalar_character(staninput_class),
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
#' @noRd
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

#' Get registered model families in MVBU
#'
#' @return Character vector of registered family names.
#' @export
get_registered_model_families <- function() {
  list_model_families()
}


#' Get the registration information for a model family
#'
#' @param family Model-family name.
#' @return A list with the registered class names.
#' @export
get_model_family_registration <- function(family) {
  .get_model_family_registration(family)
}

#' Register Stan-family extension hooks
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

#' Get Stan-family extension hooks
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
      posterior_functions[[key]] <- function(new_data,
                                             categories = NULL,
                                             .noise_treatment = noise_treatment_i,
                                             .lapse_treatment = lapse_treatment_i) {
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
#' @param decision_rule Categorization decision rule: `"sampling"` or
#'   `"argmax"`. Defaults to `"sampling"`.
#' @param category_prior Optional numeric vector of prior category
#'   probabilities summing to 1. Defaults to equal priors.
#' @param lapse_rate Numeric scalar lapse probability in `[0, 1]`. Defaults
#'   to 0.
#' @param lapse_bias Optional numeric vector of lapse category probabilities
#'   summing to 1. Defaults to uniform bias across categories.
#' @param Sigma_noise Optional perceptual/measurement noise covariance matrix
#'   or scalar variance. Defaults to `NULL`.
#' @param noise_treatment Treatment of noise: `"no_noise"`, `"sample"`, or
#'   `"marginalize"`. Defaults to `.infer_noise_treatment(Sigma_noise)`:
#'   `"marginalize"` if noise is present, `"no_noise"` otherwise.
#' @param lapse_treatment Treatment of lapses: `"no_lapses"`, `"sample"`, or
#'   `"marginalize"`. Defaults to `.infer_lapse_treatment(lapse_rate)`:
#'   `"marginalize"` if `lapse_rate > 0`, `"no_lapses"` otherwise.
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
      group = if (!is.null(metadata$group)) as.character(metadata$group) else character(0)
    )
  )
  if (is.null(out@category_likelihood_function)) {
    fam <- tryCatch(get_family(out), error = function(e) NULL)
    if (!is.null(fam)) {
      out@category_likelihood_function <- .make_family_likelihood_function(out, fam)
    }
  }
  out
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
    lik_fun <- get_category_likelihood_function(representations[[i]])
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
        Sigma_noise = get_noise(model)
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

  lapse_rate <- get_lapse_rate(model)
  lapse_bias <- get_lapse_bias(model)

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
