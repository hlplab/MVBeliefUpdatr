#' @include S7-class.R
NULL

# -------------------------
# Multi-Cue Normal-Inverse-Chi-Squared (MNIX) Family
# -------------------------

#' Multi-Cue Normal-Inverse-Chi-Squared (MNIX) Family:
#' Representations, Templates, and Models
#'
#' Constructs and manages multi-cue Normal-Inverse-\eqn{\chi^2} (MNIX) category
#' representations, templates, and ideal adaptors under early and late cue integration.
#'
#' @param category_labels Character vector of category label(s).
#' @param cue_labels Character vector of cue labels (at least 2 cue dimensions).
#' @param m Numeric vector of per-cue location hyperparameters.
#' @param kappa Positive numeric scalar or vector of precision scalings (\eqn{\kappa > 0}).
#'   For late integration, a vector of length equal to cues (or scalar recycled across cues).
#'   For early integration, a scalar reflecting composite dimension belief strength.
#' @param nu Positive numeric scalar or vector of degrees of freedom (\eqn{\nu > 2}).
#'   For late integration, a vector of length equal to cues (or scalar recycled across cues).
#'   For early integration, a scalar reflecting composite dimension degrees of freedom.
#' @param sigma2 Numeric vector of per-cue scale parameters.
#' @param weights Numeric vector or matrix of integration/attention weights.
#'   For representations and global weighted models, a numeric simplex vector summing to 1.
#'   For category-weighted models, an \eqn{M \times K} matrix whose rows each sum to 1.
#' @param weights_matrix Deprecated. Use `weights` instead.
#' @param data Data frame containing category and cue observations.
#' @param category Character string indicating the category column in `data`. Defaults to `"category"`.
#' @param cues Character vector indicating cue columns in `data` (at least 2 cues).
#' @param input_cues Deprecated. Use `cues` instead.
#' @param category_template An [MVBU_CategoryRepresentationTemplate] object
#'   containing MNIX category representations.
#' @param decision_rule Categorization decision rule: `"sampling"` or
#'   `"argmax"`. Defaults to `"sampling"`.
#' @param category_prior Optional numeric vector of prior category
#'   probabilities summing to 1. Defaults to equal priors.
#' @param lapse_rate Numeric scalar lapse probability in `[0, 1]`.
#'   Defaults to 0.
#' @param lapse_bias Optional numeric vector of lapse category probabilities
#'   summing to 1. Defaults to uniform bias `1/M`.
#' @param Sigma_noise Optional perceptual/measurement noise covariance matrix
#'   or vector.
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
#' @param ... Additional arguments passed to methods or constructors.
#'
#' @details
#' The MNIX family unifies multi-cue belief updating under both late and early
#' integration regimes:
#' \itemize{
#'   \item \strong{Late Integration}: Categories are represented by independent marginal
#'     Student-\eqn{t} predictive components along each cue dimension (\eqn{\kappa} and \eqn{\nu}
#'     have length \eqn{K}). Cue log-likelihoods are combined at decision time via Bayes rule
#'     or attention weights.
#'   \item \strong{Early Integration}: Multi-cue observations are projected to a composite
#'     dimension \eqn{z = \mathbf{w}^\top \mathbf{x}} prior to updating (\eqn{\kappa} and \eqn{\nu}
#'     are scalar length 1).
#' }
#'
#' @return An S7 object of the respective MNIX class.
#' @seealso [family-muvg], [family-niw], [family-nix]
#'
#' @name family-mnix
#' @rdname family-mnix
#' @docType class
#' @usage NULL
#' @export
MNIX_CategoryRepresentation <- S7::new_class(
  "MNIX_CategoryRepresentation",
  parent = MVBU_CategoryRepresentation,
  properties = list(
    m = S7::class_numeric,
    kappa = S7::class_numeric,
    nu = S7::class_numeric,
    sigma2 = S7::class_numeric,
    weights = S7::class_numeric
  ),
  validator = function(self) {
    n_comp <- length(self@m)

    if (n_comp < 2L) {
      return("MNIX representation must contain at least two cue components.")
    }
    if (length(self@sigma2) != n_comp || length(self@weights) != n_comp) {
      return("MNIX sigma2 and weights must have length equal to number of cue components.")
    }
    if ((length(self@kappa) != n_comp && length(self@kappa) != 1L) ||
        (length(self@nu) != n_comp && length(self@nu) != 1L)) {
      return("MNIX kappa and nu must have length equal to cues (late) or length 1 (early).")
    }
    label_information <- self@metadata$label_information
    if (length(label_information$cue) != n_comp) {
      return("MNIX cue_labels length must match number of cue components.")
    }
    if (any(self@kappa <= 0)) {
      return("MNIX kappa entries must be > 0.")
    }
    if (any(self@nu <= 0)) {
      return("MNIX nu entries must be > 0.")
    }
    if (any(self@sigma2 <= 0)) {
      return("MNIX sigma2 entries must be > 0.")
    }
    if (any(self@weights < 0) || any(self@weights > 1)) {
      return("MNIX weights entries must be in [0, 1].")
    }
    if (abs(sum(self@weights) - 1) > MVBU_PROB_TOL) {
      return("MNIX weights entries must sum to 1.")
    }

    NULL
  }
)

#' @rdname family-mnix
#' @export
new_mnix_category_representation <- function(
  category_labels,
  cue_labels,
  m,
  kappa,
  nu,
  sigma2,
  weights = NULL,
  metadata = list()
) {
  n_comp <- length(m)
  if (n_comp < 2L) {
    .stop("m must contain at least two cue elements.")
  }
  if (!is.numeric(m) || !is.numeric(kappa) ||
    !is.numeric(nu) || !is.numeric(sigma2)) {
    .stop("m, kappa, nu, and sigma2 must be numeric.")
  }
  if (length(kappa) != 1L && length(kappa) != n_comp) {
    .stop("kappa must be a positive scalar (early integration) or length-K vector (late integration).")
  }
  if (length(nu) != 1L && length(nu) != n_comp) {
    .stop("nu must be a positive scalar (early integration) or length-K vector (late integration).")
  }
  if (length(sigma2) != n_comp) {
    .stop("sigma2 must have length equal to length(m).")
  }
  if (is.null(weights)) {
    weights <- rep(1 / n_comp, n_comp)
  }
  if (!is.numeric(weights) || length(weights) != n_comp) {
    .stop("weights must be a numeric vector of length length(m).")
  }

  metadata <- as.list(metadata)
  if (is.null(metadata$parameters)) {
    metadata$parameters <- list(
      m = as.numeric(m),
      kappa = as.numeric(kappa),
      nu = as.numeric(nu),
      sigma2 = as.numeric(sigma2),
      weights = as.numeric(weights)
    )
  }
  rep <- MNIX_CategoryRepresentation(
    category_likelihood_function = NULL,
    metadata = set_labels(
      metadata,
      cue = as.character(cue_labels),
      category = as.character(category_labels),
      response_category = as.character(category_labels),
      group = if (!is.null(metadata$group)) as.character(metadata$group) else character(0)
    ),
    m = as.numeric(m),
    kappa = as.numeric(kappa),
    nu = as.numeric(nu),
    sigma2 = as.numeric(sigma2),
    weights = as.numeric(weights)
  )
  rep@category_likelihood_function <- .make_family_likelihood_function(rep)
  rep
}

#' @rdname family-mnix
#' @export
new_mnix_category_representation_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = 3,
  ...
) {
  d <- length(cues)
  if (length(kappa) == 1L) {
    kappa <- rep(kappa, d)
  }
  if (length(nu) == 1L) {
    nu <- rep(nu, d)
  }

  .assert_true(
    length(kappa) == d && is.numeric(kappa) &&
      all(!is.na(kappa)) && all(kappa > 0),
    msg = paste0(
      "kappa must be a positive numeric scalar or vector of ",
      "length length(cues)."
    )
  )
  .assert_true(
    length(nu) == d && is.numeric(nu) &&
      all(!is.na(nu)) && all(nu > 2),
    msg = paste0(
      "nu must be a numeric scalar or vector of length ",
      "length(cues) with all values > 2."
    )
  )

  muvg <- new_muvg_category_representation_from_data(
    data,
    category = category,
    cues = cues
  )
  as_mnix_category_representation(
    muvg,
    kappa = kappa,
    nu = nu
  )
}

#' @rdname family-mnix
#' @export
new_mnix_category_representation_template_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = 3,
  ...
) {
  category_labels <- sort(unique(as.character(data[[category]])))
  representations <- lapply(category_labels, function(label) {
    new_mnix_category_representation_from_data(
      data[data[[category]] == label, , drop = FALSE],
      category = category,
      cues = cues,
      kappa = kappa,
      nu = nu,
      ...
    )
  })
  names(representations) <- category_labels
  new_category_representation_template(representations)
}

# -------------------------
# Base Model Hierarchy
# -------------------------

#' @rdname family-mnix
#' @docType class
#' @usage NULL
#' @export
MNIX_IdealAdaptor <- S7::new_class(
  "MNIX_IdealAdaptor",
  parent = MVBU_CognitiveModel
)

#' @rdname family-mnix
#' @docType class
#' @usage NULL
#' @export
MNIX_LateIdealAdaptor <- S7::new_class(
  "MNIX_LateIdealAdaptor",
  parent = MNIX_IdealAdaptor
)

#' @rdname family-mnix
#' @docType class
#' @usage NULL
#' @export
MNIX_EarlyIdealAdaptor <- S7::new_class(
  "MNIX_EarlyIdealAdaptor",
  parent = MNIX_IdealAdaptor,
  properties = list(
    cues = S7::class_character
  ),
  validator = function(self) {
    if (length(self@cues) < 2L) {
      return("Early integrator cues must contain at least 2 cue dimensions.")
    }
    NULL
  }
)

# -------------------------
# Late Integration Models (Suffix: MNIX_IdealAdaptor_{Integrator})
# -------------------------

#' @rdname family-mnix
#' @docType class
#' @usage NULL
#' @export
MNIX_IdealAdaptor_IdealIntegrator <- S7::new_class(
  "MNIX_IdealAdaptor_IdealIntegrator",
  parent = MNIX_LateIdealAdaptor
)

#' @rdname family-mnix
#' @docType class
#' @usage NULL
#' @export
MNIX_IdealAdaptor_WeightedIntegrator <- S7::new_class(
  "MNIX_IdealAdaptor_WeightedIntegrator",
  parent = MNIX_LateIdealAdaptor,
  properties = list(
    weights = S7::class_numeric
  ),
  validator = function(self) {
    w <- self@weights
    if (length(w) < 2L) {
      return("MNIX weighted integrator weights must contain at least 2 cue elements.")
    }
    if (any(w < 0) || any(w > 1)) {
      return("MNIX weighted integrator weights must be in [0, 1].")
    }
    if (abs(sum(w) - 1) > MVBU_PROB_TOL) {
      return("MNIX weighted integrator weights must sum to 1.")
    }
    NULL
  }
)

#' @rdname family-mnix
#' @docType class
#' @usage NULL
#' @export
MNIX_IdealAdaptor_CategoryWeightedIntegrator <- S7::new_class(
  "MNIX_IdealAdaptor_CategoryWeightedIntegrator",
  parent = MNIX_LateIdealAdaptor,
  properties = list(
    weights = S7::class_any
  ),
  validator = function(self) {
    wm <- self@weights
    if (!is.matrix(wm) && !is.data.frame(wm)) {
      return("weights must be a matrix or data frame.")
    }
    wm <- as.matrix(wm)
    if (ncol(wm) < 2L) {
      return("weights must contain at least 2 cue columns.")
    }
    if (any(wm < 0) || any(wm > 1)) {
      return("weights entries must be in [0, 1].")
    }
    rsums <- rowSums(wm)
    if (any(abs(rsums - 1) > MVBU_PROB_TOL)) {
      return("weights rows must each sum to 1.")
    }
    NULL
  }
)


#' @rdname family-mnix
#' @docType class
#' @usage NULL
#' @export
MNIX_IdealAdaptor_InformativityWeightedIntegrator <- S7::new_class(
  "MNIX_IdealAdaptor_InformativityWeightedIntegrator",
  parent = MNIX_LateIdealAdaptor
)

# -------------------------
# Early Integration Models (Prefix: {Integrator}_MNIX_IdealAdaptor)
# -------------------------

#' @rdname family-mnix
#' @docType class
#' @usage NULL
#' @export
IdealIntegrator_MNIX_IdealAdaptor <- S7::new_class(
  "IdealIntegrator_MNIX_IdealAdaptor",
  parent = MNIX_EarlyIdealAdaptor
)

#' @rdname family-mnix
#' @docType class
#' @usage NULL
#' @export
WeightedIntegrator_MNIX_IdealAdaptor <- S7::new_class(
  "WeightedIntegrator_MNIX_IdealAdaptor",
  parent = MNIX_EarlyIdealAdaptor,
  properties = list(
    weights = S7::class_numeric
  ),
  validator = function(self) {
    w <- self@weights
    if (length(w) < 2L) {
      return("weights must contain at least 2 cue elements.")
    }
    if (length(self@cues) > 0L && length(w) != length(self@cues)) {
      return("weights length must match cues length.")
    }
    if (any(w < 0) || any(w > 1)) {
      return("weights must be in [0, 1].")
    }
    if (abs(sum(w) - 1) > MVBU_PROB_TOL) {
      return("weights must sum to 1.")
    }
    NULL
  }
)

#' @rdname family-mnix
#' @docType class
#' @usage NULL
#' @export
CategoryWeightedIntegrator_MNIX_IdealAdaptor <- S7::new_class(
  "CategoryWeightedIntegrator_MNIX_IdealAdaptor",
  parent = MNIX_EarlyIdealAdaptor,
  properties = list(
    weights = S7::class_any
  ),
  validator = function(self) {
    wm <- self@weights
    if (!is.matrix(wm) && !is.data.frame(wm)) {
      return("weights must be a matrix or data frame.")
    }
    wm <- as.matrix(wm)
    if (ncol(wm) < 2L) {
      return("weights must contain at least 2 cue columns.")
    }
    if (length(self@cues) > 0L && ncol(wm) != length(self@cues)) {
      return("weights ncol must match cues length.")
    }
    if (any(wm < 0) || any(wm > 1)) {
      return("weights entries must be in [0, 1].")
    }
    rsums <- rowSums(wm)
    if (any(abs(rsums - 1) > MVBU_PROB_TOL)) {
      return("weights rows must each sum to 1.")
    }
    NULL
  }
)


#' @rdname family-mnix
#' @docType class
#' @usage NULL
#' @export
InformativityWeightedIntegrator_MNIX_IdealAdaptor <- S7::new_class(
  "InformativityWeightedIntegrator_MNIX_IdealAdaptor",
  parent = MNIX_EarlyIdealAdaptor
)

# -------------------------
# Late Model Constructors
# -------------------------

#' @rdname family-mnix
#' @export
new_mnix_ideal_adaptor_ideal_integrator <- function(
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
  .new_family_cognitive_model(
    category_template = category_template,
    category_representation_class = MNIX_CategoryRepresentation,
    model_class = MNIX_IdealAdaptor_IdealIntegrator,
    family_label = "MNIX",
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

#' @rdname family-mnix
#' @export
new_mnix_ideal_adaptor_weighted_integrator <- function(
  category_template = NULL,
  weights = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
  metadata = list()
) {
  if (is.null(category_template)) {
    .stop("category_template must be supplied.")
  }
  first_rep <- category_template@representations[[1L]]
  n_cues <- length(get_cue_labels(first_rep))
  if (is.null(weights)) {
    weights <- rep(1 / n_cues, n_cues)
  }
  weights <- as.numeric(weights)
  if (length(weights) != n_cues) {
    .stop("weights length (", length(weights), ") must match number of cues (", n_cues, ").")
  }

  .new_family_cognitive_model(
    category_template = category_template,
    category_representation_class = MNIX_CategoryRepresentation,
    model_class = MNIX_IdealAdaptor_WeightedIntegrator,
    family_label = "MNIX",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment,
    metadata = metadata,
    weights = weights
  )
}

#' @rdname family-mnix
#' @export
new_mnix_ideal_adaptor_category_weighted_integrator <- function(
  category_template = NULL,
  weights = NULL,
  weights_matrix = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
  metadata = list()
) {
  if (!is.null(weights_matrix)) {
    lifecycle::deprecate_warn("0.1.0", "new_mnix_ideal_adaptor_category_weighted_integrator(weights_matrix)", "new_mnix_ideal_adaptor_category_weighted_integrator(weights)")
    if (is.null(weights)) weights <- weights_matrix
  }
  if (is.null(category_template)) {
    .stop("category_template must be supplied.")
  }
  reps <- category_template@representations
  n_cats <- length(reps)
  n_cues <- length(get_cue_labels(reps[[1L]]))
  if (is.null(weights)) {
    weights <- matrix(1 / n_cues, nrow = n_cats, ncol = n_cues)
    rownames(weights) <- names(reps)
  }
  wm <- as.matrix(weights)
  if (nrow(wm) != n_cats || ncol(wm) != n_cues) {
    .stop(
      "weights dimensions (", nrow(wm), "x", ncol(wm),
      ") must match categories (", n_cats, ") and cues (", n_cues, ")."
    )
  }

  .new_family_cognitive_model(
    category_template = category_template,
    category_representation_class = MNIX_CategoryRepresentation,
    model_class = MNIX_IdealAdaptor_CategoryWeightedIntegrator,
    family_label = "MNIX",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment,
    metadata = metadata,
    weights = wm
  )
}


#' @rdname family-mnix
#' @export
new_mnix_ideal_adaptor_informativity_weighted_integrator <- function(
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
  .new_family_cognitive_model(
    category_template = category_template,
    category_representation_class = MNIX_CategoryRepresentation,
    model_class = MNIX_IdealAdaptor_InformativityWeightedIntegrator,
    family_label = "MNIX",
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

#' @rdname family-mnix
#' @export
new_mnix_ideal_adaptor <- function(
  category_template = NULL,
  weights = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
  metadata = list()
) {
  if (!is.null(weights)) {
    return(new_mnix_ideal_adaptor_weighted_integrator(
      category_template = category_template,
      weights = weights,
      decision_rule = decision_rule,
      category_prior = category_prior,
      lapse_rate = lapse_rate,
      lapse_bias = lapse_bias,
      Sigma_noise = Sigma_noise,
      noise_treatment = noise_treatment,
      lapse_treatment = lapse_treatment,
      metadata = metadata
    ))
  }
  new_mnix_ideal_adaptor_ideal_integrator(
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

#' @rdname family-mnix
#' @export
new_mnix_ideal_adaptor_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = 3,
  weights = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
  ...
) {
  template <- new_mnix_category_representation_template_from_data(
    data,
    category = category,
    cues = cues,
    kappa = kappa,
    nu = nu,
    ...
  )
  new_mnix_ideal_adaptor(
    category_template = template,
    weights = weights,
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment
  )
}

# -------------------------
# Early Model Constructors
# -------------------------

#' @rdname family-mnix
#' @export
new_ideal_integrator_mnix_ideal_adaptor <- function(
  category_template = NULL,
  cues = NULL,
  input_cues = NULL,
  weights = NULL,
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
  if (!is.null(input_cues)) {
    lifecycle::deprecate_warn("0.1.0", "new_ideal_integrator_mnix_ideal_adaptor(input_cues)", "new_ideal_integrator_mnix_ideal_adaptor(cues)")
    if (is.null(cues)) cues <- input_cues
  }
  if (!is.null(weights)) {
    lifecycle::deprecate_warn("0.1.0", "new_ideal_integrator_mnix_ideal_adaptor(weights)", "new_weighted_integrator_mnix_ideal_adaptor(weights)")
    return(new_weighted_integrator_mnix_ideal_adaptor(
      category_template = category_template,
      weights = weights,
      cues = cues,
      decision_rule = decision_rule,
      category_prior = category_prior,
      lapse_rate = lapse_rate,
      lapse_bias = lapse_bias,
      Sigma_noise = Sigma_noise,
      noise_treatment = noise_treatment,
      lapse_treatment = lapse_treatment,
      metadata = metadata,
      ...
    ))
  }
  if (is.null(category_template)) {
    .stop("category_template must be supplied.")
  }
  if (is.null(cues)) {
    first_rep <- category_template@representations[[1L]]
    cues <- get_cue_labels(first_rep)
  }
  cues <- as.character(cues)

  .new_family_cognitive_model(
    category_template = category_template,
    category_representation_class = MNIX_CategoryRepresentation,
    model_class = IdealIntegrator_MNIX_IdealAdaptor,
    family_label = "MNIX",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment,
    metadata = metadata,
    cues = cues
  )
}

#' @rdname family-mnix
#' @export
new_ideal_integrator_mnix_ideal_adaptor_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = 3,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate)
) {
  .assert_true(
    length(cues) >= 2L,
    msg = "Ideal integrator requires at least two cues."
  )
  cue_vars <- vapply(cues, function(col) stats::var(data[[col]], na.rm = TRUE), numeric(1L))
  if (any(cue_vars <= 0 | is.na(cue_vars))) {
    .stop("Cue variances must be positive and non-NA to compute precision weights.")
  }
  precisions <- 1 / cue_vars
  weights <- as.numeric(precisions / sum(precisions))
  names(weights) <- cues

  z_proj <- as.numeric(as.matrix(data[, cues, drop = FALSE]) %*% weights)
  category_labels <- sort(unique(as.character(data[[category]])))
  representations <- lapply(category_labels, function(label) {
    sub_data <- data[data[[category]] == label, , drop = FALSE]
    sub_z <- z_proj[data[[category]] == label]
    m_val <- vapply(cues, function(c_col) mean(sub_data[[c_col]], na.rm = TRUE), numeric(1L))
    s2_val <- vapply(cues, function(c_col) stats::var(sub_data[[c_col]], na.rm = TRUE), numeric(1L))
    new_mnix_category_representation(
      category_labels = label,
      cue_labels = cues,
      m = m_val,
      kappa = if (length(kappa) == 1L) kappa else kappa[1L],
      nu = if (length(nu) == 1L) nu else nu[1L],
      sigma2 = s2_val,
      weights = weights
    )
  })
  names(representations) <- category_labels
  template <- new_category_representation_template(representations)

  new_ideal_integrator_mnix_ideal_adaptor(
    category_template = template,
    cues = cues,
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment
  )
}

#' @rdname family-mnix
#' @export
new_weighted_integrator_mnix_ideal_adaptor <- function(
  category_template = NULL,
  weights = NULL,
  cues = NULL,
  input_cues = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
  metadata = list()
) {
  if (!is.null(input_cues)) {
    lifecycle::deprecate_warn("0.1.0", "new_weighted_integrator_mnix_ideal_adaptor(input_cues)", "new_weighted_integrator_mnix_ideal_adaptor(cues)")
    if (is.null(cues)) cues <- input_cues
  }
  if (is.null(category_template)) {
    .stop("category_template must be supplied.")
  }
  first_rep <- category_template@representations[[1L]]
  if (is.null(cues)) {
    cues <- get_cue_labels(first_rep)
  }
  cues <- as.character(cues)
  n_cues <- length(cues)
  if (is.null(weights)) {
    weights <- rep(1 / n_cues, n_cues)
  }
  weights <- as.numeric(weights)
  if (length(weights) != n_cues) {
    .stop("weights length (", length(weights), ") must match number of cues (", n_cues, ").")
  }

  .new_family_cognitive_model(
    category_template = category_template,
    category_representation_class = MNIX_CategoryRepresentation,
    model_class = WeightedIntegrator_MNIX_IdealAdaptor,
    family_label = "MNIX",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment,
    metadata = metadata,
    weights = weights,
    cues = cues
  )
}

#' @rdname family-mnix
#' @export
new_weighted_integrator_mnix_ideal_adaptor_from_data <- function(
  data,
  category = "category",
  cues,
  weights = NULL,
  kappa = nu,
  nu = 3,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate)
) {
  n_cues <- length(cues)
  if (is.null(weights)) {
    weights <- rep(1 / n_cues, n_cues)
  }
  weights <- as.numeric(weights)

  category_labels <- sort(unique(as.character(data[[category]])))
  representations <- lapply(category_labels, function(label) {
    sub_data <- data[data[[category]] == label, , drop = FALSE]
    m_val <- vapply(cues, function(c_col) mean(sub_data[[c_col]], na.rm = TRUE), numeric(1L))
    s2_val <- vapply(cues, function(c_col) stats::var(sub_data[[c_col]], na.rm = TRUE), numeric(1L))
    new_mnix_category_representation(
      category_labels = label,
      cue_labels = cues,
      m = m_val,
      kappa = if (length(kappa) == 1L) kappa else kappa[1L],
      nu = if (length(nu) == 1L) nu else nu[1L],
      sigma2 = s2_val,
      weights = weights
    )
  })
  names(representations) <- category_labels
  template <- new_category_representation_template(representations)

  new_weighted_integrator_mnix_ideal_adaptor(
    category_template = template,
    weights = weights,
    cues = cues,
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment
  )
}

#' @rdname family-mnix
#' @export
new_category_weighted_integrator_mnix_ideal_adaptor <- function(
  category_template = NULL,
  weights = NULL,
  weights_matrix = NULL,
  cues = NULL,
  input_cues = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
  metadata = list()
) {
  if (!is.null(weights_matrix)) {
    lifecycle::deprecate_warn("0.1.0", "new_category_weighted_integrator_mnix_ideal_adaptor(weights_matrix)", "new_category_weighted_integrator_mnix_ideal_adaptor(weights)")
    if (is.null(weights)) weights <- weights_matrix
  }
  if (!is.null(input_cues)) {
    lifecycle::deprecate_warn("0.1.0", "new_category_weighted_integrator_mnix_ideal_adaptor(input_cues)", "new_category_weighted_integrator_mnix_ideal_adaptor(cues)")
    if (is.null(cues)) cues <- input_cues
  }
  if (is.null(category_template)) {
    .stop("category_template must be supplied.")
  }
  reps <- category_template@representations
  n_cats <- length(reps)
  if (is.null(cues)) {
    cues <- get_cue_labels(reps[[1L]])
  }
  cues <- as.character(cues)
  n_cues <- length(cues)
  if (is.null(weights)) {
    weights <- matrix(1 / n_cues, nrow = n_cats, ncol = n_cues)
    rownames(weights) <- names(reps)
  }
  wm <- as.matrix(weights)
  if (nrow(wm) != n_cats || ncol(wm) != n_cues) {
    .stop("weights dimensions must match categories (", n_cats, ") and cues (", n_cues, ").")
  }

  .new_family_cognitive_model(
    category_template = category_template,
    category_representation_class = MNIX_CategoryRepresentation,
    model_class = CategoryWeightedIntegrator_MNIX_IdealAdaptor,
    family_label = "MNIX",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment,
    metadata = metadata,
    weights = wm,
    cues = cues
  )
}


#' @rdname family-mnix
#' @export
new_informativity_weighted_integrator_mnix_ideal_adaptor <- function(
  category_template = NULL,
  cues = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
  metadata = list()
) {
  if (is.null(category_template)) {
    .stop("category_template must be supplied.")
  }
  if (is.null(cues)) {
    cues <- get_cue_labels(category_template@representations[[1L]])
  }
  cues <- as.character(cues)

  .new_family_cognitive_model(
    category_template = category_template,
    category_representation_class = MNIX_CategoryRepresentation,
    model_class = InformativityWeightedIntegrator_MNIX_IdealAdaptor,
    family_label = "MNIX",
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment,
    metadata = metadata,
    cues = cues
  )
}

# -------------------------
# Deprecated
# -------------------------

#' @rdname family-mnix
#' @usage NULL
#' @export
MNIX_IdealAdaptor_AdaptiveWeightedIntegrator <- S7::new_class(
  "MNIX_IdealAdaptor_AdaptiveWeightedIntegrator",
  parent = MNIX_IdealAdaptor_InformativityWeightedIntegrator
)

#' @rdname family-mnix
#' @usage NULL
#' @export
MNIX_IdealAdaptor_CategoryAdaptiveWeightedIntegrator <- S7::new_class(
  "MNIX_IdealAdaptor_CategoryAdaptiveWeightedIntegrator",
  parent = MNIX_IdealAdaptor
)

#' @rdname family-mnix
#' @export
new_mnix_ideal_adaptor_adaptive_weighted_integrator <- function(...) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "new_mnix_ideal_adaptor_adaptive_weighted_integrator()",
    "new_mnix_ideal_adaptor_informativity_weighted_integrator()"
  )
  new_mnix_ideal_adaptor_informativity_weighted_integrator(...)
}

#' @rdname family-mnix
#' @export
new_mnix_ideal_adaptor_category_adaptive_weighted_integrator <- function(...) {
  lifecycle::deprecate_stop(
    "0.1.0",
    "new_mnix_ideal_adaptor_category_adaptive_weighted_integrator()"
  )
}
