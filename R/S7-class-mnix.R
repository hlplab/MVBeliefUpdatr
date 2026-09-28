#' @include S7-class.R
NULL

# -------------------------
# Multi-Cue Normal-Inverse-Chi-Squared (MNIX) Family
# -------------------------

#' Multi-Cue Normal-Inverse-Chi-Squared (MNIX) Family:
#' Representations, Templates, and Models
#'
#' Constructs and manages multi-cue Normal-Inverse-\eqn{\chi^2} (MNIX) category
#' representations, templates, and ideal adaptors.
#'
#' @param category_labels Character vector of category label(s).
#' @param cue_labels Character vector of cue labels (at least 2 cue dimensions).
#' @param m Numeric vector of per-cue location hyperparameters.
#' @param kappa Positive numeric scalar or vector of per-cue
#'   precision scalings (\eqn{\kappa > 0}). If a scalar is provided, it is
#'   recycled across all cue components.
#' @param nu Positive numeric scalar or vector of per-cue degrees of
#'   freedom (\eqn{\nu > 2}). If a scalar is provided, it is recycled across
#'   all cue components.
#' @param sigma2 Numeric vector of per-cue scale parameters.
#' @param weights Optional numeric vector of component weights
#'   summing to 1.
#' @param data Data frame containing category and cue observations.
#' @param category Character string indicating the category column in `data`. Defaults to `"category"`.
#' @param cues Character vector indicating cue columns in `data`.
#' @param category_template An [MVBU_CategoryRepresentationTemplate] object
#'   containing MNIX category representations.
#' @param decision_rule Categorization decision rule: `"sampling"` or
#'   `"argmax"`. Defaults to `"sampling"`.
#' @param category_prior Optional numeric vector of prior category
#'   probabilities summing to 1. Defaults to equal priors.
#' @param lapse_rate Numeric scalar lapse probability in `[0, 1]`. Defaults
#'   to 0.
#' @param lapse_bias Optional numeric vector of lapse category probabilities
#'   summing to 1.
#' @param Sigma_noise Optional perceptual/measurement noise covariance matrix
#'   or vector.
#' @param noise_treatment Treatment of noise: `"no_noise"`, `"sample"`, or
#'   `"marginalize"`.
#' @param lapse_treatment Treatment of lapses: `"no_lapses"`, `"sample"`, or
#'   `"marginalize"`.
#' @param metadata Optional list of metadata.
#' @param ... Additional arguments passed to methods or constructors.
#'
#' @details
#' Under independent cue integration, the category likelihood is the product
#' of marginal Student-\eqn{t} predictive components across cues:
#' \deqn{p(\mathbf{x} \mid c) = \prod_{i=1}^D \frac{1}{s_i}
#' t_{\nu_i}\!\left(\frac{x_i - m_i}{s_i}\right)}
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

    if (n_comp < 2) {
      return("MNIX representation must contain at least two cue components.")
    }
    if (
      length(self@kappa) != n_comp ||
        length(self@nu) != n_comp ||
        length(self@sigma2) != n_comp ||
        length(self@weights) != n_comp
    ) {
      return("MNIX parameter vectors must all have equal length.")
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
  if (n_comp < 2) {
    .stop("m must contain at least two cue elements.")
  }
  if (!is.numeric(m) || !is.numeric(kappa) ||
    !is.numeric(nu) || !is.numeric(sigma2)) {
    .stop(
      "m, kappa, nu, and sigma2 must be numeric."
    )
  }
  if (is.null(weights)) {
    weights <- rep(1 / n_comp, n_comp)
  }
  if (!is.numeric(weights)) {
    .stop("weights must be numeric.")
  }

  MNIX_CategoryRepresentation(
    category_likelihood_function = {
      m0 <- as.numeric(m)
      kappa0 <- as.numeric(kappa)
      sigma20 <- as.numeric(sigma2)
      nu0 <- as.numeric(nu)
      pred_var <- sigma20 * (kappa0 + 1) / kappa0
      w0 <- as.numeric(weights)
      function(x,
               log = FALSE,
               noise_treatment = "no_noise",
               Sigma_noise = NULL) {
        x <- .as_observation_matrix(x, d = length(m0), arg_name = "x")
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
          diag(Sigma_noise)
        } else {
          rep(0, length(m0))
        }
        cue_logdens <- sapply(seq_along(m0), function(i) {
          scale_i <- sqrt(pred_var[i] + noise_variance[i])
          stats::dt(
            (x[, i] - m0[i]) / scale_i,
            df = nu0[i],
            log = TRUE
          ) - log(scale_i)
        })
        if (is.vector(cue_logdens)) {
          cue_logdens <- matrix(cue_logdens, ncol = length(m0))
        }
        total_logdens <- rowSums(cue_logdens)
        if (isTRUE(log)) {
          total_logdens
        } else {
          exp(total_logdens)
        }
      }
    },
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

#' @rdname family-mnix
#' @docType class
#' @usage NULL
#' @export
MNIX_IdealAdaptor <- S7::new_class(
  "MNIX_IdealAdaptor",
  parent = MVBU_CognitiveModel
)

#' @rdname family-mnix
#' @export
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
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = "no_noise",
  lapse_treatment = "no_lapses",
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
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment
  )
}
