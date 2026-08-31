#' @include S7-core-classes.R
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
#' @name family-mnix
#' @rdname family-mnix
#' @param category_labels Character vector of category label(s).
#' @param cue_labels Character vector of cue labels (at least 2 cue dimensions).
#' @param component_m Numeric vector of per-cue location hyperparameters.
#' @param component_kappa Numeric vector of per-cue precision scalings.
#' @param component_nu Numeric vector of per-cue degrees of freedom.
#' @param component_sigma2 Numeric vector of per-cue scale parameters.
#' @param component_weights Optional numeric vector of component weights
#'   summing to 1.
#' @param data Data frame containing category and cue observations.
#' @param category_var Bare symbol or character string indicating the category
#'   column in `data`.
#' @param cue_vars Bare tidyselect specification or character vector of cue
#'   columns in `data`.
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
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' Under independent cue integration, the category likelihood is the product
#' of marginal Student-\eqn{t} predictive components across cues:
#' \deqn{p(\mathbf{x} \mid c) = \prod_{i=1}^D \frac{1}{s_i}
#' t_{\nu_i}\!\left(\frac{x_i - m_i}{s_i}\right)}
#'
#' @return An S7 object of the respective MNIX class.
#' @seealso [family-muvg], [family-niw], [family-nix]
#' @export
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

    if (n_comp < 2) {
      return("MNIX representation must contain at least two cue components.")
    }
    if (
      length(self@component_kappa) != n_comp ||
      length(self@component_nu) != n_comp ||
      length(self@component_sigma2) != n_comp ||
      length(self@component_weights) != n_comp
    ) {
      return("MNIX component parameter vectors must all have equal length.")
    }
    label_information <- .mvbu_label_information(self@metadata)
    if (length(label_information$cue) != n_comp) {
      return("MNIX cue_labels length must match number of cue components.")
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

#' @rdname family-mnix
#' @export
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
  if (n_comp < 2) {
    .stop("component_m must contain at least two cue elements.")
  }
  if (!is.numeric(component_m) || !is.numeric(component_kappa) ||
      !is.numeric(component_nu) || !is.numeric(component_sigma2)) {
    .stop(
      "component_m, component_kappa, component_nu, and",
      " component_sigma2 must be numeric."
    )
  }
  if (is.null(component_weights)) {
    component_weights <- rep(1 / n_comp, n_comp)
  }
  if (!is.numeric(component_weights)) {
    .stop("component_weights must be numeric.")
  }

  MNIX_CategoryRepresentation(
    category_likelihood_function = {
      m0 <- as.numeric(component_m)
      kappa0 <- as.numeric(component_kappa)
      sigma20 <- as.numeric(component_sigma2)
      nu0 <- as.numeric(component_nu)
      pred_var <- sigma20 * (kappa0 + 1) / kappa0
      w0 <- as.numeric(component_weights)
      function(
        x,
        log = FALSE,
        noise_treatment = "no_noise",
        Sigma_noise = NULL
      ) {
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
        component_logdens <- sapply(seq_along(m0), function(i) {
          scale_i <- sqrt(pred_var[i] + noise_variance[i])
          stats::dt(
            (x[, i] - m0[i]) / scale_i,
            df = nu0[i],
            log = TRUE
          ) - log(scale_i)
        })
        if (is.vector(component_logdens)) {
          component_logdens <- matrix(component_logdens, ncol = length(m0))
        }
        total_logdens <- rowSums(component_logdens)
        if (isTRUE(log)) {
          total_logdens
        } else {
          exp(total_logdens)
        }
      }
    },
    metadata = .mvbu_label_metadata(
      as.character(category_labels),
      as.character(cue_labels),
      metadata
    ),
    component_m = as.numeric(component_m),
    component_kappa = as.numeric(component_kappa),
    component_nu = as.numeric(component_nu),
    component_sigma2 = as.numeric(component_sigma2),
    component_weights = as.numeric(component_weights)
  )
}

#' @rdname family-mnix
#' @export
new_mnix_category_representation_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = 3
) {
  .assert_non_NA_scalar_numeric(
    kappa,
    msg = "kappa must be a non-NA scalar numeric value."
  )
  .assert_non_NA_scalar_numeric(
    nu,
    msg = "nu must be a non-NA scalar numeric value."
  )
  .assert_true(
    nu > 2,
    msg = "nu must be greater than 2 for each MNIX component."
  )

  muvg <- new_muvg_category_representation_from_data(
    data,
    category = category,
    cues = cues
  )
  as_mnix_category_representation(muvg, kappa = kappa, nu = nu)
}

#' @rdname family-mnix
#' @export
new_mnix_category_representation_template_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = 3,
  verbose = FALSE
) {
  category_labels <- sort(unique(as.character(data[[category]])))
  representations <- lapply(category_labels, function(label) {
    new_mnix_category_representation_from_data(
      data[data[[category]] == label, , drop = FALSE],
      category = category,
      cues = cues,
      kappa = kappa,
      nu = nu
    )
  })
  names(representations) <- category_labels
  if (verbose) {
    message(
      "Constructed a MNIX category-representation template with ",
      length(category_labels), " categories and ",
      length(cues), " cue(s)."
    )
  }
  new_category_representation_template(representations)
}

#' @rdname family-mnix
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
  verbose = FALSE
) {
  template <- new_mnix_category_representation_template_from_data(
    data,
    category = category,
    cues = cues,
    kappa = kappa,
    nu = nu,
    verbose = verbose
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
