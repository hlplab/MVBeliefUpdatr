#' @include S7-class.R S7-class-mnix.R
NULL

# -------------------------
# Normal-Inverse-Chi-Squared (NIX) Family
# -------------------------

#' Normal-Inverse-Chi-Squared (NIX) Family: Representations,
#' Templates, and Models
#'
#' Constructs and manages Normal-Inverse-\eqn{\chi^2} (NIX) category
#' representations, templates, and ideal adaptor cognitive models
#' for 1D cues under belief updating.
#'
#' @param category_labels Character vector of category label(s).
#' @param cue_labels Character vector of cue label(s). For NIX, exactly one
#'   cue dimension is allowed.
#' @param m Numeric scalar location hyperparameter.
#' @param kappa Positive numeric scalar precision scaling (\eqn{\kappa > 0}).
#' @param nu Positive numeric scalar degrees of freedom (\eqn{\nu > 0}).
#' @param sigma2 Positive numeric scalar scale parameter (\eqn{\sigma^2 > 0}).
#' @param data Data frame containing category and cue observations.
#' @param category Character string indicating the category column in `data`. Defaults to `"category"`.
#' @param cues Character string indicating the cue column in `data`. Defaults to `"cue"`.
#' @param category_template An [MVBU_CategoryRepresentationTemplate] object
#'   containing NIX category representations.
#' @param decision_rule Categorization decision rule: `"sampling"` or
#'   `"argmax"`. Defaults to `"sampling"`.
#' @param category_prior Optional numeric vector of prior category
#'   probabilities summing to 1. Defaults to equal priors.
#' @param lapse_rate Numeric scalar lapse probability in `[0, 1]`. Defaults
#'   to 0.
#' @param lapse_bias Optional numeric vector of lapse category probabilities
#'   summing to 1. Defaults to uniform bias `1/M`.
#' @param Sigma_noise Optional perceptual/measurement noise covariance matrix
#'   or scalar variance.
#' @param noise_treatment Treatment of noise: `"no_noise"`, `"sample"`, or
#'   `"marginalize"`. Defaults to `.infer_noise_treatment(Sigma_noise)`:
#'   `"marginalize"` if noise is present, `"no_noise"` otherwise.
#' @param lapse_treatment Treatment of lapses: `"no_lapses"`, `"sample"`, or
#'   `"marginalize"`. Defaults to `.infer_lapse_treatment(lapse_rate)`:
#'   `"marginalize"` if `lapse_rate > 0`, `"no_lapses"` otherwise.
#' @param metadata Optional list of metadata. Defaults to `list()`.
#' @param ... Additional arguments passed to methods or constructors.
#'
#' @details
#' The Normal-Inverse-\eqn{\chi^2} (`NIX`) family models adaptor-style
#' uncertainty over univariate Gaussian category structure. The predictive
#' likelihood is a Student-\eqn{t} distribution:
#' With \eqn{s_c = \sqrt{\sigma_c^2(\kappa_c+1)/\kappa_c}},
#' \deqn{p(x \mid c) =
#'   \frac{1}{s_c} \, t_{\nu_c}\!\left(\frac{x-m_c}{s_c}\right)}
#'
#' @return An S7 object of the respective NIX class.
#' @seealso [family-uvg], [family-niw], [family-mnix]
#'
#' @name family-nix
#' @rdname family-nix
#' @docType class
#' @usage NULL
#' @export
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
    label_information <- self@metadata$label_information
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

#' @rdname family-nix
#' @export
new_nix_category_representation <- function(
  category_labels,
  cue_labels,
  m,
  kappa,
  nu,
  sigma2,
  metadata = list()
) {
  metadata <- as.list(metadata)
  if (is.null(metadata$parameters)) {
    metadata$parameters <- list(
      m = as.numeric(m),
      kappa = as.numeric(kappa),
      nu = as.numeric(nu),
      sigma2 = as.numeric(sigma2)
    )
  }
  rep <- NIX_CategoryRepresentation(
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
    sigma2 = as.numeric(sigma2)
  )
  rep@category_likelihood_function <- .make_family_likelihood_function(rep)
  rep
}

#' @rdname family-nix
#' @export
new_nix_category_representation_from_data <- function(
  data,
  category = "category",
  cues = "cue",
  kappa = nu,
  nu = 3,
  ...
) {
  .assert_true(
    length(cues) == 1L,
    msg = "NIX representations are strictly 1D and require length(cues) == 1."
  )
  .assert_true(
    length(kappa) == 1L && is.numeric(kappa) && !is.na(kappa) && kappa > 0,
    msg = "kappa must be a single positive numeric value."
  )
  .assert_true(
    length(nu) == 1L && is.numeric(nu) && !is.na(nu) && nu > 2,
    msg = "nu must be a single numeric value > 2."
  )

  uvg <- new_uvg_category_representation_from_data(
    data,
    category = category,
    cue = cues
  )
  as_nix_category_representation(
    uvg,
    kappa = kappa,
    nu = nu
  )
}

#' @rdname family-nix
#' @export
new_nix_category_representation_template_from_data <- function(
  data,
  category = "category",
  cues = "cue",
  kappa = nu,
  nu = 3,
  ...
) {
  category_labels <- sort(unique(as.character(data[[category]])))
  representations <- lapply(category_labels, function(label) {
    new_nix_category_representation_from_data(
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

#' @rdname family-nix
#' @docType class
#' @usage NULL
#' @export
NIX_IdealAdaptor <- S7::new_class(
  "NIX_IdealAdaptor",
  parent = MVBU_CognitiveModel
)

#' @rdname family-nix
#' @export
new_nix_ideal_adaptor <- function(
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
    category_representation_class = NIX_CategoryRepresentation,
    model_class = NIX_IdealAdaptor,
    family_label = "NIX",
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

#' @rdname family-nix
#' @export
new_nix_ideal_adaptor_from_data <- function(
  data,
  category = "category",
  cues = "cue",
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
  template <- new_nix_category_representation_template_from_data(
    data,
    category = category,
    cues = cues,
    kappa = kappa,
    nu = nu
  )
  new_nix_ideal_adaptor(
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
