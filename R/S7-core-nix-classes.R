#' @include S7-core-classes.R
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
#' @name family-nix
#' @rdname family-nix
#' @param category_labels Character vector of category label(s).
#' @param cue_labels Character vector of cue label(s). For NIX, exactly one
#'   cue dimension is allowed.
#' @param m Numeric scalar location hyperparameter.
#' @param kappa Positive numeric scalar precision scaling (\eqn{\kappa > 0}).
#' @param nu Positive numeric scalar degrees of freedom (\eqn{\nu > 0}).
#' @param sigma2 Positive numeric scalar scale parameter (\eqn{\sigma^2 > 0}).
#' @param data Data frame containing category and cue observations.
#' @param category_var Bare symbol or character string indicating the category
#'   column in `data`.
#' @param cue_var Bare symbol or character string indicating the cue column
#'   in `data`.
#' @param category_template An [MVBU_CategoryRepresentationTemplate] object
#'   containing NIX category representations.
#' @param decision_rule Categorization decision rule: `"sampling"` or
#'   `"argmax"`. Defaults to `"sampling"`.
#' @param category_prior Optional numeric vector of prior category
#'   probabilities summing to 1. Defaults to equal priors.
#' @param lapse_rate Numeric scalar lapse probability in `[0, 1]`. Defaults
#'   to 0.
#' @param lapse_bias Optional numeric vector of lapse category probabilities
#'   summing to 1.
#' @param Sigma_noise Optional perceptual/measurement noise covariance matrix
#'   or scalar variance.
#' @param noise_treatment Treatment of noise: `"no_noise"`, `"sample"`, or
#'   `"marginalize"`.
#' @param lapse_treatment Treatment of lapses: `"no_lapses"`, `"sample"`, or
#'   `"marginalize"`.
#' @param metadata Optional list of metadata.
#' @param ... Additional arguments passed to methods.
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
  NIX_CategoryRepresentation(
    category_likelihood_function = {
      m0 <- as.numeric(m)
      kappa0 <- as.numeric(kappa)
      nu0 <- as.numeric(nu)
      sigma20 <- as.numeric(sigma2)
      function(
        x,
        log = FALSE,
        noise_treatment = "no_noise",
        Sigma_noise = NULL
      ) {
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
        scale_eff <- sqrt((sigma20 * (kappa0 + 1) / kappa0) + noise_variance)
        z <- (x[, 1] - m0) / scale_eff
        if (isTRUE(log)) {
          stats::dt(z, df = nu0, log = TRUE) - log(scale_eff)
        } else {
          stats::dt(z, df = nu0) / scale_eff
        }
      }
    },
    metadata = .mvbu_label_metadata(
      as.character(category_labels),
      as.character(cue_labels),
      metadata
    ),
    m = as.numeric(m),
    kappa = as.numeric(kappa),
    nu = as.numeric(nu),
    sigma2 = as.numeric(sigma2)
  )
}

#' @rdname family-nix
#' @export
new_nix_category_representation_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = 3
) {
  .assert_numeric_scalar(
    kappa,
    msg = "kappa must be a non-NA scalar numeric value."
  )
  .assert_numeric_scalar(
    nu,
    msg = "nu must be a non-NA scalar numeric value."
  )
  .assert_true(
    nu > 2,
    msg = "nu must be greater than 2 for a univariate NIX representation."
  )

  uvg <- new_uvg_category_representation_from_data(
    data,
    category = category,
    cues = cues
  )
  as_nix_category_representation(uvg, kappa = kappa, nu = nu)
}

#' @rdname family-nix
#' @export
new_nix_category_representation_template_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = 3,
  verbose = FALSE
) {
  .assert_true(
    length(cues) == 1L,
    msg = "NIX templates require exactly one cue."
  )
  category_labels <- sort(unique(as.character(data[[category]])))
  representations <- lapply(category_labels, function(label) {
    new_nix_category_representation_from_data(
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
      "Constructed a NIX category-representation template with ",
      length(category_labels), " categories."
    )
  }
  new_category_representation_template(representations)
}

#' @rdname family-nix
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
  template <- new_nix_category_representation_template_from_data(
    data,
    category = category,
    cues = cues,
    kappa = kappa,
    nu = nu,
    verbose = verbose
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
