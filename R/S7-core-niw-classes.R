#' @include S7-core-classes.R
NULL

# -------------------------
# Normal-Inverse-Wishart (NIW) Family
# -------------------------

#' Normal-Inverse-Wishart (NIW) Family: Representations, Templates, and Models
#'
#' Constructs and manages Normal-Inverse-Wishart (NIW) category
#' representations, templates, and ideal adaptor cognitive models.
#'
#' @name family-niw
#' @rdname family-niw
#' @param category_labels Character vector of category label(s).
#' @param cue_labels Character vector of cue labels.
#' @param m Numeric vector of location hyperparameters.
#' @param kappa Positive numeric scalar precision scaling (\eqn{\kappa > 0}).
#' @param nu Numeric scalar degrees of freedom (\eqn{\nu > D - 1}).
#' @param S Numeric square symmetric scale matrix.
#' @param data Data frame containing category and cue observations.
#' @param category_var Bare symbol or character string indicating the category
#'   column in `data`.
#' @param cue_vars Bare tidyselect specification or character vector of cue
#'   columns in `data`.
#' @param category_template An [MVBU_CategoryRepresentationTemplate] object
#'   containing NIW category representations.
#' @param decision_rule Categorization decision rule: `"sampling"` or
#'   `"argmax"`. Defaults to `"sampling"`.
#' @param category_prior Optional numeric vector of prior category
#'   probabilities summing to 1. Defaults to equal priors.
#' @param lapse_rate Numeric scalar lapse probability in `[0, 1]`. Defaults
#'   to 0.
#' @param lapse_bias Optional numeric vector of lapse category probabilities
#'   summing to 1.
#' @param Sigma_noise Optional perceptual/measurement noise covariance matrix.
#' @param noise_treatment Treatment of noise: `"no_noise"`, `"sample"`, or
#'   `"marginalize"`.
#' @param lapse_treatment Treatment of lapses: `"no_lapses"`, `"sample"`, or
#'   `"marginalize"`.
#' @param metadata Optional list of metadata.
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' The Normal-Inverse-Wishart (`NIW`) family induces a multivariate
#' Student-\eqn{t} posterior predictive distribution:
#' For cue dimensionality \eqn{D}, \eqn{\nu_t = \nu - D + 1} and
#' \eqn{\mathbf{\Sigma}_t = \frac{\kappa + 1}{\kappa \nu_t} \mathbf{S}}:
#' \deqn{p(\mathbf{x} \mid c) = t_{\nu_t}(\mathbf{x}; \mathbf{m},
#' \mathbf{\Sigma}_t)}
#'
#' @return An S7 object of the respective NIW class.
#' @seealso [family-mvg], [family-nix], [family-mnix]
#' @export
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
      return(
        "NIW nu must be a scalar greater than cue dimensionality minus one."
      )
    }
    if (
      !is.matrix(self@S) ||
      !is.numeric(self@S) ||
      nrow(self@S) != d ||
      ncol(self@S) != d
    ) {
      return(
        "NIW S must be a numeric square matrix with cue dimensionality."
      )
    }
    if (!isTRUE(all.equal(self@S, t(self@S), tolerance = MVBU_PROB_TOL))) {
      return("NIW S must be symmetric.")
    }

    NULL
  }
)

#' @rdname family-niw
#' @export
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
      function(
        x,
        log = FALSE,
        noise_treatment = "no_noise",
        Sigma_noise = NULL
      ) {
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(
            n = nrow(x),
            mean = rep(0, ncol(x)),
            sigma = Sigma_noise
          )
        }
        if (!is.null(Sigma_noise) &&
            (identical(noise_treatment, "sample") ||
             identical(noise_treatment, "marginalize"))) {
          Sigma_eff <- Sigma0 + Sigma_noise
        } else {
          Sigma_eff <- Sigma0
        }
        .dmvt_density(x, mean = m0, Sigma = Sigma_eff, df = df0, log = log)
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
    S = as.matrix(S)
  )
}

#' @rdname family-niw
#' @export
new_niw_category_representation_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = length(cues) + 2
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
    nu > length(cues) + 1,
    msg = paste0(
      "nu must be larger than dimensionality of cues + 1 (>",
      length(cues) + 1, ")."
    )
  )

  mvg <- new_mvg_category_representation_from_data(
    data,
    category = category,
    cues = cues
  )
  as_niw_category_representation(mvg, kappa = kappa, nu = nu)
}

#' @rdname family-niw
#' @export
new_niw_category_representation_template_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = length(cues) + 2,
  verbose = FALSE
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
    nu > length(cues) + 1,
    msg = paste0(
      "nu must be larger than dimensionality of cues + 1 (>",
      length(cues) + 1, ")."
    )
  )

  mvg_template <- new_mvg_category_representation_template_from_data(
    data,
    category = category,
    cues = cues,
    verbose = FALSE
  )

  representations <- lapply(
    names(mvg_template@representations),
    function(label) {
      new_niw_category_representation_from_data(
        data[data[[category]] == label, , drop = FALSE],
        category = category,
        cues = cues,
        kappa = kappa,
        nu = nu
      )
    }
  )
  names(representations) <- names(mvg_template@representations)

  if (verbose) {
    message(
      "S is set so that the expected category covariance matrix Sigma ",
      "matches the category covariance in the sample (given nu). ",
      "It might be safer to fit an Inverse-Wishart distribution to the ",
      "entire set of covariance matrices."
    )
  }

  new_category_representation_template(representations)
}

#' @rdname family-niw
#' @export
NIW_IdealAdaptor <- S7::new_class(
  "NIW_IdealAdaptor",
  parent = MVBU_CognitiveModel
)

#' @rdname family-niw
#' @export
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
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment,
    metadata = metadata
  )
}

#' @rdname family-niw
#' @export
new_niw_ideal_adaptor_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = length(cues) + 2,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = "no_noise",
  lapse_treatment = "no_lapses",
  verbose = FALSE
) {
  template <- new_niw_category_representation_template_from_data(
    data,
    category = category,
    cues = cues,
    kappa = kappa,
    nu = nu,
    verbose = verbose
  )
  new_niw_ideal_adaptor(
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
