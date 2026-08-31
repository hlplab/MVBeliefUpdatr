#' @include S7-core-classes.R
NULL

# -------------------------
# Univariate Gaussian (UVG) Family
# -------------------------

#' Univariate Gaussian (UVG) Family: Representations, Templates, and Models
#'
#' Constructs and manages univariate Gaussian category representations,
#' templates, and ideal observer cognitive models.
#'
#' @name family-uvg
#' @rdname family-uvg
#' @param category_labels Character vector of category label(s).
#' @param cue_labels Character vector of cue label(s). For UVG, exactly one cue
#'   dimension is allowed.
#' @param mu Numeric scalar mean.
#' @param sigma2 Positive numeric scalar variance (\eqn{\sigma^2 > 0}).
#' @param data Data frame containing category and cue observations.
#' @param category_var Bare symbol or character string indicating the category
#'   column in `data`.
#' @param cue_var Bare symbol or character string indicating the cue column
#'   in `data`.
#' @param category_template An [MVBU_CategoryRepresentationTemplate] object
#'   containing UVG category representations.
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
#' The univariate Gaussian (`UVG`) family models category structure over a
#' single continuous cue dimension:
#' \deqn{p(x \mid c) = \mathcal{N}(x; \mu_c, \sigma_c^2)}
#'
#' Category representations are constructed with
#' [new_uvg_category_representation] or
#' [new_uvg_category_representation_from_data]. Cognitive ideal observer
#' models are constructed with [new_uvg_ideal_observer] or
#' [new_uvg_ideal_observer_from_data].
#'
#' @return An S7 object of the respective UVG class.
#' @seealso [family-nix], [family-mvg], [family-muvg], [family-exemplar]
#' @export
UVG_CategoryRepresentation <- S7::new_class(
  "UVG_CategoryRepresentation",
  parent = MVBU_CategoryRepresentation,
  properties = list(
    mu = S7::class_numeric,
    sigma2 = S7::class_numeric
  ),
  validator = function(self) {
    label_information <- .mvbu_label_information(self@metadata)
    if (length(label_information$cue) != 1) {
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

#' @rdname family-uvg
#' @export
UVG_IdealObserver <- S7::new_class(
  "UVG_IdealObserver",
  parent = MVBU_CognitiveModel
)

#' @rdname family-uvg
#' @export
new_uvg_category_representation <- function(
  category_labels,
  cue_labels,
  mu,
  sigma2,
  metadata = list()
) {
  UVG_CategoryRepresentation(
    category_likelihood_function = {
      mu0 <- as.numeric(mu)
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
        stats::dnorm(
          x[, 1],
          mean = mu0,
          sd = sqrt(sigma20 + noise_variance),
          log = log
        )
      }
    },
    metadata = .mvbu_label_metadata(
      as.character(category_labels),
      as.character(cue_labels),
      metadata
    ),
    mu = as.numeric(mu),
    sigma2 = as.numeric(sigma2)
  )
}

#' @rdname family-uvg
#' @export
new_uvg_ideal_observer <- function(
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
    category_representation_class = UVG_CategoryRepresentation,
    model_class = UVG_IdealObserver,
    family_label = "UVG",
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
