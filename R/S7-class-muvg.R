#' @include S7-class.R
NULL

# -------------------------
# Multi-Cue Univariate Gaussian (MUVG) Family
# -------------------------

#' Multi-Cue Univariate Gaussian (MUVG) Family:
#' Representations, Templates, and Models
#'
#' Constructs and manages multi-cue univariate Gaussian integration
#' representations, templates, and ideal observers.
#'
#' @param category_labels Character vector of category label(s).
#' @param cue_labels Character vector of cue labels matching the length of
#'   cue components.
#' @param mu Numeric vector of per-cue means.
#' @param sigma2 Numeric vector of per-cue positive variances.
#' @param weights Optional numeric vector of integration weights
#'   in `[0, 1]` summing to 1. If `NULL`, defaults to precision weights
#'   \eqn{w_i \propto 1/\sigma_i^2}.
#' @param data Data frame containing category and cue observations.
#' @param category Character string indicating the category column in `data`. Defaults to `"category"`.
#' @param cues Character vector indicating cue columns in `data`.
#' @param category_template An [MVBU_CategoryRepresentationTemplate] object
#'   containing MUVG category representations.
#' @param decision_rule Categorization decision rule: `"sampling"`,
#'   `"criterion"`, or `"proportional"`. Defaults to `"sampling"`.
#' @param category_prior Optional numeric vector of prior category
#'   probabilities summing to 1. Defaults to equal priors.
#' @param lapse_rate Numeric scalar lapse probability in `[0, 1]`. Defaults
#'   to 0.
#' @param lapse_bias Optional numeric vector of lapse category probabilities
#'   summing to 1.
#' @param Sigma_noise Optional perceptual/measurement noise covariance matrix
#'   or vector.
#' @param noise_treatment Treatment of perceptual noise during evaluation:
#'   `"no_noise"`, `"sample"`, or `"marginalize"`.
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
#' @param lapse_treatment Treatment of attentional lapses during evaluation:
#'   `"no_lapses"`, `"sample"`, or `"marginalize"`.
#'   \itemize{
#'     \item \code{"no_lapses"}: Lapses are ignored; categorization reflects
#'       pure normative evidence evaluation (even if \code{lapse_rate > 0}).
#'     \item \code{"sample"}: Simulates stochastic attentional lapses on
#'       individual observations. With probability \code{lapse_rate}, a lapse
#'       occurs: the posterior is replaced by \code{lapse_bias}.
#'     \item \code{"marginalize"}: Evaluates categorization in the limit across
#'       infinitely many encounters (expected response probabilities). The
#'       category posterior is the mixture
#'       \eqn{(1 - \lambda) p(c \mid x) + \lambda b_c}.
#'   }
#'   Defaults to `.infer_lapse_treatment(lapse_rate)`: `"marginalize"` if
#'   `lapse_rate > 0`, `"no_lapses"` otherwise.
#' @param metadata Optional list of metadata.
#'
#' @details
#' Under conditional cue independence and integrated cue
#' \eqn{z = \sum_i w_i x_i}, the category likelihood is:
#' \deqn{p(z \mid c) = \mathcal{N}\!\left(z;\; \sum_i w_i \mu_i,\;
#' \sum_i w_i^2 \sigma_i^2\right)}
#'
#' @return An S7 object of the respective MUVG class.
#' @seealso [family-mnix], [family-mvg], [family-uvg]
#'
#' @name family-muvg
#' @rdname family-muvg
#' @docType class
#' @usage NULL
#' @export
#' @examples
#' # MUVG: Multivariate independent Gaussian (diagonal covariance)
#' mu_b <- c(VOT = 10, f0 = 200)
#' sigma2_b <- c(VOT = 25, f0 = 36)
#' rep_b <- new_muvg_category_representation(
#'   category_labels = "/b/", cue_labels = c("VOT", "f0"),
#'   mu = mu_b, sigma2 = sigma2_b
#' )
#' mu_d <- c(VOT = 50, f0 = 240)
#' sigma2_d <- c(VOT = 25, f0 = 36)
#' rep_d <- new_muvg_category_representation(
#'   category_labels = "/d/", cue_labels = c("VOT", "f0"),
#'   mu = mu_d, sigma2 = sigma2_d
#' )
#' tmpl <- new_category_representation_template(list(rep_b, rep_d))
#' m_muvg <- new_muvg_ideal_observer(category_template = tmpl)
#' m_muvg
MUVG_CategoryRepresentation <- S7::new_class(
  "MUVG_CategoryRepresentation",
  parent = MVBU_CategoryRepresentation,
  properties = list(
    mu = S7::class_numeric,
    sigma2 = S7::class_numeric,
    weights = S7::class_numeric
  ),
  validator = function(self) {
    n_comp <- length(self@mu)

    if (n_comp < 1) {
      return("MUVG representation must contain at least one cue component.")
    }
    if (
      length(self@sigma2) != n_comp ||
      length(self@weights) != n_comp
    ) {
      return("MUVG parameter vectors must all have equal length.")
    }
    label_information <- self@metadata$label_information
    if (length(label_information$cue) != n_comp) {
      return("MUVG cue_labels length must match number of cue components.")
    }
    if (any(is.na(self@mu))) {
      return("MUVG mu entries must be numeric values.")
    }
    if (any(is.na(self@sigma2))) {
      return("MUVG sigma2 entries must be numeric values.")
    }
    if (any(self@sigma2 <= 0)) {
      return("MUVG sigma2 entries must be positive.")
    }
    if (any(self@weights < 0) || any(self@weights > 1)) {
      return("MUVG weights entries must be in [0, 1].")
    }
    if (abs(sum(self@weights) - 1) > MVBU_PROB_TOL) {
      return("MUVG weights entries must sum to 1.")
    }

    NULL
  }
)

#' @rdname family-muvg
#' @export
new_muvg_category_representation <- function(
  category_labels,
  cue_labels,
  mu,
  sigma2,
  weights = NULL,
  metadata = list()
) {
  n_comp <- length(mu)
  if (n_comp < 1) {
    .stop("mu must contain at least one element.")
  }
  if (!is.numeric(mu) || !is.numeric(sigma2)) {
    .stop("mu and sigma2 must be numeric.")
  }
  if (is.null(weights)) {
    precision <- 1 / as.numeric(sigma2)
    weights <- precision / sum(precision)
  }
  if (!is.numeric(weights)) {
    .stop("weights must be numeric.")
  }

  # Rationale for retaining typed class constructor:
  # Concrete typed constructors provide strict property validation, IDE tab completion,
  # and fast S7 method dispatch, while obtaining the likelihood closure uniformly
  # from the registered model family density function.
  metadata <- as.list(metadata)
  if (is.null(metadata$parameters)) {
    metadata$parameters <- list(
      mu = as.numeric(mu),
      sigma2 = as.numeric(sigma2),
      weights = as.numeric(weights)
    )
  }
  rep <- MUVG_CategoryRepresentation(
    category_likelihood_function = NULL,
    metadata = set_labels(
      metadata,
      cue = as.character(cue_labels),
      category = as.character(category_labels),
      response_category = as.character(category_labels),
      group = if (!is.null(metadata$group)) as.character(metadata$group) else character(0)
    ),
    mu = as.numeric(mu),
    sigma2 = as.numeric(sigma2),
    weights = as.numeric(weights)
  )
  rep@category_likelihood_function <- .make_family_likelihood_function(rep)
  rep
}

#' @rdname family-muvg
#' @export
new_muvg_category_representation_from_data <- function(
  data,
  category = "category",
  cues
) {
  .assert_data_frame_like(data)
  .assert_character_scalar(
    category,
    msg = "category must be a non-empty scalar character value."
  )
  .assert_true(
    is.character(cues) && length(cues) > 0,
    msg = "cues must be a non-empty character vector."
  )
  .assert_data_contains_cols(data, category)
  .assert_data_contains_cols(data, cues)
  category_labels <- unique(as.character(data[[category]]))
  .assert_true(
    length(category_labels) == 1L,
    msg = "data must contain exactly one category."
  )

  cue_values <- as.matrix(data[, cues, drop = FALSE])
  new_muvg_category_representation(
    category_labels = category_labels,
    cue_labels = cues,
    mu = .colMeans(cue_values),
    sigma2 = diag(.cov(cue_values))
  )
}

#' @rdname family-muvg
#' @export
new_muvg_category_representation_template_from_data <- function(
  data,
  category = "category",
  cues
) {
  category_labels <- sort(unique(as.character(data[[category]])))
  representations <- lapply(category_labels, function(label) {
    new_muvg_category_representation_from_data(
      data[data[[category]] == label, , drop = FALSE],
      category = category,
      cues = cues
    )
  })
  names(representations) <- category_labels
  new_category_representation_template(representations)
}

#' @rdname family-muvg
#' @docType class
#' @usage NULL
#' @export
MUVG_IdealObserver <- S7::new_class(
  "MUVG_IdealObserver",
  parent = MVBU_CognitiveModel
)

#' @rdname family-muvg
#' @export
new_muvg_ideal_observer <- function(
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
    category_representation_class = MUVG_CategoryRepresentation,
    model_class = MUVG_IdealObserver,
    family_label = "MUVG",
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

#' @rdname family-muvg
#' @export
new_muvg_ideal_observer_from_data <- function(
  data,
  category = "category",
  cues,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate)
) {
  template <- new_muvg_category_representation_template_from_data(
    data,
    category = category,
    cues = cues
  )
  new_muvg_ideal_observer(
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
