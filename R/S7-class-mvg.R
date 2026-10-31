#' @include S7-class.R
NULL

# -------------------------
# Multivariate Gaussian (MVG) Family
# -------------------------

#' Multivariate Gaussian (MVG) Family: Representations, Templates, and Models
#'
#' Constructs and manages multivariate Gaussian category representations,
#' templates, and ideal observer models.
#'
#' @param category_labels Character vector of category label(s).
#' @param cue_labels Character vector of cue labels.
#' @param mu Numeric vector of category means.
#' @param Sigma Numeric square symmetric positive-definite covariance matrix.
#' @param data Data frame containing category and cue observations.
#' @param category Character string indicating the category column in `data`. Defaults to `"category"`.
#' @param cues Character vector indicating cue columns in `data`.
#' @param category_template An [MVBU_CategoryRepresentationTemplate] object
#'   containing MVG category representations.
#' @param decision_rule Categorization decision rule: `"sampling"`,
#'   `"criterion"`, or `"proportional"`. Defaults to `"sampling"`.
#' @param category_prior Optional numeric vector of prior category
#'   probabilities summing to 1. Defaults to equal priors.
#' @param lapse_rate Numeric scalar lapse probability in `[0, 1]`. Defaults
#'   to 0.
#' @param lapse_bias Optional numeric vector of lapse category probabilities
#'   summing to 1.
#' @param Sigma_noise Optional perceptual/measurement noise covariance matrix.
#' @param noise_treatment Treatment of perceptual noise during evaluation:
#'   `"no_noise"`, `"sample"`, or `"marginalize"`.
#'   \itemize{
#'     \item \code{"no_noise"}: No perceptual noise is applied to input cues,
#'       and no noise is assumed during categorization.
#'     \item \code{"sample"}: Simulates sensory noise on a specific observation.
#'       Noise drawn from \eqn{\mathcal{N}(0, \Sigma_{\text{noise}})} is added
#'       to input cues, and the observer accounts for their own noise
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
#' @references
#' \insertRef{xie2023}{MVBeliefUpdatr}
#'
#' @details
#' The multivariate Gaussian (`MVG`) family models category structure over
#' \eqn{D}-dimensional continuous cue spaces with correlated cue dimensions:
#' \deqn{p(\mathbf{x} \mid c) = \mathcal{N}(\mathbf{x}; \boldsymbol{\mu}_c,
#' \mathbf{\Sigma}_c)}
#'
#' @return An S7 object of the respective MVG class.
#' @seealso [family-niw], [family-uvg], [family-exemplar]
#'
#' @name family-mvg
#' @rdname family-mvg
#' @docType class
#' @usage NULL
#' @export
#' @examples
#' # Create 2D MVG representations
#' mu_b <- c(VOT = 10, f0 = 200)
#' Sigma_b <- matrix(c(25, 5, 5, 36), nrow = 2,
#'   dimnames = list(names(mu_b), names(mu_b)))
#' rep_b <- new_mvg_category_representation(
#'   category_labels = "/b/", cue_labels = c("VOT", "f0"),
#'   mu = mu_b, Sigma = Sigma_b
#' )
#' mu_d <- c(VOT = 50, f0 = 240)
#' Sigma_d <- matrix(c(25, -2, -2, 36), nrow = 2,
#'   dimnames = list(names(mu_d), names(mu_d)))
#' rep_d <- new_mvg_category_representation(
#'   category_labels = "/d/", cue_labels = c("VOT", "f0"),
#'   mu = mu_d, Sigma = Sigma_d
#' )
#' tmpl <- new_category_representation_template(list(rep_b, rep_d))
#' m_mvg <- new_mvg_ideal_observer(category_template = tmpl)
#' m_mvg
MVG_CategoryRepresentation <- S7::new_class(
  "MVG_CategoryRepresentation",
  parent = MVBU_CategoryRepresentation,
  properties = list(
    mu = S7::class_numeric,
    Sigma = S7::class_any
  ),
  validator = function(self) {
    label_information <- self@metadata$label_information
    d <- length(label_information$cue)

    if (length(self@mu) != d) {
      return("MVG mu must have length equal to cue dimensionality.")
    }
    if (
      !is.matrix(self@Sigma) ||
        !is.numeric(self@Sigma) ||
        nrow(self@Sigma) != d ||
        ncol(self@Sigma) != d
    ) {
      return(
        "MVG Sigma must be a numeric square matrix with cue dimensionality."
      )
    }
    if (
      !isTRUE(
        all.equal(self@Sigma, t(self@Sigma), tolerance = MVBU_PROB_TOL)
      )
    ) {
      return("MVG Sigma must be symmetric.")
    }

    NULL
  }
)

#' @rdname family-mvg
#' @export
new_mvg_category_representation <- function(
  category_labels,
  cue_labels,
  mu,
  Sigma,
  metadata = list()
) {
  # Rationale for retaining typed class constructor:
  # Concrete typed constructors provide strict property validation, IDE tab completion,
  # and fast S7 method dispatch, while obtaining the likelihood closure uniformly
  # from the registered model family density function.
  metadata <- as.list(metadata)
  if (is.null(metadata$parameters)) {
    metadata$parameters <- list(
      mu = as.numeric(mu),
      Sigma = as.matrix(Sigma)
    )
  }
  rep <- MVG_CategoryRepresentation(
    category_likelihood_function = NULL,
    metadata = set_labels(
      metadata,
      cue = as.character(cue_labels),
      category = as.character(category_labels),
      response_category = as.character(category_labels),
      group = if (!is.null(metadata$group)) as.character(metadata$group) else character(0)
    ),
    mu = as.numeric(mu),
    Sigma = as.matrix(Sigma)
  )
  rep@category_likelihood_function <- .make_family_likelihood_function(rep)
  rep
}

#' @rdname family-mvg
#' @export
new_mvg_category_representation_from_data <- function(
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
  new_mvg_category_representation(
    category_labels = category_labels,
    cue_labels = cues,
    mu = .colMeans(cue_values),
    Sigma = .cov(cue_values)
  )
}

#' @rdname family-mvg
#' @export
new_mvg_category_representation_template_from_data <- function(
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

  category_labels <- sort(unique(as.character(data[[category]])))
  representations <- vector("list", length(category_labels))
  names(representations) <- category_labels

  for (label in category_labels) {
    representations[[label]] <- new_mvg_category_representation_from_data(
      data[data[[category]] == label, , drop = FALSE],
      category = category,
      cues = cues
    )
  }

  new_category_representation_template(representations)
}

#' @rdname family-mvg
#' @docType class
#' @usage NULL
#' @export
MVG_IdealObserver <- S7::new_class(
  "MVG_IdealObserver",
  parent = MVBU_CognitiveModel
)

#' @rdname family-mvg
#' @export
new_mvg_ideal_observer <- function(
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
    category_representation_class = MVG_CategoryRepresentation,
    model_class = MVG_IdealObserver,
    family_label = "MVG",
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

#' @rdname family-mvg
#' @export
new_mvg_ideal_observer_from_data <- function(
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
  template <- new_mvg_category_representation_template_from_data(
    data,
    category = category,
    cues = cues
  )
  new_mvg_ideal_observer(
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
