#' @include S7-core-classes.R
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
#' @name family-muvg
#' @rdname family-muvg
#' @param category_labels Character vector of category label(s).
#' @param cue_labels Character vector of cue labels matching the length of
#'   cue components.
#' @param component_mu Numeric vector of per-cue means.
#' @param component_sigma2 Numeric vector of per-cue positive variances.
#' @param component_weights Optional numeric vector of integration weights
#'   in `[0, 1]` summing to 1. If `NULL`, defaults to precision weights
#'   \eqn{w_i \propto 1/\sigma_i^2}.
#' @param data Data frame containing category and cue observations.
#' @param category_var Bare symbol or character string indicating the category
#'   column in `data`.
#' @param cue_vars Bare tidyselect specification or character vector of cue
#'   columns in `data`.
#' @param category_template An [MVBU_CategoryRepresentationTemplate] object
#'   containing MUVG category representations.
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
#' Under conditional cue independence and integrated cue
#' \eqn{z = \sum_i w_i x_i}, the category likelihood is:
#' \deqn{p(z \mid c) = \mathcal{N}\!\left(z;\; \sum_i w_i \mu_i,\;
#' \sum_i w_i^2 \sigma_i^2\right)}
#'
#' @return An S7 object of the respective MUVG class.
#' @seealso [family-mnix], [family-mvg], [family-uvg]
#' @export
MUVG_CategoryRepresentation <- S7::new_class(
  "MUVG_CategoryRepresentation",
  parent = MVBU_CategoryRepresentation,
  properties = list(
    component_mu = S7::class_numeric,
    component_sigma2 = S7::class_numeric,
    component_weights = S7::class_numeric
  ),
  validator = function(self) {
    n_comp <- length(self@component_mu)

    if (n_comp < 1) {
      return("MUVG representation must contain at least one cue component.")
    }
    if (
      length(self@component_sigma2) != n_comp ||
      length(self@component_weights) != n_comp
    ) {
      return("MUVG component parameter vectors must all have equal length.")
    }
    label_information <- .mvbu_label_information(self@metadata)
    if (length(label_information$cue) != n_comp) {
      return("MUVG cue_labels length must match number of cue components.")
    }
    if (any(is.na(self@component_mu))) {
      return("MUVG component_mu entries must be numeric values.")
    }
    if (any(is.na(self@component_sigma2))) {
      return("MUVG component_sigma2 entries must be numeric values.")
    }
    if (any(self@component_sigma2 <= 0)) {
      return("MUVG component_sigma2 entries must be positive.")
    }
    if (any(self@component_weights < 0) || any(self@component_weights > 1)) {
      return("MUVG component_weights entries must be in [0, 1].")
    }
    if (abs(sum(self@component_weights) - 1) > MVBU_PROB_TOL) {
      return("MUVG component_weights entries must sum to 1.")
    }

    NULL
  }
)

#' @rdname family-muvg
#' @export
new_muvg_category_representation <- function(
  category_labels,
  cue_labels,
  component_mu,
  component_sigma2,
  component_weights = NULL,
  metadata = list()
) {
  n_comp <- length(component_mu)
  if (n_comp < 1) {
    .stop("component_mu must contain at least one element.")
  }
  if (!is.numeric(component_mu) || !is.numeric(component_sigma2)) {
    .stop("component_mu and component_sigma2 must be numeric.")
  }
  if (is.null(component_weights)) {
    precision <- 1 / as.numeric(component_sigma2)
    component_weights <- precision / sum(precision)
  }
  if (!is.numeric(component_weights)) {
    .stop("component_weights must be numeric.")
  }

  MUVG_CategoryRepresentation(
    category_likelihood_function = {
      mu0 <- as.numeric(component_mu)
      sigma20 <- as.numeric(component_sigma2)
      w0 <- as.numeric(component_weights)
      z_mean <- sum(w0 * mu0)
      z_var <- sum((w0^2) * sigma20)
      function(
        x,
        log = FALSE,
        noise_treatment = "no_noise",
        Sigma_noise = NULL
      ) {
        x <- .as_observation_matrix(x, d = length(w0), arg_name = "x")
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
          sum((w0^2) * diag(Sigma_noise))
        } else {
          0
        }
        z <- as.numeric(x %*% w0)
        stats::dnorm(
          z,
          mean = z_mean,
          sd = sqrt(z_var + noise_variance),
          log = log
        )
      }
    },
    metadata = .mvbu_label_metadata(
      as.character(category_labels),
      as.character(cue_labels),
      metadata
    ),
    component_mu = as.numeric(component_mu),
    component_sigma2 = as.numeric(component_sigma2),
    component_weights = as.numeric(component_weights)
  )
}

#' @rdname family-muvg
#' @export
new_muvg_category_representation_from_data <- function(
  data,
  category = "category",
  cues
) {
  .assert_data_frame_like(data)
  .assert_non_NA_scalar_character(
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
    component_mu = .colMeans(cue_values),
    component_sigma2 = diag(.cov(cue_values))
  )
}

#' @rdname family-muvg
#' @export
new_muvg_category_representation_template_from_data <- function(
  data,
  category = "category",
  cues,
  verbose = FALSE
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
  if (verbose) {
    message(
      "Constructed a MUVG category-representation template with ",
      length(category_labels), " categories and ",
      length(cues), " cue(s)."
    )
  }
  new_category_representation_template(representations)
}

#' @rdname family-muvg
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
  noise_treatment = "no_noise",
  lapse_treatment = "no_lapses",
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
  noise_treatment = "no_noise",
  lapse_treatment = "no_lapses",
  verbose = FALSE
) {
  template <- new_muvg_category_representation_template_from_data(
    data,
    category = category,
    cues = cues,
    verbose = verbose
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
