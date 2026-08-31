#' @include S7-core-classes.R
NULL

# -------------------------
# Multivariate Gaussian (MVG) Family
# -------------------------

#' Multivariate Gaussian (MVG) Family: Representations, Templates, and Models
#'
#' Constructs and manages multivariate Gaussian category representations,
#' templates, and ideal observer models.
#'
#' @name family-mvg
#' @rdname family-mvg
#' @param category_labels Character vector of category label(s).
#' @param cue_labels Character vector of cue labels.
#' @param mu Numeric vector of category means.
#' @param Sigma Numeric square symmetric positive-definite covariance matrix.
#' @param data Data frame containing category and cue observations.
#' @param category_var Bare symbol or character string indicating the category
#'   column in `data`.
#' @param cue_vars Bare tidyselect specification or character vector of cue
#'   columns in `data`.
#' @param category_template An [MVBU_CategoryRepresentationTemplate] object
#'   containing MVG category representations.
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
#' The multivariate Gaussian (`MVG`) family models category structure over
#' \eqn{D}-dimensional continuous cue spaces with correlated cue dimensions:
#' \deqn{p(\mathbf{x} \mid c) = \mathcal{N}(\mathbf{x}; \boldsymbol{\mu}_c,
#' \mathbf{\Sigma}_c)}
#'
#' @return An S7 object of the respective MVG class.
#' @seealso [family-niw], [family-uvg], [family-exemplar]
#' @export
MVG_CategoryRepresentation <- S7::new_class(
  "MVG_CategoryRepresentation",
  parent = MVBU_CategoryRepresentation,
  properties = list(
    mu = S7::class_numeric,
    Sigma = S7::class_any
  ),
  validator = function(self) {
    label_information <- .mvbu_label_information(self@metadata)
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
  MVG_CategoryRepresentation(
    category_likelihood_function = {
      mu0 <- as.numeric(mu)
      Sigma0 <- as.matrix(Sigma)
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
        .dmvnorm_density(x, mean = mu0, Sigma = Sigma_eff, log = log)
      }
    },
    metadata = .mvbu_label_metadata(
      as.character(category_labels),
      as.character(cue_labels),
      metadata
    ),
    mu = as.numeric(mu),
    Sigma = as.matrix(Sigma)
  )
}

#' @rdname family-mvg
#' @export
new_mvg_category_representation_from_data <- function(
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
  cues,
  verbose = FALSE
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

  if (verbose) {
    message(
      "Constructed an MVG category-representation template with ",
      length(category_labels), " categories and ",
      length(cues), " cue(s)."
    )
  }

  new_category_representation_template(representations)
}

#' @rdname family-mvg
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
  noise_treatment = "no_noise",
  lapse_treatment = "no_lapses",
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
  noise_treatment = "no_noise",
  lapse_treatment = "no_lapses",
  verbose = FALSE
) {
  template <- new_mvg_category_representation_template_from_data(
    data,
    category = category,
    cues = cues,
    verbose = verbose
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
