#' @include S7-class.R
NULL

# -------------------------
# Normal-Inverse-Wishart (NIW) Family
# -------------------------

#' Normal-Inverse-Wishart (NIW) Family: Representations, Templates, and Models
#'
#' Constructs and manages Normal-Inverse-Wishart (NIW) category
#' representations, templates, and ideal adaptor cognitive models.
#'
#' @param category_labels Character vector of category label(s).
#' @param cue_labels Character vector of cue labels.
#' @param m Numeric vector of location hyperparameters.
#' @param kappa Positive numeric scalar precision scaling (\eqn{\kappa > 0}).
#' @param nu Numeric scalar degrees of freedom (\eqn{\nu > D - 1}).
#' @param S Numeric square symmetric scale matrix.
#' @param data Data frame containing category and cue observations.
#' @param category Character string indicating the category column in `data`. Defaults to `"category"`.
#' @param cues Character vector indicating cue columns in `data`.
#' @param category_template An [MVBU_CategoryRepresentationTemplate] object
#'   containing NIW category representations.
#' @param decision_rule Categorization decision rule: `"sampling"`,
#'   `"criterion"`, or `"proportional"`. Defaults to `"sampling"`.
#' @param category_prior Optional numeric vector of prior category
#'   probabilities summing to 1. Defaults to equal priors.
#' @param lapse_rate Numeric scalar lapse probability in `[0, 1]`. Defaults
#'   to 0.
#' @param lapse_bias Optional numeric vector of lapse category probabilities
#'   summing to 1.
#' @param Sigma_noise Optional perceptual/measurement noise covariance matrix.
#' @param noise_treatment Treatment of perceptual noise during evaluation and
#'   updating: `"no_noise"`, `"sample"`, or `"marginalize"`.
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
#' @param lapse_treatment Treatment of attentional lapses during evaluation and
#'   updating: `"no_lapses"`, `"sample"`, or `"marginalize"`.
#'   \itemize{
#'     \item \code{"no_lapses"}: Lapses are ignored; categorization and
#'       updating reflect pure normative evidence evaluation (even if
#'       \code{lapse_rate > 0}).
#'     \item \code{"sample"}: Simulates stochastic attentional lapses on
#'       individual observations. With probability \code{lapse_rate}, a lapse
#'       occurs: during categorization, the posterior is replaced by
#'       \code{lapse_bias}; during updating, the lapsing observation is treated
#'       as uninformative.
#'     \item \code{"marginalize"}: Evaluates categorization in the limit across
#'       infinitely many encounters (expected response probabilities). The
#'       category posterior is the mixture
#'       \eqn{(1 - \lambda) p(c \mid x) + \lambda b_c}, and evidence during
#'       updating is downweighted by \eqn{1 - \lambda}.
#'   }
#'   Defaults to `.infer_lapse_treatment(lapse_rate)`: `"marginalize"` if
#'   `lapse_rate > 0`, `"no_lapses"` otherwise.
#' @param metadata Optional list of metadata.
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
#'
#' @name family-niw
#' @rdname family-niw
#' @docType class
#' @usage NULL
#' @export
#' @examples
#' # Normal-Inverse-Wishart (NIW) representation for multidimensional cues
#' m_b <- c(VOT = 10, f0 = 200)
#' S_b <- matrix(c(25, 5, 5, 36), nrow = 2,
#'   dimnames = list(names(m_b), names(m_b)))
#' rep_b <- new_niw_category_representation(
#'   category_labels = "/b/", cue_labels = c("VOT", "f0"),
#'   m = m_b, kappa = 5, nu = 10, S = S_b
#' )
#' m_d <- c(VOT = 50, f0 = 240)
#' S_d <- matrix(c(25, -2, -2, 36), nrow = 2,
#'   dimnames = list(names(m_d), names(m_d)))
#' rep_d <- new_niw_category_representation(
#'   category_labels = "/d/", cue_labels = c("VOT", "f0"),
#'   m = m_d, kappa = 5, nu = 10, S = S_d
#' )
#' tmpl <- new_category_representation_template(list(rep_b, rep_d))
#' m_niw <- new_niw_ideal_adaptor(category_template = tmpl)
#' m_niw
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
    d <- length(self@metadata$label_information$cue)

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
  # Rationale for retaining typed class constructor:
  # Concrete typed constructors provide strict property validation, IDE tab completion,
  # and fast S7 method dispatch, while obtaining the likelihood closure uniformly
  # from the registered model family density function.
  metadata <- as.list(metadata)
  if (is.null(metadata$parameters)) {
    metadata$parameters <- list(
      m = as.numeric(m),
      kappa = as.numeric(kappa),
      nu = as.numeric(nu),
      S = as.matrix(S)
    )
  }
  rep <- NIW_CategoryRepresentation(
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
    S = as.matrix(S)
  )
  rep@category_likelihood_function <- .make_family_likelihood_function(rep)
  rep
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

  mvg_template <- new_mvg_category_representation_template_from_data(
    data,
    category = category,
    cues = cues
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

  new_category_representation_template(representations)
}

#' @rdname family-niw
#' @docType class
#' @usage NULL
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
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate),
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
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = .infer_lapse_treatment(lapse_rate)
) {
  template <- new_niw_category_representation_template_from_data(
    data,
    category = category,
    cues = cues,
    kappa = kappa,
    nu = nu
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
