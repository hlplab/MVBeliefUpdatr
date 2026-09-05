#' @include S7-core-classes.R
NULL

# -------------------------
# Exemplar Model Family
# -------------------------

#' Exemplar Model Family: Representations, Templates, and Models
#'
#' Constructs and manages non-parametric exemplar-based category
#' representations, templates, and cognitive models,
#' providing a direct mathematical bridge between Kernel Density
#' Estimation (KDE) and the Generalized Context Model (GCM;
#' \insertCite{nosofsky1986,ashby1995}{MVBeliefUpdatr}).
#'
#' @name family-exemplar
#' @rdname family-exemplar
#' @param category_labels Character vector of category label(s).
#' @param cue_labels Character vector of cue dimension label(s).
#' @param exemplars Numeric matrix or data frame of exemplars (rows are
#'   exemplars; columns correspond to `cue_labels`).
#' @param exemplar_weights Optional numeric vector of exemplar weights in
#'   `[0, 1]` summing to 1. Defaults to equal weights (`1 / n_exemplars`).
#' @param c Optional positive numeric scalar specifying the sensitivity /
#'   inverse kernel bandwidth scale. If `NULL` (default), `c` is estimated
#'   adaptively via Silverman's (1986) multivariate rule of thumb:
#'   \deqn{h^2 = \left(\frac{4}{(d + 2) n}\right)^{\frac{2}{d + 4}}, \quad
#'   c = \frac{1}{2 h^2}}
#' @param data Data frame containing category and cue observations.
#' @param category_var Bare symbol or character string indicating the category
#'   column in `data`.
#' @param cue_vars Bare tidyselect specification or character vector of cue
#'   columns in `data`.
#' @param category_template An [MVBU_CategoryRepresentationTemplate] object
#'   containing Exemplar category representations.
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
#' ## Relation to the Generalized Context Model (GCM)
#' In the Generalized Context Model
#' \insertCite{nosofsky1986,kruschke1992}{MVBeliefUpdatr},
#' psychological similarity between an observed stimulus \eqn{\mathbf{x}} and
#' a stored exemplar \eqn{\mathbf{e}_i} is given by:
#' \deqn{\operatorname{sim}(\mathbf{x}, \mathbf{e}_i) =
#' \exp\left(-c \cdot d(\mathbf{x}, \mathbf{e}_i)^p\right)}
#' where \eqn{c > 0} is the perceptual sensitivity parameter and
#' \eqn{d(\mathbf{x}, \mathbf{e}_i)} is the attention-weighted Minkowski
#' distance:
#' \deqn{d(\mathbf{x}, \mathbf{e}_i) =
#' \left(\sum_{k=1}^d w_k |x_k - e_{i,k}|^\tau\right)^{\frac{1}{\tau}}}
#' with attention weights \eqn{w_k} (\eqn{\sum w_k = 1}).
#'
#' Under the Luce choice rule with category response biases \eqn{b_C}, the
#' categorization probability in GCM is:
#' \deqn{P(C \mid \mathbf{x}) =
#' \frac{b_C \sum_{j \in C} \operatorname{sim}(\mathbf{x}, \mathbf{e}_{C,j})}{
#' \sum_K b_K \sum_{j \in K} \operatorname{sim}(\mathbf{x}, \mathbf{e}_{K,j})}}
#'
#' ## Density Estimation (KDE) Formulation
#' In `MVBeliefUpdatr`, category representations are framed in terms of
#' probability density estimation
#' \insertCite{ashby1995,shi2010}{MVBeliefUpdatr}:
#' \deqn{p(\mathbf{x} \mid C) =
#' \sum_{i=1}^n \pi_{C,i} \, \mathcal{N}\left(\mathbf{x} \mid
#' \mathbf{e}_{C,i}, \mathbf{\Sigma}_{\text{kernel}}\right)}
#'
#' When the kernel covariance is matched to GCM sensitivity and attention:
#' \deqn{\mathbf{\Sigma}_{\text{kernel}}^{-1} = 2c \mathbf{W} \iff
#' \mathbf{\Sigma}_{\text{kernel}} = \frac{1}{2c} \mathbf{W}^{-1}}
#' the Gaussian likelihood is directly proportional to GCM similarity:
#' \deqn{\mathcal{N}\left(\mathbf{x} \mid \mathbf{e}_i,
#' \mathbf{\Sigma}_{\text{kernel}}\right) \propto
#' \operatorname{sim}(\mathbf{x}, \mathbf{e}_i)}
#' When plugged into Bayes' rule \eqn{P(C \mid \mathbf{x}) =
#' \frac{P(C) p(\mathbf{x} \mid C)}{\sum_K P(K) p(\mathbf{x} \mid K)}},
#' the normalizing constant cancels out, resulting in an exact mathematical
#' equivalence to the GCM \insertCite{ashby1995}{MVBeliefUpdatr}.
#'
#' ## Implicit Assumptions on \eqn{\tau} and \eqn{p}
#' This density-based formulation corresponds to the GCM under two specific
#' geometric assumptions:
#' 1. **Euclidean distance metric (\eqn{\tau = 2})**: Assumes integral
#'    (holistic) perceptual dimensions where dimensions combine Euclidean-wise
#'    \insertCite{shepard1987,nosofsky1986}{MVBeliefUpdatr}.
#' 2. **Gaussian similarity gradient (\eqn{p = 2})**: Assumes a Gaussian
#'    similarity decay profile \eqn{\exp(-c \cdot d^2)} rather than an
#'    exponential profile \eqn{\exp(-c \cdot d)}, characteristic of integral
#'    dimensions and Gaussian kernel density estimators.
#'
#' In addition, `MVBeliefUpdatr` generalizes standard GCM by allowing a full
#' covariance matrix
#' \eqn{\mathbf{\Sigma}_{\text{kernel}} =
#' \frac{1}{2c}\operatorname{Cov}(\mathbf{X})},
#' accommodating correlated cue dimensions (e.g., acoustic formant cues).
#'
#' ## Bandwidth Scale \eqn{c} and Silverman's Rule
#' If `c` is not provided (`NULL`), `c` defaults to Silverman's (1986)
#' multivariate rule-of-thumb bandwidth factor:
#' \deqn{h^2 = \left(\frac{4}{(d + 2) n}\right)^{\frac{2}{d + 4}}, \quad
#' c = \frac{1}{2 h^2}}
#' which asymptotically minimizes the Mean Integrated Squared Error (MISE)
#' to the true category distribution.
#'
#' @return An S7 object of class [Exemplar_CategoryRepresentation].
#' @seealso [family-mvg], [family-uvg]
#' @references
#' \insertRef{nosofsky1986}{MVBeliefUpdatr}
#'
#' \insertRef{shepard1987}{MVBeliefUpdatr}
#'
#' \insertRef{ashby1995}{MVBeliefUpdatr}
#'
#' \insertRef{kruschke1992}{MVBeliefUpdatr}
#'
#' \insertRef{shi2010}{MVBeliefUpdatr}
#'
#' \insertRef{silverman1986}{MVBeliefUpdatr}
#' @export
Exemplar_CategoryRepresentation <- S7::new_class(
  "Exemplar_CategoryRepresentation",
  parent = MVBU_CategoryRepresentation,
  properties = list(
    exemplars = S7::class_any,
    exemplar_weights = S7::class_numeric,
    c = S7::class_numeric
  ),
  validator = function(self) {
    if (!is.matrix(self@exemplars) || !is.numeric(self@exemplars)) {
      return("Exemplar exemplars must be a numeric matrix.")
    }

    n_ex <- nrow(self@exemplars)
    d <- ncol(self@exemplars)
    if (n_ex < 1) {
      return("Exemplar exemplars must contain at least one row.")
    }
    label_information <- .mvbu_label_information(self@metadata)
    if (d != length(label_information$cue)) {
      return("Exemplar exemplar column count must match cue dimensionality.")
    }
    if (length(self@exemplar_weights) != n_ex) {
      return(
        "Exemplar exemplar_weights length must match number of exemplars."
      )
    }
    if (any(self@exemplar_weights < 0) || any(self@exemplar_weights > 1)) {
      return("Exemplar exemplar_weights entries must be in [0, 1].")
    }
    if (abs(sum(self@exemplar_weights) - 1) > MVBU_PROB_TOL) {
      return("Exemplar exemplar_weights entries must sum to 1.")
    }
    if (
      length(self@c) != 1L ||
      !is.numeric(self@c) ||
      is.na(self@c) ||
      self@c <= 0
    ) {
      return("Exemplar c must be a single positive numeric value.")
    }

    NULL
  }
)

#' @rdname family-exemplar
#' @export
new_exemplar_category_representation <- function(
  category_labels,
  cue_labels,
  exemplars,
  exemplar_weights = NULL,
  c = NULL,
  metadata = list()
) {
  exemplars <- as.matrix(exemplars)
  if (!is.numeric(exemplars)) {
    .stop("exemplars must be numeric.")
  }
  n_ex <- nrow(exemplars)
  if (is.null(exemplar_weights)) {
    exemplar_weights <- rep(1 / n_ex, n_ex)
  }
  if (!is.numeric(exemplar_weights)) {
    .stop("exemplar_weights must be numeric.")
  }

  d0 <- ncol(exemplars)
  n0 <- nrow(exemplars)

  if (!is.null(c)) {
    .assert_true(
      is.numeric(c) && length(c) == 1L && !is.na(c) && c > 0,
      msg = "c must be a single positive numeric scalar."
    )
    c_val <- as.numeric(c)
    bandwidth_scale <- 1 / (2 * c_val)
  } else {
    if (n0 > 1) {
      bandwidth_scale <- (4 / ((d0 + 2) * n0))^(2 / (d0 + 4))
    } else {
      bandwidth_scale <- 1.0
    }
    c_val <- 1 / (2 * bandwidth_scale)
  }

  Exemplar_CategoryRepresentation(
    category_likelihood_function = {
      ex0 <- exemplars
      w0 <- as.numeric(exemplar_weights)
      if (n0 > 1) {
        Sigma0 <- stats::cov(ex0) * bandwidth_scale
      } else {
        Sigma0 <- diag(1, d0) * bandwidth_scale
      }
      if (!is.matrix(Sigma0) || any(!is.finite(Sigma0))) {
        Sigma0 <- diag(1, d0) * bandwidth_scale
      }
      Sigma0 <- Sigma0 + diag(MVBU_PROB_TOL, d0)

      function(
        x,
        log = FALSE,
        noise_treatment = "no_noise",
        Sigma_noise = NULL
      ) {
        x <- .as_observation_matrix(x, d = d0, arg_name = "x")
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(
            n = nrow(x),
            mean = rep(0, ncol(x)),
            sigma = Sigma_noise
          )
        }
        if (
          !is.null(Sigma_noise) &&
          (identical(noise_treatment, "sample") ||
           identical(noise_treatment, "marginalize"))
        ) {
          Sigma_eff <- Sigma0 + Sigma_noise
        } else {
          Sigma_eff <- Sigma0
        }

        # Vectorized evaluation across all exemplars:
        # Pre-factorize Sigma_eff via Cholesky decomposition once rather than N*M times.
        # Compute squared Mahalanobis distances using BLAS matrix multiplication.
        chol_sigma <- tryCatch(chol(Sigma_eff), error = function(e) NULL)
        if (is.null(chol_sigma)) {
          chol_sigma <- chol(Sigma_eff + diag(MVBU_PROB_TOL, d0))
        }
        log_det <- 2 * sum(log(diag(chol_sigma)))

        # Transform observations and exemplars into standardized coordinate space
        x_std <- t(backsolve(chol_sigma, t(x), transpose = TRUE))
        ex_std <- t(backsolve(chol_sigma, t(ex0), transpose = TRUE))

        sq_x <- rowSums(x_std^2)
        sq_ex <- rowSums(ex_std^2)
        # ||x_std - ex_std||^2 = ||x_std||^2 + ||ex_std||^2 - 2 * x_std %*% ex_std^T
        dist_mat <- outer(sq_x, sq_ex, "+") - 2 * tcrossprod(x_std, ex_std)
        dist_mat[dist_mat < 0] <- 0

        log_dens_by_exemplar <- -0.5 * (d0 * log(2 * pi) + log_det + dist_mat)
        weighted_logdens <- sweep(log_dens_by_exemplar, 2, log(w0), "+")
        log_mix <- .logsumexp_rows(weighted_logdens)
        if (isTRUE(log)) {
          log_mix
        } else {
          exp(log_mix)
        }
      }
    },
    metadata = .mvbu_label_metadata(
      as.character(category_labels),
      as.character(cue_labels),
      metadata
    ),
    exemplars = exemplars,
    exemplar_weights = as.numeric(exemplar_weights),
    c = as.numeric(c_val)
  )
}

#' @rdname family-exemplar
#' @export
new_exemplar_category_representation_from_data <- function(
  data,
  category = "category",
  cues,
  c = NULL
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

  new_exemplar_category_representation(
    category_labels = category_labels,
    cue_labels = cues,
    exemplars = as.matrix(data[, cues, drop = FALSE]),
    c = c
  )
}

#' @rdname family-exemplar
#' @export
new_exemplar_category_representation_template_from_data <- function(
  data,
  category = "category",
  cues,
  c = NULL,
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
    representations[[label]] <- new_exemplar_category_representation_from_data(
      data[data[[category]] == label, , drop = FALSE],
      category = category,
      cues = cues,
      c = c
    )
  }

  if (verbose) {
    message(
      "Constructed an exemplar category-representation template with ",
      length(category_labels), " categories and ", length(cues), " cue(s)."
    )
  }

  new_category_representation_template(representations)
}

#' @rdname family-exemplar
#' @export
Exemplar_Model <- S7::new_class(
  "Exemplar_Model",
  parent = MVBU_CognitiveModel
)

#' @rdname family-exemplar
#' @export
new_exemplar_model <- function(
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
    category_representation_class = Exemplar_CategoryRepresentation,
    model_class = Exemplar_Model,
    family_label = "EXEMPLAR",
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

#' @rdname family-exemplar
#' @export
new_exemplar_model_from_data <- function(
  data,
  category = "category",
  cues,
  c = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = "no_noise",
  lapse_treatment = "no_lapses",
  verbose = FALSE
) {
  template <- new_exemplar_category_representation_template_from_data(
    data,
    category = category,
    cues = cues,
    c = c,
    verbose = verbose
  )
  new_exemplar_model(
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
