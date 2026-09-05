#' @include S7-family-coercion.R
#' @include S7-core-uvg-classes.R
#' @include S7-core-nix-classes.R
#' @include S7-core-muvg-classes.R
#' @include S7-core-mnix-classes.R
#' @include S7-core-mvg-classes.R
#' @include S7-core-niw-classes.R
#' @include S7-core-exemplar-classes.R
NULL

# -----------------------------------------------------------------------------
# Dynamic Object Construction from Data
# -----------------------------------------------------------------------------

#' Construct S7 Category Representations, Templates, and Models from Data
#'
#' Dynamic constructors that inspect the input `type` (or `family`) argument
#' and dispatch to the appropriate family-specific constructor from data.
#'
#' @name new-objects-from-data
#' @rdname new-objects-from-data
#' @param data A data frame containing observation cues and category labels.
#'   For individual category representations, data must contain observations
#'   from exactly one category.
#' @param type Family type string: `"UVG"`, `"NIX"`, `"MUVG"`, `"MNIX"`,
#'   `"MVG"`, `"NIW"`, or `"EXEMPLAR"`.
#' @param category Character string giving the category column name in `data`.
#'   Defaults to `"category"`.
#' @param cues Character vector of cue column names in `data`.
#' @param category_prior Optional numeric vector of prior category
#'   probabilities summing to 1.
#' @param decision_rule Categorization decision rule: `"sampling"` or
#'   `"argmax"`. Defaults to `"sampling"`.
#' @param lapse_rate Numeric scalar lapse probability in `[0, 1]`. Defaults
#'   to 0.
#' @param lapse_bias Optional numeric vector of lapse category probabilities.
#' @param Sigma_noise Optional perceptual noise covariance matrix.
#' @param noise_treatment Treatment of noise: `"no_noise"`, `"sample"`, or
#'   `"marginalize"`.
#' @param lapse_treatment Treatment of lapses: `"no_lapses"`, `"sample"`, or
#'   `"marginalize"`.
#' @param ... Additional family-specific parameters (such as `kappa`, `nu`,
#'   `m_0`, `S_0`, `c`, or `bandwidth`).
#'
#' @return An S7 representation, template, or cognitive model object.
#' @seealso [family-uvg], [family-nix], [family-muvg], [family-mnix],
#'   [family-mvg], [family-niw], [family-exemplar]
#' @export
new_category_representation_from_data <- function(
  data,
  type,
  category = "category",
  cues,
  ...
) {
  type <- toupper(type)
  constructor <- switch(
    type,
    UVG = new_uvg_category_representation_from_data,
    NIX = new_nix_category_representation_from_data,
    MUVG = new_muvg_category_representation_from_data,
    MNIX = new_mnix_category_representation_from_data,
    MVG = new_mvg_category_representation_from_data,
    NIW = new_niw_category_representation_from_data,
    EXEMPLAR = new_exemplar_category_representation_from_data,
    .stop("type must be one of UVG, NIX, MUVG, MNIX, MVG, NIW, or EXEMPLAR.")
  )
  constructor(data = data, category = category, cues = cues, ...)
}

#' @rdname new-objects-from-data
#' @export
new_category_representation_template_from_data <- function(
  data,
  type,
  category = "category",
  cues,
  ...
) {
  type <- toupper(type)
  constructor <- switch(
    type,
    UVG = new_uvg_category_representation_template_from_data,
    NIX = new_nix_category_representation_template_from_data,
    MUVG = new_muvg_category_representation_template_from_data,
    MNIX = new_mnix_category_representation_template_from_data,
    MVG = new_mvg_category_representation_template_from_data,
    NIW = new_niw_category_representation_template_from_data,
    EXEMPLAR = new_exemplar_category_representation_template_from_data,
    .stop("type must be one of UVG, NIX, MUVG, MNIX, MVG, NIW, or EXEMPLAR.")
  )
  constructor(data = data, category = category, cues = cues, ...)
}

#' @rdname new-objects-from-data
#' @export
new_model_from_data <- function(
  data,
  type,
  category = "category",
  cues,
  ...
) {
  type <- toupper(type)
  constructor <- switch(
    type,
    UVG = new_uvg_ideal_observer_from_data,
    NIX = new_nix_ideal_adaptor_from_data,
    MUVG = new_muvg_ideal_observer_from_data,
    MNIX = new_mnix_ideal_adaptor_from_data,
    MVG = new_mvg_ideal_observer_from_data,
    NIW = new_niw_ideal_adaptor_from_data,
    EXEMPLAR = new_exemplar_model_from_data,
    .stop("type must be one of UVG, NIX, MUVG, MNIX, MVG, NIW, or EXEMPLAR.")
  )
  constructor(data = data, category = category, cues = cues, ...)
}

# -----------------------------------------------------------------------------
# =============================================================================
# Aggregation of Representations, Templates, and Cognitive Models
# =============================================================================

.aggregate_weights <- function(n, weights = NULL) {
  if (is.null(weights)) {
    return(rep(1 / n, n))
  }
  .assert_true(
    is.numeric(weights) && length(weights) == n &&
      all(weights >= 0) && sum(weights) > 0,
    msg = paste0(
      "weights must be a numeric vector of positive weights ",
      "with the same length as the number of objects being aggregated."
    )
  )
  weights / sum(weights)
}

.aggregate_category_representations <- function(reps, weights = NULL) {
  .assert_true(
    is.list(reps) && length(reps) > 0,
    msg = "reps must be a non-empty list of category representation objects."
  )
  first_rep <- reps[[1]]
  .assert_true(
    S7::S7_inherits(first_rep, MVBU_CategoryRepresentation),
    msg = "All elements in reps must be S7 MVBU_CategoryRepresentation objects."
  )

  target_class <- S7::S7_class(first_rep)
  for (i in seq_along(reps)) {
    .assert_true(
      S7::S7_inherits(reps[[i]], target_class),
      msg = paste0(
        "All representations must have the same S7 class (expected ",
        target_class@name, ")."
      )
    )
  }

  n_reps <- length(reps)
  weights <- .aggregate_weights(n_reps, weights)

  if (n_reps == 1) {
    return(first_rep)
  }

  cat_labels <- get_category_labels(first_rep)
  cue_labels <- get_cue_labels(first_rep)

  for (i in seq_along(reps)) {
    .assert_true(
      identical(get_category_labels(reps[[i]]), cat_labels),
      msg = "All representations must have identical category labels."
    )
    .assert_true(
      identical(get_cue_labels(reps[[i]]), cue_labels),
      msg = "All representations must have identical cue labels."
    )
  }

  if (S7::S7_inherits(first_rep, MVG_CategoryRepresentation)) {
    mu_weighted <- Reduce(
      `+`,
      mapply(function(r, w) r@mu * w, reps, weights, SIMPLIFY = FALSE)
    )
    Sigma_weighted <- Reduce(
      `+`,
      mapply(function(r, w) r@Sigma * w, reps, weights, SIMPLIFY = FALSE)
    )
    new_mvg_category_representation(
      category_labels = cat_labels,
      cue_labels = cue_labels,
      mu = mu_weighted,
      Sigma = Sigma_weighted
    )
  } else if (S7::S7_inherits(first_rep, UVG_CategoryRepresentation)) {
    mu_weighted <- sum(
      vapply(reps, function(r) r@mu, numeric(1)) * weights
    )
    sigma2_weighted <- sum(
      vapply(reps, function(r) r@sigma2, numeric(1)) * weights
    )
    new_uvg_category_representation(
      category_labels = cat_labels,
      cue_labels = cue_labels,
      mu = mu_weighted,
      sigma2 = sigma2_weighted
    )
  } else if (S7::S7_inherits(first_rep, MUVG_CategoryRepresentation)) {
    mu_weighted <- Reduce(
      `+`,
      mapply(function(r, w) r@mu * w, reps, weights, SIMPLIFY = FALSE)
    )
    sigma2_weighted <- Reduce(
      `+`,
      mapply(function(r, w) r@sigma2 * w, reps, weights, SIMPLIFY = FALSE)
    )
    new_muvg_category_representation(
      category_labels = cat_labels,
      cue_labels = cue_labels,
      mu = mu_weighted,
      sigma2 = sigma2_weighted
    )
  } else if (S7::S7_inherits(first_rep, NIW_CategoryRepresentation)) {
    m_weighted <- Reduce(
      `+`,
      mapply(function(r, w) r@m * w, reps, weights, SIMPLIFY = FALSE)
    )
    S_weighted <- Reduce(
      `+`,
      mapply(function(r, w) r@S * w, reps, weights, SIMPLIFY = FALSE)
    )
    kappas <- vapply(reps, function(r) r@kappa, numeric(1))
    nus <- vapply(reps, function(r) r@nu, numeric(1))
    new_niw_category_representation(
      category_labels = cat_labels,
      cue_labels = cue_labels,
      m = m_weighted,
      S = S_weighted,
      kappa = sum(kappas * weights),
      nu = sum(nus * weights)
    )
  } else if (S7::S7_inherits(first_rep, NIX_CategoryRepresentation)) {
    m_weighted <- sum(
      vapply(reps, function(r) r@m, numeric(1)) * weights
    )
    sigma2_weighted <- sum(
      vapply(reps, function(r) r@sigma2, numeric(1)) * weights
    )
    kappas <- vapply(reps, function(r) r@kappa, numeric(1))
    nus <- vapply(reps, function(r) r@nu, numeric(1))
    new_nix_category_representation(
      category_labels = cat_labels,
      cue_labels = cue_labels,
      m = m_weighted,
      kappa = sum(kappas * weights),
      nu = sum(nus * weights),
      sigma2 = sigma2_weighted
    )
  } else if (S7::S7_inherits(first_rep, MNIX_CategoryRepresentation)) {
    m_weighted <- Reduce(
      `+`,
      mapply(function(r, w) r@m * w, reps, weights, SIMPLIFY = FALSE)
    )
    S_weighted <- Reduce(
      `+`,
      mapply(function(r, w) r@S * w, reps, weights, SIMPLIFY = FALSE)
    )
    kappas <- Reduce(
      `+`,
      mapply(function(r, w) r@kappa * w, reps, weights, SIMPLIFY = FALSE)
    )
    nus <- Reduce(
      `+`,
      mapply(function(r, w) r@nu * w, reps, weights, SIMPLIFY = FALSE)
    )
    new_mnix_category_representation(
      category_labels = cat_labels,
      cue_labels = cue_labels,
      m = m_weighted,
      S = S_weighted,
      kappa = kappas,
      nu = nus
    )
  } else if (S7::S7_inherits(first_rep, Exemplar_CategoryRepresentation)) {
    all_exemplars <- do.call(rbind, lapply(reps, function(r) r@exemplars))
    all_weights <- unlist(
      mapply(
        function(r, w) r@exemplar_weights * w,
        reps,
        weights,
        SIMPLIFY = FALSE
      )
    )
    all_weights <- all_weights / sum(all_weights)
    c_vals <- vapply(reps, function(r) r@c, numeric(1))
    c_weighted <- sum(c_vals * weights)
    new_exemplar_category_representation(
      category_labels = cat_labels,
      cue_labels = cue_labels,
      exemplars = all_exemplars,
      exemplar_weights = all_weights,
      c = c_weighted
    )
  } else {
    .stop(
      "Aggregation is not implemented for category representation class: ",
      target_class@name
    )
  }
}

.aggregate_category_representation_templates <- function(templates, weights = NULL) {
  .assert_true(
    is.list(templates) && length(templates) > 0,
    msg = "templates must be a non-empty list of category representation template objects."
  )
  first_template <- templates[[1]]
  .assert_true(
    S7::S7_inherits(first_template, MVBU_CategoryRepresentationTemplate),
    msg = "All elements in templates must be S7 MVBU_CategoryRepresentationTemplate objects."
  )

  for (i in seq_along(templates)) {
    .assert_true(
      S7::S7_inherits(templates[[i]], MVBU_CategoryRepresentationTemplate),
      msg = "All elements in templates must be S7 MVBU_CategoryRepresentationTemplate objects."
    )
  }

  n_templates <- length(templates)
  weights <- .aggregate_weights(n_templates, weights)

  if (n_templates == 1) {
    return(first_template)
  }

  first_reps <- first_template@representations
  categories <- names(first_reps)

  for (i in seq_along(templates)) {
    t_cats <- names(templates[[i]]@representations)
    .assert_true(
      identical(sort(t_cats), sort(categories)),
      msg = "All templates must have identical category labels."
    )
  }

  agg_reps <- list()
  for (cat in categories) {
    cat_reps <- lapply(templates, function(t) t@representations[[cat]])
    agg_reps[[cat]] <- .aggregate_category_representations(cat_reps, weights = weights)
  }

  new_category_representation_template(representations = agg_reps)
}

.aggregate_cognitive_models <- function(models, weights = NULL) {
  .assert_true(
    is.list(models) && length(models) > 0,
    msg = "models must be a non-empty list of cognitive model objects."
  )

  first_model <- models[[1]]
  .assert_true(
    S7::S7_inherits(first_model, MVBU_CognitiveModel),
    msg = "All elements in models must be S7 MVBU_CognitiveModel objects."
  )

  target_class <- S7::S7_class(first_model)
  for (i in seq_along(models)) {
    .assert_true(
      S7::S7_inherits(models[[i]], target_class),
      msg = paste0(
        "All elements in models must have the same S7 class (expected ",
        target_class@name, ")."
      )
    )
  }

  n_models <- length(models)
  weights <- .aggregate_weights(n_models, weights)

  if (n_models == 1) {
    return(first_model)
  }

  # 1. Aggregate model-level parameters
  # category_prior
  priors_list <- lapply(models, function(m) get_category_prior(m))
  first_prior <- priors_list[[1]]
  agg_category_prior <- if (!is.null(first_prior)) {
    prior_names <- names(first_prior)
    weighted_prior <- Reduce(
      `+`,
      mapply(function(p, w) p * w, priors_list, weights, SIMPLIFY = FALSE)
    )
    if (!is.null(prior_names)) names(weighted_prior) <- prior_names
    weighted_prior
  } else {
    NULL
  }

  # lapse_rate
  lapse_rates <- vapply(models, function(m) get_lapse_rate(m), numeric(1))
  agg_lapse_rate <- sum(lapse_rates * weights)

  # lapse_bias
  lapse_biases_list <- lapply(models, function(m) get_lapse_bias(m))
  agg_lapse_bias <- if (!is.null(lapse_biases_list[[1]])) {
    lb_names <- names(lapse_biases_list[[1]])
    weighted_lb <- Reduce(
      `+`,
      mapply(
        function(lb, w) lb * w,
        lapse_biases_list,
        weights,
        SIMPLIFY = FALSE
      )
    )
    if (!is.null(lb_names)) names(weighted_lb) <- lb_names
    weighted_lb
  } else {
    NULL
  }

  # Sigma_noise
  sigma_noise_list <- lapply(models, function(m) m@noise_behavior$Sigma_noise)
  agg_sigma_noise <- if (!is.null(sigma_noise_list[[1]])) {
    Reduce(
      `+`,
      mapply(
        function(sn, w) sn * w,
        sigma_noise_list,
        weights,
        SIMPLIFY = FALSE
      )
    )
  } else {
    NULL
  }

  # 2. Aggregate representations / category templates
  templates <- lapply(models, function(m) m@category_template)
  agg_template <- .aggregate_category_representation_templates(templates, weights = weights)

  # 3. Construct aggregated model
  args <- list(
    category_template = agg_template,
    category_prior = agg_category_prior,
    lapse_rate = agg_lapse_rate,
    lapse_bias = agg_lapse_bias,
    Sigma_noise = agg_sigma_noise,
    decision_rule = first_model@decision_rule,
    noise_treatment = get_noise_treatment(first_model),
    lapse_treatment = get_lapse_treatment(first_model),
    metadata = first_model@metadata
  )

  constructor <- if (S7::S7_inherits(first_model, MVG_IdealObserver)) {
    new_mvg_ideal_observer
  } else if (S7::S7_inherits(first_model, NIW_IdealAdaptor)) {
    new_niw_ideal_adaptor
  } else if (S7::S7_inherits(first_model, UVG_IdealObserver)) {
    new_uvg_ideal_observer
  } else if (S7::S7_inherits(first_model, NIX_IdealAdaptor)) {
    new_nix_ideal_adaptor
  } else if (S7::S7_inherits(first_model, MUVG_IdealObserver)) {
    new_muvg_ideal_observer
  } else if (S7::S7_inherits(first_model, MNIX_IdealAdaptor)) {
    new_mnix_ideal_adaptor
  } else if (S7::S7_inherits(first_model, Exemplar_Model)) {
    new_exemplar_model
  } else {
    .stop("Aggregation is not supported for model class: ", target_class@name)
  }

  do.call(constructor, args)
}

# -----------------------------------------------------------------------------
# S7 Methods for aggregate
# -----------------------------------------------------------------------------

S7::method(aggregate, MVBU_CategoryRepresentation) <- function(x, ..., weights = NULL) {
  dots <- list(...)
  reps <- c(list(x), dots)
  .aggregate_category_representations(reps, weights = weights)
}

S7::method(aggregate, MVBU_CategoryRepresentationTemplate) <- function(x, ..., weights = NULL) {
  dots <- list(...)
  templates <- c(list(x), dots)
  .aggregate_category_representation_templates(templates, weights = weights)
}

S7::method(aggregate, MVBU_CognitiveModel) <- function(x, ..., weights = NULL) {
  dots <- list(...)
  models <- c(list(x), dots)
  .aggregate_cognitive_models(models, weights = weights)
}

S7::method(aggregate, S7::class_list) <- function(x, weights = NULL, ...) {
  .assert_true(length(x) > 0, msg = "Cannot aggregate an empty list.")
  first_obj <- x[[1]]
  if (S7::S7_inherits(first_obj, MVBU_CategoryRepresentation)) {
    .aggregate_category_representations(x, weights = weights)
  } else if (S7::S7_inherits(first_obj, MVBU_CategoryRepresentationTemplate)) {
    .aggregate_category_representation_templates(x, weights = weights)
  } else if (S7::S7_inherits(first_obj, MVBU_CognitiveModel)) {
    .aggregate_cognitive_models(x, weights = weights)
  } else {
    .stop(
      "List elements must be S7 category representations, templates, or cognitive models."
    )
  }
}

# -----------------------------------------------------------------------------
# aggregate_models public function
# -----------------------------------------------------------------------------

#' Aggregate Multiple Category Representations, Templates, or Cognitive Models
#'
#' @description
#' Aggregates multiple S7 category representations, category representation
#' templates, or cognitive models of the same family by computing weighted
#' averages of their category distributions and model parameters.
#'
#' @name aggregate_models
#' @rdname aggregate_models
#' @aliases aggregate aggregate_models
#'
#' @param x,models An S7 category representation, template, cognitive model, or
#'   a non-empty list of such objects belonging to the same S7 class. For
#'   \code{aggregate}, additional objects can also be passed via \code{...}.
#'   For \code{aggregate_models}, a non-empty list of cognitive models.
#' @param weights Optional numeric vector of positive weights with the same
#'   length as the number of objects being aggregated. If \code{NULL} (default),
#'   uniform weights \eqn{1/N} are used. Weights are automatically normalized to
#'   sum to 1.
#' @param ... Additional S7 objects of the same class (for variadic usage of
#'   \code{aggregate}) or arguments passed to methods.
#'
#' @details
#' Aggregation operates hierarchically depending on the class of the input:
#'
#' \subsection{Parameter Aggregation Functions}{
#' \tabular{lll}{
#'   \strong{Object Level} \tab \strong{Parameter / Slot} \tab \strong{Aggregation Function} \cr
#'   \emph{Category Representation (MVG)} \tab \code{mu} (mean vector) \tab Weighted arithmetic mean: \eqn{\sum_i w_i \mu_i} \cr
#'   \emph{Category Representation (MVG)} \tab \code{Sigma} (covariance matrix) \tab Weighted arithmetic mean: \eqn{\sum_i w_i \Sigma_i} \cr
#'   \emph{Category Representation (UVG)} \tab \code{mu} (mean scalar) \tab Weighted arithmetic mean: \eqn{\sum_i w_i \mu_i} \cr
#'   \emph{Category Representation (UVG)} \tab \code{sigma2} (variance scalar) \tab Weighted arithmetic mean: \eqn{\sum_i w_i \sigma^2_i} \cr
#'   \emph{Category Representation (MUVG)} \tab \code{mu} (mean vector) \tab Weighted arithmetic mean: \eqn{\sum_i w_i \mu_i} \cr
#'   \emph{Category Representation (MUVG)} \tab \code{sigma2} (variance vector) \tab Weighted arithmetic mean: \eqn{\sum_i w_i \sigma^2_i} \cr
#'   \emph{Category Representation (NIW)} \tab \code{m} (location vector) \tab Weighted arithmetic mean: \eqn{\sum_i w_i m_i} \cr
#'   \emph{Category Representation (NIW)} \tab \code{S} (scatter matrix) \tab Weighted arithmetic mean: \eqn{\sum_i w_i S_i} \cr
#'   \emph{Category Representation (NIW)} \tab \code{kappa} (degrees of freedom) \tab Weighted arithmetic mean: \eqn{\sum_i w_i \kappa_i} \cr
#'   \emph{Category Representation (NIW)} \tab \code{nu} (degrees of freedom) \tab Weighted arithmetic mean: \eqn{\sum_i w_i \nu_i} \cr
#'   \emph{Category Representation (NIX)} \tab \code{m} (location scalar) \tab Weighted arithmetic mean: \eqn{\sum_i w_i m_i} \cr
#'   \emph{Category Representation (NIX)} \tab \code{S} (scatter scalar) \tab Weighted arithmetic mean: \eqn{\sum_i w_i S_i} \cr
#'   \emph{Category Representation (NIX)} \tab \code{kappa} (degrees of freedom) \tab Weighted arithmetic mean: \eqn{\sum_i w_i \kappa_i} \cr
#'   \emph{Category Representation (NIX)} \tab \code{nu} (degrees of freedom) \tab Weighted arithmetic mean: \eqn{\sum_i w_i \nu_i} \cr
#'   \emph{Category Representation (MNIX)} \tab \code{m}, \code{S}, \code{kappa}, \code{nu} \tab Weighted arithmetic means across cues \cr
#'   \emph{Category Representation (EXEMPLAR)} \tab \code{exemplars} \tab Row-binding (pooling) of all exemplars: \code{rbind} \cr
#'   \emph{Category Representation (EXEMPLAR)} \tab \code{exemplar_weights} \tab Re-weighted by model weights \eqn{w_i W_i} and normalized \cr
#'   \emph{Category Representation (EXEMPLAR)} \tab \code{c} (sensitivity scalar) \tab Weighted arithmetic mean: \eqn{\sum_i w_i c_i} \cr
#'   \emph{Category Template} \tab \code{representations} \tab Aggregates each representation category-by-category \cr
#'   \emph{Cognitive Model} \tab \code{category_template} \tab Template aggregation across models \cr
#'   \emph{Cognitive Model} \tab \code{category_prior} \tab Weighted arithmetic mean: \eqn{\sum_i w_i P(c)_i} \cr
#'   \emph{Cognitive Model} \tab \code{lapse_rate} \tab Weighted arithmetic mean: \eqn{\sum_i w_i \lambda_i} \cr
#'   \emph{Cognitive Model} \tab \code{lapse_bias} \tab Weighted arithmetic mean: \eqn{\sum_i w_i \beta_i} \cr
#'   \emph{Cognitive Model} \tab \code{Sigma_noise} \tab Weighted arithmetic mean: \eqn{\sum_i w_i \Sigma_{\text{noise}, i}} \cr
#'   \emph{Cognitive Model} \tab \code{decision_rule} \tab Inherited from first model \cr
#'   \emph{Cognitive Model} \tab \code{noise_treatment} \tab Inherited from first model \cr
#'   \emph{Cognitive Model} \tab \code{lapse_treatment} \tab Inherited from first model \cr
#' }
#' }
#'
#' @return An S7 object of the same class as the inputs with aggregated parameters.
#' @seealso \code{\link{new_mvg_ideal_observer}},
#'   \code{\link{new_niw_ideal_adaptor}},
#'   \code{\link{new_exemplar_model}}
#' @export
aggregate_models <- function(models, weights = NULL, ...) {
  .aggregate_cognitive_models(models, weights = weights)
}
