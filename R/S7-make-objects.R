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
# Model Aggregation
# -----------------------------------------------------------------------------

#' Aggregate Multiple Cognitive Models
#'
#' Aggregates a list of S7 cognitive models of the same family by computing
#' weighted averages of their category representations and model parameters.
#'
#' @name aggregate_models
#' @rdname aggregate_models
#' @param models A non-empty list of S7 cognitive model objects belonging to the
#'   same class.
#' @param weights Optional numeric vector of positive weights with the same
#'   length as `models`. If `NULL`, uniform weights are used.
#' @param ... Additional arguments (currently unused).
#'
#' @return An S7 cognitive model of the same class as the elements in `models`,
#'   with averaged parameters.
#' @seealso \code{\link{new_mvg_ideal_observer}},
#'   \code{\link{new_niw_ideal_adaptor}}
#' @export
aggregate_models <- function(models, weights = NULL, ...) {
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
  if (is.null(weights)) {
    weights <- rep(1 / n_models, n_models)
  } else {
    .assert_true(
      is.numeric(weights) && length(weights) == n_models &&
        all(weights >= 0) && sum(weights) > 0,
      msg = paste0(
        "weights must be a numeric vector of positive weights ",
        "with the same length as models."
      )
    )
    weights <- weights / sum(weights)
  }

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
  first_reps <- get_category_representations(first_model)
  categories <- names(first_reps)

  agg_reps <- list()
  for (cat in categories) {
    cat_reps <- lapply(models, function(m) {
      get_category_representations(m)[[cat]]
    })
    first_cat_rep <- cat_reps[[1]]
    cat_rep_class <- S7::S7_class(first_cat_rep)

    cat_labels <- get_category_labels(first_cat_rep)
    cue_labels <- get_cue_labels(first_cat_rep)

    if (S7::S7_inherits(first_cat_rep, MVG_CategoryRepresentation)) {
      mu_weighted <- Reduce(
        `+`,
        mapply(function(r, w) r@mu * w, cat_reps, weights, SIMPLIFY = FALSE)
      )
      Sigma_weighted <- Reduce(
        `+`,
        mapply(function(r, w) r@Sigma * w, cat_reps, weights, SIMPLIFY = FALSE)
      )
      agg_reps[[cat]] <- new_mvg_category_representation(
        category_labels = cat_labels,
        cue_labels = cue_labels,
        mu = mu_weighted,
        Sigma = Sigma_weighted
      )
    } else if (S7::S7_inherits(first_cat_rep, NIW_CategoryRepresentation)) {
      m_weighted <- Reduce(
        `+`,
        mapply(function(r, w) r@m * w, cat_reps, weights, SIMPLIFY = FALSE)
      )
      S_weighted <- Reduce(
        `+`,
        mapply(function(r, w) r@S * w, cat_reps, weights, SIMPLIFY = FALSE)
      )
      kappas <- vapply(cat_reps, function(r) r@kappa, numeric(1))
      nus <- vapply(cat_reps, function(r) r@nu, numeric(1))
      agg_reps[[cat]] <- new_niw_category_representation(
        category_labels = cat_labels,
        cue_labels = cue_labels,
        m = m_weighted,
        S = S_weighted,
        kappa = sum(kappas * weights),
        nu = sum(nus * weights)
      )
    } else if (S7::S7_inherits(first_cat_rep, UVG_CategoryRepresentation)) {
      mu_weighted <- sum(
        vapply(cat_reps, function(r) r@mu, numeric(1)) * weights
      )
      sigma2_weighted <- sum(
        vapply(cat_reps, function(r) r@sigma2, numeric(1)) * weights
      )
      agg_reps[[cat]] <- new_uvg_category_representation(
        category_labels = cat_labels,
        cue_labels = cue_labels,
        mu = mu_weighted,
        sigma2 = sigma2_weighted
      )
    } else if (S7::S7_inherits(first_cat_rep, NIX_CategoryRepresentation)) {
      m_weighted <- sum(
        vapply(cat_reps, function(r) r@m, numeric(1)) * weights
      )
      S_weighted <- sum(
        vapply(cat_reps, function(r) r@S, numeric(1)) * weights
      )
      kappas <- vapply(cat_reps, function(r) r@kappa, numeric(1))
      nus <- vapply(cat_reps, function(r) r@nu, numeric(1))
      agg_reps[[cat]] <- new_nix_category_representation(
        category_labels = cat_labels,
        cue_labels = cue_labels,
        m = m_weighted,
        S = S_weighted,
        kappa = sum(kappas * weights),
        nu = sum(nus * weights)
      )
    } else if (S7::S7_inherits(first_cat_rep, MNIX_CategoryRepresentation)) {
      m_weighted <- Reduce(
        `+`,
        mapply(function(r, w) r@m * w, cat_reps, weights, SIMPLIFY = FALSE)
      )
      S_weighted <- Reduce(
        `+`,
        mapply(function(r, w) r@S * w, cat_reps, weights, SIMPLIFY = FALSE)
      )
      kappas <- Reduce(
        `+`,
        mapply(function(r, w) r@kappa * w, cat_reps, weights, SIMPLIFY = FALSE)
      )
      nus <- Reduce(
        `+`,
        mapply(function(r, w) r@nu * w, cat_reps, weights, SIMPLIFY = FALSE)
      )
      agg_reps[[cat]] <- new_mnix_category_representation(
        category_labels = cat_labels,
        cue_labels = cue_labels,
        m = m_weighted,
        S = S_weighted,
        kappa = kappas,
        nu = nus
      )
    } else {
      .stop(
        "Aggregation is not implemented for category representation class: ",
        cat_rep_class@name
      )
    }

  }

  agg_template <- new_category_representation_template(
    representations = agg_reps
  )

  # 3. Construct aggregated model
  args <- list(
    category_template = agg_template,
    category_prior = agg_category_prior,
    lapse_rate = agg_lapse_rate,
    lapse_bias = agg_lapse_bias,
    Sigma_noise = agg_sigma_noise,
    decision_rule = first_model@decision_rule,
    noise_treatment = first_model@noise_behavior$noise_treatment,
    lapse_treatment = first_model@lapse_behavior$lapse_treatment,
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
