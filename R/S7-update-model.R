#' @include S7-methods.R S7-generics.R
NULL

# Order of application matters: all update helper functions assume inputs
# (kappa, nu, m, S) that reflect the state prior to the current update step.

# Conjugate update of prior weight parameter kappa:
# kappa_N = kappa_0 + n_exposure
.update_NIW_category_representation_kappa <- function(kappa_0, n_exposure) {
  kappa_0 + n_exposure
}

# Conjugate update of degrees of freedom nu:
# nu_N = nu_0 + n_exposure
.update_NIW_category_representation_nu <- function(nu_0, n_exposure) {
  nu_0 + n_exposure
}

# Conjugate update of mean vector m:
# m_N = (kappa_0 / kappa_N) * m_0 + (n_exposure / kappa_N) * mean_exposure
.update_NIW_category_representation_m <- function(
  kappa_0, m_0, n_exposure, mean_exposure
) {
  denom <- kappa_0 + n_exposure
  res <- (kappa_0 / denom) * m_0 +
    (n_exposure / denom) * as.numeric(mean_exposure)
  names(res) <- names(m_0)
  res
}

# Conjugate update of scatter matrix S with between-mean scatter contribution:
# S_N = S_0 + css_exposure +
#       ((kappa_0 * n_exposure) / kappa_N) * (mean_exposure - m_0)(mean_exposure - m_0)'
.update_NIW_category_representation_S <- function(
  kappa_0, m_0, S_0, n_exposure, mean_exposure, css_exposure
) {
  diff_m <- matrix(as.numeric(mean_exposure - m_0), ncol = 1L)
  between_scatter <- ((kappa_0 * n_exposure) / (kappa_0 + n_exposure)) *
    tcrossprod(diff_m)
  S_update <- S_0 + as.matrix(css_exposure) + between_scatter
  dimnames(S_update) <- dimnames(S_0)
  S_update
}

# Update an NIW category representation from sufficient statistics
# (sample mean, centered sum of squares, and sample size).
.update_NIW_category_representation_by_sufficient_statistics <- function(
  representation, mean_exposure, css_exposure, n_exposure
) {
  if (n_exposure == 0) {
    return(representation)
  }
  if (anyNA(mean_exposure)) {
    .stop("mean_exposure cannot contain NA values when n_exposure > 0.")
  }
  if (anyNA(css_exposure)) {
    .stop("css_exposure cannot contain NA values when n_exposure > 0.")
  }

  .assert_true(
    is.numeric(mean_exposure), msg = "mean_exposure must be numeric."
  )
  .assert_true(
    is.numeric(css_exposure), msg = "css_exposure must be numeric."
  )

  D <- length(get_cue_labels(representation))
  mean_exposure_vec <- as.numeric(mean_exposure)
  .assert_true(
    length(mean_exposure_vec) == D,
    msg = "mean_exposure must contain one value per cue."
  )
  css_exposure_mat <- as.matrix(css_exposure)
  .assert_true(
    all(dim(css_exposure_mat) == c(D, D)),
    msg = "css_exposure must be a D x D matrix."
  )

  new_niw_category_representation(
    category_labels = get_category_labels(representation)[1L],
    cue_labels = get_cue_labels(representation),
    m = .update_NIW_category_representation_m(
      representation@kappa, representation@m, n_exposure, mean_exposure_vec
    ),
    kappa = .update_NIW_category_representation_kappa(
      representation@kappa, n_exposure
    ),
    nu = .update_NIW_category_representation_nu(
      representation@nu, n_exposure
    ),
    S = .update_NIW_category_representation_S(
      representation@kappa, representation@m, representation@S,
      n_exposure, mean_exposure_vec, css_exposure_mat
    ),
    metadata = representation@metadata
  )
}

# Update a univariate NIX category representation from sufficient statistics
# (sample mean, centered sum of squares, and sample size).
.update_NIX_category_representation_by_sufficient_statistics <- function(
  representation, mean_exposure, css_exposure, n_exposure
) {
  if (n_exposure == 0) {
    return(representation)
  }
  if (anyNA(mean_exposure)) {
    .stop("mean_exposure cannot contain NA values when n_exposure > 0.")
  }
  if (anyNA(css_exposure)) {
    .stop("css_exposure cannot contain NA values when n_exposure > 0.")
  }

  .assert_true(
    is.numeric(mean_exposure), msg = "mean_exposure must be numeric."
  )
  .assert_true(
    is.numeric(css_exposure), msg = "css_exposure must be numeric."
  )

  mean_exposure_scal <- as.numeric(mean_exposure)[1L]
  css_exposure_scal <- if (is.matrix(css_exposure)) {
    as.numeric(css_exposure)[1L]
  } else {
    as.numeric(css_exposure)[1L]
  }

  k0 <- representation@kappa
  nu0 <- representation@nu
  m0 <- representation@m
  sig0 <- representation@sigma2

  k_n <- k0 + n_exposure
  nu_n <- nu0 + n_exposure
  m_n <- (k0 / k_n) * m0 + (n_exposure / k_n) * mean_exposure_scal
  sig_n <- (nu0 * sig0 + css_exposure_scal +
    ((k0 * n_exposure) / k_n) * (mean_exposure_scal - m0)^2) / nu_n
  names(m_n) <- names(m0)
  names(sig_n) <- names(sig0)

  new_nix_category_representation(
    category_labels = get_category_labels(representation)[1L],
    cue_labels = get_cue_labels(representation),
    m = m_n,
    sigma2 = sig_n,
    kappa = k_n,
    nu = nu_n,
    metadata = representation@metadata
  )
}

# Update a multivariate independent MNIX category representation from
# sufficient statistics (sample mean, centered sum of squares, sample size).
.update_MNIX_category_representation_by_sufficient_statistics <- function(
  representation, mean_exposure, css_exposure, n_exposure
) {
  if (n_exposure == 0) {
    return(representation)
  }
  if (anyNA(mean_exposure)) {
    .stop("mean_exposure cannot contain NA values when n_exposure > 0.")
  }
  if (anyNA(css_exposure)) {
    .stop("css_exposure cannot contain NA values when n_exposure > 0.")
  }

  .assert_true(
    is.numeric(mean_exposure), msg = "mean_exposure must be numeric."
  )
  .assert_true(
    is.numeric(css_exposure), msg = "css_exposure must be numeric."
  )

  D <- length(get_cue_labels(representation))
  mean_exposure_vec <- as.numeric(mean_exposure)
  .assert_true(
    length(mean_exposure_vec) == D,
    msg = "mean_exposure must contain one value per cue."
  )

  css_exposure_vec <- if (is.matrix(css_exposure) &&
    nrow(css_exposure) == ncol(css_exposure) && nrow(css_exposure) == D) {
    diag(css_exposure)
  } else {
    as.numeric(css_exposure)
  }
  .assert_true(
    length(css_exposure_vec) == D,
    msg = "css_exposure must contain one value per cue or be a D x D matrix."
  )

  k0 <- representation@kappa
  nu0 <- representation@nu
  m0 <- representation@m
  sig0 <- representation@sigma2

  k_n <- k0 + n_exposure
  nu_n <- nu0 + n_exposure
  m_n <- (k0 / k_n) * m0 + (n_exposure / k_n) * mean_exposure_vec
  sig_n <- (nu0 * sig0 + css_exposure_vec +
    ((k0 * n_exposure) / k_n) * (mean_exposure_vec - m0)^2) / nu_n
  names(m_n) <- names(m0)
  names(sig_n) <- names(sig0)

  new_mnix_category_representation(
    category_labels = get_category_labels(representation)[1L],
    cue_labels = get_cue_labels(representation),
    m = m_n,
    sigma2 = sig_n,
    kappa = k_n,
    nu = nu_n,
    weights = representation@weights,
    metadata = representation@metadata
  )
}

# Distribute evidence across categories for unlabeled or semi-supervised
# observation updating according to the chosen update rule.
.distribute_evidence_across_categories <- function(
  model, observation, update_method
) {
  categories <- get_category_labels(model)
  if (update_method == "nolabel-uniform") {
    return(rep(1 / length(categories), length(categories)))
  }
  obs_mat <- matrix(as.numeric(observation), nrow = 1L)
  probabilities <- posterior(model, obs_mat, categories = categories)[1L, ]
  if (update_method == "nolabel-criterion") {
    weights <- numeric(length(categories))
    weights[which.max(probabilities)] <- 1
    return(weights)
  }
  if (update_method == "nolabel-sampling") {
    weights <- numeric(length(categories))
    weights[sample.int(length(categories), 1L, prob = probabilities)] <- 1
    return(weights)
  }
  probabilities
}

#' @rdname update_category_representation
#' @export
S7::method(
  update_category_representation,
  list(
    NIW_CategoryRepresentation,
    S7::class_any,
    S7::class_any,
    S7::class_any
  )
) <- function(
  x, n_exposure, mean_exposure, css_exposure, ...
) {
  .assert_true(
    is.numeric(n_exposure) && length(n_exposure) == 1L &&
      !is.na(n_exposure) && n_exposure >= 0,
    msg = "n_exposure must be a non-negative numeric scalar."
  )
  .update_NIW_category_representation_by_sufficient_statistics(
    representation = x, mean_exposure = mean_exposure,
    css_exposure = css_exposure, n_exposure = n_exposure
  )
}

#' @rdname update_category_representation
#' @export
S7::method(
  update_category_representation,
  list(
    NIX_CategoryRepresentation,
    S7::class_any,
    S7::class_any,
    S7::class_any
  )
) <- function(
  x, n_exposure, mean_exposure, css_exposure, ...
) {
  .assert_true(
    is.numeric(n_exposure) && length(n_exposure) == 1L &&
      !is.na(n_exposure) && n_exposure >= 0,
    msg = "n_exposure must be a non-negative numeric scalar."
  )
  .update_NIX_category_representation_by_sufficient_statistics(
    representation = x, mean_exposure = mean_exposure,
    css_exposure = css_exposure, n_exposure = n_exposure
  )
}

#' @rdname update_category_representation
#' @export
S7::method(
  update_category_representation,
  list(
    MNIX_CategoryRepresentation,
    S7::class_any,
    S7::class_any,
    S7::class_any
  )
) <- function(
  x, n_exposure, mean_exposure, css_exposure, ...
) {
  .assert_true(
    is.numeric(n_exposure) && length(n_exposure) == 1L &&
      !is.na(n_exposure) && n_exposure >= 0,
    msg = "n_exposure must be a non-negative numeric scalar."
  )
  .update_MNIX_category_representation_by_sufficient_statistics(
    representation = x, mean_exposure = mean_exposure,
    css_exposure = css_exposure, n_exposure = n_exposure
  )
}

# Reconstruct an ideal adaptor model with updated category representations.
#
# Arguments:
#   x: An ideal adaptor model object to use as the template.
#   updated_representations: List of updated category representation objects.
#
# Returns:
#   A new ideal adaptor model of the same class as `x`, retaining its settings
#   (cues, weights, decision rule, priors, lapse, noise, metadata) with updated
#   representations in its category template.
.reconstruct_ideal_adaptor <- function(x, updated_representations) {
  new_template <- new_category_representation_template(
    updated_representations,
    metadata = x@category_template@metadata
  )
  # Early Integration Models ({Integrator}_MNIX_IdealAdaptor)
  if (S7::S7_inherits(x, IdealIntegrator_MNIX_IdealAdaptor)) {
    new_ideal_integrator_mnix_ideal_adaptor(
      category_template = new_template,
      cues = x@cues,
      decision_rule = x@decision_rule,
      category_prior = x@category_prior,
      lapse_rate = get_lapse_rate(x),
      lapse_bias = get_lapse_bias(x),
      Sigma_noise = get_noise(x),
      noise_treatment = get_noise_treatment(x),
      lapse_treatment = get_lapse_treatment(x),
      metadata = x@metadata
    )
  } else if (S7::S7_inherits(x, WeightedIntegrator_MNIX_IdealAdaptor)) {
    new_weighted_integrator_mnix_ideal_adaptor(
      category_template = new_template,
      weights = x@weights,
      cues = x@cues,
      decision_rule = x@decision_rule,
      category_prior = x@category_prior,
      lapse_rate = get_lapse_rate(x),
      lapse_bias = get_lapse_bias(x),
      Sigma_noise = get_noise(x),
      noise_treatment = get_noise_treatment(x),
      lapse_treatment = get_lapse_treatment(x),
      metadata = x@metadata
    )
  } else if (S7::S7_inherits(x, CategoryWeightedIntegrator_MNIX_IdealAdaptor)) {
    new_category_weighted_integrator_mnix_ideal_adaptor(
      category_template = new_template,
      weights = x@weights,
      cues = x@cues,
      decision_rule = x@decision_rule,
      category_prior = x@category_prior,
      lapse_rate = get_lapse_rate(x),
      lapse_bias = get_lapse_bias(x),
      Sigma_noise = get_noise(x),
      noise_treatment = get_noise_treatment(x),
      lapse_treatment = get_lapse_treatment(x),
      metadata = x@metadata
    )
  } else if (
    S7::S7_inherits(x, InformativityWeightedIntegrator_MNIX_IdealAdaptor)
  ) {
    new_informativity_weighted_integrator_mnix_ideal_adaptor(
      category_template = new_template,
      cues = x@cues,
      decision_rule = x@decision_rule,
      category_prior = x@category_prior,
      lapse_rate = get_lapse_rate(x),
      lapse_bias = get_lapse_bias(x),
      Sigma_noise = get_noise(x),
      noise_treatment = get_noise_treatment(x),
      lapse_treatment = get_lapse_treatment(x),
      metadata = x@metadata
    )
  # Late Integration Models (MNIX_IdealAdaptor_{Integrator})
  } else if (S7::S7_inherits(x, MNIX_IdealAdaptor_IdealIntegrator)) {
    new_mnix_ideal_adaptor_ideal_integrator(
      category_template = new_template,
      decision_rule = x@decision_rule,
      category_prior = x@category_prior,
      lapse_rate = get_lapse_rate(x),
      lapse_bias = get_lapse_bias(x),
      Sigma_noise = get_noise(x),
      noise_treatment = get_noise_treatment(x),
      lapse_treatment = get_lapse_treatment(x),
      metadata = x@metadata
    )
  } else if (S7::S7_inherits(x, MNIX_IdealAdaptor_WeightedIntegrator) ||
             S7::S7_inherits(x, MNIX_IdealAdaptor_AdaptiveWeightedIntegrator)) {
    new_mnix_ideal_adaptor_weighted_integrator(
      category_template = new_template,
      weights = x@weights,
      decision_rule = x@decision_rule,
      category_prior = x@category_prior,
      lapse_rate = get_lapse_rate(x),
      lapse_bias = get_lapse_bias(x),
      Sigma_noise = get_noise(x),
      noise_treatment = get_noise_treatment(x),
      lapse_treatment = get_lapse_treatment(x),
      metadata = x@metadata
    )
  } else if (
    S7::S7_inherits(x, MNIX_IdealAdaptor_CategoryWeightedIntegrator) ||
    S7::S7_inherits(x, MNIX_IdealAdaptor_CategoryAdaptiveWeightedIntegrator)
  ) {
    new_mnix_ideal_adaptor_category_weighted_integrator(
      category_template = new_template,
      weights = x@weights,
      decision_rule = x@decision_rule,
      category_prior = x@category_prior,
      lapse_rate = get_lapse_rate(x),
      lapse_bias = get_lapse_bias(x),
      Sigma_noise = get_noise(x),
      noise_treatment = get_noise_treatment(x),
      lapse_treatment = get_lapse_treatment(x),
      metadata = x@metadata
    )
  } else if (
    S7::S7_inherits(x, MNIX_IdealAdaptor_InformativityWeightedIntegrator)
  ) {
    new_mnix_ideal_adaptor_informativity_weighted_integrator(
      category_template = new_template,
      decision_rule = x@decision_rule,
      category_prior = x@category_prior,
      lapse_rate = get_lapse_rate(x),
      lapse_bias = get_lapse_bias(x),
      Sigma_noise = get_noise(x),
      noise_treatment = get_noise_treatment(x),
      lapse_treatment = get_lapse_treatment(x),
      metadata = x@metadata
    )
  # Base Models
  } else if (S7::S7_inherits(x, MNIX_IdealAdaptor)) {
    new_mnix_ideal_adaptor(
      category_template = new_template,
      decision_rule = x@decision_rule,
      category_prior = x@category_prior,
      lapse_rate = get_lapse_rate(x),
      lapse_bias = get_lapse_bias(x),
      Sigma_noise = get_noise(x),
      noise_treatment = get_noise_treatment(x),
      lapse_treatment = get_lapse_treatment(x),
      metadata = x@metadata
    )
  } else if (S7::S7_inherits(x, NIX_IdealAdaptor)) {
    new_nix_ideal_adaptor(
      category_template = new_template,
      decision_rule = x@decision_rule,
      category_prior = x@category_prior,
      lapse_rate = get_lapse_rate(x),
      lapse_bias = get_lapse_bias(x),
      Sigma_noise = get_noise(x),
      noise_treatment = get_noise_treatment(x),
      lapse_treatment = get_lapse_treatment(x),
      metadata = x@metadata
    )
  } else if (S7::S7_inherits(x, NIW_IdealAdaptor)) {
    new_niw_ideal_adaptor(
      category_template = new_template,
      decision_rule = x@decision_rule,
      category_prior = x@category_prior,
      lapse_rate = get_lapse_rate(x),
      lapse_bias = get_lapse_bias(x),
      Sigma_noise = get_noise(x),
      noise_treatment = get_noise_treatment(x),
      lapse_treatment = get_lapse_treatment(x),
      metadata = x@metadata
    )
  } else {
    .stop("Unsupported ideal adaptor class: ", class(x)[1L])
  }
}

# Internal workhorse to update category representations of an ideal adaptor.
#
# Arguments:
#   x: An ideal adaptor model object.
#   observations: Observation data (data frame or matrix).
#   updating: Updating mode ("batch" or "incremental").
#   keep_history: Logical; whether to record intermediate update history.
#   lapse_treatment: Perceptual lapse treatment mode.
#   noise_treatment: Perceptual noise treatment mode.
#   update_method: Belief updating method / strategy.
#   ...: Additional arguments passed to update methods.
#
# Returns:
#   Updated model (if keep_history = FALSE) or a list of models across updates.
.update_ideal_adaptor_template <- function(
  x, observations, updating = c("batch", "incremental"), keep_history = FALSE,
  lapse_treatment = "no_lapses", noise_treatment = "no_noise",
  update_method = "label-certain", ...
) {
  updating <- match.arg(updating)
  update_methods <- c(
    "no-updating", "label-certain", "nolabel-criterion",
    "nolabel-sampling", "nolabel-proportional", "nolabel-uniform"
  )
  .assert_one_of(
    update_method, update_methods,
    msg = "update_method is not an acceptable updating method."
  )
  .assert_one_of(
    lapse_treatment, c("no_lapses", "sample", "marginalize"),
    msg = "lapse_treatment must be one of no_lapses, sample, or marginalize."
  )
  .assert_one_of(
    noise_treatment, c("no_noise", "sample", "marginalize"),
    msg = "noise_treatment must be one of no_noise, sample, or marginalize."
  )
  .assert_data_frame_like(observations)

  cue_labels <- get_cue_labels(x)
  .assert_data_contains_cols(observations, cue_labels)
  if (update_method == "label-certain") {
    .assert_data_contains_cols(observations, "category")
  }
  .assert_true(
    updating == "incremental" || update_method == "label-certain",
    msg = paste(
      "Batch updating currently supports only",
      "update_method = 'label-certain'."
    )
  )

  if (update_method == "no-updating") {
    return(if (keep_history && updating == "incremental") list(x) else x)
  }


  family_str <- if (S7::S7_inherits(x, NIW_IdealAdaptor)) {
    "NIW"
  } else if (S7::S7_inherits(x, NIX_IdealAdaptor)) {
    "NIX"
  } else if (S7::S7_inherits(x, MNIX_IdealAdaptor)) {
    "MNIX"
  } else {
    .stop("Unsupported ideal adaptor class: ", class(x)[1L])
  }

  categories <- get_category_labels(x)
  representations <- x@category_template@representations

  if (updating == "batch") {
    # Prepare exposure data
    values <- as.matrix(observations[, cue_labels, drop = FALSE])
    sig_noise <- get_noise(x)
    if (noise_treatment == "sample" && !is.null(sig_noise)) {
      values <- values + mvtnorm::rmvnorm(nrow(values), sigma = sig_noise)
    }
    l_rate <- get_lapse_rate(x)
    if (lapse_treatment == "sample" && l_rate > 0) {
      keep <- stats::runif(nrow(values)) >= l_rate
      values <- values[keep, , drop = FALSE]
      observations <- observations[keep, , drop = FALSE]
    }
    observations[, cue_labels] <- values
    lapse_weight <- if (lapse_treatment == "marginalize") 1 - l_rate else 1

    suff_stats <- get_sufficient_category_statistics(
      observations,
      cues = cue_labels,
      category = "category",
      model_family = family_str,
      simplify = FALSE
    )
    suff_by_cat <- split(suff_stats, suff_stats$category)

    updated <- lapply(seq_along(representations), function(i) {
      cat_name <- as.character(categories[i])
      st <- suff_by_cat[[cat_name]]
      if (is.null(st) || nrow(st) == 0L || st$n[[1L]] == 0L) {
        return(representations[[i]])
      }
      # Rationale for linear lapse weighting of css_exposure:
      # Although centered sum of squares (CSS) is quadratic in cue coordinates
      # (y - mean), it is an extensive sum of N observations:
      #   CSS = sum_{i=1}^N (y_i - mean)(y_i - mean)' = (N - 1) * S^2
      # Under lapse rate lambda, the effective number of category tokens is
      # N_eff = w * N, where w = 1 - lambda (lapse_weight).
      # The expected sum of squares across those N_eff genuine tokens scales
      # linearly with N_eff: E[CSS_eff] = N_eff * Sigma = w * CSS.
      # Scaling both n_exposure and css_exposure linearly by w preserves an
      # unbiased sample covariance:
      #   Sigma_eff = CSS_eff / N_eff = (w * CSS) / (w * N) = CSS / N
      # Quadratic weighting (w^2 * CSS) would incorrectly scale covariance by
      # w = 1 - lambda, artificially shrinking category variance toward 0 as
      # lapse rate increases. Lapses reduce degrees of freedom (confidence),
      # not the underlying category cue spread.
      update_category_representation(
        representations[[i]],
        n_exposure = lapse_weight * st$n[[1L]],
        mean_exposure = st$mean[[1L]],
        css_exposure = lapse_weight * st$css[[1L]]
      )
    })
    names(updated) <- names(representations)
    return(.reconstruct_ideal_adaptor(x, updated))
  }

  # Incremental updating
  D <- length(cue_labels)
  zero_css <- if (family_str == "NIW") {
    matrix(0, nrow = D, ncol = D)
  } else if (family_str == "NIX") {
    0
  } else {
    rep(0, D)
  }

  history <- list(x)
  current <- x
  for (i in seq_len(nrow(observations))) {
    category <- if ("category" %in% names(observations)) {
      observations$category[[i]]
    } else {
      NULL
    }
    observation <- as.numeric(observations[i, cue_labels])
    sig_noise <- get_noise(current)
    if (noise_treatment == "sample" && !is.null(sig_noise)) {
      observation <- observation + as.numeric(
        mvtnorm::rmvnorm(1L, sigma = sig_noise)
      )
    }
    update_weight <- 1
    l_rate <- get_lapse_rate(current)
    if (lapse_treatment == "sample" && stats::runif(1L) < l_rate) {
      update_weight <- 0
    } else if (lapse_treatment == "marginalize") {
      update_weight <- 1 - l_rate
    }

    current_cats <- get_category_labels(current)
    weights <- if (update_method == "label-certain") {
      .assert_true(
        !is.null(category),
        msg = "category is required for label-certain updating."
      )
      result <- numeric(length(current_cats))
      result[match(as.character(category), current_cats)] <- 1
      result
    } else {
      .distribute_evidence_across_categories(
        current, observation, update_method
      )
    }

    current_reps <- current@category_template@representations
    updated <- lapply(seq_along(current_reps), function(j) {
      update_category_representation(
        current_reps[[j]],
        n_exposure = update_weight * weights[j],
        mean_exposure = observation,
        css_exposure = zero_css
      )
    })
    names(updated) <- names(current_reps)
    current <- .reconstruct_ideal_adaptor(current, updated)
    if (keep_history) history[[length(history) + 1L]] <- current
  }
  if (keep_history) history else current
}

#' Update an S7 ideal adaptor's category template from observations
#'
#' `update_template()` updates the category template representations of an
#' ideal adaptor model (`NIW_IdealAdaptor`, `NIX_IdealAdaptor`, or
#' `MNIX_IdealAdaptor`) given new observations.
#'
#' @rdname update_template
#' @param x An `NIW_IdealAdaptor`, `NIX_IdealAdaptor`, or `MNIX_IdealAdaptor`
#'   object.
#' @param observations A data.frame or tibble containing cue columns and, for
#'   `update_method = "label-certain"`, a `category` column.
#' @param updating Either `"batch"` or `"incremental"`. Batch updates are
#'   currently supported only for `update_method = "label-certain"`. Default is
#'   `"batch"`.
#' @param keep_history If `TRUE` with incremental updating, return a list of
#'   all intermediate models followed by the final model. Default is `FALSE`.
#' @param lapse_treatment Treatment of attentional lapses during updating for
#'   label-uncertain models: `"no_lapses"`, `"sample"`, or `"marginalize"`.
#'   \itemize{
#'     \item \code{"no_lapses"}: Lapses are ignored; updating reflects pure
#'       normative evidence evaluation (even if \code{lapse_rate > 0}).
#'     \item \code{"sample"}: Simulates stochastic attentional lapses on
#'       individual observations. With probability \code{lapse_rate}, a lapse
#'       occurs and the observation is treated as uninformative during updating.
#'     \item \code{"marginalize"}: Evaluates updating in expectation across
#'       infinitely many encounters; evidence is downweighted by
#'       \eqn{1 - \lambda}.
#'   }
#'   Default is `"no_lapses"`.
#' @param noise_treatment Treatment of perceptual noise during categorization
#'   for label-uncertain updating: `"no_noise"`, `"sample"`, or
#'   `"marginalize"`.
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
#'   Default is `"no_noise"`.
#' @param update_method One of `"no-updating"`, `"label-certain"`,
#'   `"nolabel-criterion"`, `"nolabel-sampling"`,
#'   `"nolabel-proportional"`, or `"nolabel-uniform"`. Default is
#'   `"label-certain"`.
#' @param ... Additional arguments passed to methods.
#' @return The updated model, or a list of intermediate models when
#'   `updating = "incremental"` and `keep_history = TRUE`.
#' @seealso \code{\link{update_category_representation}}
#' @export
S7::method(update_template, list(NIW_IdealAdaptor, S7::class_any)) <- function(
  x, observations, updating = c("batch", "incremental"), keep_history = FALSE,
  lapse_treatment = "no_lapses", noise_treatment = "no_noise",
  update_method = "label-certain", ...
) {
  .update_ideal_adaptor_template(
    x = x, observations = observations, updating = updating,
    keep_history = keep_history, lapse_treatment = lapse_treatment,
    noise_treatment = noise_treatment, update_method = update_method, ...
  )
}

#' @rdname update_template
#' @export
S7::method(update_template, list(NIX_IdealAdaptor, S7::class_any)) <- function(
  x, observations, updating = c("batch", "incremental"), keep_history = FALSE,
  lapse_treatment = "no_lapses", noise_treatment = "no_noise",
  update_method = "label-certain", ...
) {
  .update_ideal_adaptor_template(
    x = x, observations = observations, updating = updating,
    keep_history = keep_history, lapse_treatment = lapse_treatment,
    noise_treatment = noise_treatment, update_method = update_method, ...
  )
}

#' @rdname update_template
#' @export
S7::method(update_template, list(MNIX_IdealAdaptor, S7::class_any)) <- function(
  x, observations, updating = c("batch", "incremental"), keep_history = FALSE,
  lapse_treatment = "no_lapses", noise_treatment = "no_noise",
  update_method = "label-certain", ...
) {
  .update_ideal_adaptor_template(
    x = x, observations = observations, updating = updating,
    keep_history = keep_history, lapse_treatment = lapse_treatment,
    noise_treatment = noise_treatment, update_method = update_method, ...
  )
}

