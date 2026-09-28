#' @include S7-methods.R S7-generics.R
NULL

# Order of application matters: all update helper functions assume inputs (kappa, nu, m, S)
# that reflect the state prior to the current update step.

.update_NIW_category_representation_kappa <- function(kappa_0, x_N) {
  kappa_0 + x_N
}

.update_NIW_category_representation_nu <- function(nu_0, x_N) {
  nu_0 + x_N
}

.update_NIW_category_representation_m <- function(kappa_0, m_0, x_N, x_mean) {
  res <- (kappa_0 / (kappa_0 + x_N)) * m_0 + (x_N / (kappa_0 + x_N)) * as.numeric(x_mean)
  names(res) <- names(m_0)
  res
}

.update_NIW_category_representation_S <- function(kappa_0, m_0, S_0, x_N, x_mean, x_SS) {
  # Conjugate updating of scatter matrix S with between-mean scatter contribution
  diff_m <- matrix(as.numeric(x_mean - m_0), ncol = 1)
  S_update <- S_0 + as.matrix(x_SS) + ((kappa_0 * x_N) / (kappa_0 + x_N)) * tcrossprod(diff_m)
  dimnames(S_update) <- dimnames(S_0)
  S_update
}

.update_NIW_category_representation_by_sufficient_statistics <- function(
  representation, x_mean, x_SS, x_N
) {
  .assert_true(S7::S7_inherits(representation, NIW_CategoryRepresentation),
    msg = "representation must be an NIW_CategoryRepresentation."
  )
  if (x_N == 0 || anyNA(x_mean)) {
    return(representation)
  }

  new_niw_category_representation(
    category_labels = get_category_labels(representation)[1],
    cue_labels = get_cue_labels(representation),
    m = .update_NIW_category_representation_m(representation@kappa, representation@m, x_N, x_mean),
    kappa = .update_NIW_category_representation_kappa(representation@kappa, x_N),
    nu = .update_NIW_category_representation_nu(representation@nu, x_N),
    S = .update_NIW_category_representation_S(representation@kappa, representation@m, representation@S, x_N, x_mean, x_SS),
    metadata = representation@metadata
  )
}

.distribute_evidence_across_categories <- function(model, observation, update_method) {
  categories <- get_category_labels(model)
  if (update_method == "nolabel-uniform") {
    return(rep(1 / length(categories), length(categories)))
  }
  probabilities <- posterior(model, matrix(as.numeric(observation), nrow = 1), categories = categories)[1, ]
  if (update_method == "nolabel-criterion") {
    weights <- numeric(length(categories))
    weights[which.max(probabilities)] <- 1
    return(weights)
  }
  if (update_method == "nolabel-sampling") {
    weights <- numeric(length(categories))
    weights[sample.int(length(categories), 1, prob = probabilities)] <- 1
    return(weights)
  }
  probabilities
}

#' @rdname update_category_representation
#' @export
S7::method(update_category_representation, list(NIW_CategoryRepresentation, S7::class_any, S7::class_any, S7::class_any)) <- function(
  x, x_N, x_mean, x_SS, ...
) {
  .assert_true(is.numeric(x_N) && length(x_N) == 1L && !is.na(x_N) && x_N >= 0, msg = "x_N must be a non-negative numeric scalar.")
  .assert_true(is.numeric(x_mean) && length(x_mean) == length(get_cue_labels(x)), msg = "x_mean must contain one value per cue.")
  x_SS <- as.matrix(x_SS)
  .assert_true(is.numeric(x_SS) && all(dim(x_SS) == length(get_cue_labels(x))), msg = "x_SS must be a numeric square matrix with one row and column per cue.")
  .update_NIW_category_representation_by_sufficient_statistics(
    representation = x, x_mean = x_mean, x_SS = x_SS, x_N = x_N
  )
}

#' Update an S7 NIW ideal adaptor's category template from observations
#'
#' `update_template()` is currently implemented for `NIW_IdealAdaptor` models.
#' Future methods will support additional model types.
#'
#' @rdname update_template
#' @param x An `NIW_IdealAdaptor` object.
#' @param observations A data.frame or tibble containing cue columns and, for
#'   `update_method = "label-certain"`, a `category` column.
#' @param updating Either `"batch"` or `"incremental"`. Batch updates are
#'   currently supported only for `update_method = "label-certain"`.
#' @param keep_history If `TRUE` with incremental updating, return a list of
#'   all intermediate models followed by the final model.
#' @param lapse_treatment Either `"no_lapses"`, `"sample"`, or `"marginalize"`.
#' @param noise_treatment Either `"no_noise"`, `"sample"`, or `"marginalize"`.
#' @param update_method One of `"no-updating"`, `"label-certain"`,
#'   `"nolabel-criterion"`, `"nolabel-sampling"`,
#'   `"nolabel-proportional"`, or `"nolabel-uniform"`.
#' @param ... Additional arguments passed to methods.
#' @return The updated model, or a list of intermediate models when
#'   `updating = "incremental"` and `keep_history = TRUE`.
#' @export
S7::method(update_template, list(NIW_IdealAdaptor, S7::class_any)) <- function(
  x, observations, updating = c("batch", "incremental"), keep_history = FALSE,
  lapse_treatment = "no_lapses", noise_treatment = "no_noise",
  update_method = "label-certain", ...
) {
  updating <- match.arg(updating)
  update_methods <- c("no-updating", "label-certain", "nolabel-criterion", "nolabel-sampling", "nolabel-proportional", "nolabel-uniform")
  .assert_one_of(update_method, update_methods, msg = "update_method is not an acceptable updating method.")
  .assert_one_of(lapse_treatment, c("no_lapses", "sample", "marginalize"), msg = "lapse_treatment must be one of no_lapses, sample, or marginalize.")
  .assert_one_of(noise_treatment, c("no_noise", "sample", "marginalize"), msg = "noise_treatment must be one of no_noise, sample, or marginalize.")
  .assert_data_frame_like(observations)

  cue_labels <- get_cue_labels(x)
  .assert_data_contains_cols(observations, cue_labels)
  if (update_method == "label-certain") .assert_data_contains_cols(observations, "category")
  .assert_true(updating == "incremental" || update_method == "label-certain", msg = "Batch updating currently supports only update_method = 'label-certain'.")

  if (update_method == "no-updating") {
    return(if (keep_history && updating == "incremental") list(x) else x)
  }

  if (updating == "batch") {
    # Prepare exposure data
    values <- as.matrix(observations[, cue_labels, drop = FALSE])
    if (noise_treatment == "sample" && !is.null(x@noise_behavior$Sigma_noise)) {
      values <- values + mvtnorm::rmvnorm(nrow(values), sigma = x@noise_behavior$Sigma_noise)
    }
    if (lapse_treatment == "sample" && x@lapse_behavior$lapse_rate > 0) {
      keep <- stats::runif(nrow(values)) >= x@lapse_behavior$lapse_rate
      values <- values[keep, , drop = FALSE]
      observations <- observations[keep, , drop = FALSE]
    }
    observations[, cue_labels] <- values
    lapse_weight <- if (lapse_treatment == "marginalize") 1 - x@lapse_behavior$lapse_rate else 1
    categories <- get_category_labels(x)
    representations <- x@category_template@representations

    suff_stats <- get_sufficient_category_statistics(
      observations,
      cues = cue_labels,
      category = "category",
      categories = categories,
      model_family = "NIW"
    )
    suff_by_cat <- split(suff_stats, suff_stats$category)

    updated <- lapply(seq_along(representations), function(i) {
      cat_name <- as.character(categories[i])
      st <- suff_by_cat[[cat_name]]
      if (is.null(st) || nrow(st) == 0L || st$x_N[[1L]] == 0L) {
        return(representations[[i]])
      }
      update_category_representation(
        representations[[i]],
        x_N = lapse_weight * st$x_N[[1L]],
        x_mean = st$x_mean[[1L]],
        x_SS = lapse_weight * st$x_css[[1L]]
      )
    })
    names(updated) <- names(representations)
    return(new_niw_ideal_adaptor(
      category_template = new_category_representation_template(updated, metadata = x@category_template@metadata),
      decision_rule = x@decision_rule, category_prior = x@category_prior,
      lapse_rate = get_lapse_rate(x), lapse_bias = get_lapse_bias(x),
      Sigma_noise = get_noise(x), noise_treatment = get_noise_treatment(x),
      lapse_treatment = get_lapse_treatment(x), metadata = x@metadata
    ))
  }

  history <- list(x)
  current <- x
  for (i in seq_len(nrow(observations))) {
    category <- if ("category" %in% names(observations)) observations$category[[i]] else NULL
    observation <- as.numeric(observations[i, cue_labels, drop = FALSE])
    if (noise_treatment == "sample" && !is.null(current@noise_behavior$Sigma_noise)) {
      observation <- observation + as.numeric(mvtnorm::rmvnorm(1, sigma = current@noise_behavior$Sigma_noise))
    }
    update_weight <- 1
    if (lapse_treatment == "sample" && stats::runif(1) < current@lapse_behavior$lapse_rate) {
      update_weight <- 0
    } else if (lapse_treatment == "marginalize") {
      update_weight <- 1 - current@lapse_behavior$lapse_rate
    }

    categories <- get_category_labels(current)
    weights <- if (update_method == "label-certain") {
      .assert_true(!is.null(category), msg = "category is required for label-certain updating.")
      result <- numeric(length(categories))
      result[match(as.character(category), categories)] <- 1
      result
    } else {
      .distribute_evidence_across_categories(current, observation, update_method)
    }
    representations <- current@category_template@representations
    updated <- lapply(seq_along(representations), function(j) {
      update_category_representation(
        representations[[j]],
        x_N = update_weight * weights[j],
        x_mean = observation,
        x_SS = matrix(0, nrow = length(observation), ncol = length(observation))
      )
    })
    names(updated) <- names(representations)
    current <- new_niw_ideal_adaptor(
      category_template = new_category_representation_template(updated, metadata = current@category_template@metadata),
      decision_rule = current@decision_rule,
      category_prior = current@category_prior,
      lapse_rate = get_lapse_rate(current),
      lapse_bias = get_lapse_bias(current),
      Sigma_noise = get_noise(current),
      noise_treatment = get_noise_treatment(current),
      lapse_treatment = get_lapse_treatment(current),
      metadata = current@metadata
    )
    if (keep_history) history[[length(history) + 1L]] <- current
  }
  if (keep_history) history else current
}

.update_NIX_category_representation_by_sufficient_statistics <- function(
  representation, x_mean, x_SS, x_N
) {
  .assert_true(S7::S7_inherits(representation, NIX_CategoryRepresentation),
    msg = "representation must be an NIX_CategoryRepresentation."
  )
  if (x_N == 0 || anyNA(x_mean)) {
    return(representation)
  }
  k0 <- representation@kappa
  nu0 <- representation@nu
  m0 <- representation@m
  sig0 <- representation@sigma2

  k_n <- k0 + x_N
  nu_n <- nu0 + x_N
  m_n <- (k0 / k_n) * m0 + (x_N / k_n) * x_mean[1]
  sig_n <- (nu0 * sig0 + x_SS[1, 1] + ((k0 * x_N) / k_n) * (x_mean[1] - m0)^2) / nu_n

  new_nix_category_representation(
    category_labels = get_category_labels(representation)[1],
    cue_labels = get_cue_labels(representation),
    m = m_n,
    sigma2 = sig_n,
    kappa = k_n,
    nu = nu_n,
    metadata = representation@metadata
  )
}

.update_MNIX_category_representation_by_sufficient_statistics <- function(
  representation, x_mean, x_SS, x_N
) {
  .assert_true(S7::S7_inherits(representation, MNIX_CategoryRepresentation),
    msg = "representation must be an MNIX_CategoryRepresentation."
  )
  if (x_N == 0 || anyNA(x_mean)) {
    return(representation)
  }
  k0 <- representation@kappa
  nu0 <- representation@nu
  m0 <- representation@m
  sig0 <- representation@sigma2

  k_n <- k0 + x_N
  nu_n <- nu0 + x_N
  m_n <- (k0 / k_n) * m0 + (x_N / k_n) * x_mean
  sig_n <- (nu0 * sig0 + diag(x_SS) + ((k0 * x_N) / k_n) * (x_mean - m0)^2) / nu_n

  new_mnix_category_representation(
    category_labels = get_category_labels(representation)[1],
    cue_labels = get_cue_labels(representation),
    m = m_n,
    sigma2 = sig_n,
    kappa = k_n,
    nu = nu_n,
    weights = representation@weights,
    metadata = representation@metadata
  )
}

#' @rdname update_category_representation
#' @export
S7::method(update_category_representation, list(NIX_CategoryRepresentation, S7::class_any, S7::class_any, S7::class_any)) <- function(
  x, x_N, x_mean, x_SS, ...
) {
  .update_NIX_category_representation_by_sufficient_statistics(x, x_mean = x_mean, x_SS = as.matrix(x_SS), x_N = x_N)
}

#' @rdname update_category_representation
#' @export
S7::method(update_category_representation, list(MNIX_CategoryRepresentation, S7::class_any, S7::class_any, S7::class_any)) <- function(
  x, x_N, x_mean, x_SS, ...
) {
  .update_MNIX_category_representation_by_sufficient_statistics(x, x_mean = x_mean, x_SS = as.matrix(x_SS), x_N = x_N)
}

#' @rdname update_template
#' @param x An `NIX_IdealAdaptor` object.
#' @param observations A data.frame or tibble containing cue columns and category label.
#' @param ... Additional options.
#' @return The updated `NIX_IdealAdaptor` model.
#' @seealso \code{\link{update_category_representation}}
#' @export
S7::method(update_template, list(NIX_IdealAdaptor, S7::class_any)) <- function(
  x, observations, updating = c("batch", "incremental"), keep_history = FALSE,
  lapse_treatment = "no_lapses", noise_treatment = "no_noise",
  update_method = "label-certain", ...
) {
  cue_labels <- get_cue_labels(x)
  .assert_data_frame_like(observations)
  .assert_data_contains_cols(observations, cue_labels)
  if (update_method == "label-certain") .assert_data_contains_cols(observations, "category")

  categories <- get_category_labels(x)
  reps <- x@category_template@representations

  suff_stats <- get_sufficient_category_statistics(
    observations,
    cues = cue_labels,
    category = "category",
    categories = categories,
    model_family = "NIX"
  )
  suff_by_cat <- split(suff_stats, suff_stats$category)

  updated <- lapply(seq_along(reps), function(i) {
    cat_name <- as.character(categories[i])
    st <- suff_by_cat[[cat_name]]
    if (is.null(st) || nrow(st) == 0L || st$x_N[[1L]] == 0L) {
      return(reps[[i]])
    }
    update_category_representation(
      reps[[i]],
      x_N = st$x_N[[1L]],
      x_mean = st$x_mean[[1L]],
      x_SS = st$x_css[[1L]]
    )
  })
  names(updated) <- names(reps)
  new_nix_ideal_adaptor(
    category_template = new_category_representation_template(updated, metadata = x@category_template@metadata),
    decision_rule = x@decision_rule, category_prior = x@category_prior,
    lapse_rate = get_lapse_rate(x), lapse_bias = get_lapse_bias(x),
    Sigma_noise = get_noise(x), noise_treatment = get_noise_treatment(x),
    lapse_treatment = get_lapse_treatment(x), metadata = x@metadata
  )
}

#' @rdname update_template
#' @param x An `MNIX_IdealAdaptor` object.
#' @param observations A data.frame or tibble containing cue columns and category label.
#' @param ... Additional options.
#' @return The updated `MNIX_IdealAdaptor` model.
#' @seealso \code{\link{update_category_representation}}
#' @export
S7::method(update_template, list(MNIX_IdealAdaptor, S7::class_any)) <- function(
  x, observations, updating = c("batch", "incremental"), keep_history = FALSE,
  lapse_treatment = "no_lapses", noise_treatment = "no_noise",
  update_method = "label-certain", ...
) {
  cue_labels <- get_cue_labels(x)
  .assert_data_frame_like(observations)
  .assert_data_contains_cols(observations, cue_labels)
  if (update_method == "label-certain") .assert_data_contains_cols(observations, "category")

  categories <- get_category_labels(x)
  reps <- x@category_template@representations

  suff_stats <- get_sufficient_category_statistics(
    observations,
    cues = cue_labels,
    category = "category",
    categories = categories,
    model_family = "MNIX"
  )
  suff_by_cat <- split(suff_stats, suff_stats$category)

  updated <- lapply(seq_along(reps), function(i) {
    cat_name <- as.character(categories[i])
    st <- suff_by_cat[[cat_name]]
    if (is.null(st) || nrow(st) == 0L || st$x_N[[1L]] == 0L) {
      return(reps[[i]])
    }
    update_category_representation(
      reps[[i]],
      x_N = st$x_N[[1L]],
      x_mean = st$x_mean[[1L]],
      x_SS = st$x_css[[1L]]
    )
  })
  names(updated) <- names(reps)
  new_mnix_ideal_adaptor(
    category_template = new_category_representation_template(updated, metadata = x@category_template@metadata),
    decision_rule = x@decision_rule, category_prior = x@category_prior,
    lapse_rate = get_lapse_rate(x), lapse_bias = get_lapse_bias(x),
    Sigma_noise = get_noise(x), noise_treatment = get_noise_treatment(x),
    lapse_treatment = get_lapse_treatment(x), metadata = x@metadata
  )
}
