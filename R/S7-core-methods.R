#' @include S7-core-classes.R
#' @include S7-generics.R
NULL

# S7 methods and method-local helpers for MVBeliefUpdatr.

# -------------------------
# Base object methods
# -------------------------

S7::method(get_model_family, MVBU_Object) <- function(x) {
  class(x)[1]
}

S7::method(get_metadata, MVBU_Object) <- function(x) {
  x@metadata
}

S7::method(construct_mvbu, MVBU_Object) <- function(x) {
  x
}

S7::method(validate_mvbu, MVBU_Object) <- function(x) {
  validate_object(x)
}

S7::method(summarize_mvbu, MVBU_Object) <- function(x) {
  list(
    class = class(x)[1],
    model_family = get_model_family(x)
  )
}

S7::method(print_mvbu, MVBU_Object) <- function(x) {
  cat("<", class(x)[1], ">\n", sep = "")
  invisible(x)
}

S7::method(plot_prep_mvbu, MVBU_Object) <- function(x) {
  .mvbu_not_implemented("plot_prep_mvbu", class(x)[1])
}

# -------------------------
# Representation accessors
# -------------------------

S7::method(get_category_likelihood_function, MVBU_CategoryRepresentation) <- function(x) {
  x@category_likelihood_function
}

S7::method(get_category_template, MVBU_CognitiveModel) <- function(x) {
  x@category_template
}

S7::method(get_category_representations, MVBU_CognitiveModel) <- function(x) {
  template <- S7::method(get_category_template, MVBU_CognitiveModel)(x)
  get_category_representations(template)
}

S7::method(get_category_likelihood_function, MVBU_CognitiveModel) <- function(x) {
  template <- S7::method(get_category_template, MVBU_CognitiveModel)(x)
  get_category_likelihood_function(template)
}

S7::method(get_category_representations, MVBU_CategoryRepresentationTemplate) <- function(x) {
  x@representations
}

S7::method(get_parameters, MVBU_Object) <- function(x) {
  .mvbu_not_implemented("get_parameters", class(x)[1])
}

S7::method(get_parameters, UVG_CategoryRepresentation) <- function(x) {
  list(
    mu = x@mu,
    sigma2 = x@sigma2
  )
}

S7::method(get_parameters, NIX_CategoryRepresentation) <- function(x) {
  list(
    m = x@m,
    kappa = x@kappa,
    nu = x@nu,
    sigma2 = x@sigma2
  )
}

S7::method(get_parameters, MUVG_CategoryRepresentation) <- function(x) {
  list(
    component_mu = x@component_mu,
    component_sigma2 = x@component_sigma2,
    component_weights = x@component_weights
  )
}

S7::method(get_parameters, MNIX_CategoryRepresentation) <- function(x) {
  list(
    component_m = x@component_m,
    component_kappa = x@component_kappa,
    component_nu = x@component_nu,
    component_sigma2 = x@component_sigma2,
    component_weights = x@component_weights
  )
}

S7::method(get_parameters, MVG_CategoryRepresentation) <- function(x) {
  list(
    mu = x@mu,
    Sigma = x@Sigma
  )
}

S7::method(get_parameters, NIW_CategoryRepresentation) <- function(x) {
  list(
    m = x@m,
    kappa = x@kappa,
    nu = x@nu,
    S = x@S
  )
}

S7::method(get_parameters, Exemplar_CategoryRepresentation) <- function(x) {
  list(
    exemplars = x@exemplars,
    exemplar_weights = x@exemplar_weights
  )
}

# -------------------------
# Model property accessors
# -------------------------

#' Normalize category-scoped values from legacy list/data-frame inputs.
#'
#' This helper resolves a single value, a vector of values, or a named vector
#' into a value vector aligned to the requested category labels. It is used by
#' the S7 compatibility methods for legacy objects so that category priors and
#' lapse biases can be read from older list/data-frame shapes without duplicating
#' the same coercion logic in each accessor.
#'
#' @deprecated This compatibility helper is only needed while legacy
#'   list/data-frame inputs are still supported. It can be removed once the
#'   S7 interface is the only supported representation.
#' @keywords internal
.mvbu_resolve_values_for_categories <- function(values, categories = NULL) {
  if (is.null(values)) {
    return(NULL)
  }

  values <- unlist(values, recursive = TRUE, use.names = TRUE)
  if (length(values) == 0) {
    return(NULL)
  }

  if (missing(categories) || is.null(categories) || length(categories) == 0) {
    return(as.numeric(values[1]))
  }

  categories <- as.character(categories)
  if (!is.null(names(values)) && any(nzchar(names(values)))) {
    matched <- values[match(categories, names(values))]
    if (!anyNA(matched)) {
      return(as.numeric(matched))
    }
  }

  if (length(values) == 1) {
    return(rep(as.numeric(values[1]), length(categories)))
  }

  if (length(values) == length(categories)) {
    return(as.numeric(values))
  }

  rep(as.numeric(values[1]), length(categories))
}

S7::method(get_model_family, MVBU_ModelDistribution) <- function(x) {
  x@model_family
}

S7::method(get_category_prior, list(MVBU_Object, S7::class_any)) <- function(x, categories) {
  .mvbu_not_implemented("get_category_prior", class(x)[1])
}

#' @deprecated Legacy compatibility method for list/data-frame inputs; remove once
#'   S7-only representations are required.
#' @keywords internal
S7::method(get_category_prior, list(S7::class_any, S7::class_any)) <- function(x, categories) {
  if (is.list(x) && !is.null(x[["prior"]])) {
    prior <- x[["prior"]]
  } else if (is.data.frame(x) && "prior" %in% names(x)) {
    prior <- x[["prior"]]
  } else if (is.list(x) && !is.null(x[["category"]]) && !is.null(x[["prior"]])) {
    prior <- x[["prior"]]
  } else {
    return(NULL)
  }

  if (is.list(x) && !is.null(x[["category"]]) && !is.null(x[["prior"]])) {
    prior_table <- x
    prior_values <- prior_table[["prior"]]
    category_values <- prior_table[["category"]]
    if (!missing(categories) && !is.null(categories)) {
      categories <- as.character(categories)
      category_values <- as.character(category_values)
      return(as.numeric(prior_values[match(categories, category_values)]))
    }
    return(as.numeric(prior_values))
  }

  .mvbu_resolve_values_for_categories(prior, categories)
}

S7::method(get_category_prior, list(MVBU_CognitiveModel, S7::class_any)) <- function(x, categories) {
  prior <- x@category_prior
  if (missing(categories) || is.null(categories)) {
    return(prior)
  }

  prior_names <- names(prior)
  if (!is.null(prior_names) && all(nzchar(prior_names))) {
    return(as.numeric(prior[match(as.character(categories), prior_names)]))
  }

  as.numeric(prior)
}

S7::method(get_lapse_rate, MVBU_Object) <- function(x) {
  .mvbu_not_implemented("get_lapse_rate", class(x)[1])
}

#' @deprecated Legacy compatibility method for list/data-frame inputs; remove once
#'   S7-only representations are required.
#' @keywords internal
S7::method(get_lapse_rate, S7::class_any) <- function(x) {
  if (is.list(x) && !is.null(x[["lapse_rate"]])) {
    lapse_rate <- x[["lapse_rate"]]
    if (is.list(lapse_rate) && length(lapse_rate) > 0) {
      lapse_rate <- lapse_rate[[1]]
    }
    if (is.null(lapse_rate)) {
      return(NULL)
    }
    return(as.numeric(lapse_rate[1]))
  }
  if (is.data.frame(x) && "lapse_rate" %in% names(x)) {
    lapse_rate <- x[["lapse_rate"]]
    if (is.list(lapse_rate) && length(lapse_rate) > 0) {
      lapse_rate <- lapse_rate[[1]]
    }
    if (is.null(lapse_rate)) {
      return(NULL)
    }
    return(as.numeric(lapse_rate[1]))
  }
  NULL
}

S7::method(get_lapse_rate, MVBU_CognitiveModel) <- function(x) {
  x@lapse_behavior$lapse_rate
}

S7::method(get_lapse_bias, list(MVBU_Object, S7::class_any)) <- function(x, categories) {
  .mvbu_not_implemented("get_lapse_bias", class(x)[1])
}

#' @deprecated Legacy compatibility method for list/data-frame inputs; remove once
#'   S7-only representations are required.
#' @keywords internal
S7::method(get_lapse_bias, list(S7::class_any, S7::class_any)) <- function(x, categories) {
  if (is.list(x) && !is.null(x[["lapse_bias"]])) {
    lapse_bias <- x[["lapse_bias"]]
  } else if (is.data.frame(x) && "lapse_bias" %in% names(x)) {
    lapse_bias <- x[["lapse_bias"]]
  } else {
    return(NULL)
  }

  .mvbu_resolve_values_for_categories(lapse_bias, categories)
}

S7::method(get_lapse_bias, list(MVBU_CognitiveModel, S7::class_any)) <- function(x, categories) {
  lapse_bias <- x@lapse_behavior$lapse_bias
  if (missing(categories) || is.null(categories)) {
    return(lapse_bias)
  }

  bias_names <- names(lapse_bias)
  if (!is.null(bias_names) && all(nzchar(bias_names))) {
    return(as.numeric(lapse_bias[match(as.character(categories), bias_names)]))
  }

  as.numeric(lapse_bias)
}

# -------------------------
# Label accessors
# -------------------------

S7::method(get_cue_labels, list(MVBU_CategoryRepresentation, S7::class_any)) <- function(x, indices) {
  cue_labels <- .mvbu_extract_label_metadata(x)$cue
  if (missing(indices) || is.null(indices)) {
    return(cue_labels)
  }
  cue_labels[indices]
}

S7::method(get_cue_labels, list(MVBU_CategoryRepresentationTemplate, S7::class_any)) <- function(x, indices) {
  cue_labels <- .mvbu_extract_label_metadata(x)$cue
  if (missing(indices) || is.null(indices)) {
    return(cue_labels)
  }
  cue_labels[indices]
}

S7::method(get_cue_labels, list(MVBU_CognitiveModel, S7::class_any)) <- function(x, indices) {
  template <- S7::method(get_category_template, MVBU_CognitiveModel)(x)
  if (missing(indices) || is.null(indices)) {
    return(S7::method(get_cue_labels, list(MVBU_CategoryRepresentationTemplate, S7::class_any))(template))
  }
  S7::method(get_cue_labels, list(MVBU_CategoryRepresentationTemplate, S7::class_any))(template, indices)
}

#' @deprecated Legacy compatibility method for list/data-frame inputs; remove once
#'   S7-only representations are required.
#' @keywords internal
S7::method(get_category_labels, list(S7::class_any, S7::class_any)) <- function(x, indices) {
  if (is.data.frame(x) && "category" %in% names(x)) {
    category_labels <- sort(unique(as.character(x[["category"]])))
  } else if (is.list(x) && !is.null(x[["category"]])) {
    category_labels <- sort(unique(as.character(x[["category"]])))
  } else {
    return(character(0))
  }

  if (missing(indices) || is.null(indices)) {
    return(category_labels)
  }
  category_labels[indices]
}

S7::method(get_category_labels, list(MVBU_CategoryRepresentation, S7::class_any)) <- function(x, indices) {
  category_labels <- sort(unique(.mvbu_extract_label_metadata(x)$category))
  if (missing(indices) || is.null(indices)) {
    return(category_labels)
  }
  category_labels[indices]
}

S7::method(get_category_labels, list(MVBU_CategoryRepresentationTemplate, S7::class_any)) <- function(x, indices) {
  category_labels <- sort(unique(.mvbu_extract_label_metadata(x)$category))
  if (missing(indices) || is.null(indices)) {
    return(category_labels)
  }
  category_labels[indices]
}

S7::method(get_category_labels, list(MVBU_CognitiveModel, S7::class_any)) <- function(x, indices) {
  template <- S7::method(get_category_template, MVBU_CognitiveModel)(x)
  if (missing(indices) || is.null(indices)) {
    return(S7::method(get_category_labels, list(MVBU_CategoryRepresentationTemplate, S7::class_any))(template))
  }
  S7::method(get_category_labels, list(MVBU_CategoryRepresentationTemplate, S7::class_any))(template, indices)
}

S7::method(get_group_labels, list(MVBU_Object, S7::class_any)) <- function(x, indices) {
  if (missing(indices) || is.null(indices)) {
    return(character(0))
  }
  character(0)[indices]
}

# NOTE: group labels currently identify model instances in combinations of
# models. Revisit later for richer grouped-model containers.
S7::method(get_group_labels, list(MVBU_ModelDistribution, S7::class_any)) <- function(x, indices) {
  group_labels <- x@group_label
  if (missing(indices) || is.null(indices)) {
    return(group_labels)
  }
  group_labels[indices]
}

# -------------------------
# Posterior and categorization methods
# -------------------------

.mvbu_category_names <- function(representations) {
  repr_names <- names(representations)
  if (!is.null(repr_names) && all(nzchar(repr_names))) {
    return(as.character(repr_names))
  }

  vapply(representations, function(r) {
    labels <- get_category_labels(r)
    if (length(labels) > 0) {
      as.character(labels[[1]])
    } else {
      ""
    }
  }, character(1))
}

.mvbu_prior_in_repr_order <- function(x, category_names) {
  prior <- as.numeric(get_category_prior(x))
  prior_names <- names(get_category_prior(x))

  if (!is.null(prior_names) && all(nzchar(prior_names)) && setequal(prior_names, category_names)) {
    prior <- prior[match(category_names, prior_names)]
  }

  prior
}

.mvbu_posterior_matrix <- function(x, new_data, categories = NULL, noise_treatment = "no_noise", lapse_treatment = "no_lapses") {
  representations <- get_category_representations(x)
  n_cat <- length(representations)
  category_names <- .mvbu_category_names(representations)
  prior <- .mvbu_prior_in_repr_order(x, category_names)

  first_d <- length(get_cue_labels(representations[[1]]))
  x_mat <- .as_observation_matrix(new_data, d = first_d, arg_name = "new_data")
  n_obs <- nrow(x_mat)

  log_lik <- matrix(NA_real_, nrow = n_obs, ncol = n_cat)
  for (j in seq_len(n_cat)) {
    rep_j <- representations[[j]]
    d_j <- length(get_cue_labels(rep_j))
    x_j <- .as_observation_matrix(new_data, d = d_j, arg_name = "new_data")
    lik_fn <- rep_j@category_likelihood_function
    log_lik[, j] <- as.numeric(lik_fn(
      x_j,
      log = TRUE,
      noise_treatment = noise_treatment,
      Sigma_noise = x@noise_behavior$Sigma_noise
    ))
  }

  log_joint <- sweep(log_lik, 2, log(prior), "+")
  log_norm <- .logsumexp_rows(log_joint)
  posterior <- exp(log_joint - log_norm)

  if (identical(lapse_treatment, "sample")) {
    lapse_rate <- as.numeric(x@lapse_behavior$lapse_rate)
    lapse_bias <- as.numeric(x@lapse_behavior$lapse_bias)
    if (length(lapse_bias) != n_cat) {
      stop("lapse_bias length must match the number of category representations.", call. = FALSE)
    }
    if (lapse_rate > 0) {
      for (i in seq_len(n_obs)) {
        if (stats::runif(1) < lapse_rate) {
          posterior[i, ] <- lapse_bias
        }
      }
    }
  } else if (identical(lapse_treatment, "marginalize")) {
    lapse_rate <- as.numeric(x@lapse_behavior$lapse_rate)
    lapse_bias <- as.numeric(x@lapse_behavior$lapse_bias)
    if (length(lapse_bias) != n_cat) {
      stop("lapse_bias length must match the number of category representations.", call. = FALSE)
    }
    posterior <- (1 - lapse_rate) * posterior + lapse_rate * matrix(lapse_bias, nrow = n_obs, ncol = n_cat, byrow = TRUE)
  }

  if (!is.null(categories)) {
    categories <- as.character(categories)
    category_idx <- match(categories, names(representations))
    posterior <- posterior[, category_idx, drop = FALSE]
    category_names <- category_names[category_idx]
  }

  colnames(posterior) <- category_names
  posterior
}

S7::method(get_category_posterior_function, list(MVBU_CognitiveModel, S7::class_any, S7::class_any)) <- function(x, noise_treatment, lapse_treatment) {
  if (missing(noise_treatment) && missing(lapse_treatment)) {
    noise_treatment <- x@noise_behavior$noise_treatment
    lapse_treatment <- x@lapse_behavior$lapse_treatment
  } else {
    noise_treatment <- if (missing(noise_treatment) || is.null(noise_treatment)) x@noise_behavior$noise_treatment else as.character(noise_treatment)
    lapse_treatment <- if (missing(lapse_treatment) || is.null(lapse_treatment)) x@lapse_behavior$lapse_treatment else as.character(lapse_treatment)
  }
  key <- paste(noise_treatment, lapse_treatment, sep = "__")

  if (is.null(x@category_posterior_functions[[key]])) {
    x@category_posterior_functions[[key]] <- function(new_data, categories = NULL) {
      .mvbu_posterior_matrix(x, new_data, categories = categories, noise_treatment = noise_treatment, lapse_treatment = lapse_treatment)
    }
  }

  x@category_posterior_functions[[key]]
}

S7::method(posterior, list(MVBU_CognitiveModel, S7::class_list, S7::class_any)) <- function(x, new_data, categories) {
  lapply(new_data, function(batch) posterior(x, batch, categories))
}

S7::method(posterior, list(MVBU_CognitiveModel, S7::class_any, S7::class_any)) <- function(x, new_data, categories) {
  pf <- S7::method(get_category_posterior_function, list(MVBU_CognitiveModel, S7::class_any, S7::class_any))(x)
  pf(new_data, categories = categories)
}

S7::method(categorize, list(MVBU_CognitiveModel, S7::class_list, S7::class_any)) <- function(x, new_data, decision_rule) {
  lapply(new_data, function(batch) categorize(x, batch, decision_rule))
}

S7::method(categorize, list(MVBU_CognitiveModel, S7::class_any, S7::class_any)) <- function(x, new_data, decision_rule) {
  if (missing(decision_rule) || is.null(decision_rule)) {
    decision_rule <- x@decision_rule
  }
  posterior_matrix <- posterior(x, new_data, categories = NULL)
  if (identical(decision_rule, "sampling")) {
    sampled <- vapply(seq_len(nrow(posterior_matrix)), function(i) {
      sample(colnames(posterior_matrix), size = 1, prob = posterior_matrix[i, ])
    }, character(1))
    chosen_idx <- match(sampled, colnames(posterior_matrix))
    probability <- posterior_matrix[cbind(seq_len(nrow(posterior_matrix)), chosen_idx)]
    return(data.frame(category = as.character(sampled), probability = as.numeric(probability), stringsAsFactors = FALSE))
  }
  if (identical(decision_rule, "criterion")) {
    best_idx <- apply(posterior_matrix, 1, which.max)
    category <- colnames(posterior_matrix)[best_idx]
    probability <- posterior_matrix[cbind(seq_len(nrow(posterior_matrix)), best_idx)]
    return(data.frame(category = as.character(category), probability = as.numeric(probability), stringsAsFactors = FALSE))
  }
  if (identical(decision_rule, "proportional")) {
    best_idx <- apply(posterior_matrix, 1, which.max)
    category <- colnames(posterior_matrix)[best_idx]
    probability <- posterior_matrix[cbind(seq_len(nrow(posterior_matrix)), best_idx)]
    return(data.frame(category = as.character(category), probability = as.numeric(probability), stringsAsFactors = FALSE))
  }
  stop("unsupported decision_rule", call. = FALSE)
}
