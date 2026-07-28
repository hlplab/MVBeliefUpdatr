#' @include S7-classes.R S7-generics.R
NULL

# S7 methods and method-local helpers for MVBeliefUpdatr.

# -------------------------
# Minimal base methods
# -------------------------

S7::method(get_model_family, MVBU_Object) <- function(x) {
  class(x)[1]
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

S7::method(categorize_mvbu, list(MVBU_CognitiveModel, S7::class_any)) <- function(x, new_data) {
  get_category_posterior(x, new_data)
}

S7::method(predict_mvbu, list(MVBU_CognitiveModel, S7::class_any)) <- function(x, new_data) {
  get_category(x, new_data)
}

S7::method(posterior_mvbu, MVBU_Object) <- function(x) {
  get_posterior(x)
}

S7::method(plot_prep_mvbu, MVBU_Object) <- function(x) {
  .mvbu_not_implemented("plot_prep_mvbu", class(x)[1])
}

S7::method(get_model_family, MVBU_ModelDistribution) <- function(x) {
  x@model_family
}

S7::method(get_category_likelihood, MVBU_CognitiveModel) <- function(x) {
  x@category_likelihood
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

S7::method(get_category_prior, MVBU_Object) <- function(x) {
  .mvbu_not_implemented("get_category_prior", class(x)[1])
}

S7::method(get_category_prior, MVBU_CognitiveModel) <- function(x) {
  x@category_prior
}

S7::method(get_cue_labels, MVBU_CategoryRepresentation) <- function(x) {
  x@cue_labels
}

S7::method(get_category_labels, MVBU_CategoryRepresentation) <- function(x) {
  x@category_labels
}

S7::method(get_category_labels, MVBU_CategoryRepresentationTemplate) <- function(x) {
  .mvbu_category_names(x@representations)
}

S7::method(get_category_labels, MVBU_CognitiveModel) <- function(x) {
  get_category_labels(get_category_likelihood(x))
}

S7::method(get_group_labels, MVBU_Object) <- function(x) {
  character(0)
}

# NOTE: group labels currently identify model instances in combinations of
# models. Revisit later for richer grouped-model containers.
S7::method(get_group_labels, MVBU_ModelDistribution) <- function(x) {
  x@group_label
}

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

.mvbu_posterior_matrix <- function(x, new_data) {
  representations <- get_category_likelihood(x)@representations
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
    lik_fn <- rep_j@category_likelihood
    log_lik[, j] <- as.numeric(lik_fn(x_j, log = TRUE))
  }

  log_joint <- sweep(log_lik, 2, log(prior), "+")
  log_norm <- .logsumexp_rows(log_joint)
  posterior <- exp(log_joint - log_norm)

  colnames(posterior) <- category_names
  posterior
}

S7::method(get_category_posterior, list(MVBU_CognitiveModel, S7::class_list)) <- function(x, new_data) {
  lapply(new_data, function(batch) get_category_posterior(x, batch))
}

S7::method(get_category_posterior, list(MVBU_CognitiveModel, S7::class_any)) <- function(x, new_data) {
  .mvbu_posterior_matrix(x, new_data)
}

S7::method(get_category, list(MVBU_CognitiveModel, S7::class_list)) <- function(x, new_data) {
  lapply(new_data, function(batch) get_category(x, batch))
}

S7::method(get_category, list(MVBU_CognitiveModel, S7::class_any)) <- function(x, new_data) {
  posterior <- get_category_posterior(x, new_data)
  best_idx <- apply(posterior, 1, which.max)
  category <- colnames(posterior)[best_idx]
  probability <- posterior[cbind(seq_len(nrow(posterior)), best_idx)]

  data.frame(
    category = as.character(category),
    probability = as.numeric(probability),
    stringsAsFactors = FALSE
  )
}

S7::method(get_category_posterior_prediction, list(MVBU_CognitiveModel, S7::class_list)) <- function(x, new_data) {
  lapply(new_data, function(batch) get_category_posterior_prediction(x, batch))
}

S7::method(get_category_posterior_prediction, list(MVBU_CognitiveModel, S7::class_any)) <- function(x, new_data) {
  list(
    category_posterior = get_category_posterior(x, new_data),
    category = get_category(x, new_data)
  )
}