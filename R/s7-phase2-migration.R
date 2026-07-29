# Phase 2 migration adapters from legacy NIW/MVG/Exemplar structures to S7.
# These adapters are intended to be used internally by MVBeliefUpdatr and are not part of the public API.
# They will be removed after migration to S7 is complete and legacy structures are no longer supported.

.normalize_phase2_scalar_character <- function(x, arg_name) {
  if (!is.character(x) || length(x) != 1 || nchar(x) == 0) {
    stop(paste0(arg_name, " must be a non-empty scalar character value."), call. = FALSE)
  }
  x
}

.normalize_phase2_family <- function(family, allowed) {
  family <- toupper(.normalize_phase2_scalar_character(family, "family"))
  if (!(family %in% allowed)) {
    stop("Unsupported family for Phase 2 migration adapter.", call. = FALSE)
  }
  family
}

.validate_legacy_table <- function(x, required, context) {
  if (!is.data.frame(x)) {
    stop("x must be a data.frame or tibble.", call. = FALSE)
  }
  missing_cols <- setdiff(required, names(x))
  if (length(missing_cols) > 0) {
    stop(paste0("x must contain ", paste(required, collapse = ", "), " columns for ", context, " conversion."), call. = FALSE)
  }
}

#' Coerce legacy MVG rows to S7 representation objects
#' @keywords internal
as_s7_mvg_representations <- function(x, category = "category") {
  category <- .normalize_phase2_scalar_character(category, "category")
  .validate_legacy_table(x, required = c(category, "mu", "Sigma"), context = "MVG")

  reps <- vector("list", length = nrow(x))
  names(reps) <- as.character(x[[category]])

  for (i in seq_len(nrow(x))) {
    mu_i <- x$mu[[i]]
    Sigma_i <- x$Sigma[[i]]
    cue_labels <- names(mu_i)
    if (is.null(cue_labels)) {
      cue_labels <- colnames(Sigma_i)
    }
    reps[[i]] <- new_mvg_category_representation(
      category_labels = as.character(x[[category]][i]),
      cue_labels = cue_labels,
      mu = mu_i,
      Sigma = Sigma_i
    )
  }

  reps
}

#' Coerce legacy NIW rows to S7 representation objects
#' @keywords internal
as_s7_niw_representations <- function(x, category = "category") {
  category <- .normalize_phase2_scalar_character(category, "category")
  required <- c(category, "m", "kappa", "nu", "S")
  .validate_legacy_table(x, required = required, context = "NIW")

  reps <- vector("list", length = nrow(x))
  names(reps) <- as.character(x[[category]])

  for (i in seq_len(nrow(x))) {
    m_i <- x$m[[i]]
    cue_labels <- names(m_i)
    if (is.null(cue_labels)) {
      cue_labels <- colnames(x$S[[i]])
    }

    reps[[i]] <- new_niw_category_representation(
      category_labels = as.character(x[[category]][i]),
      cue_labels = cue_labels,
      m = m_i,
      kappa = x$kappa[[i]],
      nu = x$nu[[i]],
      S = x$S[[i]]
    )
  }

  reps
}

#' Coerce legacy exemplar rows to S7 representation objects
#' @keywords internal
as_s7_exemplar_representations <- function(x, category = "category") {
  category <- .normalize_phase2_scalar_character(category, "category")
  .validate_legacy_table(x, required = c(category, "exemplars"), context = "exemplar")

  reps <- vector("list", length = nrow(x))
  names(reps) <- as.character(x[[category]])

  for (i in seq_len(nrow(x))) {
    ex_i <- x$exemplars[[i]]
    if (is.data.frame(ex_i)) {
      ex_i <- as.matrix(ex_i)
    }
    cue_labels <- colnames(ex_i)
    if (is.null(cue_labels)) {
      cue_labels <- paste0("cue", seq_len(ncol(ex_i)))
    }

    reps[[i]] <- new_exemplar_category_representation(
      category_labels = as.character(x[[category]][i]),
      cue_labels = cue_labels,
      exemplars = ex_i
    )
  }

  reps
}

#' Coerce legacy MUVG rows to S7 representation objects
#' @keywords internal
as_s7_muvg_representations <- function(x, category = "category") {
  category <- .normalize_phase2_scalar_character(category, "category")
  .validate_legacy_table(x, required = c(category, "component_mu", "component_sigma2"), context = "MUVG")

  reps <- vector("list", length = nrow(x))
  names(reps) <- as.character(x[[category]])

  for (i in seq_len(nrow(x))) {
    component_mu_i <- x$component_mu[[i]]
    component_sigma2_i <- x$component_sigma2[[i]]
    cue_labels <- names(component_mu_i)
    if (is.null(cue_labels)) {
      cue_labels <- names(component_sigma2_i)
    }
    if (is.null(cue_labels)) {
      cue_labels <- paste0("cue", seq_len(length(component_mu_i)))
    }

    if ("component_weights" %in% names(x)) {
      reps[[i]] <- new_muvg_category_representation(
        category_labels = as.character(x[[category]][i]),
        cue_labels = cue_labels,
        component_mu = component_mu_i,
        component_sigma2 = component_sigma2_i,
        component_weights = x$component_weights[[i]]
      )
    } else {
      reps[[i]] <- new_muvg_category_representation(
        category_labels = as.character(x[[category]][i]),
        cue_labels = cue_labels,
        component_mu = component_mu_i,
        component_sigma2 = component_sigma2_i
      )
    }
  }

  reps
}

#' Coerce legacy MNIX rows to S7 representation objects
#' @keywords internal
as_s7_mnix_representations <- function(x, category = "category") {
  category <- .normalize_phase2_scalar_character(category, "category")
  .validate_legacy_table(
    x,
    required = c(category, "component_m", "component_kappa", "component_nu", "component_sigma2"),
    context = "MNIX"
  )

  reps <- vector("list", length = nrow(x))
  names(reps) <- as.character(x[[category]])

  for (i in seq_len(nrow(x))) {
    cue_labels <- if ("cue_labels" %in% names(x)) {
      as.character(x$cue_labels[[i]])
    } else {
      "cue1"
    }

    if ("component_weights" %in% names(x)) {
      reps[[i]] <- new_mnix_category_representation(
        category_labels = as.character(x[[category]][i]),
        cue_labels = cue_labels,
        component_m = x$component_m[[i]],
        component_kappa = x$component_kappa[[i]],
        component_nu = x$component_nu[[i]],
        component_sigma2 = x$component_sigma2[[i]],
        component_weights = x$component_weights[[i]]
      )
    } else {
      reps[[i]] <- new_mnix_category_representation(
        category_labels = as.character(x[[category]][i]),
        cue_labels = cue_labels,
        component_m = x$component_m[[i]],
        component_kappa = x$component_kappa[[i]],
        component_nu = x$component_nu[[i]],
        component_sigma2 = x$component_sigma2[[i]]
      )
    }
  }

  reps
}

#' Build an S7 category-representation template from legacy family objects
#' @keywords internal
as_s7_category_representation_template <- function(x, family, category = "category") {
  family <- .normalize_phase2_family(family, allowed = c("MVG", "NIW", "EXEMPLAR", "MUVG", "MNIX"))
  category <- .normalize_phase2_scalar_character(category, "category")

  reps <- switch(
    family,
    MVG = as_s7_mvg_representations(x, category = category),
    NIW = as_s7_niw_representations(x, category = category),
    EXEMPLAR = as_s7_exemplar_representations(x, category = category),
    MUVG = as_s7_muvg_representations(x, category = category),
    MNIX = as_s7_mnix_representations(x, category = category),
    stop("Unsupported family for Phase 2 migration adapter.", call. = FALSE)
  )

  new_category_representation_template(representations = reps)
}

.legacy_model_priors <- function(x, category) {
  category <- .normalize_phase2_scalar_character(category, "category")
  if ("prior" %in% names(x)) {
    p <- as.numeric(x$prior)
  } else {
    n <- nrow(x)
    p <- rep(1 / n, n)
  }
  names(p) <- as.character(x[[category]])
  p
}

.legacy_model_lapse_rate <- function(x) {
  if ("lapse_rate" %in% names(x)) {
    as.numeric(x$lapse_rate[[1]])
  } else {
    0
  }
}

.legacy_model_lapse_bias <- function(x, category) {
  category <- .normalize_phase2_scalar_character(category, "category")
  if ("lapse_bias" %in% names(x)) {
    b <- as.numeric(x$lapse_bias)
  } else {
    n <- nrow(x)
    b <- rep(1 / n, n)
  }
  names(b) <- as.character(x[[category]])
  b
}

#' Coerce legacy MVG ideal observer-like object to S7 model object
#' @keywords internal
as_s7_mvg_ideal_observer <- function(x, category = "category", decision_rule = "sampling") {
  category <- .normalize_phase2_scalar_character(category, "category")
  decision_rule <- .normalize_phase2_scalar_character(decision_rule, "decision_rule")
  template <- as_s7_category_representation_template(x, family = "MVG", category = category)
  new_mvg_ideal_observer(
    category_template = template,
    decision_rule = decision_rule,
    category_prior = .legacy_model_priors(x, category = category),
    lapse_rate = .legacy_model_lapse_rate(x),
    lapse_bias = .legacy_model_lapse_bias(x, category = category)
  )
}

#' Coerce legacy NIW ideal adaptor-like object to S7 model object
#' @keywords internal
as_s7_niw_ideal_adaptor <- function(x, category = "category", decision_rule = "sampling") {
  category <- .normalize_phase2_scalar_character(category, "category")
  decision_rule <- .normalize_phase2_scalar_character(decision_rule, "decision_rule")
  template <- as_s7_category_representation_template(x, family = "NIW", category = category)
  new_niw_ideal_adaptor(
    category_template = template,
    decision_rule = decision_rule,
    category_prior = .legacy_model_priors(x, category = category),
    lapse_rate = .legacy_model_lapse_rate(x),
    lapse_bias = .legacy_model_lapse_bias(x, category = category)
  )
}

#' Coerce legacy exemplar model-like object to S7 model object
#' @keywords internal
as_s7_exemplar_model <- function(x, category = "category", decision_rule = "sampling") {
  category <- .normalize_phase2_scalar_character(category, "category")
  decision_rule <- .normalize_phase2_scalar_character(decision_rule, "decision_rule")
  template <- as_s7_category_representation_template(x, family = "EXEMPLAR", category = category)
  new_exemplar_model(
    category_template = template,
    decision_rule = decision_rule,
    category_prior = .legacy_model_priors(x, category = category),
    lapse_rate = .legacy_model_lapse_rate(x),
    lapse_bias = .legacy_model_lapse_bias(x, category = category)
  )
}

#' Coerce legacy MUVG ideal observer-like object to S7 model object
#' @keywords internal
as_s7_muvg_ideal_observer <- function(x, category = "category", decision_rule = "sampling") {
  category <- .normalize_phase2_scalar_character(category, "category")
  decision_rule <- .normalize_phase2_scalar_character(decision_rule, "decision_rule")
  template <- as_s7_category_representation_template(x, family = "MUVG", category = category)
  new_muvg_ideal_observer(
    category_template = template,
    decision_rule = decision_rule,
    category_prior = .legacy_model_priors(x, category = category),
    lapse_rate = .legacy_model_lapse_rate(x),
    lapse_bias = .legacy_model_lapse_bias(x, category = category)
  )
}

#' Coerce legacy MNIX ideal adaptor-like object to S7 model object
#' @keywords internal
as_s7_mnix_ideal_adaptor <- function(x, category = "category", decision_rule = "sampling") {
  category <- .normalize_phase2_scalar_character(category, "category")
  decision_rule <- .normalize_phase2_scalar_character(decision_rule, "decision_rule")
  template <- as_s7_category_representation_template(x, family = "MNIX", category = category)
  new_mnix_ideal_adaptor(
    category_template = template,
    decision_rule = decision_rule,
    category_prior = .legacy_model_priors(x, category = category),
    lapse_rate = .legacy_model_lapse_rate(x),
    lapse_bias = .legacy_model_lapse_bias(x, category = category)
  )
}

.legacy_distribution_payload <- function(x) {
  if (is.null(x)) {
    return(list())
  }

  if (is.list(x)) {
    return(x)
  }

  # Keep a non-list legacy fit object in cache so migration remains lossless.
  list(legacy_object = x)
}

#' Coerce a legacy NIW inferred/fit-like object to S7 model-distribution object
#' @keywords internal
as_s7_niw_model_distribution <- function(x = NULL, group_label = "") {
  group_label <- .normalize_phase2_scalar_character(as.character(group_label), "group_label")
  payload <- .legacy_distribution_payload(x)

  new_niw_model_distribution(
    cache = payload,
    metadata = list(
      migrated = TRUE,
      migrated_from = class(x),
      payload_names = names(payload)
    ),
    group_label = group_label
  )
}

#' Coerce a legacy MVG inferred/fit-like object to S7 model-distribution object
#' @keywords internal
as_s7_mvg_model_distribution <- function(x = NULL, group_label = "") {
  group_label <- .normalize_phase2_scalar_character(as.character(group_label), "group_label")
  payload <- .legacy_distribution_payload(x)

  new_mvg_model_distribution(
    cache = payload,
    metadata = list(
      migrated = TRUE,
      migrated_from = class(x),
      payload_names = names(payload)
    ),
    group_label = group_label
  )
}

#' Coerce a legacy exemplar inferred/fit-like object to S7 model-distribution object
#' @keywords internal
as_s7_exemplar_model_distribution <- function(x = NULL, group_label = "") {
  group_label <- .normalize_phase2_scalar_character(as.character(group_label), "group_label")
  payload <- .legacy_distribution_payload(x)

  new_exemplar_model_distribution(
    cache = payload,
    metadata = list(
      migrated = TRUE,
      migrated_from = class(x),
      payload_names = names(payload)
    ),
    group_label = group_label
  )
}

#' Coerce a legacy MUVG inferred/fit-like object to S7 model-distribution object
#' @keywords internal
as_s7_muvg_model_distribution <- function(x = NULL, group_label = "") {
  group_label <- .normalize_phase2_scalar_character(as.character(group_label), "group_label")
  payload <- .legacy_distribution_payload(x)

  new_muvg_model_distribution(
    cache = payload,
    metadata = list(
      migrated = TRUE,
      migrated_from = class(x),
      payload_names = names(payload)
    ),
    group_label = group_label
  )
}

#' Coerce a legacy MNIX inferred/fit-like object to S7 model-distribution object
#' @keywords internal
as_s7_mnix_model_distribution <- function(x = NULL, group_label = "") {
  group_label <- .normalize_phase2_scalar_character(as.character(group_label), "group_label")
  payload <- .legacy_distribution_payload(x)

  new_mnix_model_distribution(
    cache = payload,
    metadata = list(
      migrated = TRUE,
      migrated_from = class(x),
      payload_names = names(payload)
    ),
    group_label = group_label
  )
}

#' Coerce a legacy inferred/fit-like object to an S7 model-distribution object
#' @keywords internal
as_s7_model_distribution <- function(x = NULL, family, group_label = "") {
  family <- .normalize_phase2_family(family, allowed = c("NIW", "MVG", "EXEMPLAR", "MUVG", "MNIX"))
  group_label <- .normalize_phase2_scalar_character(as.character(group_label), "group_label")

  switch(
    family,
    NIW = as_s7_niw_model_distribution(x = x, group_label = group_label),
    MVG = as_s7_mvg_model_distribution(x = x, group_label = group_label),
    EXEMPLAR = as_s7_exemplar_model_distribution(x = x, group_label = group_label),
    MUVG = as_s7_muvg_model_distribution(x = x, group_label = group_label),
    MNIX = as_s7_mnix_model_distribution(x = x, group_label = group_label),
    stop("Unsupported family for Phase 2 model-distribution migration adapter.", call. = FALSE)
  )
}