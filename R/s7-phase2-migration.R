# S7 migration adapters from legacy NIW/MVG/Exemplar structures to S7.
# These adapters are intended to be used internally by MVBeliefUpdatr and are not part of the public API.
# They will be removed after migration to S7 is complete and legacy structures are no longer supported.

#' Normalize a family name for phase-2 migration helpers
#'
#' @param family Family name.
#' @param allowed Allowed family names.
#' @return A normalized family name.
#' @keywords internal
#' @noRd
.normalize_phase2_family <- function(family, allowed) {
  .assert_non_NA_scalar_character(family, msg = paste0("family must be a non-empty scalar character value."))
  family <- toupper(family)
  if (!(family %in% allowed)) {
    .stop("Unsupported family for Phase 2 migration adapter.")
  }
  family
}

#' Validate a legacy table input for phase-2 migration helpers
#'
#' @param x Input table.
#' @param required Required column names.
#' @param context Context string for error messages.
#' @return The validated input object.
#' @keywords internal
#' @noRd
.validate_legacy_table <- function(x, required, context) {
  if (!is.data.frame(x)) {
    .stop("x must be a data.frame or tibble.")
  }
  missing_cols <- setdiff(required, names(x))
  if (length(missing_cols) > 0) {
    .stop(paste0("x must contain ", paste(required, collapse = ", "), " columns for ", context, " conversion."))
  }
}

#' Coerce legacy MVG rows to S7 representation objects
#' @keywords internal
#' @noRd
.as_s7_mvg_representations <- function(x, category = "category") {
  .assert_non_NA_scalar_character(category, msg = paste0("category must be a non-empty scalar character value."))
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
#' @noRd
.as_s7_niw_representations <- function(x, category = "category") {
  .assert_non_NA_scalar_character(category, msg = paste0("category must be a non-empty scalar character value."))
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
#' @noRd
.as_s7_exemplar_representations <- function(x, category = "category") {
  .assert_non_NA_scalar_character(category, msg = paste0("category must be a non-empty scalar character value."))
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
#' @noRd
.as_s7_muvg_representations <- function(x, category = "category") {
  .assert_non_NA_scalar_character(category, msg = paste0("category must be a non-empty scalar character value."))
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
#' @noRd
.as_s7_mnix_representations <- function(x, category = "category") {
  .assert_non_NA_scalar_character(category, msg = paste0("category must be a non-empty scalar character value."))
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
      paste0("cue", seq_along(x$component_m[[i]]))
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
#' @noRd
.as_s7_category_representation_template <- function(x, family, category = "category") {
  family <- .normalize_phase2_family(family, allowed = c("MVG", "NIW", "EXEMPLAR", "MUVG", "MNIX"))
  .assert_non_NA_scalar_character(category, msg = paste0("category must be a non-empty scalar character value."))

  reps <- switch(
    family,
    MVG = .as_s7_mvg_representations(x, category = category),
    NIW = .as_s7_niw_representations(x, category = category),
    EXEMPLAR = .as_s7_exemplar_representations(x, category = category),
    MUVG = .as_s7_muvg_representations(x, category = category),
    MNIX = .as_s7_mnix_representations(x, category = category),
    .stop("Unsupported family for Phase 2 migration adapter.")
  )

  new_category_representation_template(representations = reps)
}

#' @export
as_s7_category_representation_template <- function(x, family, category = "category") {
  .as_s7_category_representation_template(x = x, family = family, category = category)
}

#' @export
as_s7_mvg_representations <- function(x, category = "category") {
  .as_s7_mvg_representations(x = x, category = category)
}

#' @export
as_s7_niw_representations <- function(x, category = "category") {
  .as_s7_niw_representations(x = x, category = category)
}

#' @export
as_s7_exemplar_representations <- function(x, category = "category") {
  .as_s7_exemplar_representations(x = x, category = category)
}

#' @export
as_s7_muvg_representations <- function(x, category = "category") {
  .as_s7_muvg_representations(x = x, category = category)
}

#' @export
as_s7_mnix_representations <- function(x, category = "category") {
  .as_s7_mnix_representations(x = x, category = category)
}

.legacy_model_priors <- function(x, category) {
  .assert_non_NA_scalar_character(category, msg = paste0("category must be a non-empty scalar character value."))
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
  .assert_non_NA_scalar_character(category, msg = paste0("category must be a non-empty scalar character value."))
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
#' @noRd
.as_s7_mvg_ideal_observer <- function(x, category = "category", decision_rule = "sampling") {
  .assert_non_NA_scalar_character(category, msg = paste0("category must be a non-empty scalar character value."))
  .assert_non_NA_scalar_character(decision_rule, msg = paste0("decision_rule must be a non-empty scalar character value."))
  template <- .as_s7_category_representation_template(x, family = "MVG", category = category)
  new_mvg_ideal_observer(
    category_template = template,
    decision_rule = decision_rule,
    category_prior = .legacy_model_priors(x, category = category),
    lapse_rate = .legacy_model_lapse_rate(x),
    lapse_bias = .legacy_model_lapse_bias(x, category = category)
  )
}

#' @export
as_s7_mvg_ideal_observer <- function(x, category = "category", decision_rule = "sampling") {
  .as_s7_mvg_ideal_observer(x = x, category = category, decision_rule = decision_rule)
}

#' Coerce legacy NIW ideal adaptor-like object to S7 model object
#' @keywords internal
#' @noRd
.as_s7_niw_ideal_adaptor <- function(x, category = "category", decision_rule = "sampling") {
  .assert_non_NA_scalar_character(category, msg = paste0("category must be a non-empty scalar character value."))
  .assert_non_NA_scalar_character(decision_rule, msg = paste0("decision_rule must be a non-empty scalar character value."))
  template <- .as_s7_category_representation_template(x, family = "NIW", category = category)
  new_niw_ideal_adaptor(
    category_template = template,
    decision_rule = decision_rule,
    category_prior = .legacy_model_priors(x, category = category),
    lapse_rate = .legacy_model_lapse_rate(x),
    lapse_bias = .legacy_model_lapse_bias(x, category = category)
  )
}

#' @export
as_s7_niw_ideal_adaptor <- function(x, category = "category", decision_rule = "sampling") {
  .as_s7_niw_ideal_adaptor(x = x, category = category, decision_rule = decision_rule)
}

#' Coerce legacy exemplar model-like object to S7 model object
#' @keywords internal
#' @noRd
.as_s7_exemplar_model <- function(x, category = "category", decision_rule = "sampling") {
  .assert_non_NA_scalar_character(category, msg = paste0("category must be a non-empty scalar character value."))
  .assert_non_NA_scalar_character(decision_rule, msg = paste0("decision_rule must be a non-empty scalar character value."))
  template <- .as_s7_category_representation_template(x, family = "EXEMPLAR", category = category)
  new_exemplar_model(
    category_template = template,
    decision_rule = decision_rule,
    category_prior = .legacy_model_priors(x, category = category),
    lapse_rate = .legacy_model_lapse_rate(x),
    lapse_bias = .legacy_model_lapse_bias(x, category = category)
  )
}

#' @export
as_s7_exemplar_model <- function(x, category = "category", decision_rule = "sampling") {
  .as_s7_exemplar_model(x = x, category = category, decision_rule = decision_rule)
}

#' Coerce legacy MUVG ideal observer-like object to S7 model object
#' @keywords internal
#' @noRd
.as_s7_muvg_ideal_observer <- function(x, category = "category", decision_rule = "sampling") {
  .assert_non_NA_scalar_character(category, msg = paste0("category must be a non-empty scalar character value."))
  .assert_non_NA_scalar_character(decision_rule, msg = paste0("decision_rule must be a non-empty scalar character value."))
  template <- .as_s7_category_representation_template(x, family = "MUVG", category = category)
  new_muvg_ideal_observer(
    category_template = template,
    decision_rule = decision_rule,
    category_prior = .legacy_model_priors(x, category = category),
    lapse_rate = .legacy_model_lapse_rate(x),
    lapse_bias = .legacy_model_lapse_bias(x, category = category)
  )
}

#' @export
as_s7_muvg_ideal_observer <- function(x, category = "category", decision_rule = "sampling") {
  .as_s7_muvg_ideal_observer(x = x, category = category, decision_rule = decision_rule)
}

#' Coerce legacy MNIX ideal adaptor-like object to S7 model object
#' @keywords internal
#' @noRd
.as_s7_mnix_ideal_adaptor <- function(x, category = "category", decision_rule = "sampling") {
  .assert_non_NA_scalar_character(category, msg = paste0("category must be a non-empty scalar character value."))
  .assert_non_NA_scalar_character(decision_rule, msg = paste0("decision_rule must be a non-empty scalar character value."))
  template <- .as_s7_category_representation_template(x, family = "MNIX", category = category)
  new_mnix_ideal_adaptor(
    category_template = template,
    decision_rule = decision_rule,
    category_prior = .legacy_model_priors(x, category = category),
    lapse_rate = .legacy_model_lapse_rate(x),
    lapse_bias = .legacy_model_lapse_bias(x, category = category)
  )
}

#' @export
as_s7_mnix_ideal_adaptor <- function(x, category = "category", decision_rule = "sampling") {
  .as_s7_mnix_ideal_adaptor(x = x, category = category, decision_rule = decision_rule)
}

