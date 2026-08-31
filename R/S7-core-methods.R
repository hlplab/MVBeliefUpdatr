#' @include S7-core-classes.R
#' @include S7-core-uvg-classes.R
#' @include S7-core-nix-classes.R
#' @include S7-core-muvg-classes.R
#' @include S7-core-mnix-classes.R
#' @include S7-core-mvg-classes.R
#' @include S7-core-niw-classes.R
#' @include S7-core-exemplar-classes.R
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

S7::method(print, MVBU_Object) <- function(x, ...) {
  cat("<", class(x)[1], ">\n", sep = "")
  invisible(x)
}

.format_param_inline <- function(p) {
  if (is.null(p)) {
    return("NULL")
  }
  if (length(p) == 1L && is.numeric(p)) {
    return(format(as.numeric(p), digits = 4))
  }
  if (is.matrix(p)) {
    return(paste0("matrix(", nrow(p), ", ", ncol(p), ")"))
  }
  if (is.array(p) && length(dim(p)) > 2) {
    return(paste0("array(", paste(dim(p), collapse = ", "), ")"))
  }
  if (is.numeric(p)) {
    return(paste0("vector(", length(p), ")"))
  }
  if (is.character(p) && length(p) == 1L) {
    return(paste0("'", p, "'"))
  }
  class(p)[1]
}

.format_category_representation_inline <- function(r) {
  if (S7::S7_inherits(r, UVG_CategoryRepresentation)) {
    sprintf(
      "UVG(mu = %s, sigma2 = %s)",
      format(r@mu, digits = 4),
      format(r@sigma2, digits = 4)
    )
  } else if (S7::S7_inherits(r, NIX_CategoryRepresentation)) {
    sprintf(
      "NIX(kappa = %s, nu = %s, m = %s, S = %s)",
      format(r@kappa, digits = 4),
      format(r@nu, digits = 4),
      format(r@m, digits = 4),
      format(r@S, digits = 4)
    )
  } else if (S7::S7_inherits(r, MVG_CategoryRepresentation)) {
    sprintf(
      "MVG(mu = %s, Sigma = %s)",
      .format_param_inline(r@mu),
      .format_param_inline(r@Sigma)
    )
  } else if (S7::S7_inherits(r, MUVG_CategoryRepresentation)) {
    sprintf(
      "MUVG(mu = %s, sigma2 = %s)",
      .format_param_inline(r@mu),
      .format_param_inline(r@sigma2)
    )
  } else if (S7::S7_inherits(r, NIW_CategoryRepresentation)) {
    sprintf(
      "NIW(kappa = %s, nu = %s, m = %s, S = %s)",
      format(r@kappa, digits = 4),
      format(r@nu, digits = 4),
      .format_param_inline(r@m),
      .format_param_inline(r@S)
    )
  } else if (S7::S7_inherits(r, MNIX_CategoryRepresentation)) {
    sprintf(
      "MNIX(kappa = %s, nu = %s, m = %s, S = %s)",
      .format_param_inline(r@kappa),
      .format_param_inline(r@nu),
      .format_param_inline(r@m),
      .format_param_inline(r@S)
    )
  } else if (S7::S7_inherits(r, Exemplar_CategoryRepresentation)) {
    ex <- r@exemplars
    n_pts <- if (is.null(ex)) 0L else nrow(ex)
    if (n_pts == 0L) {
      "Exemplar(N = 0)"
    } else {
      n_show <- min(5L, n_pts)
      pt_strs <- vapply(seq_len(n_show), function(i) {
        row_vals <- format(as.numeric(ex[i, ]), digits = 3)
        paste0("(", paste(trimws(row_vals), collapse = ", "), ")")
      }, character(1))
      if (n_pts > 5L) {
        sprintf("Exemplar(N = %d, %s, ...)", n_pts, paste(pt_strs, collapse = ", "))
      } else {
        sprintf("Exemplar(N = %d, %s)", n_pts, paste(pt_strs, collapse = ", "))
      }
    }
  } else {
    class(r)[1]
  }
}

S7::method(print, MVBU_CategoryRepresentation) <- function(x, ...) {
  cls <- class(x)[1]
  cat("<", cls, ">\n", sep = "")
  cats <- get_category_labels(x)
  cues <- get_cue_labels(x)
  cat("  Category: ", paste(cats, collapse = ", "), "\n", sep = "")
  cat("  Cues (", length(cues), "): ", paste(cues, collapse = ", "), "\n", sep = "")

  if (S7::S7_inherits(x, UVG_CategoryRepresentation)) {
    cat("  mu: ", format(x@mu, digits = 4), "\n", sep = "")
    cat("  sigma2: ", format(x@sigma2, digits = 4), "\n", sep = "")
  } else if (S7::S7_inherits(x, NIX_CategoryRepresentation)) {
    cat("  kappa: ", format(x@kappa, digits = 4), "\n", sep = "")
    cat("  nu: ", format(x@nu, digits = 4), "\n", sep = "")
    cat("  m: ", format(x@m, digits = 4), "\n", sep = "")
    cat("  S: ", format(x@S, digits = 4), "\n", sep = "")
  } else if (S7::S7_inherits(x, MVG_CategoryRepresentation)) {
    if (length(x@mu) == 1L) {
      cat("  mu: ", format(x@mu, digits = 4), "\n", sep = "")
    } else {
      cat("  mu: c(", paste(format(x@mu, digits = 4), collapse = ", "), ")\n", sep = "")
    }
    if (length(x@Sigma) == 1L) {
      cat("  Sigma: ", format(as.numeric(x@Sigma), digits = 4), "\n", sep = "")
    } else {
      cat("  Sigma:\n")
      print(x@Sigma)
    }
  } else if (S7::S7_inherits(x, NIW_CategoryRepresentation)) {
    cat("  kappa: ", format(x@kappa, digits = 4), "\n", sep = "")
    cat("  nu: ", format(x@nu, digits = 4), "\n", sep = "")
    if (length(x@m) == 1L) {
      cat("  m: ", format(x@m, digits = 4), "\n", sep = "")
    } else {
      cat("  m: c(", paste(format(x@m, digits = 4), collapse = ", "), ")\n", sep = "")
    }
    if (length(x@S) == 1L) {
      cat("  S: ", format(as.numeric(x@S), digits = 4), "\n", sep = "")
    } else {
      cat("  S:\n")
      print(x@S)
    }
  } else if (S7::S7_inherits(x, MNIX_CategoryRepresentation)) {
    if (length(x@kappa) == 1L) {
      cat("  kappa: ", format(x@kappa, digits = 4), "\n", sep = "")
    } else {
      cat("  kappa: c(", paste(format(x@kappa, digits = 4), collapse = ", "), ")\n", sep = "")
    }
    if (length(x@nu) == 1L) {
      cat("  nu: ", format(x@nu, digits = 4), "\n", sep = "")
    } else {
      cat("  nu: c(", paste(format(x@nu, digits = 4), collapse = ", "), ")\n", sep = "")
    }
    if (length(x@m) == 1L) {
      cat("  m: ", format(x@m, digits = 4), "\n", sep = "")
    } else {
      cat("  m: c(", paste(format(x@m, digits = 4), collapse = ", "), ")\n", sep = "")
    }
    if (length(x@S) == 1L) {
      cat("  S: ", format(as.numeric(x@S), digits = 4), "\n", sep = "")
    } else {
      cat("  S: c(", paste(format(x@S, digits = 4), collapse = ", "), ")\n", sep = "")
    }
  } else if (S7::S7_inherits(x, Exemplar_CategoryRepresentation)) {
    ex <- x@exemplars
    n_pts <- if (is.null(ex)) 0L else nrow(ex)
    cat("  Exemplars (", n_pts, " points):\n", sep = "")
    if (n_pts > 0L) {
      n_show <- min(5L, n_pts)
      for (i in seq_len(n_show)) {
        row_vals <- format(as.numeric(ex[i, ]), digits = 4)
        cat("    [", i, "] (", paste(trimws(row_vals), collapse = ", "), ")\n", sep = "")
      }
      if (n_pts > 5L) {
        cat("    ...\n")
      }
    }
  }
  invisible(x)
}

S7::method(print, MVBU_CategoryRepresentationTemplate) <- function(x, ...) {
  cls <- class(x)[1]
  cats <- get_category_labels(x)
  cues <- get_cue_labels(x)
  cat("<", cls, ">\n", sep = "")
  cat("  Categories (", length(cats), "): ", paste(cats, collapse = ", "), "\n", sep = "")
  cat("  Cues (", length(cues), "): ", paste(cues, collapse = ", "), "\n", sep = "")
  cat("  Representations (", length(x@representations), "):\n", sep = "")
  for (cat_name in names(x@representations)) {
    r <- x@representations[[cat_name]]
    cat("    $", cat_name, ": ", .format_category_representation_inline(r), "\n", sep = "")
  }
  invisible(x)
}

S7::method(print, MVBU_CognitiveModel) <- function(x, ...) {
  cls <- class(x)[1]
  cats <- get_category_labels(x)
  cues <- get_cue_labels(x)
  cat("<", cls, ">\n", sep = "")
  cat("  Categories (", length(cats), "): ", paste(cats, collapse = ", "), "\n", sep = "")
  cat("  Cues (", length(cues), "): ", paste(cues, collapse = ", "), "\n", sep = "")
  cat("  Decision rule: ", x@decision_rule, "\n", sep = "")

  prior <- get_category_prior(x)
  if (!is.null(prior)) {
    if (!is.null(names(prior))) {
      prior_str <- paste(paste(names(prior), format(prior, digits = 3), sep = ": "), collapse = ", ")
    } else {
      prior_str <- paste(format(prior, digits = 3), collapse = ", ")
    }
    cat("  Category prior: ", prior_str, "\n", sep = "")
  }

  lapse_r <- get_lapse_rate(x)
  lapse_trt <- x@lapse_behavior$lapse_treatment %||% "no_lapses"
  cat("  Lapse rate: ", format(lapse_r, digits = 3), " (treatment: ", lapse_trt, ")\n", sep = "")

  noise_trt <- x@noise_behavior$noise_treatment %||% "no_noise"
  if (is.null(x@noise_behavior$Sigma_noise)) {
    cat("  Perceptual noise: none (treatment: ", noise_trt, ")\n", sep = "")
  } else {
    sig_str <- .format_param_inline(x@noise_behavior$Sigma_noise)
    cat("  Perceptual noise: ", sig_str, " (treatment: ", noise_trt, ")\n", sep = "")
  }

  reps <- get_category_representations(x)
  cat("  Category representations (", length(reps), "):\n", sep = "")
  for (cat_name in names(reps)) {
    r <- reps[[cat_name]]
    cat("    $", cat_name, ": ", .format_category_representation_inline(r), "\n", sep = "")
  }
  invisible(x)
}

S7::method(summary, MVBU_Object) <- function(object, ...) {
  print(object, ...)
  invisible(object)
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

# Helper to map S7 class name to model family
.get_family_from_class_name <- function(class_name) {
  families <- .list_model_families()
  for (fam in families) {
    reg <- .get_model_family_entry(fam)
    if (class_name %in% c(reg$category_representation, reg$cognitive_model)) {
      return(fam)
    }
  }
  clean_name <- gsub(
    paste0(
      "_(IdealObserver|IdealAdaptor|IdealAdaptorStanfit|",
      "CategoryRepresentation|CategoryRepresentationTemplate|Model)$"
    ),
    "",
    class_name
  )
  if (toupper(clean_name) == "EXEMPLAR") return("EXEMPLAR")
  toupper(clean_name)
}

S7::method(get_model_type, MVBU_Object) <- function(x) {
  .get_family_from_class_name(class(x)[1])
}

S7::method(get_representation_type, MVBU_CategoryRepresentation) <- function(x) {
  .get_family_from_class_name(class(x)[1])
}

S7::method(get_representation_type, MVBU_CategoryRepresentationTemplate) <- function(x) {
  if (length(x@representations) == 0) return("")
  get_representation_type(x@representations[[1]])
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

S7::method(get_parameters, MVBU_CategoryRepresentationTemplate) <- function(x) {
  lapply(x@representations, get_parameters)
}

S7::method(get_parameters, MVBU_CognitiveModel) <- function(x) {
  get_parameters(x@category_template)
}

S7::method(get_parameter_names, MVBU_CategoryRepresentation) <- function(x, ...) {
  setdiff(names(S7::props(x)), c("metadata", "category_likelihood_function"))
}

S7::method(get_parameter_names, MVBU_CategoryRepresentationTemplate) <- function(x, ...) {
  unique(unlist(lapply(x@representations, get_parameter_names)))
}

S7::method(get_parameter_names, MVBU_CognitiveModel) <- function(x, ...) {
  c(
    get_parameter_names(x@category_template),
    "category_prior",
    "lapse_rate",
    "lapse_bias",
    "Sigma_noise"
  )
}


# -------------------------------------------------------------
# Expected parameter extractors (for belief distributions & models)
# -------------------------------------------------------------

S7::method(get_expected_mu, NIX_CategoryRepresentation) <- function(x) {
  x@m
}

S7::method(get_expected_mu, MNIX_CategoryRepresentation) <- function(x) {
  x@component_m
}

S7::method(get_expected_mu, NIW_CategoryRepresentation) <- function(x) {
  x@m
}

S7::method(get_expected_mu, MVBU_CategoryRepresentationTemplate) <- function(x) {
  lapply(x@representations, get_expected_mu)
}

S7::method(get_expected_mu, MVBU_CognitiveModel) <- function(x) {
  get_expected_mu(x@category_template)
}

S7::method(get_expected_sigma, NIX_CategoryRepresentation) <- function(x) {
  if (x@nu > 2) x@nu * x@sigma2 / (x@nu - 2) else NA_real_
}

S7::method(get_expected_sigma, MNIX_CategoryRepresentation) <- function(x) {
  diag(ifelse(x@component_nu > 2, x@component_nu * x@component_sigma2 / (x@component_nu - 2), NA_real_), nrow = length(x@component_nu))
}

S7::method(get_expected_sigma, NIW_CategoryRepresentation) <- function(x) {
  d <- nrow(x@S)
  if (x@nu > d + 1) x@S / (x@nu - d - 1) else matrix(NA_real_, nrow = d, ncol = d)
}

S7::method(get_expected_sigma, MVBU_CategoryRepresentationTemplate) <- function(x) {
  lapply(x@representations, get_expected_sigma)
}

S7::method(get_expected_sigma, MVBU_CognitiveModel) <- function(x) {
  get_expected_sigma(x@category_template)
}

S7::method(get_expected_category_statistic, MVBU_Object) <- function(
  x,
  statistic,
  ...
) {
  stat <- tolower(statistic)
  if (stat %in% c("mu", "m", "mean")) {
    get_expected_mu(x)
  } else if (
    stat %in% c("sigma", "sigma2", "s", "cov", "covariance", "var", "variance")
  ) {
    get_expected_sigma(x)
  } else {
    .stop(
      paste0("Unknown statistic '", statistic, "'. Expected 'mu' or 'sigma'.")
    )
  }
}


# -------------------------
# Model property accessors
# -------------------------

#' Normalize category-scoped values from legacy list/data-frame inputs
#'
#' This helper resolves a single value, a vector of values, or a named vector
#' into a value vector aligned to the requested category labels. It is used by
#' the S7 compatibility methods for legacy objects so that category priors and
#' lapse biases can be read from older list/data-frame shapes without duplicating
#' the same coercion logic in each accessor.
#'
#' @description Deprecated. This compatibility helper is only needed while legacy
#'   list/data-frame inputs are still supported. It can be removed once the
#'   S7 interface is the only supported representation.
#' @keywords internal
#' @noRd
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


S7::method(get_category_prior, list(MVBU_Object, S7::class_any)) <- function(x, categories) {
  .mvbu_not_implemented("get_category_prior", class(x)[1])
}

#' @name get_category_prior
#' @title Legacy compatibility method for category priors
#' @description Deprecated. Legacy compatibility method for list/data-frame inputs; remove once
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

#' @name get_lapse_rate
#' @title Legacy compatibility method for lapse rates
#' @description Deprecated. Legacy compatibility method for list/data-frame inputs; remove once
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

#' @name get_lapse_bias
#' @title Legacy compatibility method for lapse biases
#' @description Deprecated. Legacy compatibility method for list/data-frame inputs; remove once
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

S7::method(get_noise, MVBU_Object) <- function(x) {
  .mvbu_not_implemented("get_noise", class(x)[1])
}

#' @name get_noise
#' @title Get perceptual noise from a cognitive model
#' @description Extract perceptual noise covariance from a cognitive model.
#' @keywords internal
S7::method(get_noise, MVBU_CognitiveModel) <- function(x) {
  x@noise_behavior$Sigma_noise
}

S7::method(get_noise, S7::class_any) <- function(x) {
  if (is.list(x) && !is.null(x[["Sigma_noise"]])) {
    return(x[["Sigma_noise"]])
  }
  NULL
}

# -------------------------
# Label accessors
# -------------------------

S7::method(get_cue_labels, MVBU_CategoryRepresentation) <- function(x, indices = NULL, ...) {
  cue_labels <- .mvbu_extract_label_metadata(x)$cue
  if (missing(indices) || is.null(indices)) {
    return(cue_labels)
  }
  cue_labels[indices]
}

S7::method(get_cue_labels, MVBU_CategoryRepresentationTemplate) <- function(x, indices = NULL, ...) {
  cue_labels <- .mvbu_extract_label_metadata(x)$cue
  if (missing(indices) || is.null(indices)) {
    return(cue_labels)
  }
  cue_labels[indices]
}

S7::method(get_cue_labels, MVBU_CognitiveModel) <- function(x, indices = NULL, ...) {
  template <- S7::method(get_category_template, MVBU_CognitiveModel)(x)
  if (missing(indices) || is.null(indices)) {
    return(get_cue_labels(template))
  }
  get_cue_labels(template, indices = indices)
}

#' @name get_category_labels
#' @title Legacy compatibility method for category labels
#' @description Deprecated. Legacy compatibility method for list/data-frame inputs; remove once
#'   S7-only representations are required.
#' @keywords internal
S7::method(get_category_labels, S7::class_any) <- function(x, indices = NULL, ...) {
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

S7::method(get_category_labels, MVBU_CategoryRepresentation) <- function(x, indices = NULL, ...) {
  category_labels <- sort(unique(.mvbu_extract_label_metadata(x)$category))
  if (missing(indices) || is.null(indices)) {
    return(category_labels)
  }
  category_labels[indices]
}

S7::method(get_category_labels, MVBU_CategoryRepresentationTemplate) <- function(x, indices = NULL, ...) {
  category_labels <- sort(unique(.mvbu_extract_label_metadata(x)$category))
  if (missing(indices) || is.null(indices)) {
    return(category_labels)
  }
  category_labels[indices]
}

S7::method(get_category_labels, MVBU_CognitiveModel) <- function(x, indices = NULL, ...) {
  template <- S7::method(get_category_template, MVBU_CognitiveModel)(x)
  if (missing(indices) || is.null(indices)) {
    return(get_category_labels(template))
  }
  get_category_labels(template, indices = indices)
}

S7::method(get_group_labels, MVBU_Object) <- function(x, indices = NULL, ...) {
  group_labels <- .mvbu_extract_label_metadata(x)$group
  if (missing(indices) || is.null(indices)) {
    return(group_labels)
  }
  group_labels[indices]
}


S7::method(get_labels, MVBU_Object) <- function(x, ...) {
  list(
    cue = get_cue_labels(x, ...),
    category = get_category_labels(x, ...),
    group = get_group_labels(x, ...)
  )
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

.mvbu_likelihood_matrix <- function(x, new_data, categories = NULL, log = FALSE, noise_treatment = "no_noise", Sigma_noise = NULL) {
  representations <- get_category_representations(x)
  n_cat <- length(representations)
  category_names <- .mvbu_category_names(representations)

  first_d <- length(get_cue_labels(representations[[1]]))
  n_obs <- nrow(.as_observation_matrix(new_data, d = first_d, arg_name = "new_data"))

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
      Sigma_noise = Sigma_noise
    ))
  }

  if (!is.null(categories)) {
    category_idx <- match(as.character(categories), category_names)
    log_lik <- log_lik[, category_idx, drop = FALSE]
    category_names <- category_names[category_idx]
  }

  colnames(log_lik) <- category_names
  if (log) log_lik else exp(log_lik)
}

.mvbu_posterior_matrix <- function(x, new_data, categories = NULL, noise_treatment = "no_noise", lapse_treatment = "no_lapses") {
  representations <- get_category_representations(x)
  n_cat <- length(representations)
  category_names <- .mvbu_category_names(representations)
  prior <- .mvbu_prior_in_repr_order(x, category_names)

  log_lik <- .mvbu_likelihood_matrix(
    x, new_data,
    log = TRUE,
    noise_treatment = noise_treatment,
    Sigma_noise = x@noise_behavior$Sigma_noise
  )
  n_obs <- nrow(log_lik)

  log_joint <- sweep(log_lik, 2, log(prior), "+")
  log_norm <- .logsumexp_rows(log_joint)
  posterior <- exp(log_joint - log_norm)

  if (identical(lapse_treatment, "sample")) {
    lapse_rate <- as.numeric(x@lapse_behavior$lapse_rate)
    lapse_bias <- as.numeric(x@lapse_behavior$lapse_bias)
    if (length(lapse_bias) != n_cat) {
      .stop("lapse_bias length must match the number of category representations.")
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
      .stop("lapse_bias length must match the number of category representations.")
    }
    posterior <- (1 - lapse_rate) * posterior + lapse_rate * matrix(lapse_bias, nrow = n_obs, ncol = n_cat, byrow = TRUE)
  }

  if (!is.null(categories)) {
    category_idx <- match(as.character(categories), category_names)
    posterior <- posterior[, category_idx, drop = FALSE]
    category_names <- category_names[category_idx]
  }

  colnames(posterior) <- category_names
  posterior
}

S7::method(get_category_likelihood_function, MVBU_CategoryRepresentationTemplate) <- function(x) {
  representations <- get_category_representations(x)
  function(new_data, log = FALSE, noise_treatment = "no_noise", Sigma_noise = NULL, categories = NULL) {
    .mvbu_likelihood_matrix(
      x, new_data,
      categories = categories, log = log,
      noise_treatment = noise_treatment, Sigma_noise = Sigma_noise
    )
  }
}

S7::method(likelihood, list(MVBU_CategoryRepresentation, S7::class_any, S7::class_any)) <- function(x, new_data, categories) {
  lik_fn <- x@category_likelihood_function
  d <- length(get_cue_labels(x))
  new_data <- .as_observation_matrix(new_data, d = d, arg_name = "new_data")
  as.numeric(lik_fn(new_data, log = FALSE, noise_treatment = "no_noise", Sigma_noise = NULL))
}

S7::method(likelihood, list(MVBU_CategoryRepresentationTemplate, S7::class_any, S7::class_any)) <- function(x, new_data, categories) {
  if (missing(categories)) categories <- NULL
  .mvbu_likelihood_matrix(x, new_data, categories = categories, log = FALSE)
}

S7::method(likelihood, list(MVBU_CognitiveModel, S7::class_list, S7::class_any)) <- function(x, new_data, categories) {
  if (missing(categories)) categories <- NULL
  lapply(new_data, function(batch) likelihood(x, batch, categories))
}

S7::method(likelihood, list(MVBU_CognitiveModel, S7::class_any, S7::class_any)) <- function(x, new_data, categories) {
  if (missing(categories)) categories <- NULL
  .mvbu_likelihood_matrix(
    x, new_data,
    categories = categories,
    log = FALSE,
    noise_treatment = x@noise_behavior$noise_treatment,
    Sigma_noise = x@noise_behavior$Sigma_noise
  )
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
  if (missing(categories)) categories <- NULL
  lapply(new_data, function(batch) posterior(x, batch, categories))
}

S7::method(posterior, list(MVBU_CognitiveModel, S7::class_any, S7::class_any)) <- function(x, new_data, categories) {
  if (missing(categories)) categories <- NULL
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
  categories <- colnames(posterior_matrix)
  n_obs <- nrow(posterior_matrix)

  chosen_idx <- if (identical(decision_rule, "sampling")) {
    vapply(seq_len(n_obs), function(i) sample.int(length(categories), 1L, prob = posterior_matrix[i, ]), integer(1))
  } else if (identical(decision_rule, "criterion") || identical(decision_rule, "proportional")) {
    apply(posterior_matrix, 1, which.max)
  } else {
    .stop("Invalid decision_rule: ", decision_rule, ". Must be one of 'criterion', 'sampling', or 'proportional'.")
  }

  data.frame(
    category = as.character(categories[chosen_idx]),
    probability = as.numeric(posterior_matrix[cbind(seq_len(n_obs), chosen_idx)]),
    stringsAsFactors = FALSE
  )
}

#' @rdname evaluate_model
#' @export
S7::method(evaluate_model, MVBU_CognitiveModel) <- function(
  model,
  x = NULL,
  response_category = NULL,
  method = "log_lik",
  decision_rule = if (identical(method, "accuracy")) "criterion" else "proportional",
  return_by_x = FALSE,
  ...
) {
  valid_methods <- c("log_lik", "log_lik_permutation_constant", "likelihood-up-to-constant", "accuracy")
  .assert_that(all(method %in% valid_methods),
    msg = paste0("method must be one or more of: ", paste(valid_methods, collapse = ", "))
  )
  if (is.null(decision_rule)) {
    stop("decision_rule must be specified (e.g. 'proportional', 'criterion', or 'sampling').")
  }

  if ("likelihood-up-to-constant" %in% method) {
    lifecycle::deprecate_warn(
      "0.1.0",
      "evaluate_model(method = 'likelihood-up-to-constant')",
      "evaluate_model(method = 'log_lik')",
      always = TRUE
    )
  }

  cue_labels <- get_cue_labels(model)
  d_cue <- length(cue_labels)
  mat_x <- .as_observation_matrix(x, d = d_cue, arg_name = "x")
  if (is.null(colnames(mat_x))) {
    colnames(mat_x) <- cue_labels
  }
  n_obs <- nrow(mat_x)

  if (is.data.frame(response_category)) {
    if ("category" %in% names(response_category)) {
      response_category <- response_category$category
    } else if ("response" %in% names(response_category)) {
      response_category <- response_category$response
    } else {
      response_category <- response_category[[1]]
    }
  }

  if (length(response_category) != n_obs) {
    stop("Input x and response_category must have the same number of observations.")
  }

  # Predicted posterior probabilities
  P <- posterior(model, mat_x)
  category_names <- colnames(P)
  n_cat <- ncol(P)

  # Apply decision rule
  if (identical(decision_rule, "criterion")) {
    P_dec <- matrix(0, nrow = n_obs, ncol = n_cat, dimnames = list(NULL, category_names))
    for (i in seq_len(n_obs)) {
      max_val <- max(P[i, ])
      best <- which(P[i, ] == max_val)
      P_dec[i, best] <- 1 / length(best)
    }
    P <- P_dec
  }

  resp_char <- as.character(response_category)
  resp_idx <- match(resp_char, category_names)
  if (any(is.na(resp_idx))) {
    stop("Some observed responses in response_category do not match category labels in the model.")
  }

  # For each response, get its probability under the model
  p_obs <- vapply(seq_len(n_obs), function(i) P[i, resp_idx[i]], numeric(1))

  x_key <- apply(mat_x, 1, paste, collapse = "___")
  unique_keys <- unique(x_key)
  u_idx_list <- split(seq_len(n_obs), factor(x_key, levels = unique_keys))
  n_by_x <- vapply(u_idx_list, length, integer(1))
  u_mat_x <- mat_x[match(unique_keys, x_key), , drop = FALSE]

  r <- list()

  if ("accuracy" %in% method) {
    if (return_by_x) {
      acc_by_x <- vapply(u_idx_list, function(idxs) mean(p_obs[idxs]), numeric(1))
      res_df <- tibble::as_tibble(u_mat_x)
      res_df$N <- n_by_x
      res_df$accuracy <- acc_by_x
      r[["accuracy"]] <- res_df
    } else {
      r[["accuracy"]] <- mean(p_obs)
    }
  }

  # Categorical (trial-by-trial) log likelihood
  if (any(c("log_lik", "likelihood-up-to-constant") %in% method)) {
    if (return_by_x) {
      ll_by_x <- vapply(u_idx_list, function(idxs) {
        p_sub <- p_obs[idxs]
        if (any(p_sub <= 0)) -Inf else sum(log(p_sub))
      }, numeric(1))
      res_df <- tibble::as_tibble(u_mat_x)
      res_df$N <- n_by_x
      res_df$log_lik <- ll_by_x
      if ("log_lik" %in% method) {
        r[["log_lik"]] <- res_df
      }
      if ("likelihood-up-to-constant" %in% method) {
        res_df_legacy <- res_df
        names(res_df_legacy)[names(res_df_legacy) == "log_lik"] <- "log_likelihood"
        r[["likelihood-up-to-constant"]] <- res_df_legacy
      }
    } else {
      val <- if (any(p_obs <= 0)) -Inf else sum(log(p_obs))
      if ("log_lik" %in% method) {
        r[["log_lik"]] <- val
      }
      if ("likelihood-up-to-constant" %in% method) {
        r[["likelihood-up-to-constant"]] <- val
      }
    }
  }

  # Combinatorial permutation constant
  if ("log_lik_permutation_constant" %in% method) {
    const_by_x <- vapply(u_idx_list, function(idxs) {
      n_u <- length(idxs)
      counts <- tabulate(resp_idx[idxs], nbins = n_cat)
      lfactorial(n_u) - sum(lfactorial(counts))
    }, numeric(1))

    if (return_by_x) {
      res_df <- tibble::as_tibble(u_mat_x)
      res_df$N <- n_by_x
      res_df$log_lik_permutation_constant <- const_by_x
      r[["log_lik_permutation_constant"]] <- res_df
    } else {
      r[["log_lik_permutation_constant"]] <- sum(const_by_x)
    }
  }

  if (length(r) == 1) {
    return(r[[1]])
  }
  r
}
