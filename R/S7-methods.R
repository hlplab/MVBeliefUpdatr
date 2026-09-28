#' @include S7-class.R
#' @include S7-class-uvg.R
#' @include S7-class-nix.R
#' @include S7-class-muvg.R
#' @include S7-class-mnix.R
#' @include S7-class-mvg.R
#' @include S7-class-niw.R
#' @include S7-class-exemplar.R
#' @include S7-generics.R
NULL

# S7 methods and method-local helpers for MVBeliefUpdatr.

# -------------------------
# Base object methods
# -------------------------

S7::method(get_model_family, MVBU_Object) <- function(x, ...) {
  class(x)[1]
}

S7::method(get_metadata, MVBU_Object) <- function(x, ...) {
  x@metadata
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
      "NIX(kappa = %s, nu = %s, m = %s, sigma2 = %s)",
      format(r@kappa, digits = 4),
      format(r@nu, digits = 4),
      format(r@m, digits = 4),
      format(r@sigma2, digits = 4)
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
      "MNIX(kappa = %s, nu = %s, m = %s, sigma2 = %s)",
      .format_param_inline(r@kappa),
      .format_param_inline(r@nu),
      .format_param_inline(r@m),
      .format_param_inline(r@sigma2)
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
  if (length(cats) > 0) {
    cat("  Category: ", paste(cats, collapse = ", "), "\n", sep = "")
  }
  if (length(cues) > 0) {
    cat("  Cues (", length(cues), "): ", paste(cues, collapse = ", "), "\n", sep = "")
  }
  if (S7::S7_inherits(x, UVG_CategoryRepresentation)) {
    cat("  mu: ", format(x@mu, digits = 4), "\n", sep = "")
    cat("  sigma2: ", format(x@sigma2, digits = 4), "\n", sep = "")
  } else if (S7::S7_inherits(x, MUVG_CategoryRepresentation)) {
    cat("  mu: c(", paste(format(x@mu, digits = 4), collapse = ", "), ")\n", sep = "")
    cat("  sigma2: c(", paste(format(x@sigma2, digits = 4), collapse = ", "), ")\n", sep = "")
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
    cat("  kappa: c(", paste(format(x@kappa, digits = 4), collapse = ", "), ")\n", sep = "")
    cat("  nu: c(", paste(format(x@nu, digits = 4), collapse = ", "), ")\n", sep = "")
    cat("  m: c(", paste(format(x@m, digits = 4), collapse = ", "), ")\n", sep = "")
    cat("  sigma2: c(", paste(format(x@sigma2, digits = 4), collapse = ", "), ")\n", sep = "")
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

.format_moments_subline <- function(r) {
  if (S7::S7_inherits(r, NIW_CategoryRepresentation)) {
    e_mu <- .format_param_inline(get_expected_mu(r))
    e_sig <- .format_param_inline(get_expected_sigma(r))
    m_sig <- .format_param_inline(get_marginal_sigma(r))
    sprintf("E[mu] = %s, E[Sigma] = %s, Sigma_marg = %s", e_mu, e_sig, m_sig)
  } else if (S7::S7_inherits(r, NIX_CategoryRepresentation)) {
    e_mu <- format(get_expected_mu(r), digits = 4)
    e_sig <- format(get_expected_sigma(r), digits = 4)
    m_sig <- format(get_marginal_sigma(r), digits = 4)
    sprintf("E[mu] = %s, E[sigma2] = %s, sigma2_marg = %s", e_mu, e_sig, m_sig)
  } else if (S7::S7_inherits(r, MNIX_CategoryRepresentation)) {
    e_mu <- .format_param_inline(get_expected_mu(r))
    e_sig <- .format_param_inline(get_expected_sigma(r))
    m_sig <- .format_param_inline(get_marginal_sigma(r))
    sprintf("E[mu] = %s, E[Sigma] = %s, Sigma_marg = %s", e_mu, e_sig, m_sig)
  } else {
    NULL
  }
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
    cat("    ", cat_name, ": ", .format_category_representation_inline(r), "\n", sep = "")
    mom_line <- .format_moments_subline(r)
    if (!is.null(mom_line)) {
      cat("      ", mom_line, "\n", sep = "")
    }
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
  lapse_trt <- get_lapse_treatment(x)
  cat("  Lapse rate: ", format(lapse_r, digits = 3), " (treatment: ", lapse_trt, ")\n", sep = "")

  noise_trt <- get_noise_treatment(x)
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
    cat("    ", cat_name, ": ", .format_category_representation_inline(r), "\n", sep = "")
    mom_line <- .format_moments_subline(r)
    if (!is.null(mom_line)) {
      cat("      ", mom_line, "\n", sep = "")
    }
  }
  invisible(x)
}

.build_moments_table <- function(cats, cues, mu_list, sig_list) {
  rows <- list()
  cat_names <- as.character(cats)
  for (k in seq_along(cat_names)) {
    c_name <- cat_names[k]
    mu_c <- if (!is.null(names(mu_list)) && c_name %in% names(mu_list)) {
      mu_list[[c_name]]
    } else if (k <= length(mu_list)) {
      mu_list[[k]]
    } else {
      NULL
    }
    sig_c <- if (!is.null(names(sig_list)) && c_name %in% names(sig_list)) {
      sig_list[[c_name]]
    } else if (k <= length(sig_list)) {
      sig_list[[k]]
    } else {
      NULL
    }
    if (is.null(mu_c) || is.null(sig_c)) next

    for (i in seq_along(cues)) {
      rows[[length(rows) + 1L]] <- data.frame(
        Category = c_name,
        Parameter = "mu",
        Cue1 = cues[i],
        Cue2 = "",
        Value = as.numeric(mu_c)[i],
        stringsAsFactors = FALSE
      )
    }

    sig_mat <- as.matrix(sig_c)
    for (i in seq_along(cues)) {
      for (j in i:length(cues)) {
        rows[[length(rows) + 1L]] <- data.frame(
          Category = c_name,
          Parameter = "Sigma",
          Cue1 = cues[i],
          Cue2 = cues[j],
          Value = sig_mat[i, j],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0L) {
    return(NULL)
  }
  do.call(rbind, rows)
}

S7::method(summary, MVBU_Object) <- function(object, ...) {
  print(object, ...)
  invisible(object)
}

S7::method(summary, NIX_CategoryRepresentation) <- function(object, ...) {
  print(object, ...)
  cat("\nExpected Moments:\n")
  cat("  mu    : ", paste(get_expected_mu(object), collapse = ", "), "\n", sep = "")
  cat("  sigma2: ", paste(get_expected_sigma(object), collapse = ", "), "\n", sep = "")
  cat("Marginal Moments:\n")
  cat("  mu    : ", paste(get_marginal_mu(object), collapse = ", "), "\n", sep = "")
  cat("  sigma2: ", paste(get_marginal_sigma(object), collapse = ", "), "\n", sep = "")
  invisible(object)
}

S7::method(summary, MNIX_CategoryRepresentation) <- function(object, ...) {
  print(object, ...)
  cat("\nExpected Moments:\n")
  cat("  mu    : ", paste(get_expected_mu(object), collapse = ", "), "\n", sep = "")
  cat("  Sigma :\n")
  print(get_expected_sigma(object))
  cat("Marginal Moments:\n")
  cat("  mu    : ", paste(get_marginal_mu(object), collapse = ", "), "\n", sep = "")
  cat("  Sigma :\n")
  print(get_marginal_sigma(object))
  invisible(object)
}

S7::method(summary, NIW_CategoryRepresentation) <- function(object, ...) {
  print(object, ...)
  cat("\nExpected Moments:\n")
  cat("  mu   : ", paste(get_expected_mu(object), collapse = ", "), "\n", sep = "")
  cat("  Sigma:\n")
  print(get_expected_sigma(object))
  cat("Marginal Moments:\n")
  cat("  mu   : ", paste(get_marginal_mu(object), collapse = ", "), "\n", sep = "")
  cat("  Sigma:\n")
  print(get_marginal_sigma(object))
  invisible(object)
}

S7::method(summary, MVBU_CategoryRepresentationTemplate) <- function(object, ...) {
  print(object, ...)
  first_cls <- if (length(object@representations) > 0L) class(object@representations[[1L]])[1] else ""
  if (grepl("NIX|MNIX|NIW", first_cls, ignore.case = TRUE)) {
    cats <- get_category_labels(object)
    cues <- get_cue_labels(object)

    exp_tab <- .build_moments_table(cats, cues, get_expected_mu(object), get_expected_sigma(object))
    if (!is.null(exp_tab) && nrow(exp_tab) > 0) {
      cat("\nExpected Category Moments:\n")
      print(exp_tab, row.names = FALSE)
    }

    marg_tab <- .build_moments_table(cats, cues, get_marginal_mu(object), get_marginal_sigma(object))
    if (!is.null(marg_tab) && nrow(marg_tab) > 0) {
      cat("\nMarginal Category Moments:\n")
      print(marg_tab, row.names = FALSE)
    }
  }
  invisible(object)
}

S7::method(summary, MVBU_CognitiveModel) <- function(object, ...) {
  print(object, ...)
  fam <- get_model_family(object)
  first_cls <- tryCatch(class(object@category_template@representations[[1L]])[1], error = function(e) "")
  if (grepl("NIX|MNIX|NIW", fam, ignore.case = TRUE) || grepl("NIX|MNIX|NIW", first_cls, ignore.case = TRUE)) {
    cats <- get_category_labels(object)
    cues <- get_cue_labels(object)

    exp_tab <- .build_moments_table(cats, cues, get_expected_mu(object), get_expected_sigma(object))
    if (!is.null(exp_tab) && nrow(exp_tab) > 0) {
      cat("\nExpected Category Moments:\n")
      print(exp_tab, row.names = FALSE)
    }

    marg_tab <- .build_moments_table(cats, cues, get_marginal_mu(object), get_marginal_sigma(object))
    if (!is.null(marg_tab) && nrow(marg_tab) > 0) {
      cat("\nMarginal Category Moments:\n")
      print(marg_tab, row.names = FALSE)
    }
  }
  invisible(object)
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

S7::method(get_category_representations, MVBU_CategoryRepresentationTemplate) <- function(x) {
  x@representations
}

# Helper to map S7 class name to model family
.get_family_from_class_name <- function(class_name) {
  class_name <- sub("^.*::", "", class_name)
  families <- list_model_families()
  for (fam in families) {
    reg <- .get_model_family_registration(fam)
    if (class_name %in% c(reg$category_representation, reg$cognitive_model)) {
      return(fam)
    }
  }
  clean_name <- gsub(
    paste0(
      "_(IdealObserver|IdealAdaptor|IdealAdaptorStanfit|IdealAdaptorStaninput|",
      "CategoryRepresentation|CategoryRepresentationTemplate|Model)$"
    ),
    "",
    class_name
  )
  if (toupper(clean_name) == "EXEMPLAR") {
    return("EXEMPLAR")
  }
  toupper(clean_name)
}

S7::method(get_model_type, MVBU_Object) <- function(x, ...) {
  .get_family_from_class_name(class(x)[1])
}

S7::method(get_representation_type, MVBU_CategoryRepresentation) <- function(
  x,
  ...
) {
  .get_family_from_class_name(class(x)[1])
}

S7::method(
  get_representation_type,
  MVBU_CategoryRepresentationTemplate
) <- function(x, ...) {
  if (length(x@representations) == 0) {
    return("")
  }
  get_representation_type(x@representations[[1]], ...)
}

S7::method(get_representation_type, MVBU_CognitiveModel) <- function(x, ...) {
  get_representation_type(x@category_template, ...)
}

S7::method(get_parameters, MVBU_Object) <- function(x, ...) {
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
    mu = x@mu,
    sigma2 = x@sigma2,
    weights = x@weights
  )
}

S7::method(get_parameters, MNIX_CategoryRepresentation) <- function(x) {
  list(
    m = x@m,
    kappa = x@kappa,
    nu = x@nu,
    sigma2 = x@sigma2,
    weights = x@weights
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

S7::method(get_parameter_names, MVBU_CategoryRepresentation) <- function(
  x,
  original_pars = FALSE,
  ...
) {
  setdiff(names(S7::props(x)), c("metadata", "category_likelihood_function"))
}

S7::method(
  get_parameter_names,
  MVBU_CategoryRepresentationTemplate
) <- function(
  x,
  original_pars = FALSE,
  ...
) {
  unique(
    unlist(
      lapply(
        x@representations,
        get_parameter_names,
        original_pars = original_pars,
        ...
      )
    )
  )
}

S7::method(get_parameter_names, MVBU_CognitiveModel) <- function(
  x,
  original_pars = FALSE,
  ...
) {
  c(
    get_parameter_names(
      x@category_template,
      original_pars = original_pars,
      ...
    ),
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
  x@m
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
  diag(
    ifelse(x@nu > 2, x@nu * x@sigma2 / (x@nu - 2), NA_real_),
    nrow = length(x@nu)
  )
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
  categories = NULL,
  groups = NULL,
  statistic = c("mu", "Sigma"),
  ...
) {
  if (missing(statistic) || length(statistic) > 1L) {
    if (missing(statistic)) statistic <- c("mu", "Sigma")
    res <- list()
    for (s in statistic) {
      res[[s]] <- get_expected_category_statistic(
        x,
        categories = categories,
        groups = groups,
        statistic = s,
        ...
      )
    }
    return(res)
  }
  stat <- tolower(statistic)
  if (stat %in% c("mu", "m", "mean")) {
    get_expected_mu(x, ...)
  } else if (
    stat %in% c("sigma", "sigma2", "s", "cov", "covariance", "var", "variance")
  ) {
    get_expected_sigma(x, ...)
  } else {
    .stop(
      paste0("Unknown statistic '", statistic, "'. Expected 'mu' or 'sigma'.")
    )
  }
}

# -------------------------------------------------------------
# Marginal parameter extractors (predictive moments of tokens)
# -------------------------------------------------------------

S7::method(get_marginal_mu, NIX_CategoryRepresentation) <- function(x) {
  x@m
}

S7::method(get_marginal_mu, MNIX_CategoryRepresentation) <- function(x) {
  x@m
}

S7::method(get_marginal_mu, NIW_CategoryRepresentation) <- function(x) {
  x@m
}

S7::method(get_marginal_mu, MVBU_CategoryRepresentationTemplate) <- function(x) {
  lapply(x@representations, get_marginal_mu)
}

S7::method(get_marginal_mu, MVBU_CognitiveModel) <- function(x) {
  get_marginal_mu(x@category_template)
}

S7::method(get_marginal_sigma, NIX_CategoryRepresentation) <- function(x) {
  if (x@kappa > 0) ((x@kappa + 1) / x@kappa) * get_expected_sigma(x) else get_expected_sigma(x)
}

S7::method(get_marginal_sigma, MNIX_CategoryRepresentation) <- function(x) {
  s_exp <- get_expected_sigma(x)
  kap_factor <- ifelse(x@kappa > 0, (x@kappa + 1) / x@kappa, 1)
  diag(kap_factor * diag(s_exp), nrow = length(kap_factor))
}

S7::method(get_marginal_sigma, NIW_CategoryRepresentation) <- function(x) {
  if (x@kappa > 0) ((x@kappa + 1) / x@kappa) * get_expected_sigma(x) else get_expected_sigma(x)
}

S7::method(get_marginal_sigma, MVBU_CategoryRepresentationTemplate) <- function(x) {
  lapply(x@representations, get_marginal_sigma)
}

S7::method(get_marginal_sigma, MVBU_CognitiveModel) <- function(x) {
  get_marginal_sigma(x@category_template)
}

S7::method(get_marginal_category_statistic, MVBU_Object) <- function(
  x,
  categories = NULL,
  groups = NULL,
  statistic = c("mu", "Sigma"),
  ...
) {
  if (missing(statistic) || length(statistic) > 1L) {
    if (missing(statistic)) statistic <- c("mu", "Sigma")
    res <- list()
    for (s in statistic) {
      res[[s]] <- get_marginal_category_statistic(
        x,
        categories = categories,
        groups = groups,
        statistic = s,
        ...
      )
    }
    return(res)
  }
  stat <- tolower(statistic)
  if (stat %in% c("mu", "m", "mean")) {
    get_marginal_mu(x, ...)
  } else if (
    stat %in% c("sigma", "sigma2", "s", "cov", "covariance", "var", "variance")
  ) {
    get_marginal_sigma(x, ...)
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


S7::method(get_category_prior, list(MVBU_Object, S7::class_any)) <- function(
  x,
  categories,
  ...
) {
  .mvbu_not_implemented("get_category_prior", class(x)[1])
}

#' @name get_category_prior
#' @title Legacy compatibility method for category priors
#' @description Deprecated. Legacy compatibility method for list/data-frame inputs; remove once
#'   S7-only representations are required.
#' @keywords internal
S7::method(get_category_prior, list(S7::class_any, S7::class_any)) <- function(
  x,
  categories,
  ...
) {
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

S7::method(
  get_category_prior,
  list(MVBU_CognitiveModel, S7::class_any)
) <- function(
  x,
  categories,
  ...
) {
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

S7::method(get_lapse_rate, MVBU_Object) <- function(x, ...) {
  .mvbu_not_implemented("get_lapse_rate", class(x)[1])
}

#' @name get_lapse_rate
#' @title Legacy compatibility method for lapse rates
#' @description Deprecated. Legacy compatibility method for list/data-frame inputs; remove once
#'   S7-only representations are required.
#' @keywords internal
S7::method(get_lapse_rate, S7::class_any) <- function(x, ...) {
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

S7::method(get_lapse_rate, MVBU_CognitiveModel) <- function(x, ...) {
  x@lapse_behavior$lapse_rate
}

S7::method(get_lapse_bias, list(MVBU_Object, S7::class_any)) <- function(
  x,
  categories,
  ...
) {
  .mvbu_not_implemented("get_lapse_bias", class(x)[1])
}

#' @name get_lapse_bias
#' @title Legacy compatibility method for lapse biases
#' @description Deprecated. Legacy compatibility method for list/data-frame inputs; remove once
#'   S7-only representations are required.
#' @keywords internal
S7::method(get_lapse_bias, list(S7::class_any, S7::class_any)) <- function(
  x,
  categories,
  ...
) {
  if (is.list(x) && !is.null(x[["lapse_bias"]])) {
    lapse_bias <- x[["lapse_bias"]]
  } else if (is.data.frame(x) && "lapse_bias" %in% names(x)) {
    lapse_bias <- x[["lapse_bias"]]
  } else {
    return(NULL)
  }

  .mvbu_resolve_values_for_categories(lapse_bias, categories)
}

S7::method(
  get_lapse_bias,
  list(MVBU_CognitiveModel, S7::class_any)
) <- function(
  x,
  categories,
  ...
) {
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

S7::method(get_noise, MVBU_Object) <- function(x, ...) {
  .mvbu_not_implemented("get_noise", class(x)[1])
}

#' @title Get perceptual noise from a cognitive model
#' @description Extract perceptual noise covariance from a cognitive model.
#' @name get_noise
#' @keywords internal
S7::method(get_noise, MVBU_CognitiveModel) <- function(x, ...) {
  x@noise_behavior$Sigma_noise
}

S7::method(get_noise, S7::class_any) <- function(x, ...) {
  if (is.list(x) && !is.null(x[["Sigma_noise"]])) {
    return(x[["Sigma_noise"]])
  }
  NULL
}

S7::method(get_lapse_treatment, MVBU_Object) <- function(x, ...) {
  .mvbu_not_implemented("get_lapse_treatment", class(x)[1])
}

#' @name get_lapse_treatment
#' @title Get lapse treatment from a cognitive model
#' @description Extract lapse treatment from a cognitive model.
#' @keywords internal
S7::method(get_lapse_treatment, MVBU_CognitiveModel) <- function(x, ...) {
  x@lapse_behavior$lapse_treatment %||% "no_lapses"
}

S7::method(get_lapse_treatment, S7::class_any) <- function(x, ...) {
  if (is.list(x) && !is.null(x[["lapse_treatment"]])) {
    return(as.character(x[["lapse_treatment"]][1]))
  }
  if (is.data.frame(x) && "lapse_treatment" %in% names(x)) {
    return(as.character(x[["lapse_treatment"]][1]))
  }
  "no_lapses"
}

S7::method(get_noise_treatment, MVBU_Object) <- function(x, ...) {
  .mvbu_not_implemented("get_noise_treatment", class(x)[1])
}

#' @name get_noise_treatment
#' @title Get perceptual noise treatment from a cognitive model
#' @description Extract perceptual noise treatment from a cognitive model.
#' @keywords internal
S7::method(get_noise_treatment, MVBU_CognitiveModel) <- function(x, ...) {
  x@noise_behavior$noise_treatment %||% "no_noise"
}

S7::method(get_noise_treatment, S7::class_any) <- function(x, ...) {
  if (is.list(x) && !is.null(x[["noise_treatment"]])) {
    return(as.character(x[["noise_treatment"]][1]))
  }
  if (is.data.frame(x) && "noise_treatment" %in% names(x)) {
    return(as.character(x[["noise_treatment"]][1]))
  }
  "no_noise"
}

# -------------------------
# -------------------------
# Label accessors and setters
# -------------------------

#' @rdname set_labels
#' @export
S7::method(set_labels, S7::class_list) <- function(
  x,
  cue = character(0),
  category = character(0),
  response_category = category,
  group = character(0),
  ...
) {
  x$label_information <- list(
    cue = as.character(cue),
    category = as.character(category),
    response_category = as.character(response_category),
    group = as.character(group)
  )
  x
}

#' @rdname set_labels
#' @export
S7::method(set_labels, MVBU_Object) <- function(
  x,
  cue = get_cue_labels(x),
  category = get_category_labels(x),
  response_category = get_response_category_labels(x),
  group = get_group_labels(x, include_prior = FALSE),
  ...
) {
  x@metadata <- set_labels(
    x@metadata,
    cue = cue,
    category = category,
    response_category = response_category,
    group = group
  )
  x
}

#' @rdname get_labels
#' @export
S7::method(get_labels, MVBU_Object) <- function(x, ...) {
  lbls <- x@metadata$label_information
  if (is.null(lbls) || !is.list(lbls)) {
    .stop("No label information found in object metadata.")
  }
  lbls
}

#' @rdname get_cue_labels
#' @export
S7::method(get_cue_labels, MVBU_Object) <- function(x, indices = NULL, ...) {
  cue_labels <- get_labels(x)$cue
  if (missing(indices) || is.null(indices)) {
    return(cue_labels)
  }
  cue_labels[indices]
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

#' @rdname get_category_labels
#' @export
S7::method(get_category_labels, MVBU_Object) <- function(x, indices = NULL, ...) {
  category_labels <- get_labels(x)$category
  if (missing(indices) || is.null(indices)) {
    return(category_labels)
  }
  category_labels[indices]
}

#' @rdname get_response_category_labels
#' @export
S7::method(get_response_category_labels, MVBU_Object) <- function(
  x,
  indices = NULL,
  ...
) {
  resp_labels <- get_labels(x)$response_category
  if (missing(indices) || is.null(indices)) {
    return(resp_labels)
  }
  resp_labels[indices]
}

#' @rdname get_group_labels
#' @export
S7::method(get_group_labels, MVBU_Object) <- function(
  x,
  indices = NULL,
  include_prior = FALSE,
  ...
) {
  group_labels <- get_labels(x)$group
  if (include_prior) {
    group_labels <- append("prior", group_labels)
  }
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

S7::method(
  get_category_likelihood_function,
  MVBU_CategoryRepresentationTemplate
) <- function(x, ...) {
  representations <- get_category_representations(x)
  function(
    new_data,
    log = FALSE,
    noise_treatment = "no_noise",
    Sigma_noise = NULL,
    categories = NULL
  ) {
    .mvbu_likelihood_matrix(
      x,
      new_data,
      categories = categories,
      log = log,
      noise_treatment = noise_treatment,
      Sigma_noise = Sigma_noise
    )
  }
}

S7::method(
  likelihood,
  list(MVBU_CategoryRepresentation, S7::class_any, S7::class_any)
) <- function(x, new_data, categories, ...) {
  lik_fn <- x@category_likelihood_function
  d <- length(get_cue_labels(x))
  new_data <- .as_observation_matrix(new_data, d = d, arg_name = "new_data")
  as.numeric(
    lik_fn(
      new_data,
      log = FALSE,
      noise_treatment = "no_noise",
      Sigma_noise = NULL
    )
  )
}

S7::method(
  likelihood,
  list(MVBU_CategoryRepresentationTemplate, S7::class_any, S7::class_any)
) <- function(x, new_data, categories, ...) {
  if (missing(categories)) categories <- NULL
  .mvbu_likelihood_matrix(x, new_data, categories = categories, log = FALSE)
}

S7::method(
  likelihood,
  list(MVBU_CognitiveModel, S7::class_list, S7::class_any)
) <- function(x, new_data, categories, ...) {
  if (missing(categories)) categories <- NULL
  lapply(new_data, function(batch) likelihood(x, batch, categories, ...))
}

S7::method(
  likelihood,
  list(MVBU_CognitiveModel, S7::class_any, S7::class_any)
) <- function(x, new_data, categories, ...) {
  if (missing(categories)) categories <- NULL
  lik_fn <- S7::method(get_category_likelihood_function, MVBU_CognitiveModel)(x)
  lik_fn(new_data, categories = categories)
}

S7::method(get_category_likelihood_function, MVBU_CognitiveModel) <- function(
  x,
  noise_treatment,
  ...
) {
  if (missing(noise_treatment) || is.null(noise_treatment)) {
    noise_treatment <- get_noise_treatment(x)
  } else {
    noise_treatment <- as.character(noise_treatment)
  }
  key <- noise_treatment

  c_list <- .get_cache(x)
  lik_funcs <- c_list$category_likelihood_functions
  if (is.null(lik_funcs[[key]])) {
    if (is.null(lik_funcs)) lik_funcs <- list()
    lik_funcs[[key]] <- function(new_data, categories = NULL) {
      .mvbu_likelihood_matrix(x, new_data, categories = categories, log = FALSE, noise_treatment = noise_treatment, Sigma_noise = x@noise_behavior$Sigma_noise)
    }
    c_list$category_likelihood_functions <- lik_funcs
    .set_cache(x, c_list)
  }

  lik_funcs[[key]]
}

S7::method(get_category_posterior_function, list(MVBU_CognitiveModel, S7::class_any, S7::class_any)) <- function(x, noise_treatment, lapse_treatment) {
  if (missing(noise_treatment) && missing(lapse_treatment)) {
    noise_treatment <- get_noise_treatment(x)
    lapse_treatment <- get_lapse_treatment(x)
  } else {
    noise_treatment <- if (missing(noise_treatment) || is.null(noise_treatment)) get_noise_treatment(x) else as.character(noise_treatment)
    lapse_treatment <- if (missing(lapse_treatment) || is.null(lapse_treatment)) get_lapse_treatment(x) else as.character(lapse_treatment)
  }
  key <- paste(noise_treatment, lapse_treatment, sep = "__")

  c_list <- .get_cache(x)
  post_funcs <- c_list$category_posterior_functions
  if (is.null(post_funcs[[key]])) {
    if (is.null(post_funcs)) post_funcs <- list()
    post_funcs[[key]] <- function(new_data, categories = NULL) {
      .mvbu_posterior_matrix(x, new_data, categories = categories, noise_treatment = noise_treatment, lapse_treatment = lapse_treatment)
    }
    c_list$category_posterior_functions <- post_funcs
    .set_cache(x, c_list)
  }

  post_funcs[[key]]
}

S7::method(
  posterior,
  list(MVBU_CognitiveModel, S7::class_list, S7::class_any)
) <- function(x, new_data, categories, ...) {
  if (missing(categories)) categories <- NULL
  lapply(new_data, function(batch) posterior(x, batch, categories, ...))
}

S7::method(
  posterior,
  list(MVBU_CognitiveModel, S7::class_any, S7::class_any)
) <- function(x, new_data, categories, ...) {
  if (missing(categories)) categories <- NULL
  pf <- S7::method(
    get_category_posterior_function,
    list(MVBU_CognitiveModel, S7::class_any, S7::class_any)
  )(x)
  pf(new_data, categories = categories)
}

S7::method(categorize, list(MVBU_CognitiveModel, S7::class_list, S7::class_any)) <- function(x, new_data, decision_rule, simplify = NULL, ...) {
  lapply(new_data, function(batch) categorize(x, batch, decision_rule = decision_rule, simplify = simplify, ...))
}

S7::method(categorize, list(MVBU_CognitiveModel, S7::class_any, S7::class_any)) <- function(x, new_data, decision_rule, simplify = NULL, ...) {
  if (missing(decision_rule) || is.null(decision_rule)) {
    decision_rule <- x@decision_rule
  }
  if (is.null(simplify)) {
    simplify <- !identical(decision_rule, "proportional")
  }

  posterior_matrix <- posterior(x, new_data, categories = NULL)
  categories <- colnames(posterior_matrix)
  n_obs <- nrow(posterior_matrix)

  if (identical(decision_rule, "sampling") || identical(decision_rule, "sample")) {
    chosen_idx <- vapply(seq_len(n_obs), function(i) {
      p_i <- posterior_matrix[i, ]
      if (all(is.na(p_i)) || sum(p_i, na.rm = TRUE) == 0) {
        sample.int(length(categories), 1L)
      } else {
        sample.int(length(categories), 1L, prob = p_i)
      }
    }, integer(1))
    resp_prob <- rep(1.0, n_obs)
  } else if (identical(decision_rule, "criterion")) {
    chosen_idx <- apply(posterior_matrix, 1, which.max)
    resp_prob <- rep(1.0, n_obs)
  } else if (identical(decision_rule, "proportional")) {
    chosen_idx <- apply(posterior_matrix, 1, which.max)
    resp_prob <- as.numeric(posterior_matrix[cbind(seq_len(n_obs), chosen_idx)])
  } else {
    .stop("Invalid decision_rule: ", decision_rule, ". Must be one of 'criterion', 'sample', or 'proportional'.")
  }

  chosen_cats <- as.character(categories[chosen_idx])

  if (isTRUE(simplify)) {
    chosen_cats
  } else {
    data.frame(
      response_category = chosen_cats,
      response_probability = resp_prob,
      stringsAsFactors = FALSE
    )
  }
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
    if ("response_category" %in% names(response_category)) {
      response_category <- response_category$response_category
    } else if (ncol(response_category) == 1L) {
      response_category <- response_category[[1]]
    } else {
      stop("data frame passed to response_category must contain a 'response_category' column.")
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

# -----------------------------------------------------------------------------
# sample_observations methods
# -----------------------------------------------------------------------------

S7::method(sample_observations, MVBU_CategoryRepresentation) <- function(
  x,
  n = 1L,
  with_replacement = TRUE,
  randomize_order = TRUE,
  ...
) {
  dots <- list(...)
  if ("Ns" %in% names(dots)) n <- dots$Ns
  if ("randomize.order" %in% names(dots)) randomize_order <- dots$randomize.order

  .assert_true(.is_scalar_count(n) && n >= 0L, msg = "n must be a non-negative whole number.")
  cues <- get_cue_labels(x)
  cat_lbl <- get_category_labels(x)

  if (n == 0L) {
    empty_df <- as.data.frame(matrix(numeric(0), nrow = 0, ncol = length(cues), dimnames = list(NULL, cues)))
    empty_df <- cbind(data.frame(category = factor(character(0), levels = cat_lbl)), empty_df)
    return(tibble::as_tibble(empty_df))
  }

  family <- .mvbu_family_of_representation(x)
  mat <- .mvbu_sample_from_representation(x, from = family, n = n, with_replacement = with_replacement)
  colnames(mat) <- cues

  df <- cbind(data.frame(category = factor(rep(cat_lbl, n), levels = cat_lbl)), as.data.frame(mat))
  if (isTRUE(randomize_order) && nrow(df) > 1L) {
    df <- df[sample.int(nrow(df)), , drop = FALSE]
  }
  tibble::as_tibble(df)
}

S7::method(sample_observations, MVBU_CategoryRepresentationTemplate) <- function(
  x,
  n = 1L,
  with_replacement = TRUE,
  randomize_order = TRUE,
  ...
) {
  dots <- list(...)
  if ("Ns" %in% names(dots)) n <- dots$Ns
  if ("randomize.order" %in% names(dots)) randomize_order <- dots$randomize.order

  reps <- x@representations
  K <- length(reps)
  cat_labels <- get_category_labels(x)
  cue_labels <- get_cue_labels(x)

  .assert_true(
    .is_scalar_count(n) || (is.numeric(n) && length(n) == K),
    msg = "n must be a non-negative whole number or a vector of counts with one element per category."
  )

  if (length(n) == 1L) {
    if (n == 0L) {
      empty_df <- as.data.frame(matrix(numeric(0), nrow = 0, ncol = length(cue_labels), dimnames = list(NULL, cue_labels)))
      empty_df <- cbind(data.frame(category = factor(character(0), levels = cat_labels)), empty_df)
      return(tibble::as_tibble(empty_df))
    }

    # Total available exemplars check if without replacement
    has_ex <- any(vapply(reps, function(r) S7::S7_inherits(r, Exemplar_CategoryRepresentation), logical(1)))
    if (!isTRUE(with_replacement) && has_ex) {
      avail_per_cat <- vapply(reps, function(r) {
        if (S7::S7_inherits(r, Exemplar_CategoryRepresentation)) nrow(r@exemplars) else Inf
      }, numeric(1))
      total_avail <- sum(avail_per_cat)
      if (n > total_avail) {
        .stop(sprintf(
          "Requested %d samples without replacement, but template only contains %d total exemplars across categories.",
          n, as.integer(total_avail)
        ))
      }
    }

    # Uniform sampling across categories (n total draws)
    cat_draws <- sample.int(K, size = n, replace = TRUE)
    n_per_cat <- as.vector(table(factor(cat_draws, levels = seq_len(K))))
  } else {
    n_per_cat <- as.integer(n)
  }

  rows <- lapply(seq_len(K), function(i) {
    n_i <- n_per_cat[i]
    if (n_i == 0L) {
      return(NULL)
    }
    rep_i <- reps[[i]]
    fam_i <- .mvbu_family_of_representation(rep_i)

    # Check without replacement condition for exemplar
    if (!isTRUE(with_replacement) && fam_i == "EXEMPLAR") {
      n_avail_i <- nrow(rep_i@exemplars)
      if (n_i > n_avail_i) {
        .stop(sprintf(
          "Random sampling without replacement allocated %d observations to category '%s', but it only contains %d exemplars. Because category assignment is a random process, re-evoking the function or increasing category exemplars may resolve this.",
          n_i, cat_labels[i], n_avail_i
        ))
      }
    }

    mat_i <- .mvbu_sample_from_representation(rep_i, from = fam_i, n = n_i, with_replacement = with_replacement)
    colnames(mat_i) <- cue_labels
    cbind(data.frame(category = cat_labels[i], stringsAsFactors = FALSE), as.data.frame(mat_i))
  })

  rows <- rows[!vapply(rows, is.null, logical(1))]
  out <- if (length(rows) == 0L) {
    empty_df <- as.data.frame(matrix(numeric(0), nrow = 0, ncol = length(cue_labels), dimnames = list(NULL, cue_labels)))
    cbind(data.frame(category = factor(character(0), levels = cat_labels)), empty_df)
  } else {
    do.call(rbind, rows)
  }
  out$category <- factor(out$category, levels = cat_labels)

  if (isTRUE(randomize_order) && nrow(out) > 1L) {
    out <- out[sample.int(nrow(out)), , drop = FALSE]
  }

  tibble::as_tibble(out)
}

S7::method(sample_observations, MVBU_CognitiveModel) <- function(
  x,
  n = 1L,
  with_replacement = TRUE,
  randomize_order = TRUE,
  ...
) {
  dots <- list(...)
  if ("Ns" %in% names(dots)) n <- dots$Ns
  if ("randomize.order" %in% names(dots)) randomize_order <- dots$randomize.order

  reps <- x@category_template@representations
  K <- length(reps)
  cat_labels <- get_category_labels(x)
  cue_labels <- get_cue_labels(x)

  .assert_true(
    .is_scalar_count(n) || (is.numeric(n) && length(n) == K),
    msg = "n must be a non-negative whole number or a vector of counts with one element per category."
  )

  if (length(n) == 1L) {
    if (n == 0L) {
      empty_df <- as.data.frame(matrix(numeric(0), nrow = 0, ncol = length(cue_labels), dimnames = list(NULL, cue_labels)))
      empty_df <- cbind(data.frame(category = factor(character(0), levels = cat_labels)), empty_df)
      return(tibble::as_tibble(empty_df))
    }

    # Total available exemplars check if without replacement
    has_ex <- any(vapply(reps, function(r) S7::S7_inherits(r, Exemplar_CategoryRepresentation), logical(1)))
    if (!isTRUE(with_replacement) && has_ex) {
      avail_per_cat <- vapply(reps, function(r) {
        if (S7::S7_inherits(r, Exemplar_CategoryRepresentation)) nrow(r@exemplars) else Inf
      }, numeric(1))
      total_avail <- sum(avail_per_cat)
      if (n > total_avail) {
        .stop(sprintf(
          "Requested %d samples without replacement, but model template only contains %d total exemplars across categories.",
          n, as.integer(total_avail)
        ))
      }
    }

    # Proportional sampling to category_prior
    cp <- x@category_prior
    if (is.null(cp) || length(cp) != K || any(is.na(cp)) || sum(cp) <= 0) {
      cp <- rep(1 / K, K)
    } else {
      cp <- cp / sum(cp)
    }

    cat_draws <- sample.int(K, size = n, replace = TRUE, prob = as.numeric(cp))
    n_per_cat <- as.vector(table(factor(cat_draws, levels = seq_len(K))))
  } else {
    n_per_cat <- as.integer(n)
  }

  rows <- lapply(seq_len(K), function(i) {
    n_i <- n_per_cat[i]
    if (n_i == 0L) {
      return(NULL)
    }
    rep_i <- reps[[i]]
    fam_i <- .mvbu_family_of_representation(rep_i)

    # Check without replacement condition for exemplar
    if (!isTRUE(with_replacement) && fam_i == "EXEMPLAR") {
      n_avail_i <- nrow(rep_i@exemplars)
      if (n_i > n_avail_i) {
        .stop(sprintf(
          "Random sampling without replacement allocated %d observations to category '%s', but it only contains %d exemplars. Because category assignment is proportional to category priors and involves a random process, re-evoking the function or increasing category exemplars may resolve this.",
          n_i, cat_labels[i], n_avail_i
        ))
      }
    }

    mat_i <- .mvbu_sample_from_representation(rep_i, from = fam_i, n = n_i, with_replacement = with_replacement)
    colnames(mat_i) <- cue_labels
    cbind(data.frame(category = cat_labels[i], stringsAsFactors = FALSE), as.data.frame(mat_i))
  })

  rows <- rows[!vapply(rows, is.null, logical(1))]
  out <- if (length(rows) == 0L) {
    empty_df <- as.data.frame(matrix(numeric(0), nrow = 0, ncol = length(cue_labels), dimnames = list(NULL, cue_labels)))
    cbind(data.frame(category = factor(character(0), levels = cat_labels)), empty_df)
  } else {
    do.call(rbind, rows)
  }
  out$category <- factor(out$category, levels = cat_labels)

  if (isTRUE(randomize_order) && nrow(out) > 1L) {
    out <- out[sample.int(nrow(out)), , drop = FALSE]
  }

  tibble::as_tibble(out)
}

#' @rdname get_sufficient_category_statistics
#' @export
S7::method(get_sufficient_category_statistics, MVBU_CognitiveModel) <- function(
  x,
  categories = NULL,
  groups = NULL,
  untransform_cues = FALSE,
  ...
) {
  reps <- x@category_template@representations
  cats <- get_category_labels(x)
  fam <- get_model_family(x)

  out <- list(category = factor(cats, levels = cats))
  out$x_N <- vapply(reps, function(r) {
    if ("n" %in% names(r)) {
      as.integer(r$n)
    } else if ("N" %in% names(r)) {
      as.integer(r$N)
    } else if ("x_N" %in% names(r)) {
      as.integer(r$x_N)
    } else {
      0L
    }
  }, integer(1L))

  out$x_mean <- lapply(reps, function(r) {
    mc <- .get_rep_mean_and_cov(r)
    mc$mu
  })

  out$x_ss <- lapply(reps, function(r) {
    mc <- .get_rep_mean_and_cov(r)
    if (fam %in% c("NIX", "UVG")) {
      mc$Sigma[1L, 1L]
    } else if (fam == "MNIX") {
      diag(mc$Sigma)
    } else {
      mc$Sigma
    }
  })
  out$x_css <- out$x_ss
  if (fam %in% c("NIW", "MVG", "MUVG")) {
    out$x_cov <- lapply(reps, function(r) {
      .get_rep_mean_and_cov(r)$Sigma
    })
  }
  res <- tibble::as_tibble(out)
  if (!is.null(categories)) {
    res <- res[res$category %in% categories, , drop = FALSE]
  }
  res
}

# -----------------------------------------------------------------------------
# add_category_representation methods
# -----------------------------------------------------------------------------

#' Update a probability vector when adding a category representation
#'
#' @param old_vec Numeric named vector of existing probabilities (priors or
#'   biases) summing to 1.
#' @param new_input Numeric scalar, vector, or `NULL`. If `NULL`, defaults to 0
#'   for the new category. If scalar, sets the value for the new category and
#'   rescales existing categories by `old_vec * (1 - new_input)`. If vector of
#'   length matching the total new categories, overwrites all probabilities.
#' @param new_name Character scalar naming the new category.
#' @param all_names Character vector of all category names in the new template.
#' @param param_name Character scalar ("category_prior" or "lapse_bias") for
#'   error messages.
#'
#' @return A numeric vector of length `length(all_names)` summing to 1 with
#'   names `all_names`.
#' @noRd
#' @keywords internal
.mvbu_update_probability_vector <- function(
  old_vec,
  new_input,
  new_name,
  all_names,
  param_name
) {
  k_new <- length(all_names)
  k_old <- length(old_vec)

  if (is.null(new_input)) {
    new_input <- 0
  }

  if (!is.numeric(new_input)) {
    .stop(param_name, " must be numeric.")
  }

  if (length(new_input) == 1L) {
    new_val <- as.numeric(new_input)
    if (is.na(new_val) || new_val < 0 || new_val > 1) {
      .stop(param_name, " for the new category must be between 0 and 1.")
    }

    if (!is.null(names(new_input)) && nzchar(names(new_input))) {
      input_name <- names(new_input)
      if (input_name %in% names(old_vec) && input_name != new_name) {
        .stop(
          param_name,
          " name matches an existing category, not the new category."
        )
      }
    }

    if (k_old == 0L) {
      out <- stats::setNames(1, new_name)
      return(out)
    }

    rescaled_old <- old_vec * (1 - new_val)
    out <- c(rescaled_old, stats::setNames(new_val, new_name))
    names(out) <- all_names
    return(out)
  }

  if (length(new_input) == k_new) {
    if (any(is.na(new_input)) || any(new_input < 0) || any(new_input > 1)) {
      .stop(param_name, " entries must be in [0, 1].")
    }
    if (abs(sum(new_input) - 1) > MVBU_PROB_TOL) {
      .stop(param_name, " entries must sum to 1.")
    }

    input_names <- names(new_input)
    if (!is.null(input_names) && all(nzchar(input_names))) {
      if (!setequal(input_names, all_names)) {
        .stop(
          param_name,
          " names must match all category names: ",
          paste(all_names, collapse = ", ")
        )
      }
      out <- as.numeric(new_input[all_names])
    } else {
      out <- as.numeric(new_input)
    }
    names(out) <- all_names
    return(out)
  }

  .stop(
    param_name,
    " must be a scalar (for the new category) or a numeric vector of length ",
    k_new, " (for all categories)."
  )
}

#' @rdname add_category_representation
#' @export
S7::method(add_category_representation, MVBU_CategoryRepresentationTemplate) <- function(
  x,
  representation,
  name = NULL,
  category_prior = NULL,
  lapse_bias = NULL,
  ...
) {
  if (!S7::S7_inherits(representation, MVBU_CategoryRepresentation)) {
    .stop("representation must be an MVBU_CategoryRepresentation.")
  }
  if (is.null(name)) {
    rep_labels_cat <- get_category_labels(representation)
    if (length(rep_labels_cat) == 1L && is.character(rep_labels_cat) && nzchar(rep_labels_cat)) {
      name <- rep_labels_cat
    }
  }

  representations <- x@representations
  representations[[length(representations) + 1L]] <- representation

  if (!is.null(name)) {
    names(representations)[length(representations)] <- name
  }

  .mvbu_validate_cue_consistency(representations)

  metadata <- as.list(x@metadata)
  rep_cues <- get_cue_labels(representation)
  template_cues <- get_cue_labels(x)
  if (length(template_cues) > 0 &&
    length(rep_cues) > 0 &&
    !identical(as.character(template_cues), as.character(rep_cues))) {
    .stop(
      paste0(
        "cue labels must be consistent across all representations ",
        "in a template."
      )
    )
  }

  metadata <- set_labels(
    metadata,
    cue = if (length(template_cues) > 0) template_cues else rep_cues,
    category = c(get_category_labels(x), get_category_labels(representation)),
    response_category = c(get_category_labels(x), get_category_labels(representation)),
    group = get_group_labels(x, include_prior = FALSE)
  )

  MVBU_CategoryRepresentationTemplate(
    representations = representations,
    metadata = metadata
  )
}

#' @rdname add_category_representation
#' @export
S7::method(add_category_representation, MVBU_CognitiveModel) <- function(
  x,
  representation,
  name = NULL,
  category_prior = NULL,
  lapse_bias = NULL,
  ...
) {
  new_template <- add_category_representation(
    x@category_template,
    representation,
    name = name,
    ...
  )
  new_k <- length(new_template@representations)

  all_names <- names(new_template@representations)
  if (is.null(all_names) || any(all_names == "")) {
    all_names <- get_category_labels(new_template)
  }
  if (is.null(all_names) || length(all_names) != new_k) {
    all_names <- paste0("cat_", seq_len(new_k))
  }
  new_name <- all_names[new_k]
  old_names <- all_names[-new_k]

  old_prior <- x@category_prior
  if (is.null(names(old_prior)) && length(old_prior) == length(old_names)) {
    names(old_prior) <- old_names
  }

  old_lapse_bias <- x@lapse_behavior$lapse_bias
  if (is.null(names(old_lapse_bias)) &&
    length(old_lapse_bias) == length(old_names)) {
    names(old_lapse_bias) <- old_names
  }

  updated_category_prior <- .mvbu_update_probability_vector(
    old_vec = old_prior,
    new_input = category_prior,
    new_name = new_name,
    all_names = all_names,
    param_name = "category_prior"
  )

  updated_lapse_bias <- .mvbu_update_probability_vector(
    old_vec = old_lapse_bias,
    new_input = lapse_bias,
    new_name = new_name,
    all_names = all_names,
    param_name = "lapse_bias"
  )

  cls <- S7::S7_class(x)
  base_model <- new_cognitive_model(
    category_template = new_template,
    decision_rule = x@decision_rule,
    category_prior = updated_category_prior,
    lapse_rate = x@lapse_behavior$lapse_rate,
    lapse_bias = updated_lapse_bias,
    Sigma_noise = x@noise_behavior$Sigma_noise,
    noise_treatment = x@noise_behavior$noise_treatment,
    lapse_treatment = x@lapse_behavior$lapse_treatment,
    metadata = x@metadata
  )

  if (identical(cls, MVBU_CognitiveModel)) {
    return(base_model)
  }

  typed_model <- tryCatch(
    {
      res <- cls(
        category_template = base_model@category_template,
        decision_rule = base_model@decision_rule,
        category_prior = base_model@category_prior,
        lapse_behavior = base_model@lapse_behavior,
        noise_behavior = base_model@noise_behavior,
        metadata = base_model@metadata
      )
      res@category_posterior_functions <- base_model@category_posterior_functions
      res
    },
    error = function(e) base_model
  )

  typed_model
}
