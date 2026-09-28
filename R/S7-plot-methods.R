# =============================================================================
# S7 Plot Methods for MVBeliefUpdatr
# =============================================================================

#' @include S7-generics.R
#' @include S7-class.R
#' @include S7-methods.R
#' @include S7-stanfit.R
#' @include S7-stanfit-methods.R
#' @include S7-plot-engine.R
NULL

# -----------------------------------------------------------------------------
# Internal Plot Helpers
# -----------------------------------------------------------------------------

#' Consistent base theme for all MVBU plots
#' @noRd
#' @keywords internal
.mvbu_theme <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(
        fill = "white",
        color = "gray85"
      ),
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),
      strip.background = ggplot2::element_rect(
        fill = "white",
        color = "gray85"
      ),
      strip.text = ggplot2::element_text(face = "bold", color = "black"),
      panel.grid.major = ggplot2::element_line(color = "gray95"),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 12, margin = ggplot2::margin(b = 14)),
      plot.subtitle = ggplot2::element_text(size = 9.5, margin = ggplot2::margin(b = 8), lineheight = 1.15)
    )
}

#' Format model family name using clean, cross-platform mathematical notation
#' @noRd
#' @keywords internal
.format_model_family_name <- function(x) {
  D <- length(get_cue_labels(x))
  if (S7::S7_inherits(x, UVG_IdealObserver)) {
    return("1D univariate Gaussian ideal observer")
  }
  if (S7::S7_inherits(x, MVG_IdealObserver)) {
    return(sprintf("%dD multivariate Gaussian ideal observer", D))
  }
  if (S7::S7_inherits(x, MUVG_IdealObserver)) {
    return(sprintf("%dD multivariate independent Gaussian ideal observer", D))
  }
  if (S7::S7_inherits(x, NIX_IdealAdaptor)) {
    return("1D N-Inv-\u03c7\u00b2 (NIX) ideal adaptor")
  }
  if (S7::S7_inherits(x, MNIX_IdealAdaptor)) {
    return(sprintf("%dD indep. N-Inv-\u03c7\u00b2 (MNIX) ideal adaptor", D))
  }
  if (S7::S7_inherits(x, NIW_IdealAdaptor)) {
    return(sprintf("%dD N-Inv-Wishart (NIW) ideal adaptor", D))
  }
  if (S7::S7_inherits(x, Exemplar_Model)) {
    return(sprintf("%dD exemplar model", D))
  }
  if (S7::S7_inherits(x, UVG_CategoryRepresentation)) {
    return("1D univariate Gaussian category")
  }
  if (S7::S7_inherits(x, MVG_CategoryRepresentation)) {
    return(sprintf("%dD multivariate Gaussian category", D))
  }
  if (S7::S7_inherits(x, NIX_CategoryRepresentation)) {
    return("1D N-Inv-\u03c7\u00b2 (NIX) category")
  }
  if (S7::S7_inherits(x, MNIX_CategoryRepresentation)) {
    return(sprintf("%dD indep. N-Inv-\u03c7\u00b2 (MNIX) category", D))
  }
  if (S7::S7_inherits(x, NIW_CategoryRepresentation)) {
    return(sprintf("%dD N-Inv-Wishart (NIW) category", D))
  }
  if (S7::S7_inherits(x, Exemplar_CategoryRepresentation)) {
    return(sprintf("%dD exemplar category", D))
  }
  if (S7::S7_inherits(x, MVBU_CategoryRepresentationTemplate)) {
    first_rep <- if (length(x@representations) > 0L) x@representations[[1L]] else NULL
    if (!is.null(first_rep)) {
      fam_rep <- .format_model_family_name(first_rep)
      fam_rep <- sub(" category$", " category template", fam_rep)
      if (!grepl("template$", fam_rep)) fam_rep <- paste0(fam_rep, " template")
      return(fam_rep)
    }
    return(sprintf("%dD category template", D))
  }
  if (S7::S7_inherits(x, MVBU_Stanfit)) {
    mtype <- gsub("_", " ", get_model_type(x))
    return(sprintf("%dD %s (Stanfit)", D, mtype))
  }
  "MVBU object"
}

#' Construct plot title and informative subtitle
#' @noRd
#' @keywords internal
.make_plot_title_and_subtitle <- function(
  x,
  plot_type = "categories",
  cues = NULL,
  ndraws = NULL,
  target_category = NULL,
  decision_rule = NULL
) {
  titles <- c(
    categories = "Category likelihood",
    categorization = "Categorization function",
    parameters = "Model parameters",
    parameter_correlations = "Parameter correlations",
    parameters_pairwise = "Pairwise parameter distributions",
    cues = "Cue distributions",
    sample = "Sample distributions",
    diagnostics = "MCMC diagnostics"
  )
  main_title <- if (plot_type %in% names(titles)) {
    titles[[plot_type]]
  } else {
    "MVBU plot"
  }

  if (plot_type %in% c("categories", "categorization")) {
    all_cues <- get_cue_labels(x)
    if (!is.null(cues) && length(cues) < length(all_cues)) {
      unplotted <- setdiff(all_cues, cues)
      if (length(unplotted) > 0L) {
        main_title <- sprintf(
          "%s (marginalizing over %s)",
          main_title,
          paste(unplotted, collapse = ", ")
        )
      }
    }
  }

  fam <- .format_model_family_name(x)

  if (S7::S7_inherits(x, MVBU_CategoryRepresentation) ||
    S7::S7_inherits(x, MVBU_CategoryRepresentationTemplate)) {
    return(list(title = main_title, subtitle = fam))
  }

  if (S7::S7_inherits(x, MVBU_CognitiveModel)) {
    if (identical(plot_type, "parameters")) {
      return(list(title = main_title, subtitle = fam))
    }

    lapse_beh <- tryCatch(x@lapse_behavior, error = function(e) list())
    lapse_trt <- lapse_beh$lapse_treatment %||% "none"
    if (identical(lapse_trt, "no_lapse") || identical(lapse_trt, "no_lapses")) {
      lapse_trt <- "none"
    }
    lapse_trt_txt <- sprintf("Lapses: %s", lapse_trt)

    noise_beh <- tryCatch(x@noise_behavior, error = function(e) list())
    sigma_n <- noise_beh$Sigma_noise
    noise_trt <- noise_beh$noise_treatment %||% "no_noise"

    noise_txt <- if (!is.null(sigma_n) &&
      !identical(noise_trt, "no_noise") &&
      !identical(noise_trt, "none")) {
      sd_n <- round(sqrt(diag(as.matrix(sigma_n))), 2)
      sprintf(
        "Noise: %s (\u03c3_noise = %s)",
        noise_trt,
        paste(sd_n, collapse = ", ")
      )
    } else {
      "Noise: none"
    }

    if (identical(plot_type, "categories")) {
      line2 <- sprintf("%s | %s", lapse_trt_txt, noise_txt)
      sub_txt <- sprintf("%s\n%s", fam, line2)
      return(list(title = main_title, subtitle = sub_txt))
    }

    if (identical(plot_type, "categorization")) {
      drule <- if (!is.null(decision_rule)) decision_rule else x@decision_rule
      line2 <- sprintf(
        "Decision rule: %s | %s | %s",
        drule,
        lapse_trt_txt,
        noise_txt
      )

      cp <- x@category_prior
      lb <- lapse_beh$lapse_bias
      l_rate <- lapse_beh$lapse_rate %||% 0

      .format_prior_vec <- function(vec) {
        if (is.null(vec) || length(vec) == 0L) {
          return("")
        }
        if (length(vec) > 1L && length(unique(round(as.numeric(vec), 4))) == 1L) {
          return(sprintf("all = %s", round(vec[1L], 2)))
        }
        paste(names(vec), round(vec, 2), sep = "=", collapse = ", ")
      }

      priors_equal <- FALSE
      if (!is.null(cp) && !is.null(lb) && length(cp) == length(lb)) {
        if (all(names(cp) == names(lb)) &&
          isTRUE(all.equal(as.numeric(cp), as.numeric(lb)))) {
          priors_equal <- TRUE
        }
      }

      if (isTRUE(priors_equal)) {
        pri_bias_txt <- sprintf(
          "Category prior & lapse bias: %s",
          .format_prior_vec(cp)
        )
      } else {
        pri_txt <- sprintf(
          "Category prior: %s",
          .format_prior_vec(cp)
        )
        bias_txt <- if (!is.null(lb)) {
          sprintf(
            "Lapse bias: %s",
            .format_prior_vec(lb)
          )
        } else {
          ""
        }
        pri_bias_txt <- paste(
          c(pri_txt, bias_txt)[nzchar(c(pri_txt, bias_txt))],
          collapse = " | "
        )
      }

      lapse_txt <- sprintf("Lapse rate (\u03bb): %s", round(l_rate, 2))
      line3 <- sprintf("%s | %s", pri_bias_txt, lapse_txt)
      sub_txt <- sprintf("%s\n%s\n%s", fam, line2, line3)
      return(list(title = main_title, subtitle = sub_txt))
    }

    return(list(title = main_title, subtitle = fam))
  }

  if (S7::S7_inherits(x, MVBU_Stanfit)) {
    draws_txt <- if (!is.null(ndraws)) {
      sprintf("ndraws = %d", ndraws)
    } else {
      "all MCMC draws"
    }
    lapse_noise_txt <- "Lapses: marginalize | Noise: none (incl. in S)"
    if (identical(plot_type, "categorization")) {
      line2 <- sprintf(
        "Decision rule: optimal | %s | %s",
        lapse_noise_txt,
        draws_txt
      )
    } else {
      line2 <- sprintf("%s | %s", lapse_noise_txt, draws_txt)
    }
    sub_txt <- sprintf("%s\n%s", fam, line2)
    return(list(title = main_title, subtitle = sub_txt))
  }

  list(title = main_title, subtitle = "")
}

#' Parse limits argument into named list per cue
#' @noRd
#' @keywords internal
.parse_cue_limits <- function(limits, cues) {
  if (is.null(limits)) {
    return(NULL)
  }
  if (is.numeric(limits) && length(limits) == 2L && length(cues) == 1L) {
    res <- list()
    res[[cues[1L]]] <- limits
    return(res)
  }
  if (is.list(limits)) {
    res <- list()
    for (i in seq_along(cues)) {
      c_name <- cues[i]
      if (c_name %in% names(limits)) {
        res[[c_name]] <- limits[[c_name]]
      } else if (i <= length(limits) && is.null(names(limits))) {
        res[[c_name]] <- limits[[i]]
      }
    }
    if (length(res) > 0L) {
      return(res)
    }
  }
  NULL
}

#' @importFrom geomtextpath geom_textcontour
NULL

#' Parse levels argument into contour and fill levels for plot_categories
#' @noRd
#' @keywords internal
.parse_category_levels <- function(levels, aes) {
  # Default levels correspond to two-tailed 1, 2, 3 sigma
  def_levels <- 2 * stats::pnorm(1:3) - 1
  if (is.null(levels)) {
    return(list(contour = def_levels, fill = def_levels))
  }
  if (is.numeric(levels)) {
    return(list(contour = levels, fill = levels))
  }
  if (is.list(levels)) {
    res <- list(contour = def_levels, fill = def_levels)
    if (!is.null(names(levels))) {
      if ("contour" %in% names(levels)) {
        res$contour <- levels[["contour"]]
      }
      if ("fill" %in% names(levels)) {
        res$fill <- levels[["fill"]]
      }
      if ("fill-discrete" %in% names(levels)) {
        res$fill <- levels[["fill-discrete"]]
      }
      return(res)
    }
    if (length(levels) >= 1L && length(aes) >= 1L) {
      res[[aes[1L]]] <- levels[[1L]]
    }
    if (length(levels) >= 2L && length(aes) >= 2L) {
      res[[aes[2L]]] <- levels[[2L]]
    }
    return(res)
  }
  list(contour = def_levels, fill = def_levels)
}

#' Parse levels argument into contour and fill levels for categorization
#' @noRd
#' @keywords internal
.parse_categorization_levels <- function(levels, aes) {
  def_levels <- c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99)
  if (is.null(levels)) {
    return(list(contour = def_levels, fill = def_levels))
  }
  if (is.numeric(levels)) {
    return(list(contour = levels, fill = levels))
  }
  if (is.list(levels)) {
    res <- list(contour = def_levels, fill = def_levels)
    if (!is.null(names(levels))) {
      if ("contour" %in% names(levels)) {
        res$contour <- levels[["contour"]]
      }
      if ("fill" %in% names(levels)) {
        res$fill <- levels[["fill"]]
      }
      if ("fill-discrete" %in% names(levels)) {
        res$fill <- levels[["fill-discrete"]]
      }
      return(res)
    }
    if (length(levels) >= 1L && length(aes) >= 1L) {
      res[[aes[1L]]] <- levels[[1L]]
    }
    if (length(levels) >= 2L && length(aes) >= 2L) {
      res[[aes[2L]]] <- levels[[2L]]
    }
    return(res)
  }
  list(contour = def_levels, fill = def_levels)
}

#' Normalize aes argument
#' @noRd
#' @keywords internal
.normalize_plot_aes <- function(aes, default_aes = "contour") {
  if (is.null(aes)) {
    return(default_aes)
  }
  aes <- as.character(aes)
  res <- character()
  for (a in aes) {
    if (a %in% c("both", "fill_gradient_contour", "fill-gradient-contour")) {
      res <- c(res, "fill-gradient", "contour")
    } else if (a %in% c("fill_discrete_contour", "fill-discrete-contour")) {
      res <- c(res, "fill-discrete", "contour")
    } else if (a %in% c("fill", "fill_gradient", "fill-gradient")) {
      res <- c(res, "fill-gradient")
    } else if (a %in% c("fill_discrete", "fill-discrete")) {
      res <- c(res, "fill-discrete")
    } else {
      res <- c(res, a)
    }
  }
  unique(res)
}

#' Parse levels argument into contour and fill levels (legacy compatibility)
#' @noRd
#' @keywords internal
.parse_plot_levels <- function(levels, aes) {
  .parse_category_levels(levels, aes)
}

#' Format raw parameter name to plotmath expression string
#' @noRd
#' @keywords internal
.format_param_plotmath <- function(x) {
  sapply(x, function(val) {
    if (grepl("__", val)) {
      parts_g <- strsplit(val, "__")[[1L]]
      param_part <- parts_g[1L]
      grp_raw <- parts_g[2L]
      grp_sub <- if (identical(grp_raw, "prior")) {
        "0"
      } else {
        sprintf("'%s'", grp_raw)
      }

      if (grepl("^kappa", param_part)) {
        cat_part <- sub("^kappa_?", "", param_part)
        sprintf("kappa[%s*','*%s]", cat_part, grp_sub)
      } else if (grepl("^nu", param_part)) {
        cat_part <- sub("^nu_?", "", param_part)
        sprintf("nu[%s*','*%s]", cat_part, grp_sub)
      } else if (grepl("^m_", param_part)) {
        rest <- sub("^m_", "", param_part)
        parts <- strsplit(rest, "_")[[1L]]
        sprintf(
          "mu[%s*','*%s*','*%s]",
          parts[1L],
          paste(parts[-1L], collapse = "_"),
          grp_sub
        )
      } else if (grepl("^tau_", param_part)) {
        rest <- sub("^tau_", "", param_part)
        parts <- strsplit(rest, "_")[[1L]]
        sprintf(
          "tau[%s*','*%s*','*%s]",
          parts[1L],
          paste(parts[-1L], collapse = "_"),
          grp_sub
        )
      } else if (grepl("^rho_", param_part)) {
        rest <- sub("^rho_", "", param_part)
        parts <- strsplit(rest, "_")[[1L]]
        sprintf("rho[%s*','*%s]", paste(parts, collapse = "*','*"), grp_sub)
      } else if (grepl("^lapse_rate", param_part)) {
        sprintf("lambda[%s]", grp_sub)
      } else {
        sprintf("'%s'", val)
      }
    } else {
      if (grepl("^kappa", val)) {
        cat_part <- sub("^kappa_?", "", val)
        if (nzchar(cat_part)) sprintf("kappa[%s]", cat_part) else "kappa"
      } else if (grepl("^nu", val)) {
        cat_part <- sub("^nu_?", "", val)
        if (nzchar(cat_part)) sprintf("nu[%s]", cat_part) else "nu"
      } else if (grepl("^m_", val)) {
        rest <- sub("^m_", "", val)
        parts <- strsplit(rest, "_")[[1L]]
        if (length(parts) >= 2L) {
          sprintf("mu[%s*','*%s]", parts[1L], paste(parts[-1L], collapse = "_"))
        } else {
          sprintf("mu[%s]", rest)
        }
      } else if (grepl("^tau_", val)) {
        rest <- sub("^tau_", "", val)
        parts <- strsplit(rest, "_")[[1L]]
        if (length(parts) >= 2L) {
          sprintf(
            "tau[%s*','*%s]",
            parts[1L],
            paste(parts[-1L], collapse = "_")
          )
        } else {
          sprintf("tau[%s]", rest)
        }
      } else if (grepl("^rho_", val)) {
        rest <- sub("^rho_", "", val)
        sprintf("rho[%s]", rest)
      } else if (grepl("^lapse_rate", val)) {
        "lambda"
      } else {
        sprintf("'%s'", val)
      }
    }
  }, USE.NAMES = FALSE)
}

#' Marginalize a category representation down to a subset of cues
#' @noRd
#' @keywords internal
.marginalize_representation_to_cues <- function(r, cues) {
  rep_cues <- get_cue_labels(r)
  cat_labels <- get_category_labels(r)

  if (length(cues) == length(rep_cues) && all(cues == rep_cues)) {
    return(r)
  }

  .validate_requested_labels(cues, rep_cues, label_type = "cue")

  if (length(cues) == 1L) {
    idx <- match(cues[1], rep_cues)

    if (S7::S7_inherits(r, UVG_CategoryRepresentation)) {
      return(r)
    }
    if (S7::S7_inherits(r, MVG_CategoryRepresentation)) {
      return(new_uvg_category_representation(
        category_labels = cat_labels,
        cue_labels = cues,
        mu = r@mu[idx],
        sigma2 = r@Sigma[idx, idx]
      ))
    }
    if (S7::S7_inherits(r, NIX_CategoryRepresentation)) {
      return(r)
    }
    if (S7::S7_inherits(r, MNIX_CategoryRepresentation)) {
      return(new_nix_category_representation(
        category_labels = cat_labels,
        cue_labels = cues,
        m = r@m[idx],
        kappa = r@kappa[idx],
        nu = r@nu[idx],
        sigma2 = r@sigma2[idx]
      ))
    }
    if (S7::S7_inherits(r, NIW_CategoryRepresentation)) {
      D_orig <- length(rep_cues)
      nu_sub <- r@nu - (D_orig - 1L)
      return(new_nix_category_representation(
        category_labels = cat_labels,
        cue_labels = cues,
        m = r@m[idx],
        kappa = r@kappa,
        nu = max(nu_sub, 3.01),
        sigma2 = r@S[idx, idx] / max(nu_sub, 3.01)
      ))
    }
    if (S7::S7_inherits(r, Exemplar_CategoryRepresentation)) {
      coords <- r@exemplars[, idx, drop = FALSE]
      return(new_exemplar_category_representation(
        category_labels = cat_labels,
        cue_labels = cues,
        exemplars = coords,
        exemplar_weights = r@exemplar_weights,
        c = r@c
      ))
    }
  }

  if (length(cues) == 2L) {
    idx <- match(cues, rep_cues)

    if (S7::S7_inherits(r, MVG_CategoryRepresentation)) {
      return(new_mvg_category_representation(
        category_labels = cat_labels,
        cue_labels = cues,
        mu = r@mu[idx],
        Sigma = r@Sigma[idx, idx]
      ))
    }
    if (S7::S7_inherits(r, MNIX_CategoryRepresentation)) {
      return(new_mnix_category_representation(
        category_labels = cat_labels,
        cue_labels = cues,
        m = r@m[idx],
        kappa = r@kappa[idx],
        nu = r@nu[idx],
        sigma2 = r@sigma2[idx]
      ))
    }
    if (S7::S7_inherits(r, NIW_CategoryRepresentation)) {
      D_orig <- length(rep_cues)
      nu_sub <- r@nu - (D_orig - 2L)
      return(new_niw_category_representation(
        category_labels = cat_labels,
        cue_labels = cues,
        m = r@m[idx],
        S = r@S[idx, idx],
        kappa = r@kappa,
        nu = max(nu_sub, 3.01)
      ))
    }
    if (S7::S7_inherits(r, Exemplar_CategoryRepresentation)) {
      coords <- r@exemplars[, idx, drop = FALSE]
      return(new_exemplar_category_representation(
        category_labels = cat_labels,
        cue_labels = cues,
        exemplars = coords,
        exemplar_weights = r@exemplar_weights,
        c = r@c
      ))
    }
  }

  if (length(cues) == 3L) {
    idx <- match(cues, rep_cues)

    if (S7::S7_inherits(r, MVG_CategoryRepresentation)) {
      return(new_mvg_category_representation(
        category_labels = cat_labels,
        cue_labels = cues,
        mu = r@mu[idx],
        Sigma = r@Sigma[idx, idx]
      ))
    }
    if (S7::S7_inherits(r, MNIX_CategoryRepresentation)) {
      return(new_mnix_category_representation(
        category_labels = cat_labels,
        cue_labels = cues,
        m = r@m[idx],
        kappa = r@kappa[idx],
        nu = r@nu[idx],
        sigma2 = r@sigma2[idx]
      ))
    }
    if (S7::S7_inherits(r, NIW_CategoryRepresentation)) {
      D_orig <- length(rep_cues)
      nu_sub <- r@nu - (D_orig - 3L)
      return(new_niw_category_representation(
        category_labels = cat_labels,
        cue_labels = cues,
        m = r@m[idx],
        S = r@S[idx, idx],
        kappa = r@kappa,
        nu = max(nu_sub, 4.01)
      ))
    }
    if (S7::S7_inherits(r, Exemplar_CategoryRepresentation)) {
      coords <- r@exemplars[, idx, drop = FALSE]
      return(new_exemplar_category_representation(
        category_labels = cat_labels,
        cue_labels = cues,
        exemplars = coords,
        exemplar_weights = r@exemplar_weights,
        c = r@c
      ))
    }
  }

  stop("Subsetting to length(cues) > 3 is not supported.")
}

#' Extract expected mean and covariance matrix from representation
#' @noRd
#' @keywords internal
.get_rep_mean_and_cov <- function(r) {
  if (S7::S7_inherits(r, MVG_CategoryRepresentation)) {
    return(list(mu = r@mu, Sigma = r@Sigma))
  }
  if (S7::S7_inherits(r, UVG_CategoryRepresentation)) {
    return(list(mu = r@mu, Sigma = matrix(r@sigma2, 1L, 1L)))
  }
  if (S7::S7_inherits(r, MUVG_CategoryRepresentation)) {
    return(list(mu = r@mu, Sigma = diag(r@sigma2, nrow = length(r@mu))))
  }
  if (S7::S7_inherits(r, NIX_CategoryRepresentation)) {
    sig2 <- r@sigma2 * (r@nu / max(r@nu - 2, 1))
    return(list(mu = r@m, Sigma = matrix(sig2, 1L, 1L)))
  }
  if (S7::S7_inherits(r, MNIX_CategoryRepresentation)) {
    sig2 <- r@sigma2 * (r@nu / max(r@nu - 2, 1))
    return(list(mu = r@m, Sigma = diag(sig2, nrow = length(r@m))))
  }
  if (S7::S7_inherits(r, NIW_CategoryRepresentation)) {
    d <- length(r@m)
    sig_mat <- r@S / max(r@nu - d - 1, 1)
    return(list(mu = r@m, Sigma = sig_mat))
  }
  if (S7::S7_inherits(r, Exemplar_CategoryRepresentation)) {
    mat <- as.matrix(r@exemplars)
    mu <- colMeans(mat)
    cov_mat <- stats::cov(mat)
    if (ncol(mat) == 1L) {
      cov_mat <- matrix(cov_mat, 1L, 1L)
    }
    return(list(mu = mu, Sigma = cov_mat))
  }
  stop("Unsupported representation type.")
}

#' Create 1D density data frame across category representations
#' @noRd
#' @keywords internal
.make_1D_category_density_df <- function(
  reps,
  cues,
  resolution = 200,
  limits = NULL
) {
  min_val <- Inf
  max_val <- -Inf
  params <- list()

  for (cat_name in names(reps)) {
    r <- reps[[cat_name]]
    mc <- .get_rep_mean_and_cov(r)
    mu_val <- mc$mu[1]
    sd_val <- sqrt(max(mc$Sigma[1, 1], 1e-6))
    params[[cat_name]] <- list(mu = mu_val, sd = sd_val, rep = r)
    min_val <- min(min_val, mu_val - 3.5 * sd_val)
    max_val <- max(max_val, mu_val + 3.5 * sd_val)
  }

  lim_spec <- .parse_cue_limits(limits, cues)
  if (!is.null(lim_spec[[cues[1L]]])) {
    min_val <- lim_spec[[cues[1L]]][1L]
    max_val <- lim_spec[[cues[1L]]][2L]
  }

  grid_x <- seq(min_val, max_val, length.out = resolution)
  grid_mat <- matrix(grid_x, ncol = 1L)
  colnames(grid_mat) <- cues[1]
  dfs <- list()

  for (cat_name in names(params)) {
    p <- params[[cat_name]]
    dens <- likelihood(p$rep, grid_mat)
    df <- tibble::tibble(
      cue = grid_x,
      density = dens,
      Category = cat_name
    )
    names(df)[1] <- cues[1]
    dfs[[length(dfs) + 1]] <- df
  }

  dplyr::bind_rows(dfs)
}

#' Create 2D category center points data frame
#' @noRd
#' @keywords internal
.make_category_centers_df <- function(reps, cues) {
  rows <- list()
  for (cat_name in names(reps)) {
    r <- reps[[cat_name]]
    mc <- .get_rep_mean_and_cov(r)
    df <- data.frame(
      x = mc$mu[1L],
      y = mc$mu[2L],
      Category = cat_name,
      stringsAsFactors = FALSE
    )
    names(df)[1:2] <- cues[1:2]
    rows[[length(rows) + 1]] <- df
  }
  dplyr::bind_rows(rows)
}

#' Create 2D exemplar sample points data frame
#' @noRd
#' @keywords internal
.make_exemplar_sample_df <- function(reps, cues, n_exemplars = 0L) {
  if (is.null(n_exemplars) || n_exemplars <= 0L) {
    return(NULL)
  }
  rows <- list()
  for (cat_name in names(reps)) {
    r <- reps[[cat_name]]
    if (S7::S7_inherits(r, Exemplar_CategoryRepresentation)) {
      mat <- as.matrix(r@exemplars)
      N <- nrow(mat)
      n_to_draw <- min(as.integer(n_exemplars), N)
      if (n_to_draw > 0L) {
        idx <- sample.int(N, n_to_draw, replace = FALSE)
        sub_mat <- mat[idx, cues, drop = FALSE]
        df <- as.data.frame(sub_mat)
        df$Category <- cat_name
        rows[[length(rows) + 1]] <- df
      }
    }
  }
  if (length(rows) == 0L) {
    return(NULL)
  }
  dplyr::bind_rows(rows)
}

#' Create 2D ellipse contour data frame across category representations
#' @noRd
#' @keywords internal
.make_2D_category_ellipse_df <- function(reps, cues, levels = c(0.5, 0.95)) {
  dfs <- list()
  sort_lvls <- sort(levels, decreasing = TRUE)

  for (cat_name in names(reps)) {
    r <- reps[[cat_name]]
    if (S7::S7_inherits(r, Exemplar_CategoryRepresentation)) {
      dens_df <- .make_2D_category_density_grid_df(
        reps = stats::setNames(list(r), cat_name),
        cues = cues,
        resolution = 60
      )
      sub_dens <- dens_df[dens_df$Category == cat_name, ]
      gx <- sort(unique(sub_dens[[cues[1L]]]))
      gy <- sort(unique(sub_dens[[cues[2L]]]))
      z_mat <- matrix(sub_dens$Density, nrow = length(gx), ncol = length(gy))
      d_vals <- sort(sub_dens$Density, decreasing = TRUE)
      c_mass <- cumsum(d_vals) / max(sum(d_vals), 1e-12)

      for (k in seq_along(sort_lvls)) {
        lvl <- sort_lvls[k]
        brk <- d_vals[which.min(abs(c_mass - lvl))]
        cl <- grDevices::contourLines(x = gx, y = gy, z = z_mat, levels = brk)
        for (j in seq_along(cl)) {
          n_pts <- length(cl[[j]]$x)
          df_p <- data.frame(
            x_val = cl[[j]]$x,
            y_val = cl[[j]]$y,
            level = lvl,
            level_label = sprintf("%d%%", round(lvl * 100)),
            alpha_val = (k / length(sort_lvls)) * 0.35,
            Category = cat_name,
            group = paste0(cat_name, ".", lvl, ".", j),
            stringsAsFactors = FALSE
          )
          names(df_p)[1:2] <- cues
          dfs[[length(dfs) + 1L]] <- df_p
        }
      }
    } else {
      mc <- .get_rep_mean_and_cov(r)
      mu_vec <- mc$mu
      sigma_mat <- mc$Sigma

      for (lvl in levels) {
        el <- ellipse::ellipse(sigma_mat, centre = mu_vec, level = lvl)
        df <- tibble::as_tibble(el)
        names(df) <- cues[1:2]
        df$level <- lvl
        df$level_label <- sprintf("%d%%", round(lvl * 100))
        df$alpha_val <- 1 - lvl
        df$Category <- cat_name
        df$group <- interaction(cat_name, lvl)
        dfs[[length(dfs) + 1L]] <- df
      }
    }
  }

  dplyr::bind_rows(dfs)
}

#' Create 2D density surface grid data frame across category representations
#' @noRd
#' @keywords internal
.make_2D_category_density_grid_df <- function(
  reps,
  cues,
  resolution = 60,
  limits = NULL
) {
  min_x <- Inf
  max_x <- -Inf
  min_y <- Inf
  max_y <- -Inf
  params <- list()

  for (cat_name in names(reps)) {
    r <- reps[[cat_name]]
    if (S7::S7_inherits(r, Exemplar_CategoryRepresentation)) {
      mat <- as.matrix(r@exemplars)
      rng1 <- diff(range(mat[, cues[1]]))
      rng2 <- diff(range(mat[, cues[2]]))
      min_x <- min(min_x, min(mat[, cues[1]]) - 0.25 * rng1)
      max_x <- max(max_x, max(mat[, cues[1]]) + 0.25 * rng1)
      min_y <- min(min_y, min(mat[, cues[2]]) - 0.25 * rng2)
      max_y <- max(max_y, max(mat[, cues[2]]) + 0.25 * rng2)
    } else {
      mc <- .get_rep_mean_and_cov(r)
      mu_vec <- mc$mu
      sd_x <- sqrt(max(mc$Sigma[1, 1], 1e-6))
      sd_y <- sqrt(max(mc$Sigma[2, 2], 1e-6))
      min_x <- min(min_x, mu_vec[1] - 3.5 * sd_x)
      max_x <- max(max_x, mu_vec[1] + 3.5 * sd_x)
      min_y <- min(min_y, mu_vec[2] - 3.5 * sd_y)
      max_y <- max(max_y, mu_vec[2] + 3.5 * sd_y)
    }
    params[[cat_name]] <- list(rep = r)
  }

  lim_spec <- .parse_cue_limits(limits, cues)
  if (!is.null(lim_spec[[cues[1L]]])) {
    min_x <- lim_spec[[cues[1L]]][1L]
    max_x <- lim_spec[[cues[1L]]][2L]
  }
  if (!is.null(lim_spec[[cues[2L]]])) {
    min_y <- lim_spec[[cues[2L]]][1L]
    max_y <- lim_spec[[cues[2L]]][2L]
  }

  gx <- seq(min_x, max_x, length.out = resolution)
  gy <- seq(min_y, max_y, length.out = resolution)
  grid_df <- expand.grid(x = gx, y = gy)
  names(grid_df) <- cues[1:2]
  grid_mat <- as.matrix(grid_df)

  dfs <- list()
  for (cat_name in names(params)) {
    r <- params[[cat_name]]$rep
    dens <- likelihood(r, grid_mat)
    df_cat <- grid_df
    df_cat$Density <- dens
    df_cat$Category <- cat_name
    dfs[[length(dfs) + 1]] <- df_cat
  }

  dplyr::bind_rows(dfs)
}

#' Overlay empirical exposure and test data onto a stanfit plot
#' @noRd
#' @keywords internal
.add_stanfit_data_layers <- function(
  p,
  x,
  show_exposure_data,
  show_test_data,
  cues
) {
  if (isTRUE(show_exposure_data)) {
    exp_data <- tryCatch(
      suppressMessages(get_exposure_data(x, original_names = FALSE)),
      error = function(e) NULL
    )
    if (!is.null(exp_data) && nrow(exp_data) > 0L && all(cues %in% names(exp_data))) {
      exp_data$Category <- exp_data$category
      if (length(cues) == 1L) {
        p <- p + geom_rug(
          data = exp_data,
          aes(x = .data[[cues[1]]], color = .data$Category),
          alpha = 0.6,
          inherit.aes = FALSE
        )
      } else {
        p <- p + geom_point(
          data = exp_data,
          aes(
            x = .data[[cues[1]]],
            y = .data[[cues[2]]],
            color = .data$Category
          ),
          alpha = 0.6,
          size = 1.5,
          inherit.aes = FALSE
        )
      }
    }
  }

  if (isTRUE(show_test_data)) {
    test_data <- tryCatch(
      suppressMessages(get_test_data(x, original_names = FALSE)),
      error = function(e) NULL
    )
    if (!is.null(test_data) && nrow(test_data) > 0L && all(cues %in% names(test_data))) {
      if (length(cues) == 1L) {
        p <- p + geom_rug(
          data = test_data,
          aes(x = .data[[cues[1]]]),
          color = "darkgray",
          alpha = 0.6,
          inherit.aes = FALSE
        )
      } else {
        p <- p + geom_point(
          data = test_data,
          aes(x = .data[[cues[1]]], y = .data[[cues[2]]]),
          color = "darkgray",
          alpha = 0.5,
          size = 1,
          inherit.aes = FALSE
        )
      }
    }
  }

  p
}

# -----------------------------------------------------------------------------
# Category Plot Rendering Helpers
# -----------------------------------------------------------------------------

#' Render 1D category density plot
#' @noRd
#' @keywords internal
.render_1D_category_plot <- function(
  reps,
  cues,
  aes = "contour",
  levels = NULL,
  limits = NULL,
  n_exemplars = 0L,
  resolution = 200,
  t_sub = list(title = "", subtitle = "")
) {
  aes <- .normalize_plot_aes(aes, default_aes = "contour")
  lvl_spec <- .parse_category_levels(levels, aes)
  lim_spec <- .parse_cue_limits(limits, cues)

  ex_df_1d <- .make_exemplar_sample_df(reps, cues[1L], n_exemplars = n_exemplars)
  if (!is.null(ex_df_1d) && nrow(ex_df_1d) > 0L) {
    samp_txt <- sprintf(" (sampling %d exemplars)", nrow(ex_df_1d))
    t_sub$subtitle <- if (is.language(t_sub$subtitle)) {
      bquote(.(t_sub$subtitle) ~ .(samp_txt))
    } else {
      paste0(t_sub$subtitle, samp_txt)
    }
  }

  df <- .make_1D_category_density_df(
    reps,
    cues,
    resolution = resolution,
    limits = limits
  )

  all_c <- names(reps)
  if (is.null(all_c) || length(all_c) == 0L) {
    all_c <- vapply(reps, get_category_labels, character(1))
  }
  c_colors <- if (length(all_c) == 0L) character(0) else scales::hue_pal()(length(all_c))
  names(c_colors) <- all_c

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[cues[1L]]],
      y = .data$density,
      color = .data$Category,
      fill = .data$Category
    )
  ) +
    ggplot2::labs(
      title = t_sub$title,
      subtitle = t_sub$subtitle,
      x = cues[1L],
      y = "Density"
    ) +
    .mvbu_theme() +
    ggplot2::scale_color_manual(values = c_colors, name = "Category") +
    ggplot2::scale_fill_manual(values = c_colors, name = "Category")

  has_fill <- any(c("fill-discrete", "fill", "fill-gradient") %in% aes)
  if (has_fill) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = 0, ymax = .data$density),
      alpha = 0.2,
      key_glyph = ggplot2::draw_key_rect
    )
  }
  if ("contour" %in% aes || !has_fill) {
    p <- p + ggplot2::geom_line(linewidth = 0.30, key_glyph = ggplot2::draw_key_rect)
  }

  p <- p + ggplot2::guides(
    color = ggplot2::guide_legend(
      title = "Category",
      override.aes = list(fill = c_colors, alpha = 0.5)
    ),
    fill = "none"
  )

  if (!is.null(ex_df_1d) && nrow(ex_df_1d) > 0L) {
    p <- p + ggplot2::geom_rug(
      data = ex_df_1d,
      ggplot2::aes(
        x = .data[[cues[1L]]],
        color = .data$Category
      ),
      sides = "b",
      length = ggplot2::unit(0.0225, "npc"),
      alpha = 1 / length(reps),
      linewidth = 0.30,
      inherit.aes = FALSE,
      show.legend = FALSE
    )
  }

  if (!is.null(lim_spec[[cues[1L]]])) {
    p <- p + ggplot2::coord_cartesian(xlim = lim_spec[[cues[1L]]])
  }
  p
}

#' Render 2D category density / contour plot
#' @noRd
#' @keywords internal
.render_2D_category_plot <- function(
  reps,
  cues,
  aes = "contour",
  levels = NULL,
  limits = NULL,
  n_exemplars = 0L,
  resolution = 100,
  t_sub = list(title = "", subtitle = "")
) {
  aes <- .normalize_plot_aes(aes, default_aes = "contour")
  lvl_spec <- .parse_category_levels(levels, aes)
  lim_spec <- .parse_cue_limits(limits, cues)

  ex_df <- .make_exemplar_sample_df(reps, cues, n_exemplars = n_exemplars)
  if (!is.null(ex_df) && nrow(ex_df) > 0L) {
    samp_txt <- sprintf(" (sampling %d exemplars)", nrow(ex_df))
    t_sub$subtitle <- if (is.language(t_sub$subtitle)) {
      bquote(.(t_sub$subtitle) ~ .(samp_txt))
    } else {
      paste0(t_sub$subtitle, samp_txt)
    }
  }

  has_ex <- any(sapply(reps, function(r) {
    S7::S7_inherits(r, Exemplar_CategoryRepresentation)
  }))

  has_contour <- "contour" %in% aes

  all_c <- names(reps)
  if (is.null(all_c) || length(all_c) == 0L) {
    all_c <- vapply(reps, get_category_labels, character(1))
  }
  c_colors <- if (length(all_c) == 0L) character(0) else scales::hue_pal()(length(all_c))
  names(c_colors) <- all_c

  p <- ggplot2::ggplot() +
    ggplot2::labs(
      title = t_sub$title,
      subtitle = t_sub$subtitle,
      x = cues[1L],
      y = cues[2L]
    ) +
    .mvbu_theme()

  dens_df <- NULL
  # 1. Fill Layer
  if ("fill-gradient" %in% aes) {
    dens_df <- .make_2D_category_density_grid_df(
      reps,
      cues,
      resolution = max(resolution, 60),
      limits = limits
    )
    dens_df_tile <- dens_df %>%
      dplyr::group_by(.data$Category) %>%
      dplyr::mutate(
        rel_dens = .data$Density / max(.data$Density, 1e-12)
      ) %>%
      dplyr::ungroup()

    p <- p + ggplot2::geom_tile(
      data = dens_df_tile,
      ggplot2::aes(
        x = .data[[cues[1L]]],
        y = .data[[cues[2L]]],
        fill = .data$Category,
        alpha = .data$rel_dens
      ),
      show.legend = c(fill = TRUE, alpha = FALSE)
    ) +
      ggplot2::scale_fill_manual(values = c_colors, name = "Category") +
      ggplot2::scale_alpha_continuous(
        range = c(0, 0.85),
        limits = c(0, 1),
        guide = "none"
      )
  } else if (any(c("fill-discrete", "fill") %in% aes)) {
    df_el_fill <- .make_2D_category_ellipse_df(
      reps,
      cues,
      levels = lvl_spec$fill
    )
    p <- p + ggplot2::geom_polygon(
      data = df_el_fill,
      ggplot2::aes(
        x = .data[[cues[1L]]],
        y = .data[[cues[2L]]],
        fill = .data$Category,
        alpha = .data$alpha_val,
        group = .data$group
      ),
      show.legend = c(fill = TRUE, alpha = FALSE)
    ) +
      ggplot2::scale_fill_manual(values = c_colors, name = "Category") +
      ggplot2::scale_alpha_identity()
  }

  # 2. Contour Layer
  if (has_contour) {
    if (has_ex) {
      if (is.null(dens_df)) {
        dens_df <- .make_2D_category_density_grid_df(
          reps,
          cues,
          resolution = max(resolution, 60),
          limits = limits
        )
      }
      for (c_name in all_c) {
        sub_dens <- dens_df[dens_df$Category == c_name, ]
        d_vals <- sort(sub_dens$Density, decreasing = TRUE)
        c_mass <- cumsum(d_vals) / max(sum(d_vals), 1e-12)
        brks <- sapply(lvl_spec$contour, function(lvl) {
          d_vals[which.min(abs(c_mass - lvl))]
        })
        brks <- unique(brks[brks > 0])
        if (length(brks) > 0L) {
          p <- p + ggplot2::geom_contour(
            data = sub_dens,
            ggplot2::aes(
              x = .data[[cues[1L]]],
              y = .data[[cues[2L]]],
              z = .data$Density,
              color = .data$Category
            ),
            breaks = brks,
            linewidth = 0.35,
            alpha = 1.0,
            key_glyph = ggplot2::draw_key_rect,
            show.legend = c(color = TRUE)
          )
        }
      }
      p <- p + ggplot2::scale_color_manual(
        values = c_colors,
        name = "Category"
      )
    } else {
      df_el_cont <- .make_2D_category_ellipse_df(
        reps,
        cues,
        levels = lvl_spec$contour
      )
      p <- p + ggplot2::geom_path(
        data = df_el_cont,
        ggplot2::aes(
          x = .data[[cues[1L]]],
          y = .data[[cues[2L]]],
          color = .data$Category,
          group = .data$group
        ),
        linewidth = 0.35,
        alpha = 1.0,
        key_glyph = ggplot2::draw_key_rect,
        show.legend = c(color = TRUE)
      ) +
        ggplot2::scale_color_manual(values = c_colors, name = "Category")
    }
  }

  p <- p + ggplot2::guides(
    color = ggplot2::guide_legend(
      title = "Category",
      override.aes = list(fill = c_colors, alpha = 0.5)
    ),
    fill = "none"
  )

  # 3. Center Mean Point (size = 1.8)
  centers_df <- .make_category_centers_df(reps, cues)
  p <- p + ggplot2::geom_point(
    data = centers_df,
    ggplot2::aes(
      x = .data[[cues[1L]]],
      y = .data[[cues[2L]]],
      color = .data$Category
    ),
    size = 1.8,
    show.legend = FALSE
  )

  # 4. Exemplar Points (size = 1.0, alpha = 1 / length(reps))
  if (!is.null(ex_df) && nrow(ex_df) > 0L) {
    p <- p + ggplot2::geom_point(
      data = ex_df,
      ggplot2::aes(
        x = .data[[cues[1L]]],
        y = .data[[cues[2L]]],
        color = .data$Category
      ),
      alpha = 1 / length(reps),
      size = 1.0,
      shape = 16,
      show.legend = FALSE
    )
  }

  if (!is.null(lim_spec)) {
    p <- p + ggplot2::coord_cartesian(
      xlim = lim_spec[[cues[1L]]],
      ylim = lim_spec[[cues[2L]]],
      expand = FALSE
    )
  } else if ("fill-gradient" %in% aes && !is.null(dens_df)) {
    p <- p + ggplot2::coord_cartesian(
      xlim = range(dens_df[[cues[1L]]]),
      ylim = range(dens_df[[cues[2L]]]),
      expand = FALSE
    )
  }

  p
}

# -----------------------------------------------------------------------------
# plot_categories methods
# -----------------------------------------------------------------------------

#' @rdname plot_categories
#' @export
S7::method(plot_categories, MVBU_CategoryRepresentation) <- function(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  resolution = 100,
  n_exemplars = 0L,
  ...
) {
  obj_cues <- get_cue_labels(x)
  if (is.null(cues)) {
    cues <- if (length(obj_cues) > 3L) obj_cues[1:3] else obj_cues
  }
  if (length(cues) > 3L) {
    .stop(
      "Cannot plot more than 3 cue dimensions simultaneously (received ",
      length(cues), " cues: ", paste(cues, collapse = ", "),
      "). Please specify 1, 2, or 3 cues; non-plotted cues will be analytically marginalized out."
    )
  }
  if (!all(cues %in% obj_cues)) {
    .stop(
      "Invalid cues: ",
      paste(setdiff(cues, obj_cues), collapse = ", "),
      ". Available cues: ",
      paste(obj_cues, collapse = ", ")
    )
  }

  x_proj <- if (length(cues) < length(obj_cues)) {
    .marginalize_representation_to_cues(x, cues)
  } else {
    x
  }
  cats <- get_category_labels(x_proj)
  cat_name <- if (length(cats) > 0L) cats[1L] else "Category"
  reps <- list()
  reps[[cat_name]] <- x_proj

  if (!is.null(categories)) {
    reps <- reps[names(reps) %in% categories]
  }

  t_sub <- .make_plot_title_and_subtitle(x, "categories", cues = cues)

  if (length(cues) == 1L) {
    return(.render_1D_category_plot(
      reps = reps,
      cues = cues,
      aes = aes,
      levels = levels,
      limits = limits,
      n_exemplars = n_exemplars,
      resolution = resolution,
      t_sub = t_sub
    ))
  }

  if (length(cues) == 2L) {
    if (isTRUE(list(...)$interactive) || "interactive" %in% aes) {
      return(.render_2D_interactive_category_plot(
        reps = reps,
        cues = cues,
        aes = aes,
        levels = levels,
        limits = limits,
        n_exemplars = n_exemplars,
        resolution = resolution,
        t_sub = t_sub,
        ...
      ))
    }
    return(.render_2D_category_plot(
      reps = reps,
      cues = cues,
      aes = aes,
      levels = levels,
      limits = limits,
      n_exemplars = n_exemplars,
      resolution = resolution,
      t_sub = t_sub
    ))
  }

  if (length(cues) == 3L) {
    if (isTRUE(list(...)$interactive) || "interactive" %in% aes) {
      return(.render_3D_interactive_category_plot(
        reps = reps,
        cues = cues,
        levels = levels,
        n_exemplars = n_exemplars,
        resolution = resolution,
        t_sub = t_sub,
        ...
      ))
    }
    return(.render_3D_sliced_category_plot(
      reps = reps,
      cues = cues,
      aes = aes,
      levels = levels,
      limits = limits,
      resolution = resolution,
      t_sub = t_sub,
      ...
    ))
  }

  stop("Plotting categories is currently supported for 1, 2, or 3 cues.")
}

#' @rdname plot_categories
#' @export
S7::method(plot_categories, MVBU_CategoryRepresentationTemplate) <- function(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  resolution = 100,
  n_exemplars = 0L,
  ...
) {
  obj_cues <- get_cue_labels(x)
  if (is.null(cues)) {
    cues <- if (length(obj_cues) > 3L) obj_cues[1:3] else obj_cues
  }
  if (length(cues) > 3L) {
    .stop(
      "Cannot plot more than 3 cue dimensions simultaneously (received ",
      length(cues), " cues: ", paste(cues, collapse = ", "),
      "). Please specify 1, 2, or 3 cues; non-plotted cues will be analytically marginalized out."
    )
  }
  if (!all(cues %in% obj_cues)) {
    .stop(
      "Invalid cues: ",
      paste(setdiff(cues, obj_cues), collapse = ", "),
      ". Available cues: ",
      paste(obj_cues, collapse = ", ")
    )
  }

  reps <- x@representations
  if (!is.null(categories)) {
    reps <- reps[names(reps) %in% categories]
  }
  if (length(reps) == 0L) {
    stop("No category representations match the specified categories.")
  }

  reps_proj <- if (length(cues) < length(obj_cues)) {
    lapply(reps, function(r) {
      .marginalize_representation_to_cues(r, cues)
    })
  } else {
    reps
  }

  t_sub <- .make_plot_title_and_subtitle(x, "categories", cues = cues)

  if (length(cues) == 3L) {
    if (isTRUE(list(...)$interactive) || "interactive" %in% aes) {
      return(.render_3D_interactive_category_plot(
        reps = reps_proj,
        cues = cues,
        levels = levels,
        n_exemplars = n_exemplars,
        resolution = resolution,
        t_sub = t_sub,
        ...
      ))
    }
    return(.render_3D_sliced_category_plot(
      reps = reps_proj,
      cues = cues,
      aes = aes,
      levels = levels,
      limits = limits,
      resolution = resolution,
      t_sub = t_sub,
      ...
    ))
  }

  if (length(cues) == 1L) {
    return(.render_1D_category_plot(
      reps = reps_proj,
      cues = cues,
      aes = aes,
      levels = levels,
      limits = limits,
      n_exemplars = n_exemplars,
      resolution = resolution,
      t_sub = t_sub
    ))
  }

  if (length(cues) == 2L) {
    if (isTRUE(list(...)$interactive) || "interactive" %in% aes) {
      return(.render_2D_interactive_category_plot(
        reps = reps_proj,
        cues = cues,
        aes = aes,
        levels = levels,
        limits = limits,
        n_exemplars = n_exemplars,
        resolution = resolution,
        t_sub = t_sub,
        ...
      ))
    }
    return(.render_2D_category_plot(
      reps = reps_proj,
      cues = cues,
      aes = aes,
      levels = levels,
      limits = limits,
      n_exemplars = n_exemplars,
      resolution = resolution,
      t_sub = t_sub
    ))
  }

  stop("Plotting categories is currently supported for 1, 2, or 3 cues.")
}

#' @rdname plot_categories
#' @export
S7::method(plot_categories, MVBU_CognitiveModel) <- function(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  resolution = 100,
  n_exemplars = 0L,
  ...
) {
  dots <- list(...)
  parallel <- if (!is.null(dots$parallel)) isTRUE(dots$parallel) else FALSE
  n_cores <- dots$n_cores
  noise_treatment <- dots$noise_treatment
  obj_cues <- get_cue_labels(x)
  if (is.null(cues)) {
    cues <- if (length(obj_cues) > 3L) obj_cues[1:3] else obj_cues
  }
  if (length(cues) > 3L) {
    .stop(
      "Cannot plot more than 3 cue dimensions simultaneously (received ",
      length(cues), " cues: ", paste(cues, collapse = ", "),
      "). Please specify 1, 2, or 3 cues; non-plotted cues will be analytically marginalized out."
    )
  }
  if (!all(cues %in% obj_cues)) {
    .stop(
      "Invalid cues: ",
      paste(setdiff(cues, obj_cues), collapse = ", "),
      ". Available cues: ",
      paste(obj_cues, collapse = ", ")
    )
  }

  tmpl <- x@category_template
  eff_noise_trt <- noise_treatment %||% x@noise_behavior$noise_treatment
  if (identical(eff_noise_trt, "marginalize") && !is.null(x@noise_behavior$Sigma_noise)) {
    reps <- tmpl@representations
    sn <- as.matrix(x@noise_behavior$Sigma_noise)
    reps <- lapply(reps, function(r) {
      if (S7::S7_inherits(r, UVG_CategoryRepresentation)) {
        r@sigma2 <- r@sigma2 + sn[1, 1]
      } else if (S7::S7_inherits(r, MVG_CategoryRepresentation)) {
        r@Sigma <- r@Sigma + sn
      } else if (S7::S7_inherits(r, NIX_CategoryRepresentation)) {
        r@sigma2 <- r@sigma2 + sn[1, 1]
      } else if (S7::S7_inherits(r, MNIX_CategoryRepresentation)) {
        r@sigma2 <- r@sigma2 + diag(sn)
      } else if (S7::S7_inherits(r, NIW_CategoryRepresentation)) {
        D <- ncol(sn)
        denom <- r@nu - D - 1
        if (denom > 0) {
          exp_Sigma <- r@S / denom
          r@S <- (exp_Sigma + sn) * denom
        }
      }
      r
    })
    tmpl@representations <- reps
  }

  p <- plot_categories(
    tmpl,
    cues = cues,
    categories = categories,
    aes = aes,
    levels = levels,
    limits = limits,
    n_exemplars = n_exemplars,
    resolution = resolution,
    parallel = parallel,
    n_cores = n_cores,
    ...
  )
  t_sub <- .make_plot_title_and_subtitle(x, "categories", cues = cues)
  if (inherits(p, "plotly")) {
    return(p)
  }
  p + ggplot2::labs(title = t_sub$title, subtitle = t_sub$subtitle)
}

#' @rdname plot_categories
#' @export
S7::method(plot_categories, MVBU_Stanfit) <- function(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  resolution = 100,
  n_exemplars = 0L,
  ...
) {
  dots <- list(...)
  uncertainty_treatment <- if ("uncertainty_treatment" %in% names(dots)) {
    match.arg(dots$uncertainty_treatment, c("marginalize", "sample"))
  } else {
    "marginalize"
  }
  ndraws <- if ("ndraws" %in% names(dots)) dots$ndraws else 20L
  show_exposure_data <- if ("show_exposure_data" %in% names(dots)) {
    dots$show_exposure_data
  } else {
    FALSE
  }
  show_test_data <- if ("show_test_data" %in% names(dots)) {
    dots$show_test_data
  } else {
    FALSE
  }

  obj_cues <- get_cue_labels(x)
  if (is.null(cues)) {
    cues <- if (length(obj_cues) > 3L) obj_cues[1:3] else obj_cues
  }
  if (length(cues) > 3L) {
    .stop(
      "Cannot plot more than 3 cue dimensions simultaneously (received ",
      length(cues), " cues: ", paste(cues, collapse = ", "),
      "). Please specify 1, 2, or 3 cues; non-plotted cues will be analytically marginalized out."
    )
  }
  if (!all(cues %in% obj_cues)) {
    .stop(
      "Invalid cues: ",
      paste(setdiff(cues, obj_cues), collapse = ", "),
      ". Available cues: ",
      paste(obj_cues, collapse = ", ")
    )
  }

  if (is.null(categories)) {
    categories <- get_category_labels(x)
  }
  if (is.null(groups)) {
    groups <- get_group_labels(x, include_prior = FALSE)
  }

  # In Bayesian category plotting for Stan models, posterior uncertainty is
  # represented by drawing posterior parameters and calculating draw-wise
  # summaries and/or draw trajectories
  draws_args <- dots
  draws_args[c(
    "ndraws",
    "uncertainty_treatment",
    "show_exposure_data",
    "show_test_data"
  )] <- NULL
  d_raw <- do.call(
    get_draws,
    c(
      list(
        x,
        categories = categories,
        groups = groups,
        ndraws = ndraws,
        summarize = FALSE
      ),
      draws_args
    )
  )
  if ("Sigma_marg" %in% names(d_raw)) {
    d_raw$Sigma <- d_raw$Sigma_marg
  } else if (!"Sigma" %in% names(d_raw)) {
    if ("S" %in% names(d_raw) && "nu" %in% names(d_raw)) {
      sig_conv <- get_expected_Sigma_from_S(d_raw$S, d_raw$nu)
      d_raw$Sigma <- if (is.list(sig_conv)) sig_conv else as.list(sig_conv)
    } else if ("sigma2" %in% names(d_raw)) {
      # For univariate NIX or MNIX draws
      d_raw$Sigma <- lapply(d_raw$sigma2, function(s2) {
        if (is.matrix(s2)) {
          s2
        } else if (length(s2) == 1L) {
          matrix(s2, 1L, 1L)
        } else {
          diag(as.numeric(s2), nrow = length(s2))
        }
      })
    }
  }
  d_sum <- .summarise_category_parameter_draws(d_raw)

  if (is.null(aes)) {
    aes <- if (length(cues) == 1L) "contour" else c("fill", "contour")
  }
  lvl_spec <- .parse_plot_levels(levels, aes)
  lim_spec <- .parse_cue_limits(limits, cues)
  t_sub <- .make_plot_title_and_subtitle(
    x,
    "categories",
    cues = cues,
    ndraws = ndraws
  )

  if (length(cues) == 1L) {
    dfs <- list()
    cue_idx <- match(cues[1], obj_cues)
    for (i in seq_len(nrow(d_sum))) {
      grp <- d_sum$group[i]
      cat_name <- d_sum$category[i]
      mu_sub <- d_sum$mu.mean[[i]][cue_idx]
      sigma_mat <- as.matrix(d_sum$Sigma.mean[[i]])
      sigma_sub <- sqrt(max(sigma_mat[cue_idx, cue_idx], 1e-6))
      min_x <- if (!is.null(lim_spec[[cues[1L]]])) {
        lim_spec[[cues[1L]]][1L]
      } else {
        mu_sub - 3.5 * sigma_sub
      }
      max_x <- if (!is.null(lim_spec[[cues[1L]]])) {
        lim_spec[[cues[1L]]][2L]
      } else {
        mu_sub + 3.5 * sigma_sub
      }
      grid_x <- seq(min_x, max_x, length.out = resolution)
      df <- tibble::tibble(
        cue = grid_x,
        density = stats::dnorm(grid_x, mean = mu_sub, sd = sigma_sub),
        Category = cat_name,
        Group = grp
      )
      names(df)[1] <- cues[1]
      dfs[[length(dfs) + 1]] <- df
    }
    df_all <- dplyr::bind_rows(dfs)

    p <- ggplot(
      df_all,
      aes(
        x = .data[[cues[1]]],
        y = .data$density,
        color = .data$Category,
        fill = .data$Category
      )
    ) +
      labs(
        title = t_sub$title,
        subtitle = t_sub$subtitle,
        x = cues[1],
        y = "Density"
      ) +
      .mvbu_theme()

    if (identical(uncertainty_treatment, "sample") && !is.null(d_raw)) {
      dfs_samples <- list()
      for (i in seq_len(nrow(d_raw))) {
        grp <- d_raw$group[i]
        cat_name <- d_raw$category[i]
        mu_sub <- d_raw$m[[i]][cue_idx]
        sigma_mat <- as.matrix(d_raw$Sigma[[i]])
        sigma_sub <- sqrt(max(sigma_mat[cue_idx, cue_idx], 1e-6))
        min_x <- min(df_all[[cues[1]]])
        max_x <- max(df_all[[cues[1]]])
        grid_x <- seq(min_x, max_x, length.out = resolution)
        df_s <- tibble::tibble(
          cue = grid_x,
          density = stats::dnorm(grid_x, mean = mu_sub, sd = sigma_sub),
          Category = cat_name,
          Group = grp,
          .draw = d_raw$.draw[i]
        )
        names(df_s)[1] <- cues[1]
        dfs_samples[[length(dfs_samples) + 1]] <- df_s
      }
      df_samples_all <- dplyr::bind_rows(dfs_samples)
      p <- p + geom_line(
        data = df_samples_all,
        aes(
          x = .data[[cues[1]]],
          y = .data$density,
          color = .data$Category,
          group = interaction(.data$Category, .data$Group, .data$.draw)
        ),
        alpha = 0.20,
        linewidth = 0.35,
        show.legend = FALSE
      )
    }

    if ("fill" %in% aes && !identical(uncertainty_treatment, "sample")) {
      p <- p + geom_ribbon(
        aes(ymin = 0, ymax = .data$density),
        alpha = 0.2
      )
    }
    if ("contour" %in% aes) {
      p <- p + geom_line(linewidth = 1)
    }

    if (length(unique(df_all$Group)) > 1L) {
      p <- p + facet_wrap(~Group)
    }

    if (!is.null(lim_spec[[cues[1L]]])) {
      p <- p + coord_cartesian(xlim = lim_spec[[cues[1L]]])
    }

    p <- .add_stanfit_data_layers(
      p,
      x,
      show_exposure_data,
      show_test_data,
      cues
    )
    return(p)
  }

  if (length(cues) == 2L) {
    cue_idx <- match(cues[1:2], obj_cues)
    dfs_fill <- list()
    dfs_cont <- list()
    centers <- list()

    for (i in seq_len(nrow(d_sum))) {
      grp <- d_sum$group[i]
      cat_name <- d_sum$category[i]
      mu_sub <- d_sum$mu.mean[[i]][cue_idx]
      sigma_mat <- as.matrix(d_sum$Sigma.mean[[i]])[cue_idx, cue_idx]

      centers[[length(centers) + 1]] <- data.frame(
        x = mu_sub[1L],
        y = mu_sub[2L],
        Category = cat_name,
        Group = grp,
        stringsAsFactors = FALSE
      )

      for (lvl in lvl_spec$fill) {
        el <- ellipse::ellipse(sigma_mat, centre = mu_sub, level = lvl)
        df <- tibble::as_tibble(el)
        names(df) <- cues[1:2]
        df$level <- lvl
        df$level_label <- sprintf("%d%%", round(lvl * 100))
        df$alpha_val <- 1 - lvl
        df$Category <- cat_name
        df$Group <- grp
        dfs_fill[[length(dfs_fill) + 1]] <- df
      }

      for (lvl in lvl_spec$contour) {
        el <- ellipse::ellipse(sigma_mat, centre = mu_sub, level = lvl)
        df <- tibble::as_tibble(el)
        names(df) <- cues[1:2]
        df$level <- lvl
        df$level_label <- sprintf("%d%%", round(lvl * 100))
        df$alpha_val <- 1 - lvl
        df$Category <- cat_name
        df$Group <- grp
        dfs_cont[[length(dfs_cont) + 1]] <- df
      }
    }

    p <- ggplot() +
      labs(
        title = t_sub$title,
        subtitle = t_sub$subtitle,
        x = cues[1],
        y = cues[2]
      ) +
      .mvbu_theme()

    if (identical(uncertainty_treatment, "sample") && !is.null(d_raw)) {
      dfs_sample_ellipses <- list()
      top_lvl <- lvl_spec$contour[1L]
      for (i in seq_len(nrow(d_raw))) {
        grp <- d_raw$group[i]
        cat_name <- d_raw$category[i]
        mu_sub <- d_raw$m[[i]][cue_idx]
        sigma_mat <- as.matrix(d_raw$Sigma[[i]])[cue_idx, cue_idx]
        el <- ellipse::ellipse(sigma_mat, centre = mu_sub, level = top_lvl)
        df_s <- tibble::as_tibble(el)
        names(df_s) <- cues[1:2]
        df_s$Category <- cat_name
        df_s$Group <- grp
        df_s$.draw <- d_raw$.draw[i]
        dfs_sample_ellipses[[length(dfs_sample_ellipses) + 1]] <- df_s
      }
      df_sample_ell_all <- dplyr::bind_rows(dfs_sample_ellipses)
      p <- p + geom_path(
        data = df_sample_ell_all,
        aes(
          x = .data[[cues[1]]],
          y = .data[[cues[2]]],
          color = .data$Category,
          group = interaction(.data$Group, .data$Category, .data$.draw)
        ),
        alpha = 0.20,
        linewidth = 0.35,
        show.legend = FALSE
      )
    }

    if ("fill" %in% aes && length(dfs_fill) > 0L && !identical(uncertainty_treatment, "sample")) {
      df_fill <- dplyr::bind_rows(dfs_fill)
      p <- p + geom_polygon(
        data = df_fill,
        aes(
          x = .data[[cues[1]]],
          y = .data[[cues[2]]],
          fill = .data$Category,
          alpha = .data$alpha_val,
          group = interaction(.data$Group, .data$Category, .data$level)
        ),
        show.legend = c(fill = TRUE, alpha = FALSE)
      )
    }

    if ("contour" %in% aes && length(dfs_cont) > 0L) {
      df_cont <- dplyr::bind_rows(dfs_cont)
      p <- p + geom_path(
        data = df_cont,
        aes(
          x = .data[[cues[1]]],
          y = .data[[cues[2]]],
          color = .data$Category,
          alpha = .data$alpha_val,
          group = interaction(.data$Group, .data$Category, .data$level)
        ),
        linewidth = 1,
        show.legend = c(color = TRUE, alpha = FALSE)
      )
    }

    if (("fill" %in% aes && length(dfs_fill) > 0L) ||
      ("contour" %in% aes && length(dfs_cont) > 0L)) {
      p <- p + scale_alpha_identity()
    }

    # Category centers (size = 1.8)
    if (length(centers) > 0L) {
      centers_df <- dplyr::bind_rows(centers)
      names(centers_df)[1:2] <- cues[1:2]
      p <- p + geom_point(
        data = centers_df,
        aes(
          x = .data[[cues[1]]],
          y = .data[[cues[2]]],
          color = .data$Category
        ),
        size = 1.8,
        show.legend = FALSE
      )
    }

    if (length(unique(d_sum$group)) > 1L) {
      p <- p + facet_wrap(~Group)
    }

    if (!is.null(lim_spec)) {
      p <- p + coord_cartesian(
        xlim = lim_spec[[cues[1L]]],
        ylim = lim_spec[[cues[2L]]]
      )
    }

    p <- .add_stanfit_data_layers(
      p,
      x,
      show_exposure_data,
      show_test_data,
      cues
    )
    return(p)
  }

  stop("Plotting categories is currently supported for 1 or 2 cues.")
}

# -----------------------------------------------------------------------------
# plot_categorization_functions methods
# -----------------------------------------------------------------------------

#' Create a projected cognitive model on a subset of cues
#' @noRd
#' @keywords internal
.make_projected_cognitive_model <- function(x, cues, decision_rule = NULL, noise_treatment = NULL, lapse_treatment = NULL) {
  decision_rule <- decision_rule %||% tryCatch(x@decision_rule, error = function(e) "proportional")
  tpl <- x@category_template
  reps_proj <- lapply(tpl@representations, function(r) {
    .marginalize_representation_to_cues(r, cues)
  })
  first_rep <- reps_proj[[1L]]
  new_tpl <- new_category_representation_template(
    representations = reps_proj
  )

  noise_beh <- tryCatch(x@noise_behavior, error = function(e) list())
  lapse_beh <- tryCatch(x@lapse_behavior, error = function(e) list())

  noise_trt <- noise_treatment %||% noise_beh$noise_treatment %||% "no_noise"
  lapse_trt <- lapse_treatment %||% lapse_beh$lapse_treatment %||% "no_lapses"

  sig_noise <- noise_beh$Sigma_noise
  if (!is.null(sig_noise) && is.matrix(sig_noise)) {
    orig_cues <- get_cue_labels(x)
    idx <- match(cues, orig_cues)
    sig_noise <- sig_noise[idx, idx, drop = FALSE]
  }

  if (S7::S7_inherits(first_rep, UVG_CategoryRepresentation)) {
    return(new_uvg_ideal_observer(
      category_template = new_tpl,
      category_prior = x@category_prior,
      decision_rule = decision_rule,
      lapse_rate = lapse_beh$lapse_rate %||% 0,
      lapse_bias = lapse_beh$lapse_bias %||% (1 / length(reps_proj)),
      lapse_treatment = lapse_trt,
      Sigma_noise = sig_noise,
      noise_treatment = noise_trt
    ))
  }
  if (S7::S7_inherits(first_rep, MVG_CategoryRepresentation)) {
    return(new_mvg_ideal_observer(
      category_template = new_tpl,
      category_prior = x@category_prior,
      decision_rule = decision_rule,
      lapse_rate = lapse_beh$lapse_rate %||% 0,
      lapse_bias = lapse_beh$lapse_bias %||% (1 / length(reps_proj)),
      lapse_treatment = lapse_trt,
      Sigma_noise = sig_noise,
      noise_treatment = noise_trt
    ))
  }
  if (S7::S7_inherits(first_rep, MUVG_CategoryRepresentation)) {
    return(new_muvg_ideal_observer(
      category_template = new_tpl,
      category_prior = x@category_prior,
      decision_rule = decision_rule,
      lapse_rate = lapse_beh$lapse_rate %||% 0,
      lapse_bias = lapse_beh$lapse_bias %||% (1 / length(reps_proj)),
      lapse_treatment = lapse_trt,
      Sigma_noise = sig_noise,
      noise_treatment = noise_trt
    ))
  }
  if (S7::S7_inherits(first_rep, NIX_CategoryRepresentation)) {
    return(new_nix_ideal_adaptor(
      category_template = new_tpl,
      category_prior = x@category_prior,
      decision_rule = decision_rule,
      lapse_rate = lapse_beh$lapse_rate %||% 0,
      lapse_bias = lapse_beh$lapse_bias %||% (1 / length(reps_proj)),
      lapse_treatment = lapse_trt,
      Sigma_noise = sig_noise,
      noise_treatment = noise_trt
    ))
  }
  if (S7::S7_inherits(first_rep, MNIX_CategoryRepresentation)) {
    return(new_mnix_ideal_adaptor(
      category_template = new_tpl,
      category_prior = x@category_prior,
      decision_rule = decision_rule,
      lapse_rate = lapse_beh$lapse_rate %||% 0,
      lapse_bias = lapse_beh$lapse_bias %||% (1 / length(reps_proj)),
      lapse_treatment = lapse_trt,
      Sigma_noise = sig_noise,
      noise_treatment = noise_trt
    ))
  }
  if (S7::S7_inherits(first_rep, NIW_CategoryRepresentation)) {
    return(new_niw_ideal_adaptor(
      category_template = new_tpl,
      category_prior = x@category_prior,
      decision_rule = decision_rule,
      lapse_rate = lapse_beh$lapse_rate %||% 0,
      lapse_bias = lapse_beh$lapse_bias %||% (1 / length(reps_proj)),
      lapse_treatment = lapse_trt,
      Sigma_noise = sig_noise,
      noise_treatment = noise_trt
    ))
  }
  if (S7::S7_inherits(first_rep, Exemplar_CategoryRepresentation)) {
    return(new_exemplar_model(
      category_template = new_tpl,
      category_prior = x@category_prior,
      decision_rule = decision_rule,
      lapse_rate = lapse_beh$lapse_rate %||% 0,
      lapse_bias = lapse_beh$lapse_bias %||% (1 / length(reps_proj)),
      lapse_treatment = lapse_trt
    ))
  }
  stop("Unsupported representation type for projection.")
}

#' @rdname plot_categorization_functions
#' @export
S7::method(plot_categorization_functions, MVBU_CognitiveModel) <- function(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  decision_rule = "proportional",
  resolution = 100,
  noise_treatment = NULL,
  lapse_treatment = NULL,
  ...
) {
  if (is.null(categories)) {
    categories <- get_category_labels(x)[1L]
  }

  obj_cues <- get_cue_labels(x)
  if (is.null(cues)) {
    cues <- if (length(obj_cues) > 3L) obj_cues[1:3] else obj_cues
  }
  if (length(cues) > 3L) {
    .stop(
      "Cannot plot more than 3 cue dimensions simultaneously (received ",
      length(cues), " cues: ", paste(cues, collapse = ", "),
      "). Please specify 1, 2, or 3 cues; non-plotted cues will be analytically marginalized out."
    )
  }
  if (!all(cues %in% obj_cues)) {
    .stop(
      "Invalid cues: ",
      paste(setdiff(cues, obj_cues), collapse = ", "),
      ". Available cues: ",
      paste(obj_cues, collapse = ", ")
    )
  }

  all_cats <- get_category_labels(x)
  if (!is.null(categories) && length(categories) > 0L) {
    if (!all(categories %in% all_cats)) {
      .stop(
        "Invalid categories: ",
        paste(setdiff(categories, all_cats), collapse = ", "),
        ". Available categories: ",
        paste(all_cats, collapse = ", ")
      )
    }
  }

  if (length(cues) == 3L) {
    mod_proj <- if (length(obj_cues) > 3L) {
      .make_projected_cognitive_model(
        x, cues,
        decision_rule = decision_rule,
        noise_treatment = noise_treatment,
        lapse_treatment = lapse_treatment
      )
    } else {
      x
    }
    t_sub <- .make_plot_title_and_subtitle(
      mod_proj,
      "categorization",
      cues = cues,
      target_category = if (!is.null(categories)) categories[1L] else NULL,
      decision_rule = decision_rule
    )
    return(.render_3D_sliced_categorization_plot(
      model = mod_proj,
      cues = cues,
      categories = categories,
      aes = aes,
      decision_rule = decision_rule,
      noise_treatment = noise_treatment,
      lapse_treatment = lapse_treatment,
      limits = limits,
      resolution = resolution,
      t_sub = t_sub,
      ...
    ))
  }

  all_cats <- get_category_labels(x)
  if (is.null(categories) || length(categories) == 0L) {
    categories <- all_cats[1L]
  }
  target_cat <- categories[1L]

  is_interactive <- isTRUE(list(...)$interactive) || (!is.null(aes) && "interactive" %in% aes)
  default_aes <- if (length(cues) == 1L) {
    "contour"
  } else if (length(cues) == 2L && is_interactive) {
    "fill-discrete"
  } else {
    "contour"
  }
  aes <- .normalize_plot_aes(aes, default_aes = default_aes)
  if (is_interactive && !any(c("fill", "fill-discrete", "fill-gradient", "contour") %in% aes)) {
    aes <- c(aes, default_aes)
  }
  lvl_spec <- .parse_categorization_levels(levels, aes)
  lim_spec <- .parse_cue_limits(limits, cues)
  t_sub <- .make_plot_title_and_subtitle(
    x,
    "categorization",
    cues = cues,
    target_category = target_cat,
    decision_rule = decision_rule
  )

  cat_colors <- scales::hue_pal()(length(all_cats))
  names(cat_colors) <- all_cats

  if (length(cues) == 1L) {
    mod_proj <- .make_projected_cognitive_model(
      x, cues,
      decision_rule = decision_rule,
      noise_treatment = noise_treatment,
      lapse_treatment = lapse_treatment
    )
    reps_proj <- mod_proj@category_template@representations

    dens_df <- .make_1D_category_density_df(
      reps_proj,
      cues,
      resolution = resolution,
      limits = limits
    )
    grid_vals <- if (!is.null(lim_spec[[cues[1L]]])) {
      seq(
        lim_spec[[cues[1L]]][1L],
        lim_spec[[cues[1L]]][2L],
        length.out = resolution
      )
    } else {
      unique(dens_df[[cues[1L]]])
    }
    test_mat <- matrix(grid_vals, ncol = 1L)
    colnames(test_mat) <- cues[1L]

    pf <- get_category_posterior_function(
      mod_proj,
      noise_treatment = noise_treatment %||% get_noise_treatment(mod_proj),
      lapse_treatment = lapse_treatment %||% get_lapse_treatment(mod_proj)
    )
    post_raw <- pf(test_mat, categories = all_cats)

    eff_l_trt <- lapse_treatment %||% get_lapse_treatment(mod_proj)
    resp_mat <- .apply_decision_rule_to_posteriors(
      post_raw,
      decision_rule = decision_rule,
      lapse_rate = if (identical(eff_l_trt, "no_lapses")) 0 else mod_proj@lapse_behavior$lapse_rate,
      lapse_bias = if (identical(eff_l_trt, "no_lapses")) NULL else mod_proj@lapse_behavior$lapse_bias
    )
    colnames(resp_mat) <- all_cats

    post_mat <- as.data.frame(resp_mat)
    post_mat[[cues[1L]]] <- grid_vals

    long_df <- tidyr::pivot_longer(
      post_mat,
      cols = dplyr::all_of(categories),
      names_to = "Category",
      values_to = "Probability"
    )

    show_leg <- length(categories) > 1L

    p <- ggplot(
      long_df,
      aes(
        x = .data[[cues[1L]]],
        y = .data$Probability,
        color = .data$Category
      )
    ) +
      geom_line(linewidth = 0.5, show.legend = show_leg) +
      scale_color_manual(values = cat_colors, name = "Category") +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
      labs(
        title = t_sub$title,
        subtitle = t_sub$subtitle,
        x = cues[1L],
        y = if (length(categories) == 1L) {
          sprintf("P(%s | %s)", target_cat, cues[1L])
        } else {
          "Probability"
        }
      ) +
      .mvbu_theme()

    if (!is.null(lim_spec[[cues[1L]]])) {
      p <- p + coord_cartesian(
        xlim = lim_spec[[cues[1L]]],
        ylim = c(0, 1),
        expand = FALSE
      )
    } else {
      p <- p + coord_cartesian(ylim = c(0, 1), expand = FALSE)
    }

    return(p)
  }

  if (length(cues) == 2L) {
    mod_proj <- .make_projected_cognitive_model(
      x, cues,
      decision_rule = decision_rule,
      noise_treatment = noise_treatment,
      lapse_treatment = lapse_treatment
    )
    reps_proj <- mod_proj@category_template@representations

    has_ex <- any(sapply(reps_proj, function(r) {
      S7::S7_inherits(r, Exemplar_CategoryRepresentation)
    }))

    if (!is.null(lim_spec[[cues[1L]]]) && !is.null(lim_spec[[cues[2L]]])) {
      gx <- seq(
        lim_spec[[cues[1L]]][1L],
        lim_spec[[cues[1L]]][2L],
        length.out = min(resolution, 60)
      )
      gy <- seq(
        lim_spec[[cues[2L]]][1L],
        lim_spec[[cues[2L]]][2L],
        length.out = min(resolution, 60)
      )
    } else if (has_ex) {
      d_df <- .make_2D_category_density_grid_df(
        reps_proj,
        cues,
        resolution = 60,
        limits = limits
      )
      gx <- unique(d_df[[cues[1L]]])
      gy <- unique(d_df[[cues[2L]]])
    } else {
      el_df <- .make_2D_category_ellipse_df(reps_proj, cues, levels = 0.95)
      gx <- seq(
        min(el_df[[cues[1L]]]),
        max(el_df[[cues[1L]]]),
        length.out = min(resolution, 60)
      )
      gy <- seq(
        min(el_df[[cues[2L]]]),
        max(el_df[[cues[2L]]]),
        length.out = min(resolution, 60)
      )
    }

    test_grid <- expand.grid(x = gx, y = gy)
    names(test_grid) <- cues[1:2]
    test_mat <- as.matrix(test_grid)

    pf <- get_category_posterior_function(
      mod_proj,
      noise_treatment = noise_treatment %||% get_noise_treatment(mod_proj),
      lapse_treatment = lapse_treatment %||% get_lapse_treatment(mod_proj)
    )
    post_raw <- pf(test_mat, categories = all_cats)

    eff_l_trt <- lapse_treatment %||% get_lapse_treatment(mod_proj)
    resp_mat <- .apply_decision_rule_to_posteriors(
      post_raw,
      decision_rule = decision_rule,
      lapse_rate = if (identical(eff_l_trt, "no_lapses")) 0 else mod_proj@lapse_behavior$lapse_rate,
      lapse_bias = if (identical(eff_l_trt, "no_lapses")) NULL else mod_proj@lapse_behavior$lapse_bias
    )
    colnames(resp_mat) <- all_cats

    if (isTRUE(list(...)$interactive) || "interactive" %in% aes) {
      return(.render_2D_interactive_categorization_plot(
        resp_mat = resp_mat,
        gx = gx,
        gy = gy,
        cues = cues,
        categories = categories,
        cat_colors = cat_colors,
        aes = aes,
        t_sub = t_sub,
        ...
      ))
    }

    dfs_cats <- list()
    for (cat_name in categories) {
      df_c <- test_grid
      df_c$Probability <- resp_mat[, cat_name]
      df_c$Category <- cat_name
      dfs_cats[[length(dfs_cats) + 1L]] <- df_c
    }
    df_all_cats <- dplyr::bind_rows(dfs_cats)

    p <- ggplot(
      df_all_cats,
      aes(
        x = .data[[cues[1L]]],
        y = .data[[cues[2L]]]
      )
    ) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(
        title = t_sub$title,
        subtitle = t_sub$subtitle,
        x = cues[1L],
        y = cues[2L]
      ) +
      .mvbu_theme()

    has_fill <- any(c("fill-gradient", "fill-discrete") %in% aes)
    has_contour <- "contour" %in% aes

    if ("fill-discrete" %in% aes) {
      df_first <- df_all_cats[df_all_cats$Category == target_cat, ]
      breaks_all <- c(-Inf, sort(unique(lvl_spec$fill)), Inf)
      n_bands <- length(breaks_all) - 1L
      band_colors <- scales::alpha(
        cat_colors[target_cat],
        seq(0.05, 0.85, length.out = n_bands)
      )

      p <- p + ggplot2::geom_contour_filled(
        data = df_first,
        aes(
          x = .data[[cues[1L]]],
          y = .data[[cues[2L]]],
          z = .data$Probability
        ),
        breaks = breaks_all,
        show.legend = FALSE
      ) +
        ggplot2::scale_fill_manual(values = band_colors, guide = "none")
    } else if ("fill-gradient" %in% aes) {
      p <- p + geom_tile(
        data = df_all_cats,
        aes(
          x = .data[[cues[1L]]],
          y = .data[[cues[2L]]],
          fill = .data$Category,
          alpha = .data$Probability
        ),
        show.legend = c(fill = (length(categories) > 1L), alpha = TRUE)
      ) +
        scale_fill_manual(values = cat_colors, name = "Category") +
        scale_alpha_continuous(
          range = c(0, 0.85),
          limits = c(0, 1),
          name = "Response\nprobability",
          breaks = c(0, 0.25, 0.50, 0.75, 1.0),
          labels = c("0.0", "0.25", "0.50", "0.75", "1.0")
        ) +
        guides(
          fill = guide_legend(order = 1),
          alpha = guide_legend(order = 2, reverse = TRUE)
        )
    }

    if (has_contour) {
      n_dec <- max(
        2L,
        nchar(sub("^[^.]*\\.?", "", as.character(lvl_spec$contour)))
      )
      fmt_str <- paste0("%.", n_dec, "f")

      contour_alpha <- if (length(categories) >= 2L) 0.5 else 1.0

      if (has_fill) {
        df_first <- df_all_cats[df_all_cats$Category == target_cat, ]
        p <- p + geomtextpath::geom_textcontour(
          data = df_first,
          aes(
            x = .data[[cues[1L]]],
            y = .data[[cues[2L]]],
            z = .data$Probability,
            label = after_stat(sprintf(fmt_str, level))
          ),
          breaks = lvl_spec$contour,
          color = "darkgray",
          linewidth = 0.35,
          alpha = contour_alpha,
          size = 2.8,
          vjust = 0.5,
          show.legend = FALSE
        )
      } else {
        show_leg <- length(categories) > 1L
        p <- p + geomtextpath::geom_textcontour(
          data = df_all_cats,
          aes(
            x = .data[[cues[1L]]],
            y = .data[[cues[2L]]],
            z = .data$Probability,
            color = .data$Category,
            label = after_stat(sprintf(fmt_str, level))
          ),
          breaks = lvl_spec$contour,
          linewidth = 0.35,
          alpha = contour_alpha,
          size = 2.8,
          vjust = 0.5,
          show.legend = show_leg
        ) +
          scale_color_manual(values = cat_colors, name = "Category")
      }
    }

    if (!is.null(lim_spec)) {
      p <- p + coord_cartesian(
        xlim = lim_spec[[cues[1L]]],
        ylim = lim_spec[[cues[2L]]],
        expand = FALSE
      )
    } else {
      p <- p + coord_cartesian(expand = FALSE)
    }

    return(p)
  }

  stop("Plotting categorization function is supported for 1, 2, or 3 cues.")
}

#' @rdname plot_categorization_functions
#' @export
S7::method(plot_categorization_functions, MVBU_Stanfit) <- function(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  decision_rule = "proportional",
  resolution = 100,
  noise_treatment = NULL,
  lapse_treatment = NULL,
  ...
) {
  dots <- list(...)
  uncertainty_treatment <- if ("uncertainty_treatment" %in% names(dots)) {
    match.arg(dots$uncertainty_treatment, c("marginalize", "sample"))
  } else {
    "marginalize"
  }
  ndraws <- if ("ndraws" %in% names(dots)) dots$ndraws else 20L
  show_exposure_data <- if ("show_exposure_data" %in% names(dots)) {
    dots$show_exposure_data
  } else {
    FALSE
  }
  show_test_data <- if ("show_test_data" %in% names(dots)) {
    dots$show_test_data
  } else {
    FALSE
  }
  if (is.null(noise_treatment)) {
    noise_treatment <- "none"
  } else {
    noise_treatment <- match.arg(noise_treatment, c("none", "marginalize"))
  }
  if (is.null(lapse_treatment)) {
    lapse_treatment <- "none"
  } else {
    lapse_treatment <- match.arg(lapse_treatment, c("none", "marginalize"))
  }

  obj_cues <- get_cue_labels(x)
  if (is.null(cues)) {
    cues <- if (length(obj_cues) > 2L) obj_cues[1:2] else obj_cues
  }
  if (length(cues) > 3L) {
    .stop(
      "Cannot plot more than 3 cue dimensions simultaneously (received ",
      length(cues), " cues: ", paste(cues, collapse = ", "),
      "). Please specify 1, 2, or 3 cues; non-plotted cues will be analytically marginalized out."
    )
  }
  if (length(cues) == 3L) {
    .stop("Plotting categorization function for MVBU_Stanfit is currently supported for 1 or 2 cues.")
  }
  if (!all(cues %in% obj_cues)) {
    .stop(
      "Invalid cues: ",
      paste(setdiff(cues, obj_cues), collapse = ", "),
      ". Available cues: ",
      paste(obj_cues, collapse = ", ")
    )
  }

  all_cats <- get_category_labels(x)
  if (is.null(categories) || length(categories) == 0L) {
    categories <- all_cats[1L]
  }
  target_cat <- categories[1L]

  if (is.null(groups)) {
    groups <- get_group_labels(x, include_prior = FALSE)
  }

  aes <- .normalize_plot_aes(
    aes,
    default_aes = if (length(cues) == 1L) "contour" else c("fill-gradient", "contour")
  )
  lvl_spec <- .parse_categorization_levels(levels, aes)

  # Build limits
  if (is.null(limits)) {
    limits <- .compute_default_cue_limits(x, cues, n_sds = 3.0)
  }

  # Ensure posterior latents and closures are cached
  c_list <- .get_cache(x)
  if (is.null(c_list$draws) || is.null(c_list$summary)) {
    x <- add_posterior_latents(x)
    c_list <- .get_cache(x)
  }

  if (length(cues) == 1L) {
    gx <- seq(limits[[cues[1L]]][1L], limits[[cues[1L]]][2L], length.out = resolution)
    grid_df <- data.frame(gx)
    names(grid_df) <- cues[1L]

    # Determine draw IDs to evaluate
    draw_params <- c_list$draws$parameters
    draw_col <- if (".draw" %in% names(draw_params)) ".draw" else if ("draw" %in% names(draw_params)) "draw" else ".iteration"
    all_draw_ids <- sort(unique(draw_params[[draw_col]]))
    sample_draws <- if (!is.null(ndraws) && ndraws < length(all_draw_ids)) {
      sample(all_draw_ids, ndraws)
    } else {
      all_draw_ids
    }

    res_list <- list()
    sample_res_list <- list()

    for (grp in groups) {
      post_arr <- c_list$draws$posterior_function(
        grid_df,
        categories = all_cats,
        groups = grp,
        draws = sample_draws,
        noise_treatment = noise_treatment,
        lapse_treatment = lapse_treatment
      )
      # post_arr has dim [N_obs, K, S_draws]
      post_mat <- post_arr[, target_cat, , drop = FALSE]
      dim(post_mat) <- c(nrow(grid_df), length(sample_draws))

      lvl <- if (!is.null(levels) && length(levels) > 0L) levels[1L] else 0.95
      p_lower <- (1 - lvl) / 2
      p_upper <- 1 - p_lower

      mean_vals <- rowMeans(post_mat)
      q_lower_vals <- apply(post_mat, 1L, stats::quantile, probs = p_lower, na.rm = TRUE)
      q_upper_vals <- apply(post_mat, 1L, stats::quantile, probs = p_upper, na.rm = TRUE)

      res_list[[length(res_list) + 1L]] <- data.frame(
        group = grp,
        x = grid_df[[cues[1L]]],
        mean_post = mean_vals,
        q_lower = q_lower_vals,
        q_upper = q_upper_vals,
        stringsAsFactors = FALSE
      )

      if (identical(uncertainty_treatment, "sample")) {
        for (s_idx in seq_along(sample_draws)) {
          sample_res_list[[length(sample_res_list) + 1L]] <- data.frame(
            group = grp,
            .draw = sample_draws[s_idx],
            x = grid_df[[cues[1L]]],
            posterior = post_mat[, s_idx],
            stringsAsFactors = FALSE
          )
        }
      }
    }

    sum_res <- dplyr::bind_rows(res_list)
    names(sum_res)[2L] <- cues[1L]

    p <- ggplot2::ggplot(sum_res, ggplot2::aes(x = .data[[cues[1L]]]))

    if (identical(uncertainty_treatment, "marginalize")) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$q_lower, ymax = .data$q_upper),
        alpha = 0.25,
        fill = scales::hue_pal()(1)[1]
      )
    }

    if (identical(uncertainty_treatment, "sample") && length(sample_res_list) > 0L) {
      samp_res <- dplyr::bind_rows(sample_res_list)
      names(samp_res)[3L] <- cues[1L]
      p <- p + ggplot2::geom_line(
        data = samp_res,
        ggplot2::aes(y = .data$posterior, group = .data$.draw),
        color = scales::hue_pal()(1)[1],
        alpha = 0.20,
        linewidth = 0.35,
        show.legend = FALSE
      )
    }

    p <- p +
      ggplot2::geom_line(ggplot2::aes(y = .data$mean_post), color = scales::hue_pal()(1)[1], linewidth = 0.8) +
      ggplot2::facet_wrap(~group) +
      ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      .mvbu_theme()
  } else if (length(cues) == 2L) {
    gx <- seq(limits[[cues[1L]]][1L], limits[[cues[1L]]][2L], length.out = resolution)
    gy <- seq(limits[[cues[2L]]][1L], limits[[cues[2L]]][2L], length.out = resolution)
    grid_mat <- expand.grid(x = gx, y = gy)
    names(grid_mat) <- cues[1:2]

    draw_params <- c_list$draws$parameters
    draw_col <- if (".draw" %in% names(draw_params)) ".draw" else if ("draw" %in% names(draw_params)) "draw" else ".iteration"
    all_draw_ids <- sort(unique(draw_params[[draw_col]]))
    sample_draws <- if (!is.null(ndraws) && ndraws < length(all_draw_ids)) {
      sample(all_draw_ids, ndraws)
    } else {
      all_draw_ids
    }

    res_list <- list()

    for (grp in groups) {
      post_arr <- c_list$draws$posterior_function(
        grid_mat,
        categories = all_cats,
        groups = grp,
        draws = sample_draws,
        noise_treatment = noise_treatment,
        lapse_treatment = lapse_treatment
      )
      post_mat <- post_arr[, target_cat, , drop = FALSE]
      dim(post_mat) <- c(nrow(grid_mat), length(sample_draws))

      grid_grp <- grid_mat
      grid_grp$group <- grp
      grid_grp$posterior <- rowMeans(post_mat)
      grid_grp$Category <- target_cat
      res_list[[length(res_list) + 1L]] <- grid_grp
    }
    all_res <- dplyr::bind_rows(res_list)

    p <- ggplot2::ggplot(all_res, ggplot2::aes(x = .data[[cues[1L]]], y = .data[[cues[2L]]])) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::coord_cartesian(expand = FALSE)

    if ("fill-gradient" %in% aes) {
      p <- p + ggplot2::geom_tile(
        data = all_res,
        ggplot2::aes(fill = .data$posterior),
        show.legend = TRUE
      ) +
        ggplot2::scale_fill_gradient(
          low = "white",
          high = scales::hue_pal()(1)[1],
          limits = c(0, 1),
          name = sprintf("P(%s)", target_cat)
        )
    }

    if ("contour" %in% aes) {
      contour_alpha <- if (length(all_cats) >= 2L) 0.5 else 1.0
      if (requireNamespace("geomtextpath", quietly = TRUE)) {
        p <- p + geomtextpath::geom_textcontour(
          data = all_res,
          ggplot2::aes(z = .data$posterior),
          breaks = c(0.25, 0.5, 0.75),
          size = 3,
          linewidth = 0.5,
          color = "gray30",
          alpha = contour_alpha
        )
      } else {
        p <- p + ggplot2::geom_contour(
          data = all_res,
          ggplot2::aes(z = .data$posterior),
          breaks = c(0.25, 0.5, 0.75),
          linewidth = 0.5,
          color = "gray30",
          alpha = contour_alpha
        )
      }
    }

    p <- p +
      ggplot2::facet_wrap(~group) +
      .mvbu_theme()
  }

  t_sub <- .make_plot_title_and_subtitle(
    x,
    "categorization",
    cues = cues,
    ndraws = if (identical(uncertainty_treatment, "sample")) ndraws else NULL,
    target_category = target_cat
  )

  p <- p + ggplot2::labs(
    title = t_sub$title,
    subtitle = t_sub$subtitle,
    x = cues[1L],
    y = if (length(cues) == 1L) sprintf("P(%s)", target_cat) else cues[2L]
  )

  p
}

# -----------------------------------------------------------------------------
# plot_parameters methods
# -----------------------------------------------------------------------------

#' @rdname plot_parameters
#' @export
S7::method(plot_parameters, MVBU_CognitiveModel) <- function(
  x,
  pars = NULL,
  categories = NULL,
  groups = NULL,
  ndraws = 100,
  combine_into_single_plot = TRUE,
  index_panels = FALSE,
  ...
) {
  cues <- get_cue_labels(x)
  cats <- get_category_labels(x)
  reps <- x@category_template@representations
  t_sub <- .make_plot_title_and_subtitle(x, "parameters")

  # 1. Category location (mu)
  loc_rows <- list()
  for (cat_name in names(reps)) {
    r <- reps[[cat_name]]
    mc <- .get_rep_mean_and_cov(r)
    mu_v <- mc$mu
    for (k in seq_along(cues)) {
      loc_rows[[length(loc_rows) + 1]] <- data.frame(
        Category = cat_name,
        Cue = cues[k],
        Parameter = "\u03bc",
        Value = mu_v[k],
        stringsAsFactors = FALSE
      )
    }
  }
  loc_df <- dplyr::bind_rows(loc_rows)
  if (nrow(loc_df) > 0L) {
    loc_df$Cue <- factor(loc_df$Cue, levels = cues)
    loc_df$Category <- factor(loc_df$Category, levels = cats)
  }

  # 2. Category scale (tau / SD and model noise)
  scale_rows <- list()
  for (cat_name in names(reps)) {
    r <- reps[[cat_name]]
    mc <- .get_rep_mean_and_cov(r)
    sd_v <- sqrt(pmax(diag(as.matrix(mc$Sigma)), 1e-6))
    for (k in seq_along(cues)) {
      scale_rows[[length(scale_rows) + 1]] <- data.frame(
        Category = cat_name,
        Cue = cues[k],
        Parameter = "\u03c4",
        Value = sd_v[k],
        stringsAsFactors = FALSE
      )
    }
  }
  scale_df <- dplyr::bind_rows(scale_rows)
  if (nrow(scale_df) > 0L) {
    scale_df$Cue <- factor(scale_df$Cue, levels = cues)
    scale_df$Category <- factor(scale_df$Category, levels = cats)
  }

  noise_beh <- tryCatch(x@noise_behavior, error = function(e) list())
  sigma_n <- noise_beh$Sigma_noise
  noise_rows <- list()
  if (!is.null(sigma_n) &&
    !identical(noise_beh$noise_treatment, "no_noise") &&
    !identical(noise_beh$noise_treatment, "none")) {
    sd_n <- sqrt(pmax(diag(as.matrix(sigma_n)), 1e-6))
    for (k in seq_along(cues)) {
      noise_rows[[length(noise_rows) + 1]] <- data.frame(
        Cue = cues[k],
        Value = sd_n[k],
        stringsAsFactors = FALSE
      )
    }
  }
  noise_df <- dplyr::bind_rows(noise_rows)

  # 3. Category cue correlations (rho)
  corr_rows <- list()
  if (length(cues) >= 2L) {
    for (cat_name in names(reps)) {
      r <- reps[[cat_name]]
      mc <- .get_rep_mean_and_cov(r)
      c_mat <- stats::cov2cor(as.matrix(mc$Sigma))
      for (i in 1:(length(cues) - 1L)) {
        for (j in (i + 1L):length(cues)) {
          corr_rows[[length(corr_rows) + 1]] <- data.frame(
            Category = cat_name,
            Pair = paste0(cues[i], "-", cues[j]),
            Parameter = "\u03c1",
            Value = c_mat[i, j],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  corr_df <- dplyr::bind_rows(corr_rows)
  if (nrow(corr_df) > 0L) {
    corr_df$Category <- factor(corr_df$Category, levels = cats)
  }

  # 4. Confidence / sample size (kappa, nu)
  counts_rows <- list()
  for (cat_name in names(reps)) {
    r <- reps[[cat_name]]
    if ("kappa" %in% names(r)) {
      counts_rows[[length(counts_rows) + 1]] <- data.frame(
        Category = cat_name,
        Parameter = "\u03ba (location)",
        Value = r@kappa,
        stringsAsFactors = FALSE
      )
    }
    if ("nu" %in% names(r)) {
      counts_rows[[length(counts_rows) + 1]] <- data.frame(
        Category = cat_name,
        Parameter = "\u03bd (scale)",
        Value = r@nu,
        stringsAsFactors = FALSE
      )
    }
    if ("exemplars" %in% names(r)) {
      counts_rows[[length(counts_rows) + 1]] <- data.frame(
        Category = cat_name,
        Parameter = "Exemplars",
        Value = nrow(as.matrix(r@exemplars)),
        stringsAsFactors = FALSE
      )
    }
  }
  counts_df <- dplyr::bind_rows(counts_rows)
  if (nrow(counts_df) > 0L) {
    counts_df$Category <- factor(counts_df$Category, levels = cats)
  }

  # 5. Decision-making (priors, lapse bias, lapse rate)
  probs_rows <- list()
  cp <- x@category_prior
  if (!is.null(cp)) {
    for (cat_name in names(cp)) {
      probs_rows[[length(probs_rows) + 1]] <- data.frame(
        Category = cat_name,
        Parameter = "Prior P(C)",
        Value = cp[cat_name],
        stringsAsFactors = FALSE
      )
    }
  }
  lapse_beh <- tryCatch(x@lapse_behavior, error = function(e) list())
  lb <- lapse_beh$lapse_bias
  if (!is.null(lb)) {
    for (cat_name in names(lb)) {
      probs_rows[[length(probs_rows) + 1]] <- data.frame(
        Category = cat_name,
        Parameter = "Lapse bias",
        Value = lb[cat_name],
        stringsAsFactors = FALSE
      )
    }
  }
  probs_df <- dplyr::bind_rows(probs_rows)
  if (nrow(probs_df) > 0L) {
    probs_df$Category <- factor(probs_df$Category, levels = cats)
  }

  l_rate <- lapse_beh$lapse_rate
  lapse_row <- NULL
  if (!is.null(l_rate) && l_rate > 0) {
    lapse_row <- data.frame(
      Parameter = "Lapse rate (\u03bb)",
      Value = l_rate,
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(pars)) {
    p_pat <- paste(pars, collapse = "|")
    if (nrow(loc_df) > 0L) {
      loc_df <- loc_df[grepl(p_pat, loc_df$Parameter), ]
    }
    if (nrow(scale_df) > 0L) {
      scale_df <- scale_df[grepl(p_pat, scale_df$Parameter), ]
    }
    if (nrow(corr_df) > 0L) {
      corr_df <- corr_df[grepl(p_pat, corr_df$Parameter), ]
    }
    if (nrow(counts_df) > 0L) {
      counts_df <- counts_df[grepl(p_pat, counts_df$Parameter), ]
    }
    if (nrow(probs_df) > 0L) {
      probs_df <- probs_df[grepl(p_pat, probs_df$Parameter), ]
    }
  }

  panels <- list()

  # Panel 1: Location (keep only color, no fill legend)
  if (nrow(loc_df) > 0L) {
    panels$location <- ggplot(
      loc_df,
      aes(
        x = .data$Cue,
        y = .data$Value,
        color = .data$Category
      )
    ) +
      geom_point(
        position = position_dodge(width = 0.4),
        size = 3
      ) +
      guides(fill = "none") +
      labs(title = "Category location", x = NULL, y = "Location (\u03bc)") +
      .mvbu_theme()
  }

  # Panel 2: Scale (keep only color, no fill legend)
  if (nrow(scale_df) > 0L) {
    p_sc <- ggplot() +
      geom_point(
        data = scale_df,
        aes(
          x = .data$Cue,
          y = .data$Value,
          color = .data$Category
        ),
        position = position_dodge(width = 0.4),
        size = 3
      ) +
      guides(fill = "none") +
      scale_y_log10() +
      labs(title = "Category scale", x = NULL, y = "Scale (\u03c4)") +
      .mvbu_theme()

    if (nrow(noise_df) > 0L) {
      p_sc <- p_sc + geom_point(
        data = noise_df,
        aes(x = .data$Cue, y = .data$Value),
        color = "darkgray",
        shape = 17,
        size = 3,
        show.legend = FALSE
      )
    }
    panels$scale <- p_sc
  }

  # Panel 3: Correlations
  if (nrow(corr_df) > 0L) {
    panels$correlations <- ggplot(
      corr_df,
      aes(
        x = .data$Pair,
        y = .data$Value,
        color = .data$Category
      )
    ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      geom_point(
        position = position_dodge(width = 0.4),
        size = 3
      ) +
      guides(fill = "none") +
      ylim(-1, 1) +
      labs(
        title = "Category cue correlations",
        x = NULL,
        y = "Cue correlation (\u03c1)"
      ) +
      .mvbu_theme()
  }

  has_pts <- (nrow(loc_df) > 0L || nrow(scale_df) > 0L || nrow(corr_df) > 0L)
  col_leg <- !has_pts

  # Panel 4: Confidence (no fill legend, color-matched columns)
  if (nrow(counts_df) > 0L) {
    panels$confidence <- ggplot(
      counts_df,
      aes(
        x = .data$Parameter,
        y = .data$Value,
        color = .data$Category,
        fill = .data$Category
      )
    ) +
      geom_col(
        position = position_dodge(width = 0.7),
        width = 0.6,
        show.legend = col_leg
      ) +
      scale_fill_discrete(guide = "none") +
      guides(fill = "none") +
      scale_y_log10() +
      labs(
        title = "Category confidence",
        x = NULL,
        y = "Count"
      ) +
      .mvbu_theme()
  }

  # Panel 5: Decision-making (no fill legend)
  if (nrow(probs_df) > 0L || !is.null(lapse_row)) {
    p_dec <- ggplot() +
      ylim(0, 1) +
      labs(
        title = "Decision-making",
        x = NULL,
        y = "Probability"
      ) +
      .mvbu_theme() +
      theme(
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
      )

    if (nrow(probs_df) > 0L) {
      p_dec <- p_dec + geom_col(
        data = probs_df,
        aes(
          x = .data$Parameter,
          y = .data$Value,
          color = .data$Category,
          fill = .data$Category
        ),
        position = position_dodge(width = 0.7),
        width = 0.6,
        show.legend = col_leg
      ) +
        scale_fill_discrete(guide = "none") +
        guides(fill = "none")
    }

    if (!is.null(lapse_row)) {
      p_dec <- p_dec + geom_col(
        data = lapse_row,
        aes(x = .data$Parameter, y = .data$Value),
        fill = "darkgray",
        color = "darkgray",
        width = 0.4,
        show.legend = FALSE
      )
    }
    panels$decision <- p_dec
  }

  if (length(panels) == 0L) {
    stop("No matching parameters found to plot.")
  }

  if (!isTRUE(combine_into_single_plot)) {
    return(panels)
  }

  patchwork::wrap_plots(
    panels,
    ncol = if (length(panels) > 2L) 2L else 1L
  ) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = t_sub$title,
      subtitle = t_sub$subtitle
    ) &
    guides(fill = "none") &
    theme(legend.position = "right")
}

#' @rdname plot_parameters
#' @export
S7::method(plot_parameters, MVBU_Stanfit) <- function(
  x,
  pars = NULL,
  categories = NULL,
  groups = NULL,
  ndraws = 100,
  combine_into_single_plot = TRUE,
  index_panels = FALSE,
  ...
) {
  avail_cats <- get_category_labels(x)
  if (is.null(categories)) {
    categories <- avail_cats
  }
  avail_grps <- get_group_labels(x, include_prior = TRUE)
  if (is.null(groups)) {
    groups <- avail_grps
  }

  draws <- get_draws(
    x,
    categories = categories,
    groups = groups,
    ndraws = ndraws,
    summarize = FALSE,
    ...
  )
  if (is.null(draws) || nrow(draws) == 0L) {
    stop("No draws found in stanfit object.")
  }

  cues <- get_cue_labels(x)
  t_sub <- .make_plot_title_and_subtitle(
    x,
    "parameters",
    ndraws = ndraws
  )

  if (!"tau" %in% names(draws) && "S" %in% names(draws)) {
    draws$tau <- purrr::map(draws$S, function(s) {
      mat <- as.matrix(s)
      sqrt(pmax(diag(mat), 1e-6))
    })
  }

  panels <- list()

  # 1. Location (Means) with parsed facet expressions: mu[cue]
  if ("m" %in% names(draws)) {
    m_df_list <- list()
    for (i in seq_len(nrow(draws))) {
      m_vec <- draws$m[[i]]
      for (k in seq_along(cues)) {
        m_df_list[[length(m_df_list) + 1]] <- data.frame(
          Group = draws$group[i],
          Category = draws$category[i],
          Cue = sprintf("mu[%s]", cues[k]),
          Mean = m_vec[k],
          stringsAsFactors = FALSE
        )
      }
    }
    m_df <- dplyr::bind_rows(m_df_list)
    m_df$Cue <- factor(
      m_df$Cue,
      levels = sprintf("mu[%s]", cues)
    )
    m_df$Category <- factor(m_df$Category, levels = categories)

    cat_cols <- scales::hue_pal()(length(categories))
    names(cat_cols) <- categories

    panels$location <- ggplot(
      m_df,
      aes(
        x = .data$Mean,
        fill = .data$Category,
        color = .data$Category
      )
    ) +
      geom_density(
        alpha = 0.5,
        linewidth = 0.30
      ) +
      scale_fill_manual(values = cat_cols, name = "Category") +
      scale_color_manual(values = cat_cols, name = "Category") +
      guides(
        color = guide_legend(
          title = "Category",
          override.aes = list(
            fill = cat_cols,
            alpha = 0.5
          )
        ),
        fill = "none"
      ) +
      facet_grid(Group ~ Cue, scales = "free", labeller = label_parsed) +
      scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
      scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
      labs(
        title = "Category location",
        x = "Location (\u03bc)",
        y = "Density"
      ) +
      .mvbu_theme() +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = 8),
        axis.text.x = ggplot2::element_text(size = 8)
      )
  }

  # 2. Scale (Tau) with parsed facet expressions: tau[cue]
  if ("tau" %in% names(draws)) {
    tau_df_list <- list()
    for (i in seq_len(nrow(draws))) {
      tau_vec <- draws$tau[[i]]
      for (k in seq_along(cues)) {
        tau_df_list[[length(tau_df_list) + 1]] <- data.frame(
          Group = draws$group[i],
          Category = draws$category[i],
          Cue = sprintf("tau[%s]", cues[k]),
          Tau = tau_vec[k],
          stringsAsFactors = FALSE
        )
      }
    }
    tau_df <- dplyr::bind_rows(tau_df_list)
    tau_df$Cue <- factor(
      tau_df$Cue,
      levels = sprintf("tau[%s]", cues)
    )
    tau_df$Category <- factor(tau_df$Category, levels = categories)

    panels$scale <- ggplot(
      tau_df,
      aes(
        x = .data$Tau,
        fill = .data$Category,
        color = .data$Category
      )
    ) +
      geom_density(
        alpha = 0.5,
        linewidth = 0.30
      ) +
      scale_fill_manual(values = cat_cols, name = "Category") +
      scale_color_manual(values = cat_cols, name = "Category") +
      guides(
        color = guide_legend(
          title = "Category",
          override.aes = list(
            fill = cat_cols,
            alpha = 0.5
          )
        ),
        fill = "none"
      ) +
      scale_x_log10(breaks = scales::breaks_log(n = 3)) +
      scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
      facet_grid(Group ~ Cue, scales = "free", labeller = label_parsed) +
      labs(
        title = "Category scale",
        x = "Scale (\u03c4)",
        y = "Density"
      ) +
      .mvbu_theme() +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = 8),
        axis.text.x = ggplot2::element_text(size = 8)
      )
  }

  # 3. Confidence (Kappa, Nu) with parsed expressions: kappa, nu
  if ("kappa" %in% names(draws) || "nu" %in% names(draws)) {
    count_df <- draws %>%
      dplyr::select(dplyr::any_of(c("group", "category", "kappa", "nu"))) %>%
      tidyr::pivot_longer(
        cols = dplyr::any_of(c("kappa", "nu")),
        names_to = "Parameter",
        values_to = "Value"
      )
    count_df$Category <- factor(count_df$category, levels = categories)
    count_df$Group <- count_df$group
    count_df$Parameter <- factor(
      count_df$Parameter,
      levels = c("kappa", "nu")
    )

    panels$confidence <- ggplot(
      count_df,
      aes(
        x = .data$Value,
        fill = .data$Category,
        color = .data$Category
      )
    ) +
      geom_density(
        alpha = 0.5,
        linewidth = 0.30
      ) +
      scale_fill_manual(values = cat_cols, name = "Category") +
      scale_color_manual(values = cat_cols, name = "Category") +
      guides(
        color = guide_legend(
          title = "Category",
          override.aes = list(
            fill = cat_cols,
            alpha = 0.5
          )
        ),
        fill = "none"
      ) +
      scale_x_log10(breaks = scales::breaks_log(n = 3)) +
      scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
      facet_grid(Group ~ Parameter, scales = "free", labeller = label_parsed) +
      labs(
        title = "Category confidence",
        x = "Count",
        y = "Density"
      ) +
      .mvbu_theme() +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = 8),
        axis.text.x = ggplot2::element_text(size = 8)
      )
  }

  # 4. Decision-making (Lapse rate: lambda)
  if ("lapse_rate" %in% names(draws)) {
    panels$decision <- ggplot(
      draws,
      aes(x = .data$lapse_rate)
    ) +
      geom_density(
        fill = "darkgray",
        color = "darkgray",
        alpha = 0.5,
        linewidth = 0.30
      ) +
      scale_fill_discrete(guide = "none") +
      guides(fill = "none", color = "none") +
      scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
      scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
      facet_wrap(~group) +
      labs(
        title = "Decision-making",
        x = "Lapse rate (\u03bb)",
        y = "Density"
      ) +
      .mvbu_theme() +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = 8),
        axis.text.x = ggplot2::element_text(size = 8)
      )
  }

  if (length(panels) == 0L) {
    stop("No parameter distributions found to plot.")
  }

  if (!isTRUE(combine_into_single_plot)) {
    return(panels)
  }

  patchwork::wrap_plots(
    panels,
    ncol = if (length(panels) > 2L) 2L else 1L
  ) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = t_sub$title,
      subtitle = t_sub$subtitle
    ) &
    theme(legend.position = "right")
}

# -----------------------------------------------------------------------------
# plot_parameter_correlations methods
# -----------------------------------------------------------------------------

#' @rdname plot_parameters
#' @export
S7::method(plot_parameter_correlations, MVBU_Stanfit) <- function(
  x,
  pars = NULL,
  categories = NULL,
  groups = NULL,
  cues = NULL,
  ndraws = 100,
  ...
) {
  obj_cues <- get_cue_labels(x)
  if (is.null(cues)) {
    cues <- obj_cues
  }
  avail_cats <- get_category_labels(x)
  categories <- .validate_requested_labels(
    categories,
    avail_cats,
    label_type = "category"
  )
  avail_grps <- get_group_labels(x, include_prior = TRUE)
  if (is.null(groups)) {
    groups <- if ("prior" %in% avail_grps) "prior" else avail_grps[1L]
  }
  groups <- .validate_requested_labels(
    groups,
    avail_grps,
    label_type = "group"
  )

  draws <- get_draws(
    x,
    groups = groups,
    categories = categories,
    ndraws = ndraws,
    summarize = FALSE,
    ...
  )

  if (is.null(draws) || nrow(draws) == 0L) {
    stop("No draws found to compute parameter correlations.")
  }

  unique_draws <- unique(draws$.draw)
  df_mat <- list()

  for (d in unique_draws) {
    sub_d <- draws[draws$.draw == d, ]
    row_vals <- list()
    for (c_name in categories) {
      c_clean <- gsub("[^A-Za-z0-9]", "", as.character(c_name))
      for (grp in groups) {
        grp_sub_d <- sub_d[
          sub_d$category == c_name & sub_d$group == grp,
        ]
        if (nrow(grp_sub_d) > 0L) {
          if ("kappa" %in% names(grp_sub_d)) {
            row_vals[[paste0("kappa_", c_clean, "__", grp)]] <-
              grp_sub_d$kappa[1L]
          }
          if ("nu" %in% names(grp_sub_d)) {
            row_vals[[paste0("nu_", c_clean, "__", grp)]] <-
              grp_sub_d$nu[1L]
          }
          if ("m" %in% names(grp_sub_d)) {
            m_v <- grp_sub_d$m[[1L]]
            for (k in seq_along(cues)) {
              row_vals[[paste0("m_", c_clean, "_", cues[k], "__", grp)]] <-
                m_v[k]
            }
          }
          if ("S" %in% names(grp_sub_d)) {
            s_mat <- as.matrix(grp_sub_d$S[[1L]])
            sd_v <- sqrt(pmax(diag(s_mat), 1e-6))
            for (k in seq_along(cues)) {
              row_vals[[paste0("tau_", c_clean, "_", cues[k], "__", grp)]] <-
                sd_v[k]
            }
            if (length(cues) >= 2L) {
              c_mat_d <- stats::cov2cor(s_mat)
              for (ci in 1:(length(cues) - 1L)) {
                for (cj in (ci + 1L):length(cues)) {
                  row_vals[[paste0(
                    "rho_",
                    c_clean,
                    "_",
                    cues[ci],
                    "_",
                    cues[cj],
                    "__",
                    grp
                  )]] <- c_mat_d[ci, cj]
                }
              }
            }
          }
        }
      }
    }
    if ("lapse_rate" %in% names(sub_d)) {
      for (grp in groups) {
        grp_sub_d <- sub_d[sub_d$group == grp, ]
        if (nrow(grp_sub_d) > 0L) {
          row_vals[[paste0("lapse_rate__", grp)]] <-
            grp_sub_d$lapse_rate[1L]
        }
      }
    }
    df_mat[[length(df_mat) + 1L]] <- as.data.frame(
      row_vals,
      stringsAsFactors = FALSE
    )
  }

  draws_subset <- dplyr::bind_rows(df_mat)
  param_cols <- names(draws_subset)

  if (!is.null(pars)) {
    param_cols <- param_cols[grepl(paste(pars, collapse = "|"), param_cols)]
    draws_subset <- draws_subset[, param_cols, drop = FALSE]
  }

  if (length(param_cols) < 2L) {
    stop("At least two parameter columns required for correlation plot.")
  }

  c_mat <- stats::cor(draws_subset, use = "pairwise.complete.obs")

  corr_df <- as.data.frame(as.table(c_mat))
  names(corr_df) <- c("Var1", "Var2", "Correlation")

  expr_labels <- .format_param_plotmath(param_cols)
  names(expr_labels) <- param_cols

  corr_df$Var1_expr <- expr_labels[as.character(corr_df$Var1)]
  corr_df$Var2_expr <- expr_labels[as.character(corr_df$Var2)]

  corr_df$Var1_expr <- factor(corr_df$Var1_expr, levels = expr_labels)
  corr_df$Var2_expr <- factor(corr_df$Var2_expr, levels = rev(expr_labels))

  t_sub <- .make_plot_title_and_subtitle(
    x,
    "parameter_correlations",
    ndraws = ndraws
  )

  # Category diagonal outlines
  cat_boxes <- list()
  N_cols <- length(param_cols)
  for (cat_name in categories) {
    c_clean <- gsub("[^A-Za-z0-9]", "", as.character(cat_name))
    cat_idx <- which(
      grepl(paste0("_", c_clean, "(_|__)"), param_cols) |
        grepl(paste0("^", c_clean, "_"), param_cols)
    )
    if (length(cat_idx) > 0L) {
      min_idx <- min(cat_idx)
      max_idx <- max(cat_idx)
      cat_boxes[[length(cat_boxes) + 1L]] <- data.frame(
        xmin = min_idx - 0.5,
        xmax = max_idx + 0.5,
        ymin = (N_cols - max_idx + 1L) - 0.5,
        ymax = (N_cols - min_idx + 1L) + 0.5,
        stringsAsFactors = FALSE
      )
    }
  }

  p <- ggplot(
    corr_df,
    aes(x = .data$Var1_expr, y = .data$Var2_expr, fill = .data$Correlation)
  ) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limit = c(-1, 1),
      name = "Pearson \u03c1"
    ) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_discrete(labels = scales::parse_format()) +
    labs(
      title = t_sub$title,
      subtitle = t_sub$subtitle,
      x = NULL,
      y = NULL
    ) +
    .mvbu_theme() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )

  if (length(cat_boxes) > 0L) {
    cat_boxes_df <- dplyr::bind_rows(cat_boxes)
    p <- p + geom_rect(
      data = cat_boxes_df,
      aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin = .data$ymin,
        ymax = .data$ymax
      ),
      fill = NA,
      color = "gray50",
      linewidth = 0.8,
      inherit.aes = FALSE
    )
  }

  p
}

# -----------------------------------------------------------------------------
# plot_parameters_pairwise methods
# -----------------------------------------------------------------------------

#' @rdname plot_parameters
#' @export
S7::method(plot_parameters_pairwise, MVBU_Stanfit) <- function(
  x,
  pars = NULL,
  categories = NULL,
  groups = NULL,
  cues = NULL,
  ndraws = 100,
  ...
) {
  obj_cues <- get_cue_labels(x)
  if (is.null(cues)) {
    cues <- obj_cues
  }
  avail_cats <- get_category_labels(x)
  categories <- .validate_requested_labels(
    categories,
    avail_cats,
    label_type = "category"
  )
  avail_grps <- get_group_labels(x, include_prior = TRUE)
  if (is.null(groups)) {
    groups <- if ("prior" %in% avail_grps) "prior" else avail_grps[1L]
  }
  groups <- .validate_requested_labels(
    groups,
    avail_grps,
    label_type = "group"
  )

  draws <- get_draws(
    x,
    groups = groups,
    categories = categories,
    ndraws = ndraws,
    summarize = FALSE,
    ...
  )

  if (is.null(draws) || nrow(draws) == 0L) {
    stop("No draws found to compute parameter distributions.")
  }

  # Build one row per (group, draw) containing all category parameters
  unique_grps <- unique(draws$group)
  df_mat <- list()

  for (grp in unique_grps) {
    grp_draws <- draws[draws$group == grp, ]
    unique_d <- unique(grp_draws$.draw)
    for (d in unique_d) {
      sub_d <- grp_draws[grp_draws$.draw == d, ]
      row_vals <- list(Group = as.character(grp))
      for (i in seq_len(nrow(sub_d))) {
        c_name <- as.character(sub_d$category[i])
        c_clean <- gsub("[^A-Za-z0-9]", "", c_name)
        if ("kappa" %in% names(sub_d)) {
          row_vals[[paste0("kappa_", c_clean)]] <- sub_d$kappa[i]
        }
        if ("nu" %in% names(sub_d)) {
          row_vals[[paste0("nu_", c_clean)]] <- sub_d$nu[i]
        }
        if ("m" %in% names(sub_d)) {
          m_v <- sub_d$m[[i]]
          for (k in seq_along(cues)) {
            row_vals[[paste0("m_", c_clean, "_", cues[k])]] <- m_v[k]
          }
        }
        if ("S" %in% names(sub_d)) {
          s_mat <- as.matrix(sub_d$S[[i]])
          sd_v <- sqrt(pmax(diag(s_mat), 1e-6))
          for (k in seq_along(cues)) {
            row_vals[[paste0("tau_", c_clean, "_", cues[k])]] <- sd_v[k]
          }
          if (length(cues) >= 2L) {
            c_mat_d <- stats::cov2cor(s_mat)
            for (ci in 1:(length(cues) - 1L)) {
              for (cj in (ci + 1L):length(cues)) {
                row_vals[[paste0(
                  "rho_",
                  c_clean,
                  "_",
                  cues[ci],
                  "_",
                  cues[cj]
                )]] <- c_mat_d[ci, cj]
              }
            }
          }
        }
      }
      if ("lapse_rate" %in% names(sub_d)) {
        row_vals[["lapse_rate"]] <- sub_d$lapse_rate[1L]
      }
      df_mat[[length(df_mat) + 1L]] <- as.data.frame(
        row_vals,
        stringsAsFactors = FALSE
      )
    }
  }

  draws_subset <- dplyr::bind_rows(df_mat)
  param_cols <- setdiff(names(draws_subset), "Group")

  if (!is.null(pars)) {
    param_cols <- param_cols[grepl(paste(pars, collapse = "|"), param_cols)]
    draws_subset <- draws_subset[, c("Group", param_cols), drop = FALSE]
  }

  if (length(param_cols) < 2L) {
    stop("At least two parameter columns required for pairwise distributions.")
  }

  expr_labels <- .format_param_plotmath(param_cols)
  names(expr_labels) <- param_cols

  # Build pairwise data frame
  # Rows: Var2 (i), Cols: Var1 (j)
  lower_list <- list()
  upper_list <- list()
  for (i in seq_along(param_cols)) {
    for (j in seq_along(param_cols)) {
      p1 <- param_cols[j]
      p2 <- param_cols[i]
      if (i > j) {
        lower_list[[length(lower_list) + 1L]] <- data.frame(
          Var1_expr = unname(expr_labels[p1]),
          Var2_expr = unname(expr_labels[p2]),
          Val1 = draws_subset[[p1]],
          Val2 = draws_subset[[p2]],
          Group = draws_subset$Group,
          stringsAsFactors = FALSE
        )
      } else if (i < j) {
        upper_list[[length(upper_list) + 1L]] <- data.frame(
          Var1_expr = unname(expr_labels[p1]),
          Var2_expr = unname(expr_labels[p2]),
          Val1 = draws_subset[[param_cols[j]]],
          Val2 = draws_subset[[param_cols[i]]],
          Group = draws_subset$Group,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  lower_df <- dplyr::bind_rows(lower_list)
  upper_df <- dplyr::bind_rows(upper_list)

  # Filter groups with positive variance for 2D density contours
  upper_df_valid <- upper_df %>%
    dplyr::group_by(.data$Var1_expr, .data$Var2_expr, .data$Group) %>%
    dplyr::filter(
      stats::sd(.data$Val1) > 1e-6 & stats::sd(.data$Val2) > 1e-6
    ) %>%
    dplyr::ungroup()

  # Diag df with scaled densities matching each panel y range
  diag_list <- list()
  for (i in seq_along(param_cols)) {
    p_name <- param_cols[i]
    p_expr <- unname(expr_labels[p_name])
    sub_p <- draws_subset[, c("Group", p_name)]
    names(sub_p)[2L] <- "Val"
    p_min <- min(sub_p$Val, na.rm = TRUE)
    p_max <- max(sub_p$Val, na.rm = TRUE)
    p_rng <- if (p_max > p_min) p_max - p_min else 1
    for (grp in unique(sub_p$Group)) {
      v <- sub_p$Val[sub_p$Group == grp]
      v <- v[!is.na(v)]
      if (length(v) >= 2L && stats::sd(v) > 1e-9) {
        d_est <- stats::density(v, n = 100)
        scaled_y <- p_min + (d_est$y / max(d_est$y, 1e-12)) * (0.85 * p_rng)
        poly_x <- c(d_est$x[1L], d_est$x, d_est$x[length(d_est$x)])
        poly_y <- c(p_min, scaled_y, p_min)
        diag_list[[length(diag_list) + 1L]] <- data.frame(
          Var1_expr = p_expr,
          Var2_expr = p_expr,
          Val1 = poly_x,
          Val2 = poly_y,
          Group = grp,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  diag_df <- dplyr::bind_rows(diag_list)

  if (nrow(lower_df) > 0L) {
    lower_df$Var1_expr <- factor(lower_df$Var1_expr, levels = expr_labels)
    lower_df$Var2_expr <- factor(lower_df$Var2_expr, levels = expr_labels)
    lower_df$Group <- factor(
      lower_df$Group,
      levels = unique(draws_subset$Group)
    )
  }
  if (nrow(upper_df_valid) > 0L) {
    upper_df_valid$Var1_expr <- factor(
      upper_df_valid$Var1_expr,
      levels = expr_labels
    )
    upper_df_valid$Var2_expr <- factor(
      upper_df_valid$Var2_expr,
      levels = expr_labels
    )
    upper_df_valid$Group <- factor(
      upper_df_valid$Group,
      levels = unique(draws_subset$Group)
    )
  }
  if (nrow(diag_df) > 0L) {
    diag_df$Var1_expr <- factor(diag_df$Var1_expr, levels = expr_labels)
    diag_df$Var2_expr <- factor(diag_df$Var2_expr, levels = expr_labels)
    diag_df$Group <- factor(diag_df$Group, levels = unique(draws_subset$Group))
  }

  # Category diagonal bounding box outlines
  cat_borders_list <- list()
  for (c_name in categories) {
    c_clean <- gsub("[^A-Za-z0-9]", "", c_name)
    c_pars <- param_cols[
      grepl(paste0("_(", c_clean, ")_"), param_cols) |
        grepl(paste0("_(", c_clean, ")$"), param_cols) |
        grepl(paste0("^(kappa|nu)_", c_clean), param_cols) |
        grepl(paste0("^(m|tau|rho)_", c_clean, "_"), param_cols)
    ]
    if (length(c_pars) > 0L) {
      for (p1 in c_pars) {
        for (p2 in c_pars) {
          cat_borders_list[[length(cat_borders_list) + 1L]] <- data.frame(
            Var1_expr = unname(expr_labels[p1]),
            Var2_expr = unname(expr_labels[p2]),
            xmin = -Inf,
            xmax = Inf,
            ymin = -Inf,
            ymax = Inf,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  if (length(cat_borders_list) > 0L) {
    cat_borders_df <- dplyr::bind_rows(cat_borders_list)
    cat_borders_df$Var1_expr <- factor(
      cat_borders_df$Var1_expr,
      levels = expr_labels
    )
    cat_borders_df$Var2_expr <- factor(
      cat_borders_df$Var2_expr,
      levels = expr_labels
    )
  } else {
    cat_borders_df <- NULL
  }

  t_sub <- .make_plot_title_and_subtitle(
    x,
    "parameters_pairwise",
    ndraws = ndraws
  )

  show_grp_legend <- length(unique(draws_subset$Group)) > 1L

  grp_levels <- unique(draws_subset$Group)
  n_non_prior <- sum(grp_levels != "prior")
  pal_colors <- if (n_non_prior > 0L) {
    scales::hue_pal()(n_non_prior)
  } else {
    character(0)
  }
  grp_colors <- character(length(grp_levels))
  names(grp_colors) <- grp_levels
  non_prior_idx <- 1L
  for (g in grp_levels) {
    if (identical(g, "prior")) {
      grp_colors[g] <- "black"
    } else {
      grp_colors[g] <- pal_colors[non_prior_idx]
      non_prior_idx <- non_prior_idx + 1L
    }
  }

  q_breaks <- function(lims) {
    b <- c(
      lims[1L] + 0.050 * (lims[2L] - lims[1L]),
      lims[1L] + 0.500 * (lims[2L] - lims[1L]),
      lims[1L] + 0.950 * (lims[2L] - lims[1L])
    )
    round(b, 1L)
  }

  p <- ggplot()

  if (!is.null(cat_borders_df) && nrow(cat_borders_df) > 0L) {
    p <- p + geom_rect(
      data = cat_borders_df,
      aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin = .data$ymin,
        ymax = .data$ymax
      ),
      color = "gray50",
      fill = NA,
      linewidth = 0.4,
      inherit.aes = FALSE
    )
  }

  if (nrow(lower_df) > 0L) {
    p <- p +
      geom_point(
        data = lower_df,
        aes(
          x = .data$Val1,
          y = .data$Val2,
          color = .data$Group
        ),
        alpha = 0.30,
        size = 0.5,
        show.legend = show_grp_legend
      ) +
      geom_smooth(
        data = lower_df,
        aes(
          x = .data$Val1,
          y = .data$Val2,
          color = .data$Group
        ),
        method = "loess",
        formula = y ~ x,
        se = FALSE,
        linewidth = 0.35,
        show.legend = FALSE
      )
  }

  if (nrow(diag_df) > 0L) {
    p <- p +
      geom_polygon(
        data = diag_df,
        aes(
          x = .data$Val1,
          y = .data$Val2,
          color = .data$Group,
          fill = .data$Group
        ),
        linewidth = 0.25,
        alpha = 0.25,
        show.legend = FALSE
      )
  }

  if (nrow(upper_df_valid) > 0L) {
    p <- p +
      geom_density_2d(
        data = upper_df_valid,
        aes(
          x = .data$Val1,
          y = .data$Val2,
          color = .data$Group
        ),
        bins = 4,
        contour_var = "ndensity",
        linewidth = 0.30,
        show.legend = FALSE
      )
  }

  p <- p +
    facet_grid(
      Var2_expr ~ Var1_expr,
      scales = "free",
      labeller = label_parsed
    ) +
    scale_x_continuous(breaks = q_breaks) +
    scale_y_continuous(breaks = q_breaks) +
    labs(
      title = t_sub$title,
      subtitle = t_sub$subtitle,
      x = NULL,
      y = NULL
    ) +
    .mvbu_theme() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6.5),
      axis.text.y = element_text(size = 6.5),
      panel.spacing = unit(2, "pt"),
      strip.background = element_rect(fill = "white", color = "gray80"),
      strip.text = element_text(size = 8)
    )

  if (show_grp_legend) {
    p <- p +
      scale_color_manual(values = grp_colors, name = "Group") +
      scale_fill_manual(values = grp_colors, name = "Group")
  } else {
    p <- p +
      scale_color_manual(values = grp_colors, guide = "none") +
      scale_fill_manual(values = grp_colors, guide = "none")
  }

  p
}

# -----------------------------------------------------------------------------
# plot_sample methods
# -----------------------------------------------------------------------------

#' @rdname plot_sample
#' @export
S7::method(plot_sample, MVBU_Object) <- function(
  x,
  sample = c("exposure", "test"),
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  densities = c("gaussian", "kernel"),
  n_samples = NULL,
  levels = NULL,
  limits = NULL,
  resolution = 100,
  ...
) {
  sample <- match.arg(sample, choices = c("exposure", "test"), several.ok = TRUE)
  densities <- match.arg(densities)

  obj_cues <- get_cue_labels(x)
  if (is.null(cues)) {
    cues <- if (length(obj_cues) > 3L) obj_cues[1:3] else obj_cues
  }
  if (length(cues) > 3L) {
    .stop(
      "Cannot plot more than 3 cue dimensions simultaneously (received ",
      length(cues), " cues: ", paste(cues, collapse = ", "),
      "). Please specify 1-3 cues; non-plotted cues will be analytically marginalized out."
    )
  }
  if (!all(cues %in% obj_cues)) {
    .stop(
      "Invalid cues: ",
      paste(setdiff(cues, obj_cues), collapse = ", "),
      ". Available cues: ",
      paste(obj_cues, collapse = ", ")
    )
  }

  exp_data <- if ("exposure" %in% sample) {
    tryCatch(
      get_exposure_data(
        x,
        groups = groups,
        categories = categories,
        n_samples = n_samples,
        original_names = FALSE
      ),
      error = function(e) {
        # If groups or categories were specified and failed, rethrow
        if (!is.null(groups) || !is.null(categories)) {
          stop(e)
        }
        NULL
      }
    )
  } else {
    NULL
  }

  test_data <- if ("test" %in% sample) {
    tryCatch(
      get_test_data(
        x,
        groups = groups,
        response_categories = categories,
        n_samples = n_samples,
        original_names = FALSE
      ),
      error = function(e) {
        # If groups or categories were specified and failed, rethrow
        if (!is.null(groups) || !is.null(categories)) {
          stop(e)
        }
        NULL
      }
    )
  } else {
    NULL
  }

  if ((is.null(exp_data) || nrow(exp_data) == 0L) && (is.null(test_data) || nrow(test_data) == 0L)) {
    .stop("Object does not contain any exposure or test data to plot.")
  }

  # Standardize Category & Group columns
  if (!is.null(exp_data) && nrow(exp_data) > 0L) {
    exp_data$Category <- as.factor(exp_data$category)
    exp_data$Group <- as.factor(exp_data$group_unique)
  }
  if (!is.null(test_data) && nrow(test_data) > 0L) {
    test_data$Category <- as.factor(test_data$response_category)
    test_data$Group <- as.factor(test_data$group_unique)
  }

  # Parse aes
  if (is.list(aes)) {
    exp_aes <- if (!is.null(aes$exposure)) aes$exposure else c("fill-gradient", "contour")
    test_aes <- if (!is.null(aes$test)) aes$test else "points"
  } else if (is.character(aes)) {
    exp_aes <- aes
    test_aes <- aes
  } else {
    exp_aes <- if (length(cues) == 1L) "fill-discrete" else c("fill-gradient", "contour")
    test_aes <- "points"
  }

  t_sub <- .make_plot_title_and_subtitle(x, "sample", cues = cues)

  all_grps <- unique(c(
    if (!is.null(exp_data) && "Group" %in% names(exp_data)) as.character(exp_data$Group),
    if (!is.null(test_data) && "Group" %in% names(test_data)) as.character(test_data$Group)
  ))

  if (length(all_grps) == 0L) {
    all_grps <- "1"
  }

  p_list <- lapply(all_grps, function(g) {
    exp_g <- if (!is.null(exp_data) && "Group" %in% names(exp_data)) {
      exp_data[exp_data$Group == g, , drop = FALSE]
    } else {
      exp_data
    }
    test_g <- if (!is.null(test_data) && "Group" %in% names(test_data)) {
      test_data[test_data$Group == g, , drop = FALSE]
    } else {
      test_data
    }

    # 1. Base category representation plot using plot_categories
    p_g <- if (!is.null(exp_g) && nrow(exp_g) > 0L && "category" %in% names(exp_g)) {
      exp_clean <- stats::na.omit(exp_g[, c(cues, "category"), drop = FALSE])
      tpl <- if (identical(densities, "gaussian")) {
        new_mvg_category_representation_template_from_data(
          data = exp_clean,
          category = "category",
          cues = cues
        )
      } else {
        new_exemplar_category_representation_template_from_data(
          data = exp_clean,
          category = "category",
          cues = cues
        )
      }
      plot_categories(
        tpl,
        cues = cues,
        categories = categories,
        aes = exp_aes,
        levels = levels,
        limits = limits,
        resolution = resolution,
        ...
      )
    } else if (!is.null(test_g) && nrow(test_g) > 0L && any(c("contour", "fill-gradient", "fill-discrete") %in% test_aes)) {
      cat_col <- if ("Category" %in% names(test_g)) "Category" else if ("category" %in% names(test_g)) "category" else NULL
      test_clean <- if (!is.null(cat_col)) {
        test_g[, c(cues, cat_col), drop = FALSE]
      } else {
        cbind(test_g[, cues, drop = FALSE], category = "test")
      }
      names(test_clean)[ncol(test_clean)] <- "category"
      test_clean <- stats::na.omit(test_clean)
      tpl <- if (identical(densities, "gaussian")) {
        new_mvg_category_representation_template_from_data(
          data = test_clean,
          category = "category",
          cues = cues
        )
      } else {
        new_exemplar_category_representation_template_from_data(
          data = test_clean,
          category = "category",
          cues = cues
        )
      }
      plot_categories(
        tpl,
        cues = cues,
        categories = categories,
        aes = test_aes,
        levels = levels,
        limits = limits,
        resolution = resolution,
        ...
      )
    } else {
      ggplot2::ggplot() +
        .mvbu_theme()
    }

    # 2. Add exposure points if requested in aes and density is not exemplar (exemplar plot_categories already handles exemplars)
    if (!is.null(exp_g) && nrow(exp_g) > 0L && "points" %in% exp_aes && !identical(densities, "kernel") && inherits(p_g, "ggplot")) {
      if (length(cues) == 1L) {
        p_g <- p_g + ggplot2::geom_rug(
          data = exp_g,
          ggplot2::aes(x = .data[[cues[1L]]], color = .data$Category),
          alpha = 0.6,
          inherit.aes = FALSE
        )
      } else if (length(cues) == 2L) {
        p_g <- p_g + ggplot2::geom_point(
          data = exp_g,
          ggplot2::aes(x = .data[[cues[1L]]], y = .data[[cues[2L]]], color = .data$Category),
          alpha = 0.6,
          inherit.aes = FALSE
        )
      }
    }

    # 3. Add test points / layers
    if (!is.null(test_g) && nrow(test_g) > 0L && inherits(p_g, "ggplot")) {
      if ("points" %in% test_aes) {
        if (length(cues) == 1L) {
          p_g <- p_g + ggplot2::geom_rug(
            data = test_g,
            ggplot2::aes(x = .data[[cues[1L]]]),
            color = "black",
            alpha = 0.7,
            inherit.aes = FALSE
          )
        } else if (length(cues) == 2L) {
          p_g <- p_g + ggplot2::geom_point(
            data = test_g,
            ggplot2::aes(x = .data[[cues[1L]]], y = .data[[cues[2L]]]),
            color = "black",
            alpha = 0.7,
            size = 1.5,
            inherit.aes = FALSE
          )
        }
      }

      # If both exposure and test are present and test density was requested
      if (!is.null(exp_g) && nrow(exp_g) > 0L && any(c("contour", "fill-gradient", "fill-discrete") %in% test_aes)) {
        if (length(cues) == 1L) {
          vals <- test_g[[cues[1L]]][!is.na(test_g[[cues[1L]]])]
          if (length(vals) >= 2L) {
            grid <- seq(min(vals), max(vals), length.out = resolution)
            df_test_dens <- tibble::tibble(
              cue = grid,
              density = stats::dnorm(grid, mean(vals), stats::sd(vals))
            )
            names(df_test_dens)[1L] <- cues[1L]
            p_g <- p_g + ggplot2::geom_line(
              data = df_test_dens,
              ggplot2::aes(x = .data[[cues[1L]]], y = .data$density),
              color = "black",
              linetype = "dashed",
              inherit.aes = FALSE
            )
          }
        } else if (length(cues) == 2L) {
          sub_test <- stats::na.omit(test_g[, cues[1:2], drop = FALSE])
          if (nrow(sub_test) >= 3L) {
            X_test <- as.matrix(sub_test)
            mu_t <- colMeans(X_test)
            Sigma_t <- stats::cov(X_test)
            if (all(is.finite(Sigma_t)) && det(Sigma_t) > 1e-9) {
              dfs_test_cont <- list()
              lvl_spec <- .parse_category_levels(levels, test_aes)
              for (lvl in lvl_spec$contour) {
                el <- ellipse::ellipse(Sigma_t, centre = mu_t, level = lvl)
                df_el <- tibble::as_tibble(el)
                names(df_el) <- cues[1:2]
                df_el$level <- lvl
                dfs_test_cont[[length(dfs_test_cont) + 1L]] <- df_el
              }
              test_cont_data <- dplyr::bind_rows(dfs_test_cont)
              p_g <- p_g + ggplot2::geom_path(
                data = test_cont_data,
                ggplot2::aes(x = .data[[cues[1L]]], y = .data[[cues[2L]]], group = .data$level),
                color = "black",
                linetype = "dashed",
                inherit.aes = FALSE
              )
            }
          }
        }
      }
    }

    # If axis labels are missing (e.g. from empty ggplot), add them
    if (inherits(p_g, "ggplot") && is.null(p_g$labels$x)) {
      p_g <- p_g + ggplot2::labs(
        x = cues[1L],
        y = if (length(cues) == 1L) "Density" else if (length(cues) >= 2L) cues[2L] else NULL
      )
    }

    if (!is.null(limits) && inherits(p_g, "ggplot") && is.null(p_g$coordinates$limits$x)) {
      lim_spec <- .parse_cue_limits(limits, cues)
      if (length(cues) == 1L && !is.null(lim_spec[[cues[1L]]])) {
        p_g <- p_g + ggplot2::coord_cartesian(xlim = lim_spec[[cues[1L]]])
      } else if (length(cues) >= 2L) {
        p_g <- p_g + ggplot2::coord_cartesian(
          xlim = lim_spec[[cues[1L]]],
          ylim = lim_spec[[cues[2L]]]
        )
      }
    }

    if (length(all_grps) > 1L && inherits(p_g, "ggplot")) {
      p_g <- p_g + ggplot2::ggtitle(g)
    }

    p_g
  })

  if (length(all_grps) == 1L) {
    p <- p_list[[1L]]
    if (inherits(p, "ggplot")) {
      p <- p + ggplot2::labs(title = t_sub$title, subtitle = t_sub$subtitle)
    }
    return(p)
  }

  patchwork::wrap_plots(p_list, guides = "collect") +
    patchwork::plot_annotation(title = t_sub$title, subtitle = t_sub$subtitle)
}


# -----------------------------------------------------------------------------
# plot_diagnostics methods
# -----------------------------------------------------------------------------

#' @rdname plot_diagnostics
#' @export
S7::method(plot_diagnostics, MVBU_Stanfit) <- function(x, ...) {
  fit <- get_stanfit(x)
  if (is.null(fit)) {
    stop("No stanfit found in MVBU_Stanfit object.")
  }
  rstan::stan_rhat(fit, ...)
}

# -----------------------------------------------------------------------------
# Base S7 plot() dispatch
# -----------------------------------------------------------------------------

#' @export
S7::method(plot, MVBU_CategoryRepresentation) <- function(x, ...) {
  plot_categories(x, ...)
}

#' @export
S7::method(plot, MVBU_CategoryRepresentationTemplate) <- function(x, ...) {
  plot_categories(x, ...)
}

#' @export
S7::method(plot, MVBU_CognitiveModel) <- function(x, ...) {
  plot_categories(x, ...)
}

#' @export
S7::method(plot, MVBU_Stanfit) <- function(x, ...) {
  plot_categories(x, ...)
}

# -----------------------------------------------------------------------------
# plot_model_updates methods
# -----------------------------------------------------------------------------

#' @rdname plot_model_updates
#' @export
S7::method(plot_model_updates, S7::class_list) <- function(
  x,
  what = c("categories", "categorization_functions"),
  groups = NULL,
  step_labels = NULL,
  cues = NULL,
  categories = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  resolution = 100,
  ncol = NULL,
  ...
) {
  what <- match.arg(what)
  .assert_true(length(x) >= 1L, msg = "x must contain at least one model or template object.")

  if (is.null(step_labels)) {
    if (!is.null(names(x)) && any(nzchar(names(x)))) {
      step_labels <- names(x)
      empty_idx <- which(!nzchar(step_labels))
      if (length(empty_idx) > 0) {
        step_labels[empty_idx] <- sprintf("Step %d", empty_idx - 1L)
      }
    } else {
      step_labels <- sprintf("Step %d", seq_along(x) - 1L)
    }
  } else {
    .assert_true(length(step_labels) == length(x),
      msg = "length of step_labels must match length of model list."
    )
  }

  step_factors <- factor(step_labels, levels = unique(step_labels))

  plots_list <- list()
  dfs <- list()

  for (i in seq_along(x)) {
    mod <- x[[i]]
    st_name <- as.character(step_factors[i])

    if (what == "categories") {
      p_i <- plot_categories(
        mod,
        cues = cues, categories = categories, aes = aes,
        levels = levels, limits = limits, resolution = resolution, ...
      )
      df_i <- p_i$data
      if (!is.null(df_i) && nrow(df_i) > 0) {
        df_i$.update_step <- st_name
        dfs[[length(dfs) + 1L]] <- df_i
      }
      plots_list[[i]] <- p_i
    } else {
      p_i <- plot_categorization_functions(
        mod,
        cues = cues, categories = categories, aes = aes,
        levels = levels, limits = limits, resolution = resolution, ...
      )
      df_i <- p_i$data
      if (!is.null(df_i) && nrow(df_i) > 0) {
        df_i$.update_step <- st_name
        dfs[[length(dfs) + 1L]] <- df_i
      }
      plots_list[[i]] <- p_i
    }
  }

  if (length(dfs) == 0L) {
    stop("No plot data could be constructed for the provided model list.")
  }

  combined_df <- dplyr::bind_rows(dfs)
  combined_df$.update_step <- factor(combined_df$.update_step, levels = levels(step_factors))

  base_plot <- plots_list[[1L]]
  base_plot$data <- combined_df

  p_faceted <- base_plot + ggplot2::facet_wrap(~.update_step, ncol = ncol) +
    ggplot2::labs(
      title = "Model Belief Updates Across Stages",
      subtitle = sprintf("Comparison of %d update steps", length(x))
    )

  p_faceted
}

#' @rdname plot_model_updates
#' @export
S7::method(plot_model_updates, MVBU_Stanfit) <- function(
  x,
  what = c("categories", "categorization_functions"),
  groups = NULL,
  step_labels = NULL,
  cues = NULL,
  categories = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  resolution = 100,
  ncol = NULL,
  ...
) {
  what <- match.arg(what)
  avail_grps <- get_group_labels(x, include_prior = FALSE)
  if (is.null(groups)) {
    groups <- if (length(avail_grps) > 0) avail_grps[1L] else "group1"
  }
  if (is.null(categories)) {
    categories <- get_category_labels(x)
  }

  d_prior <- get_draws(x, categories = categories, groups = "prior", summarize = TRUE, ...)
  d_post <- get_draws(x, categories = categories, groups = groups[1L], summarize = TRUE, ...)

  cats <- get_category_labels(x)
  if (!is.null(categories)) cats <- intersect(cats, categories)
  cues_list <- get_cue_labels(x)
  if (is.null(cues)) cues <- cues_list

  prior_reps <- list()
  post_reps <- list()

  for (cat_name in cats) {
    pr_sub <- d_prior[d_prior$category == cat_name, ]
    if (nrow(pr_sub) > 0) {
      m_vec <- pr_sub$m[[1L]]
      Sigma_mat <- pr_sub$S[[1L]]
      if (length(m_vec) == 1L) {
        prior_reps[[cat_name]] <- new_uvg_category_representation(
          category_labels = cat_name, cue_labels = cues_list,
          mu = as.numeric(m_vec), sigma2 = as.numeric(Sigma_mat)
        )
      } else {
        prior_reps[[cat_name]] <- new_mvg_category_representation(
          category_labels = cat_name, cue_labels = cues_list,
          mu = as.vector(m_vec), Sigma = as.matrix(Sigma_mat)
        )
      }
    }

    po_sub <- d_post[d_post$category == cat_name, ]
    if (nrow(po_sub) > 0) {
      m_vec <- po_sub$m[[1L]]
      Sigma_mat <- po_sub$S[[1L]]
      if (length(m_vec) == 1L) {
        post_reps[[cat_name]] <- new_uvg_category_representation(
          category_labels = cat_name, cue_labels = cues_list,
          mu = as.numeric(m_vec), sigma2 = as.numeric(Sigma_mat)
        )
      } else {
        post_reps[[cat_name]] <- new_mvg_category_representation(
          category_labels = cat_name, cue_labels = cues_list,
          mu = as.vector(m_vec), Sigma = as.matrix(Sigma_mat)
        )
      }
    }
  }

  prior_tpl <- new_category_representation_template(prior_reps)
  post_tpl <- new_category_representation_template(post_reps)

  is_uvg <- S7::S7_inherits(prior_reps[[1L]], UVG_CategoryRepresentation)
  if (is_uvg) {
    prior_mod <- new_uvg_ideal_observer(prior_tpl)
    post_mod <- new_uvg_ideal_observer(post_tpl)
  } else {
    prior_mod <- new_mvg_ideal_observer(prior_tpl)
    post_mod <- new_mvg_ideal_observer(post_tpl)
  }

  mod_list <- list(
    "Prior (Expected)" = prior_mod,
    "Posterior (Expected)" = post_mod
  )

  if (!is.null(step_labels)) {
    names(mod_list) <- step_labels[seq_along(mod_list)]
  }

  plot_model_updates(
    mod_list,
    what = what, cues = cues, categories = categories,
    aes = aes, levels = levels, limits = limits, resolution = resolution,
    ncol = ncol, ...
  )
}



# Clean Plotly legend so only one legend entry per category group appears
# @noRd
.clean_plotly_legend <- function(p) {
  if (!inherits(p, "plotly")) {
    return(p)
  }
  p_built <- plotly::plotly_build(p)
  seen_groups <- list()
  if (!is.null(p_built$x$data)) {
    for (i in seq_along(p_built$x$data)) {
      lg <- as.character(p_built$x$data[[i]]$legendgroup %||% p_built$x$data[[i]]$name %||% i)
      raw_name <- p_built$x$data[[i]]$name %||% lg
      cat_name <- sub("<br />.*$", "", raw_name)
      cat_name <- sub("^\\((.*),.*\\)$", "\\1", cat_name)
      if (is.null(seen_groups[[lg]])) {
        seen_groups[[lg]] <- TRUE
        p_built$x$data[[i]]$showlegend <- TRUE
        p_built$x$data[[i]]$name <- cat_name
      } else {
        p_built$x$data[[i]]$showlegend <- FALSE
      }
    }
  }
  if (!is.null(p_built$x$frames)) {
    for (f in seq_along(p_built$x$frames)) {
      if (!is.null(p_built$x$frames[[f]]$data)) {
        seen_f <- list()
        for (i in seq_along(p_built$x$frames[[f]]$data)) {
          lg <- as.character(p_built$x$frames[[f]]$data[[i]]$legendgroup %||% p_built$x$frames[[f]]$data[[i]]$name %||% i)
          raw_name <- p_built$x$frames[[f]]$data[[i]]$name %||% lg
          cat_name <- sub("<br />.*$", "", raw_name)
          cat_name <- sub("^\\((.*),.*\\)$", "\\1", cat_name)
          if (is.null(seen_f[[lg]])) {
            seen_f[[lg]] <- TRUE
            p_built$x$frames[[f]]$data[[i]]$showlegend <- TRUE
            p_built$x$frames[[f]]$data[[i]]$name <- cat_name
          } else {
            p_built$x$frames[[f]]$data[[i]]$showlegend <- FALSE
          }
        }
      }
    }
  }
  p_built
}

.build_3d_model_list_ellipsoids <- function(model_list, cues, categories = NULL, level = 0.95, n_pts = 40L) {
  models <- model_list@models
  labels <- model_list@model_labels
  theta <- seq(0, 2 * pi, length.out = n_pts)
  r_val <- sqrt(stats::qchisq(level, df = 3))
  c45 <- cos(pi / 4)
  s45 <- sin(pi / 4)

  # Coordinate rings on unit sphere
  r_eq <- rbind(cos(theta), sin(theta), 0)
  r_m1 <- rbind(cos(theta), 0, sin(theta))
  r_m2 <- rbind(0, cos(theta), sin(theta))
  r_p1 <- rbind(c45 * cos(theta), c45 * sin(theta), s45)
  r_p2 <- rbind(c45 * cos(theta), c45 * sin(theta), -s45)
  rings <- list(r_eq, r_m1, r_m2, r_p1, r_p2)

  all_rows <- list()
  for (i in seq_along(models)) {
    m <- models[[i]]
    step_lbl <- labels[i]
    reps <- if (S7::S7_inherits(m, MVBU_CognitiveModel)) {
      m@category_template@representations
    } else if (S7::S7_inherits(m, MVBU_CategoryRepresentationTemplate)) {
      m@representations
    } else {
      list()
    }
    if (length(reps) == 0L) next

    target_cats <- if (!is.null(categories) && length(categories) > 0L) {
      intersect(names(reps), categories)
    } else {
      names(reps)
    }

    for (cat_name in target_cats) {
      r <- reps[[cat_name]]
      mc <- tryCatch(.get_rep_mean_and_cov(r), error = function(e) NULL)
      if (is.null(mc)) next
      idx_c <- match(cues, get_cue_labels(r))
      mu <- mc$mu[idx_c]
      Sigma <- mc$Sigma[idx_c, idx_c]
      chol_S <- tryCatch(t(chol(Sigma)), error = function(e) diag(sqrt(pmax(diag(Sigma), 1e-6))))

      pts_rings <- list()
      for (rg in rings) {
        tr <- chol_S %*% (rg * r_val) + mu
        tr_na <- cbind(tr, c(NA_real_, NA_real_, NA_real_))
        pts_rings[[length(pts_rings) + 1L]] <- tr_na
      }
      cat_pts <- do.call(cbind, pts_rings)
      df_cat <- data.frame(
        x = cat_pts[1, ],
        y = cat_pts[2, ],
        z = cat_pts[3, ],
        category = cat_name,
        step_label = factor(step_lbl, levels = labels)
      )
      names(df_cat)[1:3] <- cues
      all_rows[[length(all_rows) + 1L]] <- df_cat
    }
  }
  if (length(all_rows) == 0L) {
    return(NULL)
  }
  do.call(rbind, all_rows)
}

.build_3d_model_list_mesh_plot <- function(
  model_list,
  cues,
  categories = NULL,
  levels = NULL,
  opacity = 0.6,
  resolution = 20L,
  ...
) {
  rlang::check_installed("plotly", reason = "for interactive 3D category surface plots.")
  models <- model_list@models
  labels <- model_list@model_labels
  all_cats <- get_category_labels(models[[1L]])
  if (is.null(categories) || length(categories) == 0L) {
    categories <- all_cats
  } else {
    categories <- intersect(all_cats, categories)
  }
  if (length(categories) == 0L) {
    .stop("None of the requested categories were found in the models.")
  }

  if (is.null(levels) || length(levels) == 0L) {
    levels <- 0.95
  }
  lvl <- levels[1L]
  r_val <- sqrt(stats::qchisq(lvl, df = 3))

  n_theta <- as.integer(max(resolution, 18L))
  n_phi <- as.integer(max(resolution / 2, 9L))
  theta <- seq(0, 2 * pi, length.out = n_theta)
  phi <- seq(0, pi, length.out = n_phi)
  xs <- outer(cos(theta), sin(phi))
  ys <- outer(sin(theta), sin(phi))
  zs <- outer(rep(1, n_theta), cos(phi))
  unit_sphere <- rbind(as.vector(xs), as.vector(ys), as.vector(zs))

  # Triangulate sphere (0-based indices for mesh3d in plotly)
  i_idx <- integer()
  j_idx <- integer()
  k_idx <- integer()
  for (p_i in 0:(n_phi - 2)) {
    for (t_i in 0:(n_theta - 2)) {
      p1 <- p_i * n_theta + t_i
      p2 <- p_i * n_theta + (t_i + 1)
      p3 <- (p_i + 1) * n_theta + t_i
      p4 <- (p_i + 1) * n_theta + (t_i + 1)
      i_idx <- c(i_idx, p1, p2)
      j_idx <- c(j_idx, p2, p4)
      k_idx <- c(k_idx, p3, p3)
    }
  }

  c_colors <- scales::hue_pal()(length(categories))
  names(c_colors) <- categories

  # Helper to compute category trace data for a given model checkpoint
  get_cat_traces <- function(m) {
    reps <- if (S7::S7_inherits(m, MVBU_CognitiveModel)) {
      m@category_template@representations
    } else if (S7::S7_inherits(m, MVBU_CategoryRepresentationTemplate)) {
      m@representations
    } else {
      list()
    }
    tr_list <- list()
    for (cat_name in categories) {
      r <- reps[[cat_name]]
      if (is.null(r)) next
      mc <- tryCatch(.get_rep_mean_and_cov(r), error = function(e) NULL)
      if (is.null(mc)) next
      idx_c <- match(cues, get_cue_labels(r))
      mu <- mc$mu[idx_c]
      Sigma <- mc$Sigma[idx_c, idx_c]
      chol_S <- tryCatch(
        t(chol(Sigma)),
        error = function(e) diag(sqrt(pmax(diag(Sigma), 1e-6)))
      )
      pts <- chol_S %*% (unit_sphere * r_val) + mu
      col_k <- c_colors[cat_name]

      # Category center marker trace
      tr_marker <- list(
        x = I(mu[1L]),
        y = I(mu[2L]),
        z = I(mu[3L]),
        type = "scatter3d",
        mode = "markers",
        marker = list(color = col_k, size = 6),
        name = cat_name,
        legendgroup = cat_name,
        showlegend = TRUE,
        hoverinfo = "text",
        text = paste0(cat_name, " (center)")
      )
      # Category ellipsoid mesh surface trace
      tr_mesh <- list(
        x = as.vector(pts[1, ]),
        y = as.vector(pts[2, ]),
        z = as.vector(pts[3, ]),
        i = i_idx,
        j = j_idx,
        k = k_idx,
        type = "mesh3d",
        name = cat_name,
        legendgroup = cat_name,
        opacity = opacity,
        colorscale = list(c(0, col_k), c(1, col_k)),
        intensity = rep(0, ncol(pts)),
        showscale = FALSE,
        hoverinfo = "name",
        showlegend = FALSE
      )
      tr_list <- c(tr_list, list(tr_marker, tr_mesh))
    }
    tr_list
  }

  # Build initial base plot from first model checkpoint
  p <- plotly::plot_ly()
  base_traces <- get_cat_traces(models[[1L]])
  for (tr in base_traces) {
    p <- plotly::add_trace(
      p,
      x = tr$x, y = tr$y, z = tr$z,
      i = tr$i, j = tr$j, k = tr$k,
      type = tr$type, mode = tr$mode,
      marker = tr$marker, name = tr$name,
      legendgroup = tr$legendgroup, showlegend = tr$showlegend,
      opacity = tr$opacity, colorscale = tr$colorscale,
      intensity = tr$intensity, showscale = tr$showscale,
      hoverinfo = tr$hoverinfo, text = tr$text
    )
  }

  p <- p %>% plotly::layout(
    scene = list(
      xaxis = list(title = cues[1L]),
      yaxis = list(title = cues[2L]),
      zaxis = list(title = cues[3L]),
      aspectmode = "cube"
    )
  )

  p <- plotly::plotly_build(p)

  # Explicitly attach animation frames to avoid plotly mesh3d reordering bugs
  n_traces <- length(base_traces)
  trace_indices <- as.list(seq(0L, n_traces - 1L))

  frames <- vector("list", length(models))
  for (m_i in seq_along(models)) {
    m_traces <- get_cat_traces(models[[m_i]])
    frames[[m_i]] <- list(
      name = labels[m_i],
      data = m_traces,
      traces = trace_indices
    )
  }
  p$x$frames <- frames

  steps <- lapply(seq_along(labels), function(m_i) {
    list(
      label = labels[m_i],
      value = labels[m_i],
      method = "animate",
      args = list(
        list(labels[m_i]),
        list(
          mode = "immediate",
          frame = list(duration = 500, redraw = TRUE),
          transition = list(duration = 0)
        )
      )
    )
  })

  p$x$layout$sliders <- list(
    structure(
      list(
        active = 0,
        currentvalue = list(prefix = "Model: "),
        steps = steps,
        visible = TRUE,
        pad = list(t = 40)
      ),
      class = "aniSlider"
    )
  )

  p$x$layout$updatemenus <- list(
    structure(
      list(
        type = "buttons",
        direction = "right",
        showactive = FALSE,
        y = 0, x = 0,
        yanchor = "top", xanchor = "right",
        pad = list(t = 60, r = 5),
        buttons = list(
          list(
            label = "Play",
            method = "animate",
            args = list(
              NULL,
              list(
                fromcurrent = TRUE,
                mode = "immediate",
                transition = list(duration = 0),
                frame = list(duration = 500, redraw = TRUE)
              )
            )
          )
        )
      ),
      class = "aniButton"
    )
  )

  .clean_plotly_legend(p)
}

#' @rdname plot_categories
#' @export
S7::method(plot_categories, MVBU_ModelList) <- function(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  resolution = 100,
  n_exemplars = 0L,
  ...
) {
  dots <- list(...)
  mode <- if ("mode" %in% names(dots)) {
    match.arg(dots$mode, c("facet", "animate"))
  } else {
    "facet"
  }
  interactive <- if ("interactive" %in% names(dots)) dots$interactive else FALSE
  models <- x@models
  labels <- x@model_labels
  n_mod <- length(models)
  if (n_mod == 0L) .stop("Model list is empty.")

  first_cues <- get_cue_labels(models[[1L]])
  if (is.null(cues)) {
    cues <- if (length(first_cues) > 3L) first_cues[1:3] else first_cues
  }
  d <- length(cues)

  all_cats <- get_category_labels(models[[1L]])
  if (!is.null(categories)) {
    categories <- intersect(all_cats, categories)
    if (length(categories) == 0L) {
      .stop("None of the requested categories match model category labels.")
    }
  }

  if (isFALSE(interactive)) {
    if (identical(mode, "facet")) {
      p_list <- lapply(seq_along(models), function(i) {
        p_i <- plot_categories(models[[i]], interactive = FALSE, cues = cues, categories = categories, groups = groups, aes = aes, levels = levels, ...) + ggplot2::ggtitle(labels[i])
        if (i > 1L) {
          p_i <- p_i + ggplot2::labs(subtitle = NULL, x = NULL, y = NULL)
        }
        p_i
      })
      patchwork::wrap_plots(p_list, guides = "collect")
    } else {
      rlang::check_installed("gganimate", reason = "for animating ggplot model updates.")
      df_list <- lapply(seq_along(models), function(i) {
        p <- plot_categories(models[[i]], interactive = FALSE, cues = cues, categories = categories, groups = groups, aes = aes, levels = levels, ...)
        df <- if (!is.null(p$data) && is.data.frame(p$data) && nrow(p$data) > 0L) {
          p$data
        } else if (length(p$layers) > 0L && !is.null(p$layers[[1]]$data) && is.data.frame(p$layers[[1]]$data)) {
          p$layers[[1]]$data
        } else {
          NULL
        }
        if (!is.null(df) && is.data.frame(df) && nrow(df) > 0L) {
          df$step_label <- factor(labels[i], levels = labels)
          grp_base <- if ("group" %in% colnames(df)) df$group else if ("Category" %in% colnames(df)) df$Category else if ("category" %in% colnames(df)) df$category else NULL
          if (d == 1L) {
            # In 1D, preserve category grouping across steps for smooth curve transitions
            df$step_anim_group <- if (!is.null(grp_base)) grp_base else df$step_label
          } else {
            df$step_anim_group <- if (!is.null(grp_base)) interaction(grp_base, df$step_label) else df$step_label
          }
        }
        df
      })
      combined_df <- do.call(rbind, df_list)

      base_p <- plot_categories(models[[1L]], interactive = FALSE, cues = cues, categories = categories, groups = groups, aes = aes, levels = levels, ...)
      for (l_i in seq_along(base_p$layers)) {
        if (!is.null(base_p$layers[[l_i]]$data) && is.data.frame(base_p$layers[[l_i]]$data) && nrow(base_p$layers[[l_i]]$data) > 0L) {
          base_p$layers[[l_i]]$data <- combined_df
          base_p$layers[[l_i]]$mapping$group <- quote(step_anim_group)
        }
      }
      if (!is.null(base_p$data) && is.data.frame(base_p$data) && nrow(base_p$data) > 0L) {
        base_p$data <- combined_df
        base_p$mapping$group <- quote(step_anim_group)
      }

      base_title <- base_p$labels$title %||% "Category likelihood"
      anim_p <- base_p +
        ggplot2::labs(title = sprintf("%s: {closest_state}", base_title)) +
        ggplot2::theme(plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 20)))

      if (d == 1L) {
        # 1D: Smooth curve transitions across exposure steps
        gganimate::animate(
          anim_p +
            gganimate::transition_states(step_label, transition_length = 1, state_length = 2),
          ...
        )
      } else {
        # 2D+: Whole contour shift across states without vertex morphing
        gganimate::animate(
          anim_p +
            gganimate::transition_states(step_label, transition_length = 0, state_length = 2),
          ...
        )
      }
    }
  } else {
    rlang::check_installed("plotly", reason = "for interactive Plotly plots.")
    if (identical(mode, "facet")) {
      if (d > 1L) {
        .stop(
          "Interactive multi-panel faceting (`mode = 'facet', interactive = TRUE`) ",
          "is not supported for ", d, "D representations because Plotly cannot ",
          "tile multi-dimensional WebGL scenes. To view model updating ",
          "interactively with a slider bar, use `mode = 'animate', interactive = TRUE`. ",
          "To generate a multi-panel faceted layout, use `mode = 'facet', interactive = FALSE`."
        )
      }
      p_list <- lapply(seq_along(models), function(i) {
        p_i <- plot_categories(models[[i]], interactive = TRUE, cues = cues, categories = categories, groups = groups, aes = aes, levels = levels, ...)
        if (i > 1L) {
          p_i <- plotly::style(p_i, showlegend = FALSE)
        }
        p_i
      })
      plotly::subplot(p_list, nrows = ceiling(n_mod / 2), shareX = TRUE, shareY = TRUE)
    } else {
      df_list <- lapply(seq_along(models), function(i) {
        p <- plot_categories(models[[i]], interactive = FALSE, cues = cues, categories = categories, groups = groups, aes = aes, levels = levels, ...)
        df <- if (!is.null(p$data) && is.data.frame(p$data) && nrow(p$data) > 0) p$data else p$layers[[1]]$data
        if (!is.null(df) && is.data.frame(df) && nrow(df) > 0L) {
          df$step_label <- factor(labels[i], levels = labels)
        }
        df
      })
      combined_df <- do.call(rbind, df_list)
      cat_col <- if ("Category" %in% colnames(combined_df)) "Category" else if ("category" %in% colnames(combined_df)) "category" else NULL

      if (d == 1L) {
        plotly::plot_ly(
          combined_df,
          x = as.formula(paste0("~", cues[1])),
          y = ~density,
          color = if (!is.null(cat_col)) as.formula(paste0("~", cat_col)) else NULL,
          frame = ~step_label,
          type = "scatter",
          mode = "lines"
        )
      } else if (d == 2L) {
        split_formula <- if (!is.null(cat_col) && "level" %in% colnames(combined_df)) {
          as.formula(paste0("~interaction(", cat_col, ", level)"))
        } else if (!is.null(cat_col)) {
          as.formula(paste0("~", cat_col))
        } else {
          NULL
        }
        p_res <- plotly::plot_ly(
          combined_df,
          x = as.formula(paste0("~", cues[1])),
          y = as.formula(paste0("~", cues[2])),
          color = if (!is.null(cat_col)) as.formula(paste0("~", cat_col)) else NULL,
          split = split_formula,
          legendgroup = if (!is.null(cat_col)) as.formula(paste0("~", cat_col)) else NULL,
          frame = ~step_label,
          type = "scatter",
          mode = "lines"
        )
        .clean_plotly_legend(p_res)
      } else {
        is_wireframe <- "contour" %in% aes &&
          !any(c("fill", "fill-discrete", "surface") %in% aes)
        if (is_wireframe) {
          ell_df <- .build_3d_model_list_ellipsoids(
            model_list = x,
            cues = cues, categories = categories,
            level = if (!is.null(levels)) levels[1L] else 0.95
          )
          if (!is.null(ell_df) && nrow(ell_df) > 0L) {
            p_res <- plotly::plot_ly(
              ell_df,
              x = as.formula(paste0("~", cues[1])),
              y = as.formula(paste0("~", cues[2])),
              z = as.formula(paste0("~", cues[3])),
              color = ~category,
              frame = ~step_label,
              type = "scatter3d",
              mode = "lines",
              line = list(width = 3)
            )
            return(.clean_plotly_legend(p_res))
          }
        }
        # 3D Shaded Surfaces (mesh3d)
        mesh_dots <- dots
        mesh_dots[c("mode", "interactive")] <- NULL
        do.call(
          .build_3d_model_list_mesh_plot,
          c(
            list(
              model_list = x,
              cues = cues,
              categories = categories,
              levels = levels,
              aes = aes
            ),
            mesh_dots
          )
        )
      }
    }
  }
}


#' @rdname plot_categorization_functions
#' @export
S7::method(plot_categorization_functions, MVBU_ModelList) <- function(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  decision_rule = "proportional",
  resolution = 100,
  noise_treatment = NULL,
  lapse_treatment = NULL,
  ...
) {
  dots <- list(...)
  mode <- if ("mode" %in% names(dots)) {
    match.arg(dots$mode, c("facet", "animate"))
  } else {
    "facet"
  }
  interactive <- if ("interactive" %in% names(dots)) dots$interactive else FALSE
  models <- x@models
  labels <- x@model_labels
  n_mod <- length(models)
  if (n_mod == 0L) .stop("Model list is empty.")

  first_cues <- get_cue_labels(models[[1L]])
  if (is.null(cues)) {
    cues <- if (length(first_cues) > 3L) first_cues[1:3] else first_cues
  }
  d <- length(cues)

  all_cats <- get_category_labels(models[[1L]])
  if (!is.null(categories)) {
    categories <- intersect(all_cats, categories)
    if (length(categories) == 0L) {
      .stop("None of the requested categories match model category labels.")
    }
  }

  if (isFALSE(interactive)) {
    if (identical(mode, "facet")) {
      p_list <- lapply(seq_along(models), function(i) {
        p_i <- plot_categorization_functions(models[[i]], interactive = FALSE, cues = cues, categories = categories, groups = groups, ...) + ggplot2::ggtitle(labels[i])
        if (i > 1L) {
          p_i <- p_i + ggplot2::labs(subtitle = NULL, x = NULL, y = NULL)
        }
        p_i
      })
      patchwork::wrap_plots(p_list, guides = "collect")
    } else {
      rlang::check_installed("gganimate", reason = "for animating ggplot model updates.")
      df_list <- lapply(seq_along(models), function(i) {
        p <- plot_categorization_functions(models[[i]], interactive = FALSE, cues = cues, categories = categories, groups = groups, ...)
        df <- if (!is.null(p$data) && is.data.frame(p$data) && nrow(p$data) > 0L) {
          p$data
        } else if (length(p$layers) > 0L && !is.null(p$layers[[1]]$data) && is.data.frame(p$layers[[1]]$data)) {
          p$layers[[1]]$data
        } else {
          NULL
        }
        if (!is.null(df) && is.data.frame(df) && nrow(df) > 0L) {
          df$step_label <- factor(labels[i], levels = labels)
          grp_base <- if ("group" %in% colnames(df)) df$group else if ("Category" %in% colnames(df)) df$Category else if ("category" %in% colnames(df)) df$category else NULL
          df$step_anim_group <- if (!is.null(grp_base)) interaction(grp_base, df$step_label) else df$step_label
        }
        df
      })
      combined_df <- do.call(rbind, df_list)

      base_p <- plot_categorization_functions(models[[1L]], interactive = FALSE, cues = cues, categories = categories, groups = groups, ...)
      for (l_i in seq_along(base_p$layers)) {
        if (!is.null(base_p$layers[[l_i]]$data) && is.data.frame(base_p$layers[[l_i]]$data) && nrow(base_p$layers[[l_i]]$data) > 0L) {
          base_p$layers[[l_i]]$data <- combined_df
          base_p$layers[[l_i]]$mapping$group <- quote(step_anim_group)
        }
      }
      if (!is.null(base_p$data) && is.data.frame(base_p$data) && nrow(base_p$data) > 0L) {
        base_p$data <- combined_df
        base_p$mapping$group <- quote(step_anim_group)
      }
      base_title <- base_p$labels$title %||% "Categorization function"
      anim_p <- base_p +
        ggplot2::labs(title = sprintf("%s: {closest_state}", base_title)) +
        ggplot2::theme(plot.title = ggplot2::element_text(margin = ggplot2::margin(b = 28)))

      if (d == 1L) {
        gganimate::animate(
          anim_p +
            gganimate::transition_states(step_label, transition_length = 1, state_length = 2),
          ...
        )
      } else {
        # 2D+: Whole contour shift without individual contour line vertex morphing
        gganimate::animate(
          anim_p +
            gganimate::transition_states(step_label, transition_length = 0, state_length = 2),
          ...
        )
      }
    }
  } else {
    rlang::check_installed("plotly", reason = "for interactive Plotly plots.")
    if (identical(mode, "facet")) {
      if (d > 1L) {
        .stop(
          "Interactive multi-panel faceting (`mode = 'facet', interactive = TRUE`) ",
          "is not supported for ", d, "D representations because Plotly cannot ",
          "tile multi-dimensional WebGL scenes. To view model updating ",
          "interactively with a slider bar, use `mode = 'animate', interactive = TRUE`. ",
          "To generate a multi-panel faceted layout, use `mode = 'facet', interactive = FALSE`."
        )
      }
      p_list <- lapply(seq_along(models), function(i) {
        p_i <- plot_categorization_functions(models[[i]], interactive = TRUE, cues = cues, categories = categories, groups = groups, ...)
        if (i > 1L) {
          p_i <- plotly::style(p_i, showlegend = FALSE)
        }
        p_i
      })
      plotly::subplot(p_list, nrows = ceiling(n_mod / 2), shareX = TRUE, shareY = TRUE)
    } else {
      if (d > 2L) {
        .stop(
          "Interactive animated response functions (`mode = 'animate', interactive = TRUE`) ",
          "are supported for 1D and 2D cue spaces. For ", d, "D representations, use ",
          "`mode = 'facet', interactive = FALSE`."
        )
      }
      p_list <- lapply(seq_along(models), function(i) {
        plot_categorization_functions(
          models[[i]],
          interactive = FALSE,
          cues = cues,
          categories = categories,
          groups = groups,
          ...
        )
      })
      df1 <- if (!is.null(p_list[[1L]]$data) && is.data.frame(p_list[[1L]]$data) && nrow(p_list[[1L]]$data) > 0L) {
        p_list[[1L]]$data
      } else {
        p_list[[1L]]$layers[[1L]]$data
      }
      cat_col <- if ("Category" %in% colnames(df1)) "Category" else if ("category" %in% colnames(df1)) "category" else NULL
      candidate_prob_cols <- c("Probability", "posterior", "p", "density")
      prob_col <- candidate_prob_cols[candidate_prob_cols %in% colnames(df1)][1L]
      if (is.na(prob_col)) {
        prob_col <- colnames(df1)[sapply(df1, is.numeric)][1L]
      }

      if (d == 1L) {
        df_list <- lapply(seq_along(p_list), function(i) {
          df_i <- if (!is.null(p_list[[i]]$data) && is.data.frame(p_list[[i]]$data) && nrow(p_list[[i]]$data) > 0L) {
            p_list[[i]]$data
          } else {
            p_list[[i]]$layers[[1L]]$data
          }
          if (!is.null(df_i) && is.data.frame(df_i) && nrow(df_i) > 0L) {
            df_i$step_label <- factor(labels[i], levels = labels)
          }
          df_i
        })
        combined_df <- do.call(rbind, df_list)
        p_res <- plotly::plot_ly(
          combined_df,
          x = as.formula(paste0("~", cues[1L])),
          y = as.formula(paste0("~", prob_col)),
          color = if (!is.null(cat_col)) as.formula(paste0("~", cat_col)) else NULL,
          frame = ~step_label,
          type = "scatter",
          mode = "lines"
        )
        return(.clean_plotly_legend(p_res))
      } else {
        # d == 2L: 3D WebGL Surface animation across steps
        gx <- sort(unique(df1[[cues[1L]]]))
        gy <- sort(unique(df1[[cues[2L]]]))
        cats_present <- if (!is.null(cat_col)) unique(df1[[cat_col]]) else "Model"
        cat_colors <- scales::hue_pal()(length(cats_present))
        names(cat_colors) <- cats_present

        p <- plotly::plot_ly()
        for (cat_name in cats_present) {
          sub1 <- if (!is.null(cat_col)) df1[df1[[cat_col]] == cat_name, , drop = FALSE] else df1
          z1 <- t(matrix(sub1[[prob_col]], nrow = length(gx), ncol = length(gy)))
          col <- cat_colors[cat_name]
          col_rgb <- grDevices::col2rgb(col)
          p <- .add_2d_surface_trace(
            p = p,
            gx = gx,
            gy = gy,
            z_surf = z1,
            cat_name = cat_name,
            col = col,
            col_rgb = col_rgb,
            aes = "fill-discrete",
            hoverinfo = "name+z"
          )
        }

        trace_indices <- as.list(seq(0L, length(cats_present) - 1L))
        frames <- lapply(seq_along(models), function(i) {
          p_dat <- p_list[[i]]$data
          df_i <- if (!is.null(p_dat) && is.data.frame(p_dat) &&
                      nrow(p_dat) > 0L) {
            p_dat
          } else {
            p_list[[i]]$layers[[1L]]$data
          }
          frame_data <- lapply(cats_present, function(cat_name) {
            sub_i <- if (!is.null(cat_col)) {
              df_i[df_i[[cat_col]] == cat_name, , drop = FALSE]
            } else {
              df_i
            }
            z_i <- t(matrix(
              sub_i[[prob_col]],
              nrow = length(gx),
              ncol = length(gy)
            ))
            list(z = z_i)
          })
          list(name = labels[i], data = frame_data, traces = trace_indices)
        })

        p <- p %>% plotly::layout(
          scene = list(
            xaxis = list(title = cues[1L]),
            yaxis = list(title = cues[2L]),
            zaxis = list(title = "Response Probability", range = c(0, 1)),
            camera = list(eye = list(x = 1.5, y = 1.5, z = 1.2))
          ),
          legend = list(orientation = "h", x = 0.1, y = -0.1)
        )
        p <- plotly::plotly_build(p)
        p$x$frames <- frames

        slider_steps <- lapply(seq_along(labels), function(i) {
          list(
            label = labels[i],
            value = labels[i],
            method = "animate",
            args = list(
              list(labels[i]),
              list(
                mode = "immediate",
                frame = list(duration = 500, redraw = TRUE),
                transition = list(duration = 0)
              )
            )
          )
        })

        p$x$layout$sliders <- list(
          structure(
            list(
              active = 0,
              currentvalue = list(prefix = "Step: "),
              steps = slider_steps,
              visible = TRUE,
              pad = list(t = 40)
            ),
            class = "aniSlider"
          )
        )

        p$x$layout$updatemenus <- list(
          structure(
            list(
              type = "buttons",
              direction = "right",
              showactive = FALSE,
              y = 0, x = 0,
              yanchor = "top", xanchor = "right",
              pad = list(t = 60, r = 5),
              buttons = list(
                list(
                  label = "Play",
                  method = "animate",
                  args = list(
                    NULL,
                    list(
                      fromcurrent = TRUE,
                      mode = "immediate",
                      transition = list(duration = 0),
                      frame = list(duration = 500, redraw = TRUE)
                    )
                  )
                )
              )
            ),
            class = "aniButton"
          )
        )
        return(p)
      }
    }
  }
}


#' @rdname plot_categories
#' @export
S7::method(plot_categories, MVBU_StanfitPosterior) <- function(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  resolution = 100,
  n_exemplars = 0L,
  ...
) {
  dots <- list(...)
  interactive <- if ("interactive" %in% names(dots)) dots$interactive else FALSE
  # Coerce lightweight posterior to an expected model and plot
  exp_model <- .as_cognitive_model_from_posterior(x)
  plot_categories(
    exp_model,
    cues = cues,
    categories = categories,
    groups = groups,
    aes = aes,
    levels = levels,
    limits = limits,
    resolution = resolution,
    n_exemplars = n_exemplars,
    interactive = interactive,
    ...
  )
}

#' @rdname plot_categorization_functions
#' @export
S7::method(plot_categorization_functions, MVBU_StanfitPosterior) <- function(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  decision_rule = "proportional",
  resolution = 100,
  noise_treatment = NULL,
  lapse_treatment = NULL,
  ...
) {
  dots <- list(...)
  interactive <- if ("interactive" %in% names(dots)) dots$interactive else FALSE
  exp_model <- .as_cognitive_model_from_posterior(x)
  plot_categorization_functions(
    exp_model,
    cues = cues,
    categories = categories,
    groups = groups,
    aes = aes,
    levels = levels,
    limits = limits,
    decision_rule = decision_rule,
    resolution = resolution,
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment,
    interactive = interactive,
    ...
  )
}
