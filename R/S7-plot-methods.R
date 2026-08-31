# =============================================================================
# S7 Plot Methods for MVBeliefUpdatr
# =============================================================================

#' @include S7-generics.R
#' @include S7-core-classes.R
#' @include S7-core-methods.R
#' @include S7-stanfit.R
#' @include S7-stanfit-methods.R
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
      panel.grid.minor = ggplot2::element_blank()
    )
}

#' Format model family name in spelled-out form with math notation
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
    return(sprintf(
      "%dD multivariate independent Gaussian ideal observer",
      D
    ))
  }
  if (S7::S7_inherits(x, NIX_IdealAdaptor)) {
    return("1D normal-inverse-chi-squared (N\u03c7\u207b\u00b2) ideal adaptor")
  }
  if (S7::S7_inherits(x, MNIX_IdealAdaptor)) {
    return(sprintf(
      paste0(
        "%dD independent normal-inverse-chi-squared ",
        "(MN\u03c7\u207b\u00b2) ideal adaptor"
      ),
      D
    ))
  }
  if (S7::S7_inherits(x, NIW_IdealAdaptor)) {
    return(sprintf(
      "%dD normal-inverse-Wishart (NW\u207b\u00b9) ideal adaptor",
      D
    ))
  }
  if (S7::S7_inherits(x, Exemplar_Model)) {
    return(sprintf("%dD exemplar model", D))
  }
  if (S7::S7_inherits(x, UVG_CategoryRepresentation)) {
    return("1D univariate Gaussian category representation")
  }
  if (S7::S7_inherits(x, MVG_CategoryRepresentation)) {
    return(sprintf("%dD multivariate Gaussian category representation", D))
  }
  if (S7::S7_inherits(x, NIX_CategoryRepresentation)) {
    return(paste0(
      "1D normal-inverse-chi-squared ",
      "(N\u03c7\u207b\u00b2) category representation"
    ))
  }
  if (S7::S7_inherits(x, MNIX_CategoryRepresentation)) {
    return(sprintf(
      paste0(
        "%dD independent normal-inverse-chi-squared ",
        "(MN\u03c7\u207b\u00b2) representation"
      ),
      D
    ))
  }
  if (S7::S7_inherits(x, NIW_CategoryRepresentation)) {
    return(sprintf(
      "%dD normal-inverse-Wishart (NW\u207b\u00b9) category representation",
      D
    ))
  }
  if (S7::S7_inherits(x, Exemplar_CategoryRepresentation)) {
    return(sprintf("%dD exemplar category representation", D))
  }
  if (S7::S7_inherits(x, MVBU_CategoryRepresentationTemplate)) {
    return(sprintf(
      "%dD category template (%d categories)",
      D,
      length(x@representations)
    ))
  }
  if (S7::S7_inherits(x, MVBU_Stanfit)) {
    mtype <- get_model_type(x)
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
      sub_txt <- sprintf("%s\n%s | %s", fam, lapse_trt_txt, noise_txt)
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
          paste(names(cp), round(cp, 2), sep = "=", collapse = ", ")
        )
      } else {
        pri_txt <- sprintf(
          "Category prior: %s",
          paste(names(cp), round(cp, 2), sep = "=", collapse = ", ")
        )
        bias_txt <- if (!is.null(lb)) {
          sprintf(
            "Lapse bias: %s",
            paste(names(lb), round(lb, 2), sep = "=", collapse = ", ")
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
  if (is.null(limits)) return(NULL)
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
    if (length(res) > 0L) return(res)
  }
  NULL
}

#' @importFrom geomtextpath geom_textcontour
NULL

#' Parse levels argument into contour and fill levels for plot_categories
#' @noRd
#' @keywords internal
.parse_category_levels <- function(levels, aes) {
  # Default levels correspond to two-tailed 1, 2, 3, 4 sigma
  def_levels <- 2 * stats::pnorm(1:4) - 1
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
  # "fill" defaults to "fill-gradient" in all cases
  aes[aes == "fill"] <- "fill-gradient"
  unique(aes)
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

  if (length(cues) == 1L) {
    idx <- match(cues[1], rep_cues)
    if (is.na(idx)) {
      stop("Requested cue not found in representation cue labels.")
    }

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
        nu = max(nu_sub, 2.01),
        sigma2 = r@S[idx, idx]
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
    if (any(is.na(idx))) {
      stop("Requested cues not found in representation cue labels.")
    }

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

  stop("Marginalization is currently supported down to 1 or 2 cues.")
}

#' Extract expected mean and covariance matrix from representation
#' @noRd
#' @keywords internal
.get_rep_mean_and_cov <- function(r) {
  if (S7::S7_inherits(r, UVG_CategoryRepresentation)) {
    return(list(mu = r@mu, Sigma = matrix(r@sigma2, 1L, 1L)))
  }
  if (S7::S7_inherits(r, MVG_CategoryRepresentation)) {
    return(list(mu = r@mu, Sigma = r@Sigma))
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
    r <- p$rep
    dens <- if (S7::S7_inherits(r, Exemplar_CategoryRepresentation)) {
      r@category_likelihood_function(grid_mat)
    } else {
      stats::dnorm(grid_x, mean = p$mu, sd = p$sd)
    }
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
.make_exemplar_sample_df <- function(reps, cues, n_exemplars = 100) {
  rows <- list()
  for (cat_name in names(reps)) {
    r <- reps[[cat_name]]
    if (S7::S7_inherits(r, Exemplar_CategoryRepresentation)) {
      mat <- as.matrix(r@exemplars)
      N <- nrow(mat)
      idx <- if (is.null(n_exemplars) || N <= n_exemplars) {
        seq_len(N)
      } else {
        sample.int(N, n_exemplars, replace = FALSE)
      }
      sub_mat <- mat[idx, cues, drop = FALSE]
      df <- as.data.frame(sub_mat)
      df$Category <- cat_name
      rows[[length(rows) + 1]] <- df
    }
  }
  if (length(rows) == 0L) return(NULL)
  dplyr::bind_rows(rows)
}

#' Create 2D ellipse contour data frame across category representations
#' @noRd
#' @keywords internal
.make_2D_category_ellipse_df <- function(reps, cues, levels = c(0.5, 0.95)) {
  dfs <- list()

  for (cat_name in names(reps)) {
    r <- reps[[cat_name]]
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
      dfs[[length(dfs) + 1]] <- df
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
    dens <- if (S7::S7_inherits(r, Exemplar_CategoryRepresentation)) {
      r@category_likelihood_function(grid_mat)
    } else {
      mc <- .get_rep_mean_and_cov(r)
      mvtnorm::dmvnorm(grid_mat, mean = mc$mu, sigma = mc$Sigma)
    }
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
    exp_data <- tryCatch(get_exposure_data(x), error = function(e) NULL)
    if (!is.null(exp_data) && all(cues %in% names(exp_data))) {
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
    test_data <- tryCatch(get_test_data(x), error = function(e) NULL)
    if (!is.null(test_data) && all(cues %in% names(test_data))) {
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
  resolution = 200,
  t_sub = list(title = "", subtitle = "")
) {
  aes <- .normalize_plot_aes(aes, default_aes = "contour")
  lvl_spec <- .parse_category_levels(levels, aes)
  lim_spec <- .parse_cue_limits(limits, cues)

  df <- .make_1D_category_density_df(
    reps,
    cues,
    resolution = resolution,
    limits = limits
  )

  all_c <- names(reps)
  c_colors <- scales::hue_pal()(length(all_c))
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

  has_fill <- any(c("fill-gradient", "fill-discrete", "fill") %in% aes)
  if (has_fill || length(reps) == 1L) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = 0, ymax = .data$density),
      alpha = 0.2
    )
  }
  if ("contour" %in% aes || length(reps) == 1L) {
    p <- p + ggplot2::geom_line(linewidth = 0.30)
  }

  ex_df_1d <- .make_exemplar_sample_df(reps, cues[1L], n_exemplars = NULL)
  if (!is.null(ex_df_1d) && nrow(ex_df_1d) > 0L) {
    p <- p + ggplot2::geom_rug(
      data = ex_df_1d,
      ggplot2::aes(
        x = .data[[cues[1L]]],
        color = .data$Category
      ),
      sides = "b",
      length = ggplot2::unit(0.04, "npc"),
      alpha = 0.7,
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
  n_exemplars = 100,
  resolution = 100,
  t_sub = list(title = "", subtitle = "")
) {
  aes <- .normalize_plot_aes(aes, default_aes = "contour")
  lvl_spec <- .parse_category_levels(levels, aes)
  lim_spec <- .parse_cue_limits(limits, cues)

  has_ex <- any(sapply(reps, function(r) {
    S7::S7_inherits(r, Exemplar_CategoryRepresentation)
  }))

  if ("fill-discrete" %in% aes && has_ex) {
    .stop(
      "'fill-discrete' is only supported for parametric category ",
      "representations."
    )
  }

  has_fill <- any(c("fill-gradient", "fill-discrete") %in% aes)
  has_contour <- "contour" %in% aes

  all_c <- names(reps)
  c_colors <- scales::hue_pal()(length(all_c))
  names(c_colors) <- all_c

  p <- ggplot2::ggplot() +
    ggplot2::labs(
      title = t_sub$title,
      subtitle = t_sub$subtitle,
      x = cues[1L],
      y = cues[2L]
    ) +
    .mvbu_theme()

  # 1. Fill Layer
  if ("fill-gradient" %in% aes) {
    dens_df <- .make_2D_category_density_grid_df(
      reps,
      cues,
      resolution = max(resolution, 60),
      limits = limits
    )
    dens_df <- dens_df %>%
      dplyr::group_by(.data$Category) %>%
      dplyr::mutate(
        rel_dens = .data$Density / max(.data$Density, 1e-12)
      ) %>%
      dplyr::ungroup()

    p <- p + ggplot2::geom_tile(
      data = dens_df,
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
  } else if ("fill-discrete" %in% aes) {
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
        group = interaction(.data$Category, .data$level)
      ),
      show.legend = c(fill = TRUE, alpha = FALSE)
    ) +
      ggplot2::scale_fill_manual(values = c_colors, name = "Category") +
      ggplot2::scale_alpha_identity()
  }

  # 2. Contour Layer
  if (has_contour) {
    if (has_ex) {
      dens_df <- .make_2D_category_density_grid_df(
        reps,
        cues,
        resolution = max(resolution, 60),
        limits = limits
      )
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
          group = interaction(.data$Category, .data$level)
        ),
        linewidth = 0.35,
        alpha = 1.0,
        show.legend = c(color = TRUE)
      ) +
        ggplot2::scale_color_manual(values = c_colors, name = "Category")
    }
  }

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

  # 4. Exemplar Points (size = 1.0)
  ex_df <- .make_exemplar_sample_df(reps, cues, n_exemplars = n_exemplars)
  if (!is.null(ex_df)) {
    p <- p + ggplot2::geom_point(
      data = ex_df,
      ggplot2::aes(
        x = .data[[cues[1L]]],
        y = .data[[cues[2L]]],
        color = .data$Category
      ),
      alpha = 0.35,
      size = 1.0,
      shape = 16,
      show.legend = FALSE
    )
  }

  if (!is.null(lim_spec)) {
    p <- p + ggplot2::coord_cartesian(
      xlim = lim_spec[[cues[1L]]],
      ylim = lim_spec[[cues[2L]]]
    )
  }

  p
}

# -----------------------------------------------------------------------------
# plot_categories methods
# -----------------------------------------------------------------------------

#' @rdname plot_mvbu
#' @export
S7::method(plot_categories, MVBU_CategoryRepresentation) <- function(
  x,
  cues = NULL,
  categories = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  n_exemplars = 100,
  resolution = 100,
  ...
) {
  obj_cues <- get_cue_labels(x)
  if (is.null(cues)) {
    cues <- if (length(obj_cues) > 2L) obj_cues[1:2] else obj_cues
  }
  if (!all(cues %in% obj_cues)) {
    .stop(
      "Invalid cues: ",
      paste(setdiff(cues, obj_cues), collapse = ", "),
      ". Available cues: ",
      paste(obj_cues, collapse = ", ")
    )
  }

  x_proj <- .marginalize_representation_to_cues(x, cues)
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
      resolution = resolution,
      t_sub = t_sub
    ))
  }

  if (length(cues) == 2L) {
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

  stop("Plotting categories is currently supported for 1 or 2 cues.")
}

#' @rdname plot_mvbu
#' @export
S7::method(plot_categories, MVBU_CategoryRepresentationTemplate) <- function(
  x,
  cues = NULL,
  categories = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  n_exemplars = 100,
  resolution = 100,
  ...
) {
  obj_cues <- get_cue_labels(x)
  if (is.null(cues)) {
    cues <- if (length(obj_cues) > 2L) obj_cues[1:2] else obj_cues
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

  reps_proj <- lapply(reps, function(r) {
    .marginalize_representation_to_cues(r, cues)
  })

  t_sub <- .make_plot_title_and_subtitle(x, "categories", cues = cues)

  if (length(cues) == 1L) {
    return(.render_1D_category_plot(
      reps = reps_proj,
      cues = cues,
      aes = aes,
      levels = levels,
      limits = limits,
      resolution = resolution,
      t_sub = t_sub
    ))
  }

  if (length(cues) == 2L) {
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

  stop("Plotting categories is currently supported for 1 or 2 cues.")
}

#' @rdname plot_mvbu
#' @export
S7::method(plot_categories, MVBU_CognitiveModel) <- function(
  x,
  cues = NULL,
  categories = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  n_exemplars = 100,
  resolution = 100,
  ...
) {
  obj_cues <- get_cue_labels(x)
  if (is.null(cues)) {
    cues <- if (length(obj_cues) > 2L) obj_cues[1:2] else obj_cues
  }
  if (!all(cues %in% obj_cues)) {
    .stop(
      "Invalid cues: ",
      paste(setdiff(cues, obj_cues), collapse = ", "),
      ". Available cues: ",
      paste(obj_cues, collapse = ", ")
    )
  }

  p <- plot_categories(
    x@category_template,
    cues = cues,
    categories = categories,
    aes = aes,
    levels = levels,
    limits = limits,
    n_exemplars = n_exemplars,
    resolution = resolution,
    ...
  )
  t_sub <- .make_plot_title_and_subtitle(x, "categories", cues = cues)
  p + labs(title = t_sub$title, subtitle = t_sub$subtitle)
}

#' @rdname plot_mvbu
#' @export
S7::method(plot_categories, MVBU_Stanfit) <- function(
  x,
  cues = NULL,
  categories = NULL,
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  ndraws = 100,
  n_exemplars = 100,
  resolution = 100,
  show_exposure_data = FALSE,
  show_test_data = FALSE,
  ...
) {
  obj_cues <- get_cue_labels(x)
  if (is.null(cues)) {
    cues <- if (length(obj_cues) > 2L) obj_cues[1:2] else obj_cues
  }
  if (!all(cues %in% obj_cues)) {
    .stop(
      "Invalid cues: ",
      paste(setdiff(cues, obj_cues), collapse = ", "),
      ". Available cues: ",
      paste(obj_cues, collapse = ", ")
    )
  }

  avail_cats <- get_category_labels(x)
  if (is.null(categories)) {
    categories <- avail_cats
  }
  avail_grps <- get_group_labels(x, include_prior = FALSE)
  if (is.null(groups)) {
    groups <- avail_grps
  }

  d_sum <- get_draws(
    x,
    categories = categories,
    groups = groups,
    ndraws = ndraws,
    summarize = FALSE,
    wide = FALSE,
    ...
  )

  if (!"Sigma" %in% names(d_sum) &&
      "S" %in% names(d_sum) &&
      "nu" %in% names(d_sum)) {
    d_sum <- dplyr::mutate(
      d_sum,
      Sigma = get_expected_Sigma_from_S(.data$S, .data$nu)
    )
  }

  d_sum <- d_sum %>%
    dplyr::group_by(.data$group, .data$category) %>%
    dplyr::summarise(
      mu.mean = list(purrr::reduce(.data$m, `+`) / length(.data$m)),
      Sigma.mean = list(purrr::reduce(.data$Sigma, `+`) / length(.data$Sigma)),
      .groups = "drop"
    )

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

    if ("fill" %in% aes) {
      p <- p + geom_ribbon(
        aes(ymin = 0, ymax = .data$density),
        alpha = 0.2
      )
    }
    if ("contour" %in% aes) {
      p <- p + geom_line(linewidth = 1)
    }

    if (length(unique(df_all$Group)) > 1L) {
      p <- p + facet_wrap(~ Group)
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

    if ("fill" %in% aes && length(dfs_fill) > 0L) {
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
      p <- p + facet_wrap(~ Group)
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
# plot_categorization_function methods
# -----------------------------------------------------------------------------

#' Create a projected cognitive model on a subset of cues
#' @noRd
#' @keywords internal
.make_projected_cognitive_model <- function(x, cues, decision_rule) {
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
      lapse_treatment = lapse_beh$lapse_treatment %||% "no_lapses",
      Sigma_noise = sig_noise,
      noise_treatment = noise_beh$noise_treatment %||% "no_noise"
    ))
  }
  if (S7::S7_inherits(first_rep, MVG_CategoryRepresentation)) {
    return(new_mvg_ideal_observer(
      category_template = new_tpl,
      category_prior = x@category_prior,
      decision_rule = decision_rule,
      lapse_rate = lapse_beh$lapse_rate %||% 0,
      lapse_bias = lapse_beh$lapse_bias %||% (1 / length(reps_proj)),
      lapse_treatment = lapse_beh$lapse_treatment %||% "no_lapses",
      Sigma_noise = sig_noise,
      noise_treatment = noise_beh$noise_treatment %||% "no_noise"
    ))
  }
  if (S7::S7_inherits(first_rep, NIX_CategoryRepresentation)) {
    return(new_nix_ideal_adaptor(
      category_template = new_tpl,
      category_prior = x@category_prior,
      decision_rule = decision_rule,
      lapse_rate = lapse_beh$lapse_rate %||% 0,
      lapse_bias = lapse_beh$lapse_bias %||% (1 / length(reps_proj)),
      lapse_treatment = lapse_beh$lapse_treatment %||% "no_lapses",
      Sigma_noise = sig_noise,
      noise_treatment = noise_beh$noise_treatment %||% "no_noise"
    ))
  }
  if (S7::S7_inherits(first_rep, MNIX_CategoryRepresentation)) {
    return(new_mnix_ideal_adaptor(
      category_template = new_tpl,
      category_prior = x@category_prior,
      decision_rule = decision_rule,
      lapse_rate = lapse_beh$lapse_rate %||% 0,
      lapse_bias = lapse_beh$lapse_bias %||% (1 / length(reps_proj)),
      lapse_treatment = lapse_beh$lapse_treatment %||% "no_lapses",
      Sigma_noise = sig_noise,
      noise_treatment = noise_beh$noise_treatment %||% "no_noise"
    ))
  }
  if (S7::S7_inherits(first_rep, NIW_CategoryRepresentation)) {
    return(new_niw_ideal_adaptor(
      category_template = new_tpl,
      category_prior = x@category_prior,
      decision_rule = decision_rule,
      lapse_rate = lapse_beh$lapse_rate %||% 0,
      lapse_bias = lapse_beh$lapse_bias %||% (1 / length(reps_proj)),
      lapse_treatment = lapse_beh$lapse_treatment %||% "no_lapses",
      Sigma_noise = sig_noise,
      noise_treatment = noise_beh$noise_treatment %||% "no_noise"
    ))
  }
  if (S7::S7_inherits(first_rep, Exemplar_CategoryRepresentation)) {
    return(new_exemplar_model(
      category_template = new_tpl,
      category_prior = x@category_prior,
      decision_rule = decision_rule,
      lapse_rate = lapse_beh$lapse_rate %||% 0,
      lapse_bias = lapse_beh$lapse_bias %||% (1 / length(reps_proj)),
      lapse_treatment = lapse_beh$lapse_treatment %||% "no_lapses"
    ))
  }
  stop("Unsupported representation type for projection.")
}

#' @rdname plot_mvbu
#' @export
S7::method(plot_categorization_function, MVBU_CognitiveModel) <- function(
  x,
  cues = NULL,
  categories = get_category_labels(x)[1L],
  aes = NULL,
  levels = NULL,
  limits = NULL,
  decision_rule = "proportional",
  resolution = 100,
  ...
) {
  obj_cues <- get_cue_labels(x)
  if (is.null(cues)) {
    cues <- if (length(obj_cues) > 2L) obj_cues[1:2] else obj_cues
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

  aes <- .normalize_plot_aes(
    aes,
    default_aes = if (length(cues) == 1L) "contour" else "contour"
  )
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
    mod_proj <- .make_projected_cognitive_model(x, cues, decision_rule)
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

    post_mat <- as.data.frame(
      posterior(mod_proj, test_mat, categories = all_cats)
    )
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
    mod_proj <- .make_projected_cognitive_model(x, cues, decision_rule)
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

    post_mat <- posterior(mod_proj, test_mat, categories = all_cats)

    dfs_cats <- list()
    for (cat_name in categories) {
      df_c <- test_grid
      df_c$Probability <- post_mat[, cat_name]
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
          alpha = 1.0,
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

  stop("Plotting categorization function is supported for 1 or 2 cues.")
}

#' @rdname plot_mvbu
#' @export
S7::method(plot_categorization_function, MVBU_Stanfit) <- function(
  x,
  cues = NULL,
  categories = get_category_labels(x)[1L],
  groups = NULL,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  ndraws = 100,
  resolution = 100,
  show_exposure_data = FALSE,
  show_test_data = FALSE,
  ...
) {
  obj_cues <- get_cue_labels(x)
  if (is.null(cues)) {
    cues <- if (length(obj_cues) > 2L) obj_cues[1:2] else obj_cues
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

  p_cat <- plot_categories(
    x,
    categories = categories,
    groups = groups,
    cues = cues,
    aes = aes,
    levels = levels,
    limits = limits,
    ndraws = ndraws,
    resolution = resolution,
    show_exposure_data = show_exposure_data,
    show_test_data = show_test_data,
    ...
  )
  t_sub <- .make_plot_title_and_subtitle(
    x,
    "categorization",
    cues = cues,
    ndraws = ndraws,
    target_category = categories[1L]
  )
  p_cat <- p_cat +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(expand = FALSE) +
    labs(title = t_sub$title, subtitle = t_sub$subtitle)
  p_cat
}

# -----------------------------------------------------------------------------
# plot_parameters methods
# -----------------------------------------------------------------------------

#' @rdname plot_mvbu
#' @export
S7::method(plot_parameters, MVBU_CognitiveModel) <- function(
  x,
  pars = NULL,
  combine_into_single_plot = TRUE,
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

#' @rdname plot_mvbu
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
    wide = FALSE,
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
        linewidth = 0.30,
        show.legend = c(color = TRUE, fill = FALSE)
      ) +
      scale_fill_discrete(guide = "none") +
      guides(fill = "none") +
      facet_grid(Group ~ Cue, scales = "free", labeller = label_parsed) +
      labs(
        title = "Category location",
        x = "Location (\u03bc)",
        y = "Density"
      ) +
      .mvbu_theme()
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
        linewidth = 0.30,
        show.legend = c(color = TRUE, fill = FALSE)
      ) +
      scale_fill_discrete(guide = "none") +
      guides(fill = "none") +
      scale_x_log10() +
      facet_grid(Group ~ Cue, scales = "free", labeller = label_parsed) +
      labs(
        title = "Category scale",
        x = "Scale (\u03c4)",
        y = "Density"
      ) +
      .mvbu_theme()
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
        linewidth = 0.30,
        show.legend = c(color = TRUE, fill = FALSE)
      ) +
      scale_fill_discrete(guide = "none") +
      guides(fill = "none") +
      scale_x_log10() +
      facet_grid(Group ~ Parameter, scales = "free", labeller = label_parsed) +
      labs(
        title = "Category confidence",
        x = "Count",
        y = "Density"
      ) +
      .mvbu_theme()
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
      facet_wrap(~ group) +
      labs(
        title = "Decision-making",
        x = "Lapse rate (\u03bb)",
        y = "Density"
      ) +
      .mvbu_theme()
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
    guides(fill = "none") &
    theme(legend.position = "right")
}

# -----------------------------------------------------------------------------
# plot_parameter_correlations methods
# -----------------------------------------------------------------------------

#' @rdname plot_mvbu
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
  if (is.null(categories)) {
    categories <- avail_cats
  }
  if (!all(categories %in% avail_cats)) {
    .stop(
      "Invalid categories: ",
      paste(setdiff(categories, avail_cats), collapse = ", "),
      ". Available categories: ",
      paste(avail_cats, collapse = ", ")
    )
  }
  avail_grps <- get_group_labels(x, include_prior = TRUE)
  if (is.null(groups)) {
    groups <- if ("prior" %in% avail_grps) "prior" else avail_grps[1L]
  }
  if (!all(groups %in% avail_grps)) {
    .stop(
      "Invalid groups: ",
      paste(setdiff(groups, avail_grps), collapse = ", "),
      ". Available groups: ",
      paste(avail_grps, collapse = ", ")
    )
  }

  draws <- get_draws(
    x,
    groups = groups,
    categories = categories,
    ndraws = ndraws,
    summarize = FALSE,
    wide = FALSE,
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

#' @rdname plot_mvbu
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
  if (is.null(categories)) {
    categories <- avail_cats
  }
  if (!all(categories %in% avail_cats)) {
    .stop(
      "Invalid categories: ",
      paste(setdiff(categories, avail_cats), collapse = ", "),
      ". Available categories: ",
      paste(avail_cats, collapse = ", ")
    )
  }
  avail_grps <- get_group_labels(x, include_prior = TRUE)
  if (is.null(groups)) {
    groups <- if ("prior" %in% avail_grps) "prior" else avail_grps[1L]
  }
  if (!all(groups %in% avail_grps)) {
    .stop(
      "Invalid groups: ",
      paste(setdiff(groups, avail_grps), collapse = ", "),
      ". Available groups: ",
      paste(avail_grps, collapse = ", ")
    )
  }

  draws <- get_draws(
    x,
    groups = groups,
    categories = categories,
    ndraws = ndraws,
    summarize = FALSE,
    wide = FALSE,
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
# plot_cues methods
# -----------------------------------------------------------------------------

#' @rdname plot_mvbu
#' @export
S7::method(plot_cues, MVBU_Object) <- function(
  x,
  cues = NULL,
  categories = NULL,
  ...
) {
  obj_cues <- get_cue_labels(x)
  if (is.null(cues)) {
    cues <- if (length(obj_cues) > 2L) obj_cues[1:2] else obj_cues
  }
  if (!all(cues %in% obj_cues)) {
    .stop(
      "Invalid cues: ",
      paste(setdiff(cues, obj_cues), collapse = ", "),
      ". Available cues: ",
      paste(obj_cues, collapse = ", ")
    )
  }

  exp_data <- tryCatch(get_exposure_data(x), error = function(e) NULL)
  test_data <- tryCatch(get_test_data(x), error = function(e) NULL)

  if (is.null(exp_data) && is.null(test_data)) {
    stop(
      "Object does not contain empirical exposure or test data to plot cues."
    )
  }

  df <- if (!is.null(exp_data)) exp_data else test_data
  if (!is.null(categories) && "category" %in% names(df)) {
    df <- df[df$category %in% categories, , drop = FALSE]
  }
  if ("category" %in% names(df)) {
    df$Category <- df$category
  }

  t_sub <- .make_plot_title_and_subtitle(x, "cues", cues = cues)

  if (length(cues) == 1L) {
    p <- ggplot(df, aes(x = .data[[cues[1]]])) +
      geom_density(alpha = 0.5) +
      labs(
        title = t_sub$title,
        subtitle = t_sub$subtitle,
        x = cues[1],
        y = "Density"
      ) +
      .mvbu_theme()
    if ("Category" %in% names(df)) {
      p <- p + aes(fill = .data$Category, color = .data$Category)
    }
    return(p)
  }

  if (length(cues) == 2L) {
    p <- ggplot(df, aes(x = .data[[cues[1]]], y = .data[[cues[2]]])) +
      geom_point(alpha = 0.6) +
      labs(
        title = t_sub$title,
        subtitle = t_sub$subtitle,
        x = cues[1],
        y = cues[2]
      ) +
      .mvbu_theme()
    if ("Category" %in% names(df)) {
      p <- p + aes(color = .data$Category)
    }
    return(p)
  }

  stop("plot_cues currently supports 1 or 2 cues.")
}

# -----------------------------------------------------------------------------
# plot_diagnostics methods
# -----------------------------------------------------------------------------

#' @rdname plot_mvbu
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
