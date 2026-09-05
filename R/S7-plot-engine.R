# =============================================================================
# S7 Plot Engine: Shared Grid, Evaluation, Parallelization, and Layer Builders
# =============================================================================

#' @include S7-generics.R
#' @include S7-core-classes.R
#' @include S7-core-methods.R
#' @import ggplot2
#' @importFrom rlang .data !! !!! sym
NULL

# -----------------------------------------------------------------------------
# Parallel Execution Dispatcher
# -----------------------------------------------------------------------------

#' Execute a function over a list of chunks, optionally in parallel
#'
#' Standardizes parallel execution across plotting operations using the base
#' R `parallel` package.
#'
#' @param chunks A list of inputs/chunks to process.
#' @param eval_fn A function to apply to each element of `chunks`.
#' @param parallel Logical; if `TRUE`, execute in parallel. Defaults to `FALSE`.
#' @param n_cores Integer; number of cores to use when `parallel = TRUE`.
#'   Defaults to `max(1L, parallel::detectCores() - 1L)`.
#' @return A list containing results for each chunk.
#' @noRd
#' @keywords internal
.mvbu_parallel_eval <- function(
  chunks,
  eval_fn,
  parallel = FALSE,
  n_cores = NULL
) {
  if (!isTRUE(parallel) || length(chunks) <= 1L) {
    return(lapply(chunks, eval_fn))
  }

  if (is.null(n_cores) || n_cores < 1L) {
    n_cores <- max(1L, parallel::detectCores() - 1L)
  }
  n_cores <- min(as.integer(n_cores), length(chunks))
  if (n_cores <= 1L) {
    return(lapply(chunks, eval_fn))
  }

  res <- tryCatch(
    {
      is_darwin <- identical(Sys.info()[["sysname"]], "Darwin")
      if (.Platform$OS.type == "unix" && !is_darwin) {
        parallel::mclapply(chunks, eval_fn, mc.cores = n_cores)
      } else {
        # On Darwin and Windows, use socket cluster or sequential execution
        cl <- parallel::makeCluster(n_cores)
        on.exit(parallel::stopCluster(cl), add = TRUE)
        parallel::parLapply(cl, chunks, eval_fn)
      }
    },
    error = function(e) NULL
  )

  if (is.null(res) || any(sapply(res, function(r) is.null(r) || inherits(r, "try-error")))) {
    res <- lapply(chunks, eval_fn)
  }

  res
}

# -----------------------------------------------------------------------------
# Coordinate Grid Construction
# -----------------------------------------------------------------------------

#' Generate a regular n-dimensional coordinate grid over specified limits
#'
#' @param limits Named list of length-2 numeric vectors specifying min and max
#'   for each cue dimension.
#' @param n_points Integer vector or scalar specifying grid resolution along
#'   each dimension. Defaults to 100L.
#' @param cues Character vector of cue names. Defaults to `names(limits)`.
#' @return A data frame containing grid coordinates across all cue dimensions.
#' @noRd
#' @keywords internal
.mvbu_make_grid <- function(limits, n_points = 100L, cues = NULL) {
  if (is.null(cues)) {
    cues <- names(limits)
  }
  D <- length(cues)
  if (length(n_points) == 1L) {
    n_points <- rep(as.integer(n_points), D)
  }

  seq_list <- vector("list", D)
  names(seq_list) <- cues
  for (i in seq_len(D)) {
    c_name <- cues[i]
    c_lim <- limits[[c_name]]
    seq_list[[c_name]] <- seq(
      from = c_lim[1],
      to = c_lim[2],
      length.out = n_points[i]
    )
  }

  if (D == 1L) {
    df <- data.frame(seq_list[[1L]])
    names(df) <- cues[1L]
    return(df)
  }

  expand.grid(seq_list, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
}

#' Generate a sliced grid for 3-cue spaces
#'
#' @param limits Named list of limits.
#' @param cues Character vector of length 3.
#' @param slice_cue Name of the dimension to slice. Defaults to `cues[3L]`.
#' @param slice_values Numeric vector of values along `slice_cue`. If `NULL`,
#'   generates 3 quantiles (e.g. 25%, 50%, 75%).
#' @param n_points Resolution in the other two dimensions. Defaults to 80L.
#' @return A data frame containing grid points and slice labels.
#' @noRd
#' @keywords internal
.mvbu_make_sliced_grid <- function(
  limits,
  cues,
  slice_cue = NULL,
  slice_values = NULL,
  n_points = 80L
) {
  if (is.null(slice_cue)) {
    slice_cue <- cues[3L]
  }
  active_cues <- setdiff(cues, slice_cue)

  if (is.null(slice_values)) {
    s_lim <- limits[[slice_cue]]
    slice_values <- seq(from = s_lim[1], to = s_lim[2], length.out = 3L)
  }

  base_grid <- .mvbu_make_grid(
    limits = limits[active_cues],
    n_points = n_points,
    cues = active_cues
  )

  sliced_dfs <- vector("list", length(slice_values))
  for (i in seq_along(slice_values)) {
    val <- slice_values[i]
    df <- base_grid
    df[[slice_cue]] <- val
    df$slice_label <- sprintf("%s = %.2f", slice_cue, val)
    df$slice_value <- val
    sliced_dfs[[i]] <- df
  }

  do.call(rbind, sliced_dfs)
}

# -----------------------------------------------------------------------------
# Density and Likelihood Evaluation
# -----------------------------------------------------------------------------

#' Evaluate category densities across a grid
#'
#' @param x An S7 representation, template, or cognitive model object.
#' @param grid Data frame containing grid points.
#' @param cues Cue column names.
#' @param parallel Logical; execute in parallel. Defaults to `FALSE`.
#' @param n_cores Number of cores for parallel execution.
#' @return A long data frame with columns `cue1`, ..., `category`, `density`.
#' @noRd
#' @keywords internal
.mvbu_eval_category_densities <- function(
  x,
  grid,
  cues = NULL,
  parallel = FALSE,
  n_cores = NULL
) {
  if (is.null(cues)) {
    cues <- get_cue_labels(x)
  }
  mat <- as.matrix(grid[cues])
  cat_labels <- get_category_labels(x)

  # Chunking for optional parallel execution
  n_rows <- nrow(mat)
  if (isTRUE(parallel) && n_rows > 500L) {
    chunk_size <- ceiling(n_rows / (n_cores %||% 2L))
    chunk_indices <- split(
      seq_len(n_rows),
      ceiling(seq_len(n_rows) / chunk_size)
    )
    chunk_mats <- lapply(chunk_indices, function(idx) mat[idx, , drop = FALSE])

    eval_chunk <- function(m_chunk) {
      likelihood(x, m_chunk)
    }
    lik_res <- .mvbu_parallel_eval(
      chunks = chunk_mats,
      eval_fn = eval_chunk,
      parallel = parallel,
      n_cores = n_cores
    )
    lik_mat <- do.call(rbind, lik_res)
  } else {
    lik_mat <- likelihood(x, mat)
  }

  if (is.vector(lik_mat) && length(cat_labels) == 1L) {
    lik_mat <- matrix(lik_mat, ncol = 1L, dimnames = list(NULL, cat_labels))
  }

  res_list <- vector("list", length(cat_labels))
  for (j in seq_along(cat_labels)) {
    cat_nm <- cat_labels[j]
    df_cat <- grid
    df_cat$category <- factor(cat_nm, levels = cat_labels)
    df_cat$density <- lik_mat[, j]
    res_list[[j]] <- df_cat
  }

  do.call(rbind, res_list)
}

# -----------------------------------------------------------------------------
# Categorization Posterior Evaluation
# -----------------------------------------------------------------------------

#' Evaluate posterior probabilities and decision choices across a grid
#'
#' @param x An S7 cognitive model object.
#' @param grid Data frame containing grid points.
#' @param cues Cue column names.
#' @param decision_rule Decision rule override.
#' @param parallel Logical; execute in parallel.
#' @param n_cores Number of cores.
#' @return A long data frame with posterior probabilities and chosen category.
#' Apply decision rule to a posterior probability matrix
#'
#' Converts posterior beliefs P(c | x) to choice probabilities P(respond c | x)
#' under the specified decision rule ("proportional", "criterion", or "sampling").
#'
#' @param post_mat Numeric matrix of posterior probabilities (observations x categories).
#' @param decision_rule Character string: "proportional", "criterion", or "sampling".
#' @param lapse_rate Optional numeric scalar lapse rate in \code{[0, 1]}.
#' @param lapse_bias Optional numeric vector of lapse biases summing to 1.
#' @return Numeric matrix of choice probabilities with identical dimensions and column names.
#' @noRd
#' @keywords internal
.apply_decision_rule_to_posteriors <- function(
  post_mat,
  decision_rule = "proportional",
  lapse_rate = 0,
  lapse_bias = NULL
) {
  drule <- decision_rule %||% "proportional"
  if (identical(drule, "proportional") || identical(drule, "sampling")) {
    return(post_mat)
  }
  if (identical(drule, "criterion")) {
    n_rows <- nrow(post_mat)
    n_cols <- ncol(post_mat)
    win_idx <- max.col(post_mat, ties.method = "first")
    resp_mat <- matrix(0, nrow = n_rows, ncol = n_cols, dimnames = dimnames(post_mat))
    resp_mat[cbind(seq_len(n_rows), win_idx)] <- 1
    if (!is.null(lapse_rate) && lapse_rate > 0 && !is.null(lapse_bias) && length(lapse_bias) == n_cols) {
      resp_mat <- (1 - lapse_rate) * resp_mat + lapse_rate * matrix(lapse_bias, nrow = n_rows, ncol = n_cols, byrow = TRUE)
    }
    return(resp_mat)
  }
  post_mat
}

#' Evaluate categorization response probabilities over a grid
#'
#' @param x An S7 cognitive model.
#' @param grid A data frame of grid coordinates.
#' @param cues Character vector of cue names. Defaults to `get_cue_labels(x)`.
#' @param decision_rule Character; decision rule (e.g. `"proportional"`, `"criterion"`).
#' @param parallel Logical; whether to evaluate in parallel across chunks.
#' @param n_cores Integer; number of cores for parallel execution.
#' @param noise_treatment Optional character; noise treatment override.
#' @param lapse_treatment Optional character; lapse treatment override.
#' @return A tidy data frame with columns `cues`, `category`, and `posterior`.
#' @noRd
#' @keywords internal
.mvbu_eval_categorization_posteriors <- function(
  x,
  grid,
  cues = NULL,
  decision_rule = NULL,
  parallel = FALSE,
  n_cores = NULL,
  noise_treatment = NULL,
  lapse_treatment = NULL
) {
  if (is.null(cues)) {
    cues <- get_cue_labels(x)
  }
  mat <- as.matrix(grid[cues])
  cat_labels <- get_category_labels(x)
  drule <- decision_rule %||% (if (S7::S7_inherits(x, MVBU_CognitiveModel)) x@decision_rule else "proportional")

  # Obtain posterior function, respecting noise_treatment and lapse_treatment overrides
  pf <- if (S7::S7_inherits(x, MVBU_CognitiveModel)) {
    get_category_posterior_function(
      x,
      noise_treatment = noise_treatment %||% get_noise_treatment(x),
      lapse_treatment = lapse_treatment %||% get_lapse_treatment(x)
    )
  } else {
    function(m, categories = NULL) posterior(x, m, categories = categories)
  }

  compute_post <- function(m) {
    pf(m, categories = cat_labels)
  }

  n_rows <- nrow(mat)
  if (isTRUE(parallel) && n_rows > 500L) {
    chunk_size <- ceiling(n_rows / (n_cores %||% 2L))
    chunk_indices <- split(
      seq_len(n_rows),
      ceiling(seq_len(n_rows) / chunk_size)
    )
    chunk_mats <- lapply(chunk_indices, function(idx) mat[idx, , drop = FALSE])

    post_res <- .mvbu_parallel_eval(
      chunks = chunk_mats,
      eval_fn = compute_post,
      parallel = parallel,
      n_cores = n_cores
    )
    post_mat <- do.call(rbind, post_res)
  } else {
    post_mat <- compute_post(mat)
  }

  # Effective lapse rate and bias depending on lapse_treatment
  effective_lapse_rate <- if (S7::S7_inherits(x, MVBU_CognitiveModel)) {
    eff_lapse_trt <- lapse_treatment %||% get_lapse_treatment(x)
    if (identical(eff_lapse_trt, "no_lapses")) 0 else x@lapse_behavior$lapse_rate
  } else 0
  effective_lapse_bias <- if (S7::S7_inherits(x, MVBU_CognitiveModel)) {
    eff_lapse_trt <- lapse_treatment %||% get_lapse_treatment(x)
    if (identical(eff_lapse_trt, "no_lapses")) NULL else x@lapse_behavior$lapse_bias
  } else NULL

  resp_mat <- .apply_decision_rule_to_posteriors(
    post_mat,
    decision_rule = drule,
    lapse_rate = effective_lapse_rate,
    lapse_bias = effective_lapse_bias
  )

  res_list <- vector("list", length(cat_labels))
  for (j in seq_along(cat_labels)) {
    cat_nm <- cat_labels[j]
    df_cat <- grid
    df_cat$category <- factor(cat_nm, levels = cat_labels)
    df_cat$posterior <- resp_mat[, j]
    res_list[[j]] <- df_cat
  }
  dplyr::bind_rows(res_list)
}

# -----------------------------------------------------------------------------
# Geom and Layer Construction Helpers
# -----------------------------------------------------------------------------

#' Construct 1D category density layers
#'
#' @param data Long data frame with cue column, `category`, and `density`.
#' @param cue Name of cue column.
#' @param aes Aesthetic mode: `"both"`, `"fill"`, or `"contour"`.
#' @param alpha Fill alpha transparency.
#' @param line_size Line width.
#' @param linetype Line type.
#' @return A list of ggplot layers.
#' @noRd
#' @keywords internal
.mvbu_build_1d_density_layers <- function(
  data,
  cue,
  aes = "both",
  alpha = 0.25,
  line_size = 0.8,
  linetype = 1
) {
  layers <- list()
  if (aes %in% c("both", "fill")) {
    layers <- c(
      layers,
      list(
        ggplot2::geom_ribbon(
          data = data,
          mapping = ggplot2::aes(
            x = .data[[cue]],
            ymin = 0,
            ymax = .data$density,
            fill = .data$category
          ),
          alpha = alpha,
          inherit.aes = FALSE
        )
      )
    )
  }
  if (aes %in% c("both", "contour")) {
    layers <- c(
      layers,
      list(
        ggplot2::geom_line(
          data = data,
          mapping = ggplot2::aes(
            x = .data[[cue]],
            y = .data$density,
            color = .data$category
          ),
          linewidth = line_size,
          linetype = linetype,
          inherit.aes = FALSE
        )
      )
    )
  }
  layers
}

#' Construct 2D category density layers
#'
#' @param data Long data frame with cue columns, `category`, and `density`.
#' @param cues Names of cue columns (length 2).
#' @param aes Aesthetic mode: `"contour"`, `"fill"`, `"both"`,
#'   `"fill_gradient"`, or `"discrete"`.
#' @param alpha Transparency.
#' @param bins Number of contour bins.
#' @param level Ellipse/contour probability mass.
#' @return A list of ggplot layers.
#' @noRd
#' @keywords internal
.mvbu_build_2d_density_layers <- function(
  data,
  cues,
  aes = "contour",
  alpha = 0.25,
  bins = 6L,
  level = 0.95
) {
  layers <- list()
  c1 <- cues[1L]
  c2 <- cues[2L]

  if (aes %in% c("fill_gradient", "discrete", "fill", "both")) {
    layers <- c(
      layers,
      list(
        ggplot2::geom_contour_filled(
          data = data,
          mapping = ggplot2::aes(
            x = .data[[c1]],
            y = .data[[c2]],
            z = .data$density,
            fill = .data$category
          ),
          alpha = alpha,
          bins = bins,
          inherit.aes = FALSE
        )
      )
    )
  }

  if (aes %in% c("contour", "both")) {
    layers <- c(
      layers,
      list(
        ggplot2::geom_contour(
          data = data,
          mapping = ggplot2::aes(
            x = .data[[c1]],
            y = .data[[c2]],
            z = .data$density,
            color = .data$category
          ),
          bins = bins,
          inherit.aes = FALSE
        )
      )
    )
  }

  layers
}

#' Render multi-panel 2D slices for 3-cue category plots
#'
#' @param reps List of category representations.
#' @param cues Character vector of length 3.
#' @param slice_cue Dimension along which to slice (defaults to 3rd cue).
#' @param slice_values Specific values along slice_cue.
#' @param aes Aesthetic mode.
#' @param levels Contour levels.
#' @param limits Named list of limits.
#' @param resolution Grid resolution.
#' @param parallel Logical; evaluate in parallel.
#' @param n_cores Number of cores.
#' @param t_sub Plot title and subtitle list.
#' @return A ggplot object.
#' @noRd
#' @keywords internal
.render_3D_sliced_category_plot <- function(
  reps,
  cues,
  slice_cue = NULL,
  slice_values = NULL,
  slices = NULL,
  aes = "contour",
  levels = NULL,
  limits = NULL,
  resolution = 60L,
  parallel = FALSE,
  n_cores = NULL,
  t_sub = NULL,
  ...
) {
  if (is.null(aes) || length(aes) == 0L) {
    aes <- "contour"
  }
  if (is.null(slice_values) && !is.null(slices)) {
    slice_values <- slices
  }
  if (is.null(slice_cue)) {
    slice_cue <- cues[3L]
  }
  active_cues <- setdiff(cues, slice_cue)

  if (is.null(limits)) {
    limits <- list()
    for (c_name in cues) {
      min_v <- Inf
      max_v <- -Inf
      for (cat_name in names(reps)) {
        r <- reps[[cat_name]]
        mc <- .get_rep_mean_and_cov(r)
        idx <- match(c_name, get_cue_labels(r))
        mu_val <- mc$mu[idx]
        sd_val <- sqrt(max(mc$Sigma[idx, idx], 1e-6))
        min_v <- min(min_v, mu_val - 3 * sd_val)
        max_v <- max(max_v, mu_val + 3 * sd_val)
      }
      limits[[c_name]] <- c(min_v, max_v)
    }
  }

  if (is.null(slice_values)) {
    slice_mus <- sort(vapply(reps, function(r) {
      mc <- .get_rep_mean_and_cov(r)
      idx <- match(slice_cue, get_cue_labels(r))
      mc$mu[idx]
    }, numeric(1L)))
    slice_sds <- vapply(reps, function(r) {
      mc <- .get_rep_mean_and_cov(r)
      idx <- match(slice_cue, get_cue_labels(r))
      sqrt(max(mc$Sigma[idx, idx], 1e-6))
    }, numeric(1L))

    if (!is.null(slices) && is.numeric(slices)) {
      if (all(slices > 0 & slices < 1)) {
        bar_mu <- mean(slice_mus)
        bar_sd <- mean(slice_sds)
        slice_values <- stats::qnorm(slices, mean = bar_mu, sd = bar_sd)
      } else {
        slice_values <- slices
      }
    } else {
      bar_mu <- mean(slice_mus)
      bar_sd <- mean(slice_sds)
      slice_values <- bar_mu + c(-2, -1, 0, 1, 2) * bar_sd
    }
  }

  sliced_grid <- .mvbu_make_sliced_grid(
    limits = limits,
    cues = cues,
    slice_cue = slice_cue,
    slice_values = slice_values,
    n_points = resolution
  )

  temp_template <- new_category_representation_template(
    representations = reps
  )
  dens_df <- .mvbu_eval_category_densities(
    x = temp_template,
    grid = sliced_grid,
    cues = cues,
    parallel = parallel,
    n_cores = n_cores
  )

  all_c <- names(reps)
  c_colors <- scales::hue_pal()(length(all_c))
  names(c_colors) <- all_c

  p <- ggplot2::ggplot(
    dens_df,
    ggplot2::aes(
      x = .data[[active_cues[1L]]],
      y = .data[[active_cues[2L]]]
    )
  )

  # Fill layer
  has_fill <- "fill-gradient" %in% aes || "fill-discrete" %in% aes || "fill" %in% aes
  if (has_fill) {
    dens_df_rel <- dens_df %>%
      dplyr::group_by(.data$category, .data$slice_label) %>%
      dplyr::mutate(
        rel_dens = .data$density / max(.data$density, 1e-12)
      ) %>%
      dplyr::ungroup()

    p <- p + ggplot2::geom_tile(
      data = dens_df_rel,
      ggplot2::aes(
        fill = .data$category,
        alpha = .data$rel_dens
      ),
      show.legend = c(fill = TRUE, alpha = FALSE)
    ) +
      ggplot2::scale_fill_manual(values = c_colors, name = "Category") +
      ggplot2::scale_alpha_continuous(range = c(0, 0.85), limits = c(0, 1), guide = "none")
  }

  # Contour layer (contour lines corresponding to 3D central confidence levels)
  has_contour <- "contour" %in% aes
  if (has_contour) {
    if (is.null(levels)) {
      levels <- 2 * stats::pnorm(1:3) - 1 # 0.683, 0.954, 0.997
    }
    levels <- sort(levels)

    path_dfs <- list()
    label_dfs <- list()

    for (cat_name in all_c) {
      r <- reps[[cat_name]]
      is_ex <- S7::S7_inherits(r, Exemplar_CategoryRepresentation)

      if (!is_ex) {
        mc <- .get_rep_mean_and_cov(r)
        idx_c <- match(cues, get_cue_labels(r))
        Sigma <- mc$Sigma[idx_c, idx_c]
        det_S <- det(Sigma)
        c_vals <- stats::qchisq(levels, df = 3)
        z_crits <- (2 * pi)^(-1.5) * (det_S)^(-0.5) * exp(-0.5 * c_vals)
      }

      for (s_lbl in unique(dens_df$slice_label)) {
        sub_df <- dens_df[dens_df$category == cat_name & dens_df$slice_label == s_lbl, ]
        gx <- sort(unique(sub_df[[active_cues[1L]]]))
        gy <- sort(unique(sub_df[[active_cues[2L]]]))
        z_mat <- matrix(sub_df$density, nrow = length(gx), ncol = length(gy))

        if (is_ex) {
          d_vals <- sort(sub_df$density, decreasing = TRUE)
          c_mass <- cumsum(d_vals) / max(sum(d_vals), 1e-12)
          brks <- sapply(levels, function(lvl) d_vals[which.min(abs(c_mass - lvl))])
        } else {
          brks <- z_crits
        }

        for (k in seq_along(levels)) {
          lvl <- levels[k]
          brk <- brks[k]
          if (is.na(brk) || brk <= 0 || max(z_mat, na.rm = TRUE) < brk) next

          cl <- tryCatch(
            grDevices::contourLines(x = gx, y = gy, z = z_mat, levels = brk),
            error = function(e) list()
          )
          if (length(cl) == 0L) next

          # Linear line width gradation: thickest at central mass near 0, thinnest at 100%
          lw <- 0.85 - 0.60 * lvl
          lbl_text <- sprintf("%d%%", round(lvl * 100))

          for (j in seq_along(cl)) {
            line_x <- cl[[j]]$x
            line_y <- cl[[j]]$y
            path_dfs[[length(path_dfs) + 1L]] <- data.frame(
              x = line_x,
              y = line_y,
              category = cat_name,
              slice_label = s_lbl,
              group_id = paste(cat_name, s_lbl, k, j, sep = "."),
              linewidth = lw,
              stringsAsFactors = FALSE
            )
            idx_top <- which.max(line_y)
            label_dfs[[length(label_dfs) + 1L]] <- data.frame(
              x = line_x[idx_top],
              y = line_y[idx_top],
              label = lbl_text,
              category = cat_name,
              slice_label = s_lbl,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }

    if (length(path_dfs) > 0L) {
      all_paths <- do.call(rbind, path_dfs)
      names(all_paths)[1:2] <- active_cues[1:2]
      p <- p + ggplot2::geom_path(
        data = all_paths,
        ggplot2::aes(
          x = .data[[active_cues[1L]]],
          y = .data[[active_cues[2L]]],
          group = .data$group_id,
          color = .data$category,
          linewidth = .data$linewidth
        )
      ) +
        ggplot2::scale_linewidth_identity()
    }

    if (length(label_dfs) > 0L) {
      all_labels <- do.call(rbind, label_dfs)
      names(all_labels)[1:2] <- active_cues[1:2]
      p <- p + ggplot2::geom_text(
        data = all_labels,
        ggplot2::aes(
          x = .data[[active_cues[1L]]],
          y = .data[[active_cues[2L]]],
          label = .data$label,
          color = .data$category
        ),
        size = 2.4,
        vjust = -0.3,
        show.legend = FALSE
      )
    }

    p <- p + ggplot2::scale_color_manual(values = c_colors, name = "Category")
  }

  if (has_fill && has_contour) {
    p <- p + ggplot2::guides(
      color = ggplot2::guide_legend(
        title = "Category",
        override.aes = list(fill = c_colors, alpha = 0.5)
      ),
      fill = "none"
    )
  } else if (has_fill) {
    p <- p + ggplot2::guides(
      color = "none",
      fill = ggplot2::guide_legend(title = "Category")
    )
  } else {
    p <- p + ggplot2::guides(
      fill = "none",
      color = ggplot2::guide_legend(title = "Category")
    )
  }

  p <- p +
    ggplot2::facet_wrap(~slice_label) +
    .mvbu_theme() +
    ggplot2::labs(
      title = t_sub$title %||% "Category likelihood",
      subtitle = t_sub$subtitle %||% "",
      x = active_cues[1L],
      y = active_cues[2L]
    )

  p
}

#' Render interactive 2D density surface WebGL plot via Plotly
#'
#' @param reps Named list of category representations.
#' @param cues Character vector of length 2.
#' @param aes Aesthetic layer type ("contour", "fill-discrete", "fill-gradient", or "fill").
#' @param levels Central probability masses.
#' @param limits Named list of limits or vector.
#' @param n_exemplars Number of exemplars to sample for exemplar representations.
#' @param resolution Grid density (default 60L).
#' @param t_sub Plot title and subtitle list.
#' @param ... Additional arguments.
#' @return A plotly htmlwidget object.
#' @noRd
#' @keywords internal
.render_2D_interactive_category_plot <- function(
  reps,
  cues,
  aes = NULL,
  levels = NULL,
  limits = NULL,
  n_exemplars = 0L,
  resolution = 60L,
  t_sub = NULL,
  ...
) {
  rlang::check_installed("plotly", reason = "for interactive category plots")

  if (is.null(aes) || length(aes) == 0L) {
    aes <- "fill-discrete"
  }
  if ("fill" %in% aes && !"fill-gradient" %in% aes && !"fill-discrete" %in% aes) {
    aes <- unique(c(setdiff(aes, "fill"), "fill-discrete"))
  }

  all_c <- names(reps)
  c_colors <- scales::hue_pal()(length(all_c))
  names(c_colors) <- all_c

  min_x <- Inf
  max_x <- -Inf
  min_y <- Inf
  max_y <- -Inf

  for (cat_name in names(reps)) {
    r <- reps[[cat_name]]
    if (S7::S7_inherits(r, Exemplar_CategoryRepresentation)) {
      mat <- as.matrix(r@exemplars)
      rng1 <- diff(range(mat[, cues[1L]]))
      rng2 <- diff(range(mat[, cues[2L]]))
      min_x <- min(min_x, min(mat[, cues[1L]]) - 0.25 * rng1)
      max_x <- max(max_x, max(mat[, cues[1L]]) + 0.25 * rng1)
      min_y <- min(min_y, min(mat[, cues[2L]]) - 0.25 * rng2)
      max_y <- max(max_y, max(mat[, cues[2L]]) + 0.25 * rng2)
    } else {
      mc <- .get_rep_mean_and_cov(r)
      idx_c <- match(cues, get_cue_labels(r))
      mu_vec <- mc$mu[idx_c]
      sd_x <- sqrt(max(mc$Sigma[idx_c[1L], idx_c[1L]], 1e-6))
      sd_y <- sqrt(max(mc$Sigma[idx_c[2L], idx_c[2L]], 1e-6))
      min_x <- min(min_x, mu_vec[1L] - 3.5 * sd_x)
      max_x <- max(max_x, mu_vec[1L] + 3.5 * sd_x)
      min_y <- min(min_y, mu_vec[2L] - 3.5 * sd_y)
      max_y <- max(max_y, mu_vec[2L] + 3.5 * sd_y)
    }
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

  res <- as.integer(max(resolution %||% 60L, 20L))
  gx <- seq(min_x, max_x, length.out = res)
  gy <- seq(min_y, max_y, length.out = res)
  grid_df <- expand.grid(c1 = gx, c2 = gy)
  names(grid_df) <- cues[1:2]

  p <- plotly::plot_ly()

  for (cat_name in all_c) {
    r <- reps[[cat_name]]
    col <- c_colors[cat_name]
    col_rgb <- grDevices::col2rgb(col)

    dens <- likelihood(r, as.matrix(grid_df[cues]))

    z_mat <- matrix(dens, nrow = length(gx), ncol = length(gy))
    z_surf <- t(z_mat)

    if ("contour" %in% aes && !("fill-discrete" %in% aes) && !("fill-gradient" %in% aes)) {
      p <- p %>% plotly::add_surface(
        x = gx,
        y = gy,
        z = z_surf,
        name = cat_name,
        hidesurface = TRUE,
        contours = list(
          z = list(
            show = TRUE,
            usecolormap = FALSE,
            color = col,
            project = list(z = FALSE),
            width = 3
          )
        ),
        showscale = FALSE,
        hoverinfo = "name"
      )
    } else if ("fill-gradient" %in% aes) {
      c_start <- sprintf("rgba(%d,%d,%d,0.25)", col_rgb[1], col_rgb[2], col_rgb[3])
      c_end <- sprintf("rgba(%d,%d,%d,0.95)", col_rgb[1], col_rgb[2], col_rgb[3])
      p <- p %>% plotly::add_surface(
        x = gx,
        y = gy,
        z = z_surf,
        name = cat_name,
        opacity = 0.90,
        colorscale = list(c(0, c_start), c(1, c_end)),
        showscale = FALSE,
        contours = if ("contour" %in% aes) {
          list(z = list(show = TRUE, color = "#222222", width = 2))
        } else {
          list()
        },
        hoverinfo = "name"
      )
    } else {
      # "fill-discrete" (default): constant opaqueness per category
      p <- p %>% plotly::add_surface(
        x = gx,
        y = gy,
        z = z_surf,
        name = cat_name,
        opacity = 0.85,
        colorscale = list(c(0, col), c(1, col)),
        showscale = FALSE,
        contours = if ("contour" %in% aes) {
          list(z = list(show = TRUE, color = "#222222", width = 2))
        } else {
          list()
        },
        hoverinfo = "name"
      )
    }

    if (S7::S7_inherits(r, Exemplar_CategoryRepresentation)) {
      if (n_exemplars > 0L) {
        mat <- as.matrix(r@exemplars)
        n_pts <- min(as.integer(n_exemplars), nrow(mat))
        idx_pts <- sample.int(nrow(mat), n_pts, replace = FALSE)
        sub_pts <- mat[idx_pts, cues, drop = FALSE]
        pts_z <- .eval_exemplar_density_grid(r, as.data.frame(sub_pts), cues)
        p <- p %>% plotly::add_markers(
          x = sub_pts[, 1],
          y = sub_pts[, 2],
          z = pts_z,
          name = cat_name,
          legendgroup = cat_name,
          marker = list(size = 3, color = col, opacity = 0.5),
          hoverinfo = "name"
        )
      }
    }
  }

  main_title <- t_sub$title %||% "2D Category Density Surfaces"
  sub_title <- t_sub$subtitle %||% ""

  p <- p %>% plotly::layout(
    title = list(
      text = if (nzchar(sub_title)) {
        sprintf("<b>%s</b><br><span style='font-size:12px;color:gray;'>%s</span>", main_title, sub_title)
      } else {
        sprintf("<b>%s</b>", main_title)
      }
    ),
    scene = list(
      xaxis = list(title = cues[1L]),
      yaxis = list(title = cues[2L]),
      zaxis = list(title = "Probability Density"),
      camera = list(
        eye = list(x = 1.5, y = 1.5, z = 1.2)
      )
    ),
    legend = list(orientation = "h", x = 0.1, y = -0.1)
  )

  p
}

#' Render interactive 3D WebGL categorization surface via Plotly
#'
#' @param resp_mat Matrix of response probabilities (rows: grid points, columns: categories).
#' @param gx Numeric vector of x-coordinates (cue 1).
#' @param gy Numeric vector of y-coordinates (cue 2).
#' @param cues Character vector of length 2.
#' @param categories Character vector of category names to display.
#' @param cat_colors Named vector of colors for each category.
#' @param aes Aesthetic layer type ("fill-discrete", "fill-gradient", "contour").
#' @param t_sub Plot title and subtitle list.
#' @param ... Additional arguments.
#' @return A plotly htmlwidget object.
#' @noRd
#' @keywords internal
.render_2D_interactive_categorization_plot <- function(
  resp_mat,
  gx,
  gy,
  cues,
  categories,
  cat_colors,
  aes = NULL,
  t_sub = NULL,
  ...
) {
  rlang::check_installed("plotly", reason = "for interactive categorization plots")

  if (is.null(aes) || length(aes) == 0L) {
    aes <- "fill-discrete"
  }
  if ("fill" %in% aes && !"fill-gradient" %in% aes && !"fill-discrete" %in% aes) {
    aes <- unique(c(setdiff(aes, "fill"), "fill-discrete"))
  }

  p <- plotly::plot_ly()

  for (cat_name in categories) {
    if (!cat_name %in% colnames(resp_mat)) next
    col <- cat_colors[cat_name]
    col_rgb <- grDevices::col2rgb(col)
    prob_vec <- resp_mat[, cat_name]
    z_mat <- matrix(prob_vec, nrow = length(gx), ncol = length(gy))
    z_surf <- t(z_mat)

    if ("contour" %in% aes && !("fill-discrete" %in% aes) && !("fill-gradient" %in% aes)) {
      p <- p %>% plotly::add_surface(
        x = gx,
        y = gy,
        z = z_surf,
        name = cat_name,
        hidesurface = TRUE,
        contours = list(
          z = list(
            show = TRUE,
            usecolormap = FALSE,
            color = col,
            project = list(z = FALSE),
            width = 3
          )
        ),
        showscale = FALSE,
        hoverinfo = "name+z"
      )
    } else if ("fill-gradient" %in% aes) {
      c_start <- sprintf("rgba(%d,%d,%d,0.20)", col_rgb[1], col_rgb[2], col_rgb[3])
      c_end <- sprintf("rgba(%d,%d,%d,0.90)", col_rgb[1], col_rgb[2], col_rgb[3])
      p <- p %>% plotly::add_surface(
        x = gx,
        y = gy,
        z = z_surf,
        name = cat_name,
        opacity = 0.90,
        colorscale = list(c(0, c_start), c(1, c_end)),
        cmin = 0,
        cmax = 1,
        showscale = FALSE,
        contours = if ("contour" %in% aes) {
          list(z = list(show = TRUE, color = "#222222", width = 2))
        } else {
          list()
        },
        hoverinfo = "name+z"
      )
    } else {
      # "fill-discrete" (default): constant opaqueness per category
      p <- p %>% plotly::add_surface(
        x = gx,
        y = gy,
        z = z_surf,
        name = cat_name,
        opacity = 0.85,
        colorscale = list(c(0, col), c(1, col)),
        cmin = 0,
        cmax = 1,
        showscale = FALSE,
        contours = if ("contour" %in% aes) {
          list(z = list(show = TRUE, color = "#222222", width = 2))
        } else {
          list()
        },
        hoverinfo = "name+z"
      )
    }
  }

  main_title <- t_sub$title %||% "2D Categorization Function Surface"
  sub_title <- t_sub$subtitle %||% ""

  p <- p %>% plotly::layout(
    title = list(
      text = if (nzchar(sub_title)) {
        sprintf("<b>%s</b><br><span style='font-size:12px;color:gray;'>%s</span>", main_title, sub_title)
      } else {
        sprintf("<b>%s</b>", main_title)
      }
    ),
    scene = list(
      xaxis = list(title = cues[1L]),
      yaxis = list(title = cues[2L]),
      zaxis = list(title = "Response Probability", range = c(0, 1)),
      camera = list(
        eye = list(x = 1.5, y = 1.5, z = 1.2)
      )
    ),
    legend = list(orientation = "h", x = 0.1, y = -0.1)
  )

  p
}

#' Render interactive 3D WebGL category plot via Plotly
#'
#' @param reps Named list of category representations.
#' @param cues Character vector of length 3.
#' @param aes Aesthetic layer type ("scatter", "fill-discrete", "fill-gradient", "fill", or "contour").
#' @param levels Numeric vector of central probability masses (defaults to 2 sigma: ~0.954).
#' @param n_exemplars Number of exemplars to sample for exemplar representations.
#' @param opacity Base surface opacity for 3D ellipsoids (defaults to 0.35).
#' @param resolution Grid resolution for ellipsoid mesh.
#' @param t_sub Plot title and subtitle list.
#' @param ... Additional arguments.
#' @return A plotly htmlwidget object.
#' @noRd
#' @keywords internal
.render_3D_interactive_category_plot <- function(
  reps,
  cues,
  aes = NULL,
  levels = NULL,
  n_exemplars = 100L,
  opacity = 0.65,
  resolution = 25L,
  t_sub = NULL,
  ...
) {
  rlang::check_installed("plotly", reason = "for interactive 3D category plots")

  # Default levels for 3D ellipsis: 2 sigma (~0.9545)
  if (is.null(levels)) {
    levels <- 2 * stats::pnorm(2) - 1
  }
  levels <- sort(levels)

  if (is.null(aes) || length(aes) == 0L) {
    first_r <- reps[[1L]]
    aes <- if (S7::S7_inherits(first_r, Exemplar_CategoryRepresentation)) "scatter" else "fill-discrete"
  }
  if ("fill" %in% aes && !"fill-gradient" %in% aes && !"fill-discrete" %in% aes) {
    aes <- unique(c(setdiff(aes, "fill"), "fill-discrete"))
  }

  all_c <- names(reps)
  c_colors <- scales::hue_pal()(length(all_c))
  names(c_colors) <- all_c

  p <- plotly::plot_ly()

  for (cat_name in all_c) {
    r <- reps[[cat_name]]
    col <- c_colors[cat_name]

    if (S7::S7_inherits(r, Exemplar_CategoryRepresentation)) {
      mat <- as.matrix(r@exemplars)
      mu <- colMeans(mat[, cues, drop = FALSE])

      use_isosurface <- any(c("fill-discrete", "fill-gradient", "contour") %in% aes)

      if (use_isosurface) {
        # 3D KDE Isosurface evaluation (compute intensive)
        min_c <- apply(mat[, cues, drop = FALSE], 2, min)
        max_c <- apply(mat[, cues, drop = FALSE], 2, max)
        span_c <- pmax(max_c - min_c, 1e-3)
        res_3d <- 18L
        gx <- seq(min_c[1] - 0.1 * span_c[1], max_c[1] + 0.1 * span_c[1], length.out = res_3d)
        gy <- seq(min_c[2] - 0.1 * span_c[2], max_c[2] + 0.1 * span_c[2], length.out = res_3d)
        gz <- seq(min_c[3] - 0.1 * span_c[3], max_c[3] + 0.1 * span_c[3], length.out = res_3d)
        grid_3d <- expand.grid(x = gx, y = gy, z = gz)
        names(grid_3d) <- cues
        dens_3d <- .eval_exemplar_density_grid(r, grid_3d, cues)

        d_vals <- sort(dens_3d, decreasing = TRUE)
        c_mass <- cumsum(d_vals) / max(sum(d_vals), 1e-12)
        iso_val <- d_vals[which.min(abs(c_mass - levels[1L]))]

        p <- p %>% plotly::add_isosurface(
          x = grid_3d[[cues[1L]]],
          y = grid_3d[[cues[2L]]],
          z = grid_3d[[cues[3L]]],
          value = dens_3d,
          isomin = iso_val,
          isomax = max(dens_3d),
          surface = list(count = length(levels), fill = if ("contour" %in% aes) 0.1 else 0.7),
          opacity = if ("fill-gradient" %in% aes) 0.45 else 0.65,
          colorscale = list(c(0, col), c(1, col)),
          showscale = FALSE,
          name = cat_name,
          legendgroup = cat_name
        )
      } else {
        # Scatter mode (default)
        N <- nrow(mat)
        n_to_draw <- if (!is.null(n_exemplars) && n_exemplars > 0L) {
          min(as.integer(n_exemplars), N)
        } else {
          min(150L, N)
        }
        idx_s <- sample.int(N, n_to_draw, replace = FALSE)
        sub_mat <- mat[idx_s, cues, drop = FALSE]

        p <- p %>% plotly::add_markers(
          x = sub_mat[, 1],
          y = sub_mat[, 2],
          z = sub_mat[, 3],
          name = cat_name,
          legendgroup = cat_name,
          marker = list(
            size = 3,
            color = col,
            opacity = 0.5
          ),
          hoverinfo = "name"
        )
      }

      p <- p %>% plotly::add_markers(
        x = mu[1],
        y = mu[2],
        z = mu[3],
        name = sprintf("%s (mean)", cat_name),
        legendgroup = cat_name,
        showlegend = FALSE,
        marker = list(
          size = 6,
          color = col,
          symbol = "diamond"
        ),
        hoverinfo = "name"
      )
    } else {
      mc <- .get_rep_mean_and_cov(r)
      idx_c <- match(cues, get_cue_labels(r))
      mu <- mc$mu[idx_c]
      Sigma <- mc$Sigma[idx_c, idx_c]
      chol_S <- tryCatch(t(chol(Sigma)), error = function(e) diag(sqrt(pmax(diag(Sigma), 1e-6))))

      n_theta <- as.integer(max(resolution, 15L))
      n_phi <- as.integer(max(resolution, 15L))
      theta <- seq(0, 2 * pi, length.out = n_theta)
      phi <- seq(0, pi, length.out = n_phi)

      xs <- outer(cos(theta), sin(phi))
      ys <- outer(sin(theta), sin(phi))
      zs <- outer(rep(1, n_theta), cos(phi))
      unit_sphere <- rbind(as.vector(xs), as.vector(ys), as.vector(zs))

      for (k in seq_along(levels)) {
        lvl <- levels[k]
        r_val <- sqrt(stats::qchisq(lvl, df = 3))
        xyz_trans <- chol_S %*% (unit_sphere * r_val) + mu

        X <- matrix(xyz_trans[1, ], nrow = n_theta, ncol = n_phi)
        Y <- matrix(xyz_trans[2, ], nrow = n_theta, ncol = n_phi)
        Z <- matrix(xyz_trans[3, ], nrow = n_theta, ncol = n_phi)

        lvl_pct <- round(lvl * 100)

        # Opacity and surface style according to aes
        if ("contour" %in% aes && !("fill-discrete" %in% aes) && !("fill-gradient" %in% aes)) {
          p <- p %>% plotly::add_surface(
            x = X, y = Y, z = Z,
            name = cat_name,
            legendgroup = cat_name,
            showlegend = (k == 1L),
            hidesurface = TRUE,
            contours = list(
              x = list(show = TRUE, color = col, width = 2),
              y = list(show = TRUE, color = col, width = 2),
              z = list(show = TRUE, color = col, width = 2)
            ),
            showscale = FALSE,
            hoverinfo = "name"
          )
        } else if ("fill-gradient" %in% aes) {
          lvl_opacity <- opacity * (1.5 - 0.5 * (k / length(levels)))
          p <- p %>% plotly::add_surface(
            x = X, y = Y, z = Z,
            name = cat_name,
            legendgroup = cat_name,
            showlegend = (k == 1L),
            opacity = lvl_opacity,
            showscale = FALSE,
            colorscale = list(c(0, col), c(1, col)),
            contours = if ("contour" %in% aes) {
              list(z = list(show = TRUE, color = "#222222", width = 2))
            } else {
              list()
            },
            hoverinfo = "name"
          )
        } else {
          # "fill-discrete" (default): constant opaqueness per category
          p <- p %>% plotly::add_surface(
            x = X, y = Y, z = Z,
            name = cat_name,
            legendgroup = cat_name,
            showlegend = (k == 1L),
            opacity = opacity,
            showscale = FALSE,
            colorscale = list(c(0, col), c(1, col)),
            contours = if ("contour" %in% aes) {
              list(z = list(show = TRUE, color = "#222222", width = 2))
            } else {
              list()
            },
            hoverinfo = "name"
          )
        }
      }

      p <- p %>% plotly::add_markers(
        x = mu[1],
        y = mu[2],
        z = mu[3],
        name = sprintf("%s (mean)", cat_name),
        legendgroup = cat_name,
        showlegend = FALSE,
        marker = list(
          size = 6,
          color = col,
          symbol = "diamond"
        ),
        hoverinfo = "name"
      )
    }
  }

  main_title <- t_sub$title %||% "3D Categories"
  sub_title <- t_sub$subtitle %||% ""

  p <- p %>% plotly::layout(
    title = list(
      text = if (nzchar(sub_title)) {
        sprintf("<b>%s</b><br><span style='font-size:12px;color:gray;'>%s</span>", main_title, sub_title)
      } else {
        sprintf("<b>%s</b>", main_title)
      }
    ),
    scene = list(
      xaxis = list(title = cues[1L]),
      yaxis = list(title = cues[2L]),
      zaxis = list(title = cues[3L]),
      camera = list(
        eye = list(x = 1.6, y = 1.6, z = 1.3)
      )
    ),
    legend = list(orientation = "h", x = 0.1, y = -0.1)
  )

  p
}

#' Render multi-panel 2D slices for 3-cue categorization plots
#'
#' @param model An S7 cognitive model object.
#' @param cues Character vector of length 3.
#' @param categories Character vector of categories to plot.
#' @param aes Plot aesthetics.
#' @param slice_cue Dimension along which to slice (defaults to 3rd cue).
#' @param slice_values Specific values along slice_cue.
#' @param decision_rule Decision rule.
#' @param limits Named list of limits.
#' @param resolution Grid resolution.
#' @param parallel Logical; evaluate in parallel.
#' @param n_cores Number of cores.
#' @param t_sub Plot title and subtitle list.
#' @return A ggplot object.
#' @noRd
#' @keywords internal
.render_3D_sliced_categorization_plot <- function(
  model,
  cues,
  categories = NULL,
  aes = "contour",
  slice_cue = NULL,
  slice_values = NULL,
  slices = NULL,
  decision_rule = NULL,
  limits = NULL,
  resolution = 60L,
  parallel = FALSE,
  n_cores = NULL,
  t_sub = NULL,
  ...
) {
  if (is.null(aes) || length(aes) == 0L) {
    aes <- "contour"
  }
  if (is.null(slice_values) && !is.null(slices)) {
    slice_values <- slices
  }
  if (is.null(slice_cue)) {
    slice_cue <- cues[3L]
  }
  active_cues <- setdiff(cues, slice_cue)
  reps <- model@category_template@representations

  if (is.null(limits)) {
    limits <- list()
    for (c_name in cues) {
      min_v <- Inf
      max_v <- -Inf
      for (cat_name in names(reps)) {
        r <- reps[[cat_name]]
        mc <- .get_rep_mean_and_cov(r)
        idx <- match(c_name, get_cue_labels(r))
        mu_val <- mc$mu[idx]
        sd_val <- sqrt(max(mc$Sigma[idx, idx], 1e-6))
        min_v <- min(min_v, mu_val - 3 * sd_val)
        max_v <- max(max_v, mu_val + 3 * sd_val)
      }
      limits[[c_name]] <- c(min_v, max_v)
    }
  }

  if (is.null(slice_values)) {
    slice_mus <- sort(vapply(reps, function(r) {
      mc <- .get_rep_mean_and_cov(r)
      idx <- match(slice_cue, get_cue_labels(r))
      mc$mu[idx]
    }, numeric(1L)))
    if (length(slice_mus) == 3L) {
      slice_values <- unname(slice_mus)
    } else {
      s_lim <- limits[[slice_cue]]
      span <- s_lim[2] - s_lim[1]
      slice_values <- c(s_lim[1] + 0.25 * span, s_lim[1] + 0.50 * span, s_lim[1] + 0.75 * span)
    }
  }

  sliced_grid <- .mvbu_make_sliced_grid(
    limits = limits,
    cues = cues,
    slice_cue = slice_cue,
    slice_values = slice_values,
    n_points = resolution
  )

  post_df <- .mvbu_eval_categorization_posteriors(
    x = model,
    grid = sliced_grid,
    cues = cues,
    decision_rule = decision_rule,
    parallel = parallel,
    n_cores = n_cores
  )

  if (!is.null(categories) && length(categories) > 0L) {
    post_df <- post_df[post_df$category %in% categories, , drop = FALSE]
  }

  all_c <- unique(post_df$category)
  c_colors <- scales::hue_pal()(length(all_c))
  names(c_colors) <- all_c

  p <- ggplot2::ggplot(
    post_df,
    ggplot2::aes(
      x = .data[[active_cues[1L]]],
      y = .data[[active_cues[2L]]]
    )
  ) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(expand = FALSE)

  # Fill layer
  if ("fill-gradient" %in% aes) {
    p <- p + ggplot2::geom_tile(
      data = post_df,
      ggplot2::aes(
        fill = .data$category,
        alpha = .data$posterior
      ),
      show.legend = c(fill = TRUE, alpha = FALSE)
    ) +
      ggplot2::scale_fill_manual(values = c_colors, name = "Category") +
      ggplot2::scale_alpha_continuous(range = c(0, 0.85), limits = c(0, 1), guide = "none")
  }

  # Contour layer with labels
  if ("contour" %in% aes) {
    contour_alpha <- if (length(unique(post_df$category)) >= 2L) 0.5 else 1.0
    if (requireNamespace("metR", quietly = TRUE)) {
      p <- p + metR::geom_text_contour(
        data = post_df,
        ggplot2::aes(
          z = .data$posterior,
          color = .data$category
        ),
        breaks = c(0.25, 0.5, 0.75),
        stroke = 0.15,
        rotate = FALSE,
        size = 3,
        alpha = contour_alpha
      ) +
        ggplot2::geom_contour(
          data = post_df,
          ggplot2::aes(
            z = .data$posterior,
            color = .data$category
          ),
          breaks = c(0.25, 0.5, 0.75),
          linewidth = 0.4,
          alpha = contour_alpha
        )
    } else if (requireNamespace("geomtextpath", quietly = TRUE)) {
      p <- p + geomtextpath::geom_textcontour(
        data = post_df,
        ggplot2::aes(
          z = .data$posterior,
          color = .data$category
        ),
        breaks = c(0.25, 0.5, 0.75),
        size = 3,
        linewidth = 0.4,
        alpha = contour_alpha
      )
    } else {
      p <- p + ggplot2::geom_contour(
        data = post_df,
        ggplot2::aes(
          z = .data$posterior,
          color = .data$category
        ),
        breaks = c(0.25, 0.5, 0.75),
        linewidth = 0.4,
        alpha = contour_alpha
      )
    }
    p <- p + ggplot2::scale_color_manual(values = c_colors, name = "Category")
  }

  p <- p +
    ggplot2::facet_wrap(~slice_label) +
    .mvbu_theme() +
    ggplot2::labs(
      title = t_sub$title %||% "Categorization function",
      subtitle = t_sub$subtitle %||% "",
      x = active_cues[1L],
      y = active_cues[2L]
    )

  p
}
