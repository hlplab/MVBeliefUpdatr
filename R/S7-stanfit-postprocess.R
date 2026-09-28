#' @include S7-class.R
#' @include S7-stanfit.R
#' @include S7-stanfit-methods.R
#' @include internal-utils.R
#' @importFrom S7 new_generic method S7_dispatch
#' @importFrom loo loo waic loo_compare
#' @importFrom bayesplot pp_check
#' @export
bayesplot::pp_check
NULL

# -------------------------------------------------------------------------
# Internal Helper: .align_new_data_cues
# -------------------------------------------------------------------------
# Aligns columns of new_data with the model's full cue vector.
# This enables evaluating sub-cue slices (e.g. 1D slice of a 2D/3D model)
# without requiring new_data to contain all model cues.
#
# Arguments:
#   new_data: data.frame or numeric matrix of evaluation points.
#   cues: character vector of all cue names defined in the model.
#
# Returns:
#   A list containing:
#     - mat: numeric matrix with only the matched columns.
#     - c_idx: integer vector of 1-based indices indicating which model cues
#              are present in mat (used for sub-indexing mean vector m and
#              scatter matrix S).
# -------------------------------------------------------------------------
.align_new_data_cues <- function(new_data, cues) {
  if (is.data.frame(new_data) || (is.matrix(new_data) && !is.null(colnames(new_data)))) {
    avail_cues <- intersect(colnames(new_data), cues)
    if (length(avail_cues) > 0L) {
      c_idx <- match(avail_cues, cues)
      mat <- as.matrix(new_data[, avail_cues, drop = FALSE])
      return(list(mat = mat, c_idx = c_idx))
    }
  }
  mat <- as.matrix(new_data)
  if (ncol(mat) == length(cues)) {
    c_idx <- seq_along(cues)
  } else if (ncol(mat) == 1L) {
    c_idx <- 1L
  } else {
    .stop(
      "new_data dimensions do not match model cues (",
      paste(cues, collapse = ", "), ")"
    )
  }
  list(mat = mat, c_idx = c_idx)
}

# -------------------------------------------------------------------------
# Internal Helper: .evaluate_stan_category_density
# -------------------------------------------------------------------------
# Evaluates the exact analytical predictive density for category k given
# hyperparameters (m, kappa, nu, S) according to the model family:
#
#   1. NIW (Normal-Inverse-Wishart):
#      Predictive distribution is Multivariate Student-t (or univariate t for 1D slices).
#      Marginalizing out unobserved cues analytically preserves the Student-t form
#      with degrees of freedom df_t = nu - D_total + 1 and scale matrix
#      Sigma_scale = ((kappa + 1) / (kappa * df_t)) * S[c_idx, c_idx].
#
#   2. NIX (Normal-Inverse-Chi-squared):
#      Univariate Student-t predictive density with df = nu and scale
#      scale_eff = sqrt(sigma2 * (kappa + 1) / kappa).
#
#   3. MNIX (Multivariate Normal-Independent-Inverse-Chi-squared):
#      Product of independent univariate Student-t predictive densities across
#      the evaluated cues.
#
# Unsupported model families fail fast with an informative error (no fallbacks).
# -------------------------------------------------------------------------
.evaluate_stan_category_density <- function(
  x_mat,
  m_k,
  kappa_k,
  nu_k,
  S_k,
  sig_marg_k,
  c_idx,
  is_niw,
  is_nix,
  is_mnix,
  D_total
) {
  D_eval <- length(c_idx)
  m_sub <- m_k[c_idx]

  if (is_niw) {
    S_sub <- S_k[c_idx, c_idx, drop = FALSE]
    # For NIW predictive Student-t, marginalization preserves df = nu - D_total + 1
    df_t <- nu_k - D_total + 1
    if (df_t <= 0) df_t <- 1
    Sigma_scale <- ((kappa_k + 1) / (kappa_k * df_t)) * S_sub
    if (D_eval == 1L) {
      scale_eff <- sqrt(max(Sigma_scale[1L, 1L], 1e-8))
      z <- (x_mat[, 1L] - m_sub[1L]) / scale_eff
      return(stats::dt(z, df = df_t) / scale_eff)
    } else {
      Sigma_pd <- .make_pos_def(Sigma_scale)
      return(.dmvt_density(x_mat, mean = m_sub, Sigma = Sigma_pd, df = df_t, log = FALSE))
    }
  } else if (is_nix) {
    sig2 <- if (is.matrix(S_k)) S_k[1L, 1L] else as.numeric(S_k)[1L]
    scale_eff <- sqrt(max(sig2 * (kappa_k + 1) / kappa_k, 1e-8))
    z <- (x_mat[, 1L] - m_sub[1L]) / scale_eff
    return(stats::dt(z, df = max(nu_k, 1)) / scale_eff)
  } else if (is_mnix) {
    dens_prod <- rep(1, nrow(x_mat))
    for (d_i in seq_along(c_idx)) {
      idx_orig <- c_idx[d_i]
      m_d <- m_k[idx_orig]
      kap_d <- if (length(kappa_k) >= idx_orig) kappa_k[idx_orig] else kappa_k[1L]
      nu_d <- if (length(nu_k) >= idx_orig) nu_k[idx_orig] else nu_k[1L]
      s2_d <- if (is.matrix(S_k)) S_k[idx_orig, idx_orig] else S_k[idx_orig]
      scale_d <- sqrt(max(s2_d * (kap_d + 1) / kap_d, 1e-8))
      z_d <- (x_mat[, d_i] - m_d) / scale_d
      dens_prod <- dens_prod * (stats::dt(z_d, df = max(nu_d, 1)) / scale_d)
    }
    return(dens_prod)
  } else {
    .stop(
      "Predictive density evaluation is only supported for NIW, NIX, and MNIX models. ",
      "No analytical Student-t predictive density is defined for this model type."
    )
  }
}

#' Add Pre-extracted Posterior Latents and Evaluation Closures to a Stanfit Object
#'
#' Extracts and structures draw-level latent parameters (\eqn{m}, \eqn{S}, \eqn{\kappa},
#' \eqn{\nu}, \eqn{\Sigma_{\text{exp}}}, \eqn{\Sigma_{\text{marg}}}, and lapse rates)
#' and compiles unified likelihood and posterior evaluation closures. Pre-computed outputs
#' are stored strictly in \code{model@cache$draws} and \code{model@cache$summary}:
#' \itemize{
#'   \item \code{model@cache$draws}: Contains draw-level representations:
#'     \describe{
#'       \item{\code{parameters}}{A nested tibble of raw posterior draws across iterations/draws,
#'         including \eqn{m}, \eqn{S}, \eqn{\kappa}, \eqn{\nu}, \eqn{\Sigma_{\text{exp}}},
#'         \eqn{\Sigma_{\text{marg}}}, and optional \code{lapse_rate}.}
#'       \item{\code{likelihood_function}}{Closure \code{function(new_data, categories, groups, draws, noise_treatment)}
#'         evaluating draw-wise category likelihoods with model-family predictive densities
#'         (NIW: multivariate Student-t; NIX: univariate Student-t; MNIX: independent Student-t).}
#'       \item{\code{posterior_function}}{Closure \code{function(new_data, categories, groups, draws, noise_treatment, lapse_treatment)}
#'         evaluating draw-wise posterior category probabilities.}
#'     }
#'   \item \code{model@cache$summary}: Contains point estimates and across-draw summaries:
#'     \describe{
#'       \item{\code{parameters}}{Point estimates obtained by taking the mean across all posterior draws
#'         for each latent parameter (\eqn{m}, \eqn{S}, \eqn{\kappa}, \eqn{\nu}, \eqn{\Sigma_{\text{exp}}},
#'         \eqn{\Sigma_{\text{marg}}}, and \code{lapse_rate}) per group and category.}
#'       \item{\code{expected_moments}}{Distribution summary table (mean, standard deviation, and quantile
#'         credible intervals specified by \code{probs}) across posterior draws for the category mean \eqn{\mu}
#'         and expected category covariance \eqn{\Sigma_{\text{exp}} = S / (\nu - D - 1)}.}
#'       \item{\code{marginal_moments}}{Distribution summary table across posterior draws for the category
#'         mean \eqn{\mu} and marginal predictive covariance \eqn{\Sigma_{\text{marg}} = \frac{\kappa + 1}{\kappa} \Sigma_{\text{exp}}}.}
#'       \item{\code{likelihood_function}}{Closure evaluating category likelihoods from aggregated point estimates.}
#'       \item{\code{posterior_function}}{Closure evaluating category posterior probabilities from aggregated point estimates.}
#'     }
#' }
#'
#' @param model An \code{\link{MVBU_Stanfit}} object.
#' @param probs Numeric vector of probabilities for credible interval quantiles.
#'   (default: \code{c(0.025, 0.5, 0.975)})
#' @param ... Additional arguments passed to \code{\link{get_draws}}.
#'
#' @return The updated \code{\link{MVBU_Stanfit}} object with pre-extracted latents stored
#'   in \code{model@cache$draws} and \code{model@cache$summary}.
#' @seealso \code{\link{MVBU_Stanfit}}, \code{\link{get_draws}}, \code{\link{add_criterion}},
#'   \code{vignette("fitting-and-working-with-stanfit-models")}
#' @export
add_posterior_latents <- function(model, probs = c(0.025, 0.5, 0.975), ...) {
  .assert_true(S7::S7_inherits(model, MVBU_Stanfit), msg = "model must be an MVBU_Stanfit object.")

  # Extract nested draws across all groups (including prior) and categories
  draws_df <- get_draws(model, summarize = FALSE, nest = TRUE, ...)

  groups <- if (is.factor(draws_df$group)) levels(draws_df$group) else unique(draws_df$group)
  cats <- get_category_labels(model)
  cues <- get_cue_labels(model)
  D <- length(cues)

  model_type <- tryCatch(get_model_type(model), error = function(e) "")
  is_niw <- grepl("NIW", model_type, ignore.case = TRUE)
  is_nix <- grepl("NIX", model_type, ignore.case = TRUE) && !grepl("MNIX", model_type, ignore.case = TRUE)
  is_mnix <- grepl("MNIX", model_type, ignore.case = TRUE)
  model_fam <- if (is_niw) "NIW" else if (is_nix) "NIX" else if (is_mnix) "MNIX" else model_type

  # Ensure Sigma_exp and Sigma_marg are computed
  if (!"Sigma_exp" %in% names(draws_df) && "S" %in% names(draws_df) && "nu" %in% names(draws_df)) {
    draws_df$Sigma_exp <- get_expected_Sigma_from_S(draws_df$S, draws_df$nu)
  }
  if (!"Sigma_marg" %in% names(draws_df) && "S" %in% names(draws_df) && "nu" %in% names(draws_df) && "kappa" %in% names(draws_df)) {
    draws_df$Sigma_marg <- get_marginal_Sigma_from_S(draws_df$S, draws_df$nu, draws_df$kappa)
  }

  draw_col <- if (".draw" %in% names(draws_df)) ".draw" else if ("draw" %in% names(draws_df)) "draw" else ".iteration"
  draw_ids <- sort(unique(draws_df[[draw_col]]))

  # Build unified draw-wise likelihood closure
  draw_likelihood_fun <- function(new_data, categories = cats, groups = NULL, draws = draw_ids, noise_treatment = c("none", "marginalize")) {
    noise_treatment <- match.arg(noise_treatment)
    if (is.null(groups)) groups <- if ("all" %in% draws_df$group) "all" else setdiff(unique(draws_df$group), "prior")[1L]
    if (is.na(groups) || length(groups) == 0L) groups <- unique(draws_df$group)[1L]

    groups <- .validate_requested_labels(
      groups,
      unique(draws_df$group),
      label_type = "group"
    )
    categories <- .validate_requested_labels(
      categories,
      unique(draws_df$category),
      label_type = "category"
    )

    aligned <- .align_new_data_cues(new_data, cues)
    x_mat <- aligned$mat
    c_idx <- aligned$c_idx

    N_obs <- nrow(x_mat)
    K <- length(categories)
    S_draws <- length(draws)

    # 3D array: [N_obs, K, S_draws]
    out <- array(0, dim = c(N_obs, K, S_draws), dimnames = list(NULL, categories, as.character(draws)))

    sub_df <- draws_df[draws_df[[draw_col]] %in% draws & draws_df$group %in% groups & draws_df$category %in% categories, ]

    for (s_idx in seq_along(draws)) {
      d_id <- draws[s_idx]
      d_sub <- sub_df[sub_df[[draw_col]] == d_id, ]
      for (k_idx in seq_along(categories)) {
        cat_k <- categories[k_idx]
        row_k <- d_sub[d_sub$category == cat_k, ]
        if (nrow(row_k) == 0L) next
        mu_k <- as.numeric(row_k$m[[1L]])
        kappa_k <- as.numeric(row_k$kappa[[1L]])
        nu_k <- as.numeric(row_k$nu[[1L]])
        S_k <- as.matrix(row_k$S[[1L]])
        sig_marg_k <- as.matrix(row_k$Sigma_marg[[1L]])
        out[, k_idx, s_idx] <- .evaluate_stan_category_density(
          x_mat = x_mat,
          m_k = mu_k,
          kappa_k = kappa_k,
          nu_k = nu_k,
          S_k = S_k,
          sig_marg_k = sig_marg_k,
          c_idx = c_idx,
          is_niw = is_niw,
          is_nix = is_nix,
          is_mnix = is_mnix,
          D_total = D
        )
      }
    }
    out
  }

  # Build unified draw-wise posterior closure
  draw_posterior_fun <- function(new_data, categories = cats, groups = NULL, draws = draw_ids, noise_treatment = c("none", "marginalize"), lapse_treatment = c("none", "uniform")) {
    noise_treatment <- match.arg(noise_treatment)
    lapse_treatment <- match.arg(lapse_treatment)

    lik_arr <- draw_likelihood_fun(new_data, categories = categories, groups = groups, draws = draws, noise_treatment = noise_treatment)
    dim_a <- dim(lik_arr)
    N_obs <- dim_a[1L]
    K <- dim_a[2L]
    S_draws <- dim_a[3L]

    # Category prior
    prior_w <- rep(1 / K, K)
    post_arr <- lik_arr
    for (s_idx in seq_len(S_draws)) {
      mat_s <- lik_arr[, , s_idx, drop = FALSE]
      dim(mat_s) <- c(N_obs, K)
      # Multiply by prior
      mat_s <- t(t(mat_s) * prior_w)
      row_sums <- rowSums(mat_s)
      row_sums[row_sums == 0] <- 1e-12
      mat_s <- mat_s / row_sums

      if (identical(lapse_treatment, "uniform") && "lapse_rate" %in% names(draws_df)) {
        sub_d <- draws_df[draws_df[[draw_col]] == draws[s_idx], ]
        l_rate <- if (nrow(sub_d) > 0L && !is.na(sub_d$lapse_rate[1L])) sub_d$lapse_rate[1L] else 0
        mat_s <- (1 - l_rate) * mat_s + l_rate * (1 / K)
      }
      post_arr[, , s_idx] <- mat_s
    }
    post_arr
  }

  # Build summary parameter stats (mean, sd, quantiles) for ALL parameters
  cat_moments <- .summarize_moment_draws(
    draws_df,
    cues,
    type = c("expected", "marginal"),
    probs = probs,
    model_family = model_fam
  )
  exp_moments <- if (!is.null(cat_moments)) {
    cat_moments[cat_moments$Parameter %in% c("mu", "Sigma_exp"), ]
  } else {
    NULL
  }
  marg_moments <- if (!is.null(cat_moments)) {
    cat_moments[cat_moments$Parameter %in% c("mu", "Sigma_marg"), ]
  } else {
    NULL
  }

  # Aggregate parameter point estimates per group & category
  agg_rows <- list()
  grp_cat_combos <- unique(draws_df[, c("group", "category")])
  for (i in seq_len(nrow(grp_cat_combos))) {
    g_i <- grp_cat_combos$group[i]
    c_i <- grp_cat_combos$category[i]
    sub_gc <- draws_df[draws_df$group == g_i & draws_df$category == c_i, ]
    mean_m <- Reduce(`+`, sub_gc$m) / nrow(sub_gc)
    mean_S <- Reduce(`+`, sub_gc$S) / nrow(sub_gc)
    mean_sig_exp <- Reduce(`+`, sub_gc$Sigma_exp) / nrow(sub_gc)
    mean_sig_marg <- Reduce(`+`, sub_gc$Sigma_marg) / nrow(sub_gc)
    agg_rows[[length(agg_rows) + 1L]] <- tibble::tibble(
      group = g_i,
      category = c_i,
      m = list(mean_m),
      S = list(mean_S),
      kappa = mean(sub_gc$kappa),
      nu = mean(sub_gc$nu),
      Sigma_exp = list(mean_sig_exp),
      Sigma_marg = list(mean_sig_marg),
      lapse_rate = if ("lapse_rate" %in% names(sub_gc)) mean(sub_gc$lapse_rate) else NA_real_
    )
  }
  summary_parameters <- dplyr::bind_rows(agg_rows)

  # Build summary likelihood closure (evaluates from marginal moments)
  summary_likelihood_fun <- function(new_data, categories = cats, groups = NULL, noise_treatment = c("none", "marginalize")) {
    noise_treatment <- match.arg(noise_treatment)
    if (is.null(groups)) groups <- if ("all" %in% summary_parameters$group) "all" else setdiff(unique(summary_parameters$group), "prior")[1L]
    if (is.na(groups) || length(groups) == 0L) groups <- unique(summary_parameters$group)[1L]

    groups <- .validate_requested_labels(
      groups,
      unique(summary_parameters$group),
      label_type = "group"
    )
    categories <- .validate_requested_labels(
      categories,
      unique(summary_parameters$category),
      label_type = "category"
    )

    aligned <- .align_new_data_cues(new_data, cues)
    x_mat <- aligned$mat
    c_idx <- aligned$c_idx

    N_obs <- nrow(x_mat)
    K <- length(categories)
    out <- matrix(0, nrow = N_obs, ncol = K, dimnames = list(NULL, categories))

    sub_sum <- summary_parameters[summary_parameters$group %in% groups & summary_parameters$category %in% categories, ]
    for (k_idx in seq_along(categories)) {
      cat_k <- categories[k_idx]
      row_k <- sub_sum[sub_sum$category == cat_k, ]
      if (nrow(row_k) == 0L) next
      mu_k <- as.numeric(row_k$m[[1L]])
      kappa_k <- as.numeric(row_k$kappa[[1L]])
      nu_k <- as.numeric(row_k$nu[[1L]])
      S_k <- as.matrix(row_k$S[[1L]])
      sig_marg_k <- as.matrix(row_k$Sigma_marg[[1L]])
      out[, k_idx] <- .evaluate_stan_category_density(
        x_mat = x_mat,
        m_k = mu_k,
        kappa_k = kappa_k,
        nu_k = nu_k,
        S_k = S_k,
        sig_marg_k = sig_marg_k,
        c_idx = c_idx,
        is_niw = is_niw,
        is_nix = is_nix,
        is_mnix = is_mnix,
        D_total = D
      )
    }
    out
  }

  # Build summary posterior closure
  summary_posterior_fun <- function(new_data, categories = cats, groups = NULL, noise_treatment = c("none", "marginalize"), lapse_treatment = c("none", "uniform")) {
    noise_treatment <- match.arg(noise_treatment)
    lapse_treatment <- match.arg(lapse_treatment)

    lik_mat <- summary_likelihood_fun(new_data, categories = categories, groups = groups, noise_treatment = noise_treatment)
    K <- length(categories)
    prior_w <- rep(1 / K, K)
    post_mat <- t(t(lik_mat) * prior_w)
    row_sums <- rowSums(post_mat)
    row_sums[row_sums == 0] <- 1e-12
    post_mat <- post_mat / row_sums

    if (identical(lapse_treatment, "uniform")) {
      l_rate <- mean(summary_parameters$lapse_rate, na.rm = TRUE)
      if (is.finite(l_rate) && l_rate > 0) {
        post_mat <- (1 - l_rate) * post_mat + l_rate * (1 / K)
      }
    }
    post_mat
  }

  # Assemble cache$draws and cache$summary
  c_list <- .get_cache(model)
  c_list$draws <- list(
    parameters = draws_df,
    likelihood_function = draw_likelihood_fun,
    posterior_function = draw_posterior_fun
  )
  c_list$summary <- list(
    parameters = summary_parameters,
    expected_moments = exp_moments,
    marginal_moments = marg_moments,
    likelihood_function = summary_likelihood_fun,
    posterior_function = summary_posterior_fun
  )

  # Clear any legacy or redundant top-level moment fields
  c_list$category_moments <- NULL
  c_list$expected_moments <- NULL
  c_list$marginal_moments <- NULL

  .set_cache(model, c_list)
}

#' Add Information Criteria to a Stanfit Object (brms-compatible)
#'
#' Computes and stores model fit criteria (such as \code{"loo"}, \code{"waic"},
#' \code{"kfold"}, \code{"loo_subsample"}, \code{"bayes_R2"}, \code{"loo_R2"}, or
#' \code{"marglik"}) in an \code{\link{MVBU_Stanfit}} object's \code{@criteria}
#' slot, following the interface convention of \pkg{brms}.
#'
#' @param x An \code{\link{MVBU_Stanfit}} object.
#' @param criterion Character vector of criteria to add. Supported: \code{"loo"},
#'   \code{"waic"}, \code{"kfold"}, \code{"loo_subsample"}, \code{"bayes_R2"},
#'   \code{"loo_R2"}, \code{"marglik"}. Default: \code{"loo"}.
#' @param model_name Optional character string giving a name for the model.
#' @param overwrite Logical; whether to overwrite existing criteria. Default: \code{TRUE}.
#' @param ... Additional arguments passed to underlying criterion evaluation functions
#'   (such as \code{\link[loo]{loo}}, \code{\link[loo]{waic}}, \code{\link[loo]{kfold}},
#'   or \code{bridgesampling::bridge_sampler}).
#'
#' @return The updated \code{\link{MVBU_Stanfit}} object with criteria stored in
#'   \code{x@criteria}.
#' @seealso \code{\link{loo_compare.MVBU_Stanfit}}, \code{\link{loo.MVBU_Stanfit}}
#' @export
add_criterion <- function(x, criterion = "loo", model_name = NULL, overwrite = TRUE, ...) {
  .assert_true(S7::S7_inherits(x, MVBU_Stanfit), msg = "x must be an MVBU_Stanfit object.")
  .assert_true(is.character(criterion) && length(criterion) >= 1L, msg = "criterion must be a character vector.")

  criterion <- tolower(criterion)
  valid_criteria <- c("loo", "waic", "kfold", "loo_subsample", "bayes_r2", "loo_r2", "marglik")
  invalid <- setdiff(criterion, valid_criteria)
  if (length(invalid) > 0L) {
    .stop(
      "Unsupported criterion: ", paste(invalid, collapse = ", "),
      ". Supported criteria: 'loo', 'waic', 'kfold', 'loo_subsample', 'bayes_R2', 'loo_R2', 'marglik'."
    )
  }

  current_criteria <- x@criteria
  if (is.null(current_criteria)) current_criteria <- list()

  if (!is.null(model_name)) {
    x@metadata$model_name <- as.character(model_name)
  }

  stanfit <- get_stanfit(x)

  for (crit in criterion) {
    if (!overwrite && crit %in% names(current_criteria)) {
      next
    }
    if (identical(crit, "loo")) {
      current_criteria$loo <- loo.MVBU_Stanfit(x, ...)
    } else if (identical(crit, "waic")) {
      LLarray <- loo::extract_log_lik(stanfit = stanfit, parameter_name = "log_lik", merge_chains = FALSE)
      current_criteria$waic <- loo::waic(LLarray)
    } else if (identical(crit, "loo_subsample")) {
      LLarray <- loo::extract_log_lik(stanfit = stanfit, parameter_name = "log_lik", merge_chains = FALSE)
      current_criteria$loo_subsample <- loo::loo_subsample(LLarray, ...)
    } else if (identical(crit, "kfold")) {
      LLarray <- loo::extract_log_lik(stanfit = stanfit, parameter_name = "log_lik", merge_chains = FALSE)
      current_criteria$kfold <- loo::kfold(LLarray, ...)
    } else if (identical(crit, "bayes_r2")) {
      LLmat <- loo::extract_log_lik(stanfit = stanfit, parameter_name = "log_lik", merge_chains = TRUE)
      y_pred <- exp(LLmat)
      var_fit <- apply(y_pred, 1, stats::var)
      var_res <- apply(1 - y_pred, 1, stats::var)
      r2_draws <- var_fit / (var_fit + var_res)
      r2_draws <- r2_draws[!is.na(r2_draws)]
      current_criteria$bayes_R2 <- structure(
        c(
          Estimate = mean(r2_draws), Est.Error = stats::sd(r2_draws),
          Q2.5 = stats::quantile(r2_draws, 0.025, names = FALSE),
          Q97.5 = stats::quantile(r2_draws, 0.975, names = FALSE)
        ),
        class = "bayes_R2"
      )
    } else if (identical(crit, "loo_r2")) {
      LLmat <- loo::extract_log_lik(stanfit = stanfit, parameter_name = "log_lik", merge_chains = TRUE)
      y_pred <- exp(LLmat)
      var_fit <- apply(y_pred, 1, stats::var)
      var_res <- apply(1 - y_pred, 1, stats::var)
      r2_draws <- var_fit / (var_fit + var_res)
      r2_draws <- r2_draws[!is.na(r2_draws)]
      current_criteria$loo_R2 <- structure(
        c(
          Estimate = mean(r2_draws), Est.Error = stats::sd(r2_draws),
          Q2.5 = stats::quantile(r2_draws, 0.025, names = FALSE),
          Q97.5 = stats::quantile(r2_draws, 0.975, names = FALSE)
        ),
        class = "loo_R2"
      )
    } else if (identical(crit, "marglik")) {
      if (!requireNamespace("bridgesampling", quietly = TRUE)) {
        .stop("Package 'bridgesampling' is required for criterion 'marglik'. Please install it.")
      }
      current_criteria$marglik <- bridgesampling::bridge_sampler(stanfit, ...)
    }
  }

  x@criteria <- current_criteria
  x
}

#' Compare Models Fitted with Stan
#'
#' Compare fitted \code{\link{MVBU_Stanfit}} objects using Leave-One-Out Cross-Validation
#' (\code{loo}) or WAIC via \code{\link[loo]{loo_compare}}.
#'
#' @param x An \code{\link{MVBU_Stanfit}} object.
#' @param ... Additional \code{\link{MVBU_Stanfit}} objects or arguments passed to
#'   \code{\link[loo]{loo_compare}}.
#' @param criterion Character string specifying the criterion to use for model
#'   comparison (\code{"loo"} or \code{"waic"}). Default: \code{"loo"}.
#'
#' @return A matrix of class \code{compare.loo} (see \code{\link[loo]{loo_compare}}) containing model comparison metrics.
#' @seealso \code{\link{add_criterion}}, \code{\link{loo.MVBU_Stanfit}}
#' @export
loo_compare.MVBU_Stanfit <- function(x, ..., criterion = "loo") {
  # Note for future extensions:
  # LATER we might want to extend the stanfit loo approach to versions of the basic core classes
  # that are fitted to listener data (e.g., by fitting lapse rate, lapse biases, etc. through
  # bootstrap or cross-validation or alike). At that point, we might implement log-lik, AIC, BIC,
  # and similar information criteria for these frequentist fitted models.

  dots <- list(...)
  is_stanfit <- vapply(dots, function(obj) S7::S7_inherits(obj, MVBU_Stanfit), logical(1))
  models <- c(list(x), dots[is_stanfit])
  other_args <- dots[!is_stanfit]

  criterion <- tolower(criterion)

  # Ensure all models have the requested criterion added
  models <- lapply(models, function(m) {
    if (is.null(m@criteria[[criterion]])) {
      add_criterion(m, criterion = criterion)
    } else {
      m
    }
  })

  crit_list <- lapply(models, function(m) m@criteria[[criterion]])

  # Assign names if available
  model_names <- vapply(seq_along(models), function(i) {
    name <- models[[i]]@metadata$model_name
    if (!is.null(name) && nchar(name) > 0) name else paste0("model", i)
  }, character(1))
  names(crit_list) <- model_names

  do.call(loo::loo_compare, c(list(x = crit_list), other_args))
}


#' Graphical Posterior Predictive Checks for Stanfit Objects
#'
#' Performs posterior predictive checks using \pkg{bayesplot}'s \code{\link[bayesplot]{ppc_dens_overlay}}
#' or related functions on an \code{\link{MVBU_Stanfit}} model object.
#'
#' @param object An \code{\link{MVBU_Stanfit}} object.
#' @param type Character string giving the type of \pkg{bayesplot} plot to create
#'   (e.g., \code{"dens_overlay"}, \code{"hist"}, \code{"stat"}, \code{"bars"}, \code{"error_hist"}).
#'   Default: \code{"dens_overlay"}.
#' @param ndraws Positive integer; number of posterior draws to overlay. Default: \code{50}.
#' @param ... Additional arguments passed to the underlying \pkg{bayesplot} plotting function.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object returned by \pkg{bayesplot}.
#' @export
pp_check.MVBU_Stanfit <- function(object, type = "dens_overlay", ndraws = 50, ...) {
  .assert_true(S7::S7_inherits(object, MVBU_Stanfit), msg = "object must be an MVBU_Stanfit object.")

  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    .stop("Package 'bayesplot' is required for pp_check(). Please install it.")
  }

  cues <- get_cue_labels(object)
  data_df <- object@data

  # Determine response vector y (cue values or category responses)
  y <- if (!is.null(cues) && length(cues) > 0L && cues[1] %in% names(data_df)) {
    as.numeric(data_df[[cues[1]]])
  } else if ("category" %in% names(data_df)) {
    as.numeric(as.factor(data_df$category))
  } else {
    .stop("Could not determine observed response vector y from model data.")
  }

  n_obs <- length(y)
  draws_df <- get_draws(object, summarize = FALSE, nest = TRUE)
  available_draws <- if (".draw" %in% names(draws_df)) length(unique(draws_df$.draw)) else nrow(draws_df)

  ndraws <- min(as.integer(ndraws), available_draws)

  # Generate synthetic responses yrep matrix (ndraws x n_obs)
  yrep <- matrix(NA_real_, nrow = ndraws, ncol = n_obs)

  for (i in seq_len(ndraws)) {
    obs <- sample_observations(object, n = n_obs)
    if (!is.null(cues) && length(cues) > 0L && cues[1] %in% names(obs)) {
      yrep[i, ] <- as.numeric(obs[[cues[1]]])
    } else if ("category" %in% names(obs)) {
      yrep[i, ] <- as.numeric(as.factor(obs$category))
    }
  }

  ppc_func_name <- paste0("ppc_", type)
  ppc_func <- tryCatch(
    get(ppc_func_name, envir = asNamespace("bayesplot")),
    error = function(e) NULL
  )

  if (is.null(ppc_func)) {
    .stop("Unsupported pp_check type: '", type, "'. See ?bayesplot::PPC-overview for supported PPC types.")
  }

  ppc_func(y, yrep, ...)
}
