#' @include S7-core-classes.R
#' @include S7-stanfit.R
#' @include S7-stanfit-methods.R
#' @importFrom S7 new_generic method S7_dispatch
#' @importFrom loo loo waic loo_compare
#' @importFrom bayesplot pp_check
NULL

#' Add Pre-extracted Parameter Draw Matrices to a Stanfit Object
#'
#' Optional post-processing helper that pre-extracts parameter draw matrices
#' (\eqn{\mu}, \eqn{\Sigma}, \eqn{m}, \eqn{S}, log-priors) from an
#' \code{\link{MVBU_Stanfit}} object and caches them in the object's
#' \code{@cache} slot for accelerated plotting and matrix-level evaluations.
#'
#' @param model An \code{\link{MVBU_Stanfit}} object.
#' @param pars Character vector of parameters to pre-extract. Defaults to
#'   \code{c("mu", "Sigma", "m", "S", "log_prior")}.
#' @param ... Additional arguments passed to \code{\link{get_draws}}.
#'
#' @return The updated \code{\link{MVBU_Stanfit}} object with pre-extracted draw
#'   matrices stored in \code{model@cache$draw_matrices}.
#' @seealso \code{\link{get_draws}}, \code{\link{add_criterion}}
#' @export
add_parameter_draws <- function(model, pars = c("mu", "Sigma", "m", "S", "log_prior"), ...) {
  .assert_true(S7::S7_inherits(model, MVBU_Stanfit), msg = "model must be an MVBU_Stanfit object.")
  
  # Reuse get_draws core extraction logic to avoid code duplication
  draws_df <- get_draws(model, summarize = FALSE, nest = TRUE, untransform_cues = TRUE, ...)
  
  # Organize extracted draws by group/category for fast indexed slicing
  groups <- if (is.factor(draws_df$group)) levels(draws_df$group) else unique(draws_df$group)
  cats <- get_category_labels(model)
  cues <- get_cue_labels(model)
  
  draw_matrices <- list(
    draws_df = draws_df,
    groups = groups,
    categories = cats,
    cues = cues,
    pars = pars,
    n_draws = if (".draw" %in% names(draws_df)) {
      length(unique(draws_df$.draw))
    } else if ("draw" %in% names(draws_df)) {
      length(unique(draws_df$draw))
    } else if (".iteration" %in% names(draws_df)) {
      length(unique(draws_df$.iteration))
    } else {
      nrow(draws_df)
    }
  )
  
  c_list <- .get_cache(model)
  c_list$draw_matrices <- draw_matrices
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
    .stop("Unsupported criterion: ", paste(invalid, collapse = ", "),
          ". Supported criteria: 'loo', 'waic', 'kfold', 'loo_subsample', 'bayes_R2', 'loo_R2', 'marglik'.")
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
        c(Estimate = mean(r2_draws), Est.Error = stats::sd(r2_draws),
          Q2.5 = stats::quantile(r2_draws, 0.025, names = FALSE),
          Q97.5 = stats::quantile(r2_draws, 0.975, names = FALSE)),
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
        c(Estimate = mean(r2_draws), Est.Error = stats::sd(r2_draws),
          Q2.5 = stats::quantile(r2_draws, 0.025, names = FALSE),
          Q97.5 = stats::quantile(r2_draws, 0.975, names = FALSE)),
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
#' @return A matrix of class \code{\link[loo]{compare.loo}} containing model comparison metrics.
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
