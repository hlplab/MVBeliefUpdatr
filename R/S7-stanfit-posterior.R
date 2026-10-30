# =============================================================================
# Lightweight S7 Representation of Stanfit Posterior Parameter Draws
# =============================================================================

#' @include S7-generics.R
#' @include S7-class.R
#' @include S7-transform-information.R
#' @include S7-stanfit.R
#' @include S7-stanfit-input-methods.R
#' @importFrom S7 new_class new_object method class_character class_list
#' @importFrom S7 class_any new_property new_S3_class
#' @importFrom posterior as_draws_df
#' @importFrom stats median quantile sd
#' @importFrom utils head
NULL

#' S7 class for Stanfit posterior parameter draws
#'
#' Represents a posterior distribution over a cognitive model, as obtained
#' for example from an [MVBU_Stanfit] object. Holds MCMC parameter draw
#' distributions without carrying sampling chains, diagnostics, or raw model
#' source structures, and is thus more lightweight than the [MVBU_Stanfit]
#' class.
#' @param draws A list containing \code{raw} and \code{summary} posterior
#'   draws. Default: \code{list()}.
#' @param metadata A list containing auxiliary information including
#'   \code{label_information} (cues, categories, group), \code{model_family},
#'   \code{parameters}, \code{original_variable_names}, and
#'   \code{exposure_category_statistics}. Default: \code{list()}.
#' @param data A named \code{list} containing \code{$original} and
#'   \code{$transformed} data frames. Default: \code{list()}.
#' @param cache A list storing cached computations and summaries.
#'   Default: \code{list()}.
#' @param transform_information An optional [MVBU_TransformInformation] object.
#'   Default: \code{NULL}.
#'
#' @name MVBU_StanfitPosterior
#' @rdname MVBU_StanfitPosterior
#' @docType class
#' @export
MVBU_StanfitPosterior <- S7::new_class(
  "MVBU_StanfitPosterior",
  package = NULL,
  parent = MVBU_Object,
  properties = list(
    draws = S7::class_list,
    metadata = S7::class_list,
    data = S7::class_list,
    cache = S7::class_list,
    transform_information = S7::new_union(NULL, MVBU_TransformInformation)
  ),
  constructor = function(
    draws = list(),
    metadata = list(),
    data = list(),
    cache = list(),
    transform_information = NULL
  ) {
    if (!is.list(metadata)) {
      metadata <- list()
    }

    if (is.data.frame(data)) {
      data <- list(original = data, transformed = data)
    } else if (is.list(data)) {
      if (is.null(data$original) && is.null(data$transformed) &&
        length(data) == 0L) {
        data <- list(
          original = tibble::tibble(),
          transformed = tibble::tibble()
        )
      } else if (is.null(data$original) && !is.null(data$transformed)) {
        data$original <- data$transformed
      } else if (!is.null(data$original) && is.null(data$transformed)) {
        data$transformed <- data$original
      }
    }

    label_info <- if (!is.null(metadata$label_information) &&
      is.list(metadata$label_information)) {
      metadata$label_information
    } else {
      list()
    }

    metadata$label_information <- list(
      cue = if (!is.null(label_info$cue)) {
        as.character(label_info$cue)
      } else {
        character(0)
      },
      category = if (!is.null(label_info$category)) {
        as.character(label_info$category)
      } else {
        character(0)
      },
      response_category = if (!is.null(label_info$response_category)) {
        as.character(label_info$response_category)
      } else {
        character(0)
      },
      group = if (!is.null(label_info$group)) {
        as.character(label_info$group)
      } else {
        character(0)
      }
    )

    S7::new_object(
      MVBU_Object(),
      draws = as.list(draws),
      metadata = metadata,
      data = data,
      cache = as.list(cache),
      transform_information = transform_information
    )
  },
  validator = function(self) {
    if (!is.list(self@draws)) {
      return("`draws` must be a list.")
    }
    if (!is.list(self@metadata)) {
      return("`metadata` must be a list.")
    }
    if (!is.list(self@data) || is.data.frame(self@data)) {
      return("`data` must be a list containing $original and $transformed.")
    }
    if (!inherits(self@data$original, "data.frame")) {
      return("`data$original` must inherit from data.frame.")
    }
    if (!inherits(self@data$transformed, "data.frame")) {
      return("`data$transformed` must inherit from data.frame.")
    }
    if (!is.list(self@cache)) {
      return("`cache` must be a list.")
    }
    lbls <- self@metadata$label_information
    if (!is.null(lbls)) {
      if (length(lbls$cue) > 0L && !is.character(lbls$cue)) {
        return("`metadata$label_information$cue` must be a character vector.")
      }
      if (length(lbls$category) > 0L && !is.character(lbls$category)) {
        return(
          "`metadata$label_information$category` must be a character vector."
        )
      }
      if (
        length(lbls$response_category) > 0L &&
          !is.character(lbls$response_category)
      ) {
        return(paste0(
          "`metadata$label_information$response_category` must be a ",
          "character vector."
        ))
      }
      if (length(lbls$group) > 0L && !is.character(lbls$group)) {
        return("`metadata$label_information$group` must be a character vector.")
      }
    }
    NULL
  }
)

#' Coerce an MVBU_Stanfit object to a lightweight MVBU_StanfitPosterior
#'
#' @param x An [MVBU_Stanfit] object.
#' @param ... Additional arguments (currently unused).
#' @return An [MVBU_StanfitPosterior] object.
#' @rdname as_MVBU_stanfit_posterior
#' @export
S7::method(as_MVBU_stanfit_posterior, MVBU_Stanfit) <- function(x, ...) {
  meta <- get_metadata(x)
  if (is.null(meta$label_information) || !is.list(meta$label_information)) {
    meta$label_information <- list(
      cue = tryCatch(get_cue_labels(x), error = function(e) character(0)),
      category = tryCatch(
        get_category_labels(x),
        error = function(e) character(0)
      ),
      response_category = tryCatch(
        get_response_category_labels(x),
        error = function(e) character(0)
      ),
      group = tryCatch(
        get_group_labels(x, include_prior = FALSE),
        error = function(e) character(0)
      )
    )
  } else {
    if (length(meta$label_information$cue) == 0L) {
      meta$label_information$cue <- tryCatch(
        get_cue_labels(x),
        error = function(e) character(0)
      )
    }
    if (length(meta$label_information$category) == 0L) {
      meta$label_information$category <- tryCatch(
        get_category_labels(x),
        error = function(e) character(0)
      )
    }
    if (length(meta$label_information$response_category) == 0L) {
      meta$label_information$response_category <- tryCatch(
        get_response_category_labels(x),
        error = function(e) character(0)
      )
    }
    if (length(meta$label_information$group) == 0L) {
      meta$label_information$group <- tryCatch(
        get_group_labels(x, include_prior = FALSE),
        error = function(e) character(0)
      )
    }
  }
  meta$model_family <- get_model_family(x)
  meta$model_name <- meta$model_family

  sf <- x@stanfit
  draws_df <- if (!is.null(sf) && inherits(sf, "stanfit")) {
    posterior::as_draws_df(sf)
  } else if (is.list(sf) && "draws" %in% names(sf)) {
    sf$draws
  } else {
    list()
  }

  summary_df <- tryCatch(
    get_draws(x, summarize = TRUE, nest = TRUE),
    error = function(e) NULL
  )

  data <- if (!is.null(x@data)) {
    x@data
  } else {
    list(original = tibble::tibble(), transformed = tibble::tibble())
  }
  meta$original_variable_names <- get_original_variable_names(x)
  meta$exposure_category_statistics <- tryCatch(
    get_exposure_category_statistics(x, simplify = FALSE),
    error = function(e) NULL
  )
  if (is.null(meta$parameters)) {
    meta$parameters <- tryCatch(get_parameters(x), error = function(e) list())
  }

  cache <- if (!is.null(x@cache) && is.list(x@cache)) x@cache else list()

  MVBU_StanfitPosterior(
    draws = list(raw = draws_df, summary = summary_df),
    metadata = meta,
    data = data,
    cache = cache,
    transform_information = x@transform_information
  )
}

#' @rdname as_MVBU_stanfit_posterior
#' @export
S7::method(as_MVBU_stanfit_posterior, MVBU_StanfitPosterior) <- function(
  x,
  ...
) {
  x
}



# -----------------------------------------------------------------------------
# Shared Stanfit & StanfitPosterior Methods
# -----------------------------------------------------------------------------
# The following methods for MVBU_StanfitPosterior are defined jointly with
# MVBU_Stanfit in R/S7-stanfit-methods.R via S7::new_union(MVBU_Stanfit, ...):
# - get_model_family
# - get_family_name
# - get_family
# - get_transform_information
# - get_labels
# - get_cue_labels
# - get_category_labels
# - get_response_category_labels
# - get_group_labels
# - get_original_variable_names
# - get_data
# - get_exposure_data
# - get_test_data
# - get_expected_mu
# - get_expected_sigma
# - get_marginal_mu
# - get_marginal_sigma
# - evaluate_model
# - likelihood
# - posterior
# - categorize
# - get_category_likelihood_function
# - get_sufficient_category_statistics


#' @rdname get_draws
#' @export
S7::method(get_draws, MVBU_StanfitPosterior) <- function(
  fit,
  categories = NULL,
  groups = NULL,
  which = "posterior",
  ndraws = NULL,
  summarize = FALSE,
  nest = TRUE,
  seed = NULL,
  ...
) {
  if (is.null(categories)) categories <- get_category_labels(fit)
  if (is.null(groups)) groups <- get_group_labels(fit)
  if (summarize && !is.null(fit@draws$summary)) {
    return(fit@draws$summary)
  }
  raw_draws <- fit@draws$raw
  if (is.data.frame(raw_draws) || inherits(raw_draws, "draws")) {
    raw_draws
  } else if (!is.null(fit@draws$summary)) {
    fit@draws$summary
  } else {
    tibble::tibble()
  }
}

#' Helper to extract expected cognitive model from MVBU_StanfitPosterior
#' @param model An MVBU_StanfitPosterior or MVBU_Stanfit object.
#' @param group Optional group name. If NULL, auto-selected if only 1
#'   post-exposure group exists.
#' @keywords internal
#' @noRd
.as_cognitive_model_from_posterior <- function(model, group = NULL) {
  if (S7::S7_inherits(model, MVBU_Stanfit)) {
    model <- as_MVBU_stanfit_posterior(model)
  }
  pars_sum <- model@draws$summary
  if (is.null(pars_sum) || nrow(pars_sum) == 0L) {
    .stop("No parameter summary found in MVBU_StanfitPosterior object.")
  }

  avail_groups <- if (is.factor(pars_sum$group)) {
    levels(pars_sum$group)
  } else {
    unique(pars_sum$group)
  }
  post_groups <- setdiff(avail_groups, "prior")

  if (!is.null(group)) {
    group_str <- as.character(group)
    if (!group_str %in% avail_groups) {
      .stop(
        "Group '", group_str, "' not found in model posterior. ",
        "Available groups: ", paste(avail_groups, collapse = ", ")
      )
    }
    target_group <- group_str
  } else {
    if (length(post_groups) == 1L) {
      target_group <- post_groups[1]
    } else if (length(post_groups) == 0L) {
      target_group <- "prior"
    } else {
      .stop(
        "Multiple groups found in model posterior: ",
        paste(post_groups, collapse = ", "),
        ". Please specify which group to extract via 'group'."
      )
    }
  }

  pars_group <- pars_sum[pars_sum$group == target_group, ]

  cues <- get_cue_labels(model)
  cats <- get_category_labels(model)
  fam <- get_model_family(model)
  is_nix <- grepl("NIX", fam, ignore.case = TRUE) &&
    !grepl("MNIX", fam, ignore.case = TRUE)
  is_mnix <- grepl("MNIX", fam, ignore.case = TRUE)

  cat_reps <- lapply(cats, function(cat_name) {
    row_match <- which(pars_group$category == cat_name)
    if (length(row_match) == 0L) row_match <- 1L
    m_val <- pars_group$m[[row_match]]
    s_val <- pars_group$S[[row_match]]
    kappa_val <- pars_group$kappa[row_match]
    nu_val <- pars_group$nu[row_match]
    if (is_nix) {
      new_nix_category_representation(
        category_labels = cat_name,
        cue_labels = cues,
        m = as.numeric(m_val)[1],
        sigma2 = as.numeric(s_val)[1] / as.numeric(nu_val),
        kappa = as.numeric(kappa_val),
        nu = as.numeric(nu_val)
      )
    } else if (is_mnix) {
      k_comp <- length(as.vector(m_val))
      s_diag <- if (is.matrix(s_val)) diag(s_val) else as.numeric(s_val)
      new_mnix_category_representation(
        category_labels = cat_name,
        cue_labels = cues,
        m = as.vector(m_val),
        sigma2 = as.vector(s_diag) / as.numeric(nu_val),
        kappa = rep_len(as.numeric(kappa_val), k_comp),
        nu = rep_len(as.numeric(nu_val), k_comp)
      )
    } else {
      new_niw_category_representation(
        category_labels = cat_name,
        cue_labels = cues,
        m = as.vector(m_val),
        S = as.matrix(s_val),
        kappa = as.numeric(kappa_val),
        nu = as.numeric(nu_val)
      )
    }
  })
  names(cat_reps) <- cats

  lapse <- if ("lapse_rate" %in% names(pars_group)) {
    pars_group$lapse_rate[1]
  } else {
    0
  }

  if (is_nix) {
    new_nix_ideal_adaptor(
      category_template = new_category_representation_template(cat_reps),
      lapse_rate = lapse
    )
  } else if (is_mnix) {
    new_mnix_ideal_adaptor(
      category_template = new_category_representation_template(cat_reps),
      lapse_rate = lapse
    )
  } else {
    new_niw_ideal_adaptor(
      category_template = new_category_representation_template(cat_reps),
      lapse_rate = lapse
    )
  }
}


#' @export
S7::method(print, MVBU_StanfitPosterior) <- function(x, ...) {
  cat("MVBeliefUpdatr Lightweight Stanfit Posterior\n")
  cat("Model Family     :", get_model_family(x), "\n")
  cat("Group Context    :", paste(get_group_labels(x), collapse = ", "), "\n")
  cat("Cues             :", paste(get_cue_labels(x), collapse = ", "), "\n")
  cat(
    "Categories       :",
    paste(get_category_labels(x), collapse = ", "),
    "\n"
  )
  invisible(x)
}

#' @exportS3Method base::summary
summary.MVBU_StanfitPosterior <- function(object, ...) {
  print(object)
  fixed_params <- get_parameters(object, space = "original")$fixed
  fixed_summary <- .format_fixed_parameters_table(fixed_params, object = object)
  if (nrow(fixed_summary) > 0L) {
    cat("\nFixed Parameters:\n")
    print(fixed_summary)
  }
  if (!is.null(object@draws$summary)) {
    cat("\nParameter Summary:\n")
    print(object@draws$summary)
  }

  grps <- get_group_labels(object, include_prior = TRUE)
  moments_by_group <- list()
  for (g in grps) {
    cog_model <- tryCatch(
      .as_cognitive_model_from_posterior(object, group = g),
      error = function(e) NULL
    )
    if (!is.null(cog_model)) {
      moments_by_group[[g]] <- list(
        mu = get_expected_mu(cog_model),
        Sigma_exp = get_expected_sigma(cog_model),
        Sigma_marg = get_marginal_sigma(cog_model)
      )
    }
  }

  if (length(moments_by_group) > 0L) {
    cat("\nCategory moments:\n")
    if (length(moments_by_group) == 1L) {
      print(moments_by_group[[1L]])
    } else {
      print(moments_by_group)
    }
  }
  invisible(object)
}

S7::method(summary, MVBU_StanfitPosterior) <- summary.MVBU_StanfitPosterior

#' @rdname get_expected_category_statistic
#' @export
S7::method(get_expected_category_statistic, MVBU_StanfitPosterior) <- function(
  x,
  categories = NULL,
  groups = NULL,
  statistic = c("mu", "Sigma"),
  ...
) {
  if (is.null(categories)) categories <- get_category_labels(x)
  if (is.null(groups)) groups <- get_group_labels(x, include_prior = TRUE)
  .assert_that(all(statistic %in% c("mu", "Sigma")))
  avail_cats <- get_category_labels(x)
  categories <- .validate_requested_labels(
    categories,
    avail_cats,
    label_type = "category"
  )
  avail_grps <- get_group_labels(x, include_prior = TRUE)
  groups <- .validate_requested_labels(
    groups,
    avail_grps,
    label_type = "group"
  )

  cached <- x@cache$summary$expected_moments
  if (!is.null(cached) && is.data.frame(cached) &&
    "category" %in% names(cached) && "group" %in% names(cached) &&
    all(categories %in% cached$category) &&
    all(groups %in% cached$group)) {
    draws <- cached %>%
      dplyr::filter(
        .data$category %in% categories,
        .data$group %in% groups
      ) %>%
      dplyr::select(
        dplyr::all_of(c("group", "category", paste0(statistic, ".mean")))
      ) %>%
      dplyr::mutate(
        category = factor(.data$category, levels = categories),
        group = factor(.data$group, levels = groups)
      )
  } else {
    res_list <- lapply(groups, function(g) {
      cog <- .as_cognitive_model_from_posterior(x, group = g)
      cat_df <- tibble::tibble(
        group = factor(g, levels = groups),
        category = factor(categories, levels = categories)
      )
      if ("mu" %in% statistic) {
        cat_df$mu.mean <- get_expected_mu(cog)[categories]
      }
      if ("Sigma" %in% statistic) {
        cat_df$Sigma.mean <- get_expected_sigma(cog)[categories]
      }
      cat_df
    })
    draws <- dplyr::bind_rows(res_list)
  }

  if (nrow(draws) == 1 && length(statistic) == 1) {
    draws <- draws[[paste0(statistic, ".mean")]][[1]]
  }

  draws
}

#' @rdname get_marginal_category_statistic
#' @export
S7::method(get_marginal_category_statistic, MVBU_StanfitPosterior) <- function(
  x,
  categories = NULL,
  groups = NULL,
  statistic = c("mu", "Sigma"),
  ...
) {
  if (is.null(categories)) categories <- get_category_labels(x)
  if (is.null(groups)) groups <- get_group_labels(x, include_prior = TRUE)
  .assert_that(all(statistic %in% c("mu", "Sigma")))
  avail_cats <- get_category_labels(x)
  categories <- .validate_requested_labels(
    categories,
    avail_cats,
    label_type = "category"
  )
  avail_grps <- get_group_labels(x, include_prior = TRUE)
  groups <- .validate_requested_labels(
    groups,
    avail_grps,
    label_type = "group"
  )

  cached <- x@cache$summary$marginal_moments
  if (!is.null(cached) && is.data.frame(cached) &&
    "category" %in% names(cached) && "group" %in% names(cached) &&
    all(categories %in% cached$category) &&
    all(groups %in% cached$group)) {
    draws <- cached %>%
      dplyr::filter(
        .data$category %in% categories,
        .data$group %in% groups
      ) %>%
      dplyr::select(
        dplyr::all_of(c("group", "category", paste0(statistic, ".mean")))
      ) %>%
      dplyr::mutate(
        category = factor(.data$category, levels = categories),
        group = factor(.data$group, levels = groups)
      )
  } else {
    res_list <- lapply(groups, function(g) {
      cog <- .as_cognitive_model_from_posterior(x, group = g)
      cat_df <- tibble::tibble(
        group = factor(g, levels = groups),
        category = factor(categories, levels = categories)
      )
      if ("mu" %in% statistic) {
        cat_df$mu.mean <- get_marginal_mu(cog)[categories]
      }
      if ("Sigma" %in% statistic) {
        cat_df$Sigma.mean <- get_marginal_sigma(cog)[categories]
      }
      cat_df
    })
    draws <- dplyr::bind_rows(res_list)
  }

  if (nrow(draws) == 1 && length(statistic) == 1) {
    draws <- draws[[paste0(statistic, ".mean")]][[1]]
  }

  draws
}


