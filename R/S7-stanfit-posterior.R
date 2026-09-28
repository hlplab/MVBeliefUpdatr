# =============================================================================
# Lightweight S7 Representation of Stanfit Posterior Parameter Draws
# =============================================================================

#' @include S7-generics.R
#' @include S7-class.R
#' @include S7-stanfit.R
#' @importFrom S7 new_class new_object method class_character class_list class_any
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
#'   \code{original_variable_names}, \code{sufficient_statistics}, and
#'   attached data. Default: \code{list()}.
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
    metadata = S7::class_list
  ),
  constructor = function(
    draws = list(),
    metadata = list()
  ) {
    if (!is.list(metadata)) {
      metadata <- list()
    }
    label_info <- if (!is.null(metadata$label_information) &&
      is.list(metadata$label_information)) {
      metadata$label_information
    } else {
      list()
    }

    metadata$label_information <- list(
      cue = if (!is.null(label_info$cue)) as.character(label_info$cue) else character(0),
      category = if (!is.null(label_info$category)) as.character(label_info$category) else character(0),
      response_category = if (!is.null(label_info$response_category)) as.character(label_info$response_category) else character(0),
      group = if (!is.null(label_info$group)) as.character(label_info$group) else character(0)
    )

    S7::new_object(
      MVBU_Object(),
      draws = as.list(draws),
      metadata = metadata
    )
  },
  validator = function(self) {
    if (!is.list(self@draws)) {
      return("`draws` must be a list.")
    }
    if (!is.list(self@metadata)) {
      return("`metadata` must be a list.")
    }
    lbls <- self@metadata$label_information
    if (!is.null(lbls)) {
      if (length(lbls$cue) > 0L && !is.character(lbls$cue)) {
        return("`metadata$label_information$cue` must be a character vector.")
      }
      if (length(lbls$category) > 0L && !is.character(lbls$category)) {
        return("`metadata$label_information$category` must be a character vector.")
      }
      if (length(lbls$response_category) > 0L && !is.character(lbls$response_category)) {
        return("`metadata$label_information$response_category` must be a character vector.")
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
      category = tryCatch(get_category_labels(x), error = function(e) character(0)),
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

  test_data <- tryCatch(
    get_test_data(x, original_names = FALSE),
    error = function(e) tibble::tibble()
  )
  meta$test_data <- test_data
  meta$original_variable_names <- get_original_variable_names(x)
  meta$sufficient_statistics <- tryCatch(
    get_sufficient_category_statistics(x),
    error = function(e) NULL
  )

  MVBU_StanfitPosterior(
    draws = list(raw = draws_df, summary = summary_df),
    metadata = meta
  )
}

#' @rdname get_cue_labels
#' @export
S7::method(get_cue_labels, MVBU_StanfitPosterior) <- function(x, indices = NULL, ...) {
  cue_labels <- get_labels(x)$cue
  if (missing(indices) || is.null(indices)) {
    return(cue_labels)
  }
  cue_labels[indices]
}

#' @rdname get_category_labels
#' @export
S7::method(get_category_labels, MVBU_StanfitPosterior) <- function(x, indices = NULL, ...) {
  category_labels <- get_labels(x)$category
  if (missing(indices) || is.null(indices)) {
    return(category_labels)
  }
  category_labels[indices]
}

#' @rdname get_response_category_labels
#' @export
S7::method(get_response_category_labels, MVBU_StanfitPosterior) <- function(x, indices = NULL, ...) {
  resp_labels <- get_labels(x)$response_category
  if (missing(indices) || is.null(indices)) {
    return(resp_labels)
  }
  resp_labels[indices]
}

#' @rdname get_group_labels
#' @export
S7::method(get_group_labels, MVBU_StanfitPosterior) <- function(x, indices = NULL, include_prior = FALSE, ...) {
  group_labels <- get_labels(x)$group
  if (include_prior) {
    group_labels <- append("prior", group_labels)
  }
  if (missing(indices) || is.null(indices)) {
    return(group_labels)
  }
  group_labels[indices]
}

#' @rdname get_labels
#' @export
S7::method(get_labels, MVBU_StanfitPosterior) <- function(x, ...) {
  lbls <- x@metadata$label_information
  if (is.null(lbls) || !is.list(lbls)) {
    .stop("No label information found in object metadata.")
  }
  lbls
}

#' @rdname get_model_family
#' @export
S7::method(get_model_family, MVBU_StanfitPosterior) <- function(x, ...) {
  if (!is.null(x@metadata$model_family)) {
    x@metadata$model_family
  } else if (!is.null(x@metadata$model_name)) {
    x@metadata$model_name
  } else if (!is.null(x@metadata$stanmodel)) {
    x@metadata$stanmodel
  } else {
    "MVBU_StanfitPosterior"
  }
}

#' @rdname get_metadata
#' @export
S7::method(get_metadata, MVBU_StanfitPosterior) <- function(x, ...) {
  x@metadata
}

#' @rdname get_original_variable_names
#' @export
S7::method(get_original_variable_names, MVBU_StanfitPosterior) <- function(
  x,
  variable = c("group", "group_unique", "category", "response_category", "cues"),
  ...
) {
  orig <- x@metadata$original_variable_names
  if (is.null(orig)) {
    orig <- list(
      group = get_group_labels(x),
      group_unique = get_group_labels(x),
      category = get_category_labels(x),
      response_category = "response_category",
      cues = get_cue_labels(x)
    )
  }
  if (missing(variable)) {
    return(orig)
  }
  valid_vars <- c("group", "group_unique", "category", "response_category", "cues")
  selected <- match.arg(variable, valid_vars, several.ok = TRUE)
  if (length(selected) == 1L) {
    orig[[selected]]
  } else {
    orig[selected]
  }
}

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
#' @noRd
.as_cognitive_model_from_posterior <- function(model) {
  pars_sum <- model@draws$summary
  if (is.null(pars_sum) || nrow(pars_sum) == 0L) {
    .stop("No parameter summary found in MVBU_StanfitPosterior object.")
  }

  avail_groups <- if (is.factor(pars_sum$group)) levels(pars_sum$group) else unique(pars_sum$group)
  post_groups <- setdiff(avail_groups, "prior")
  target_group <- if (length(post_groups) > 0) post_groups[1] else "prior"
  pars_group <- pars_sum[pars_sum$group == target_group, ]

  cues <- get_cue_labels(model)
  cats <- get_category_labels(model)
  fam <- get_model_family(model)
  is_nix <- grepl("NIX", fam, ignore.case = TRUE) && !grepl("MNIX", fam, ignore.case = TRUE)

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

  lapse <- if ("lapse_rate" %in% names(pars_group)) pars_group$lapse_rate[1] else 0

  if (is_nix) {
    new_nix_ideal_adaptor(
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

#' @rdname evaluate_model
#' @export
S7::method(evaluate_model, MVBU_StanfitPosterior) <- function(
  model,
  x = NULL,
  response_category = NULL,
  method = "log_lik",
  decision_rule = if (identical(method, "accuracy")) "criterion" else "proportional",
  return_by_x = FALSE,
  ...
) {
  cog_model <- .as_cognitive_model_from_posterior(model)

  if (is.null(x) || is.null(response_category)) {
    test_df <- model@metadata$test_data
    if (is.null(test_df) || nrow(test_df) == 0L) {
      .stop("No test data found in the MVBU_StanfitPosterior object. Please supply x and response_category.")
    }
    cues <- get_cue_labels(model)
    if (!"response_category" %in% names(test_df)) {
      .stop("No 'response_category' column found in test data frame.")
    }
    x_mat <- as.matrix(test_df[, cues, drop = FALSE])
    response_category <- test_df[["response_category"]]
    valid_idx <- !is.na(response_category)
    if (any(!valid_idx)) {
      x_mat <- x_mat[valid_idx, , drop = FALSE]
      response_category <- response_category[valid_idx]
    }
    x <- x_mat
  }

  evaluate_model(
    cog_model,
    x = x,
    response_category = response_category,
    method = method,
    decision_rule = decision_rule,
    ...,
    return_by_x = return_by_x
  )
}

#' @rdname likelihood
#' @export
S7::method(
  likelihood,
  list(MVBU_StanfitPosterior, S7::class_any, S7::class_any)
) <- function(x, new_data, categories, ...) {
  cog_model <- .as_cognitive_model_from_posterior(x)
  likelihood(cog_model, new_data, categories, ...)
}

#' @rdname categorize
#' @export
S7::method(
  categorize,
  list(MVBU_StanfitPosterior, S7::class_any, S7::class_any)
) <- function(x, new_data, decision_rule, simplify = NULL, ...) {
  cog_model <- .as_cognitive_model_from_posterior(x)
  categorize(cog_model, new_data, decision_rule, simplify = simplify, ...)
}

#' @rdname get_category_likelihood_function
#' @export
S7::method(
  get_category_likelihood_function,
  MVBU_StanfitPosterior
) <- function(x, ...) {
  cog_model <- .as_cognitive_model_from_posterior(x)
  get_category_likelihood_function(cog_model, ...)
}

#' @export
S7::method(print, MVBU_StanfitPosterior) <- function(x, ...) {
  cat("MVBeliefUpdatr Lightweight Stanfit Posterior\n")
  cat("Model Family     :", get_model_family(x), "\n")
  cat("Group Context    :", paste(get_group_labels(x), collapse = ", "), "\n")
  cat("Cues             :", paste(get_cue_labels(x), collapse = ", "), "\n")
  cat("Categories       :", paste(get_category_labels(x), collapse = ", "), "\n")
  invisible(x)
}

#' @exportS3Method base::summary
summary.MVBU_StanfitPosterior <- function(object, ...) {
  print(object)
  if (!is.null(object@draws$summary)) {
    cat("\nParameter Summary:\n")
    print(object@draws$summary)
  }

  cog_model <- tryCatch(
    .as_cognitive_model_from_posterior(object),
    error = function(e) NULL
  )
  if (!is.null(cog_model)) {
    cat("\nCategory moments:\n")
    exp_mu <- get_expected_mu(cog_model)
    exp_sig <- get_expected_sigma(cog_model)
    marg_sig <- get_marginal_sigma(cog_model)
    print(list(mu = exp_mu, Sigma_exp = exp_sig, Sigma_marg = marg_sig))
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
  cog_model <- .as_cognitive_model_from_posterior(x)
  get_expected_category_statistic(
    cog_model,
    categories = categories,
    groups = groups,
    statistic = statistic,
    ...
  )
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
  cog_model <- .as_cognitive_model_from_posterior(x)
  get_marginal_category_statistic(
    cog_model,
    categories = categories,
    groups = groups,
    statistic = statistic,
    ...
  )
}

#' @rdname get_parameters
#' @export
S7::method(get_expected_mu, MVBU_StanfitPosterior) <- function(x, ...) {
  get_expected_category_statistic(x, statistic = "mu", ...)
}

#' @rdname get_parameters
#' @export
S7::method(get_expected_sigma, MVBU_StanfitPosterior) <- function(x, ...) {
  get_expected_category_statistic(x, statistic = "Sigma", ...)
}

#' @rdname get_parameters
#' @export
S7::method(get_marginal_mu, MVBU_StanfitPosterior) <- function(x, ...) {
  get_marginal_category_statistic(x, statistic = "mu", ...)
}

#' @rdname get_parameters
#' @export
S7::method(get_marginal_sigma, MVBU_StanfitPosterior) <- function(x, ...) {
  get_marginal_category_statistic(x, statistic = "Sigma", ...)
}

#' @rdname get_sufficient_category_statistics
#' @export
S7::method(get_sufficient_category_statistics, MVBU_StanfitPosterior) <- function(
  x,
  categories = NULL,
  groups = NULL,
  untransform_cues = FALSE,
  ...
) {
  if (is.null(categories)) categories <- get_category_labels(x)
  if (is.null(groups)) groups <- get_group_labels(x, include_prior = FALSE)
  meta <- get_metadata(x)
  if (!is.null(meta$sufficient_statistics)) {
    res <- meta$sufficient_statistics
    if ("category" %in% names(res) && !is.null(categories)) {
      res <- res[res$category %in% categories, , drop = FALSE]
    }
    if ("group" %in% names(res) && !is.null(groups)) {
      res <- res[res$group %in% groups, , drop = FALSE]
    }
    return(as.data.frame(res, stringsAsFactors = FALSE))
  }
  .stop("Sufficient exposure statistics are not available in this MVBU_StanfitPosterior object.")
}
