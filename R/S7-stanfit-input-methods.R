#' @include S7-class.R
#' @include S7-generics.R
#' @include S7-transform-information.R
#' @include S7-staninput.R
#' @include S7-stanfit-input.R
#' @include S7-stanfit.R
NULL

S7::method(
  set_staninput,
  list(S7::class_any, S7::class_any)
) <- function(x, staninput, ...) {
  .stop("x must be an IdealAdaptorStanfit object.")
}

S7::method(
  set_staninput,
  list(MVBU_Stanfit, S7::class_any)
) <- function(x, staninput, ...) {
  if (!S7::S7_inherits(staninput, MVBU_Staninput)) {
    .stop("staninput must be an MVBU_Staninput object.")
  }

  x@staninput <- staninput
  x
}


S7::method(get_original_variable_names, IdealAdaptorStanfitInput) <- function(
  x,
  variable = c("group", "group_unique", "category", "response_category", "cues"),
  ...
) {
  orig <- x@metadata$original_variable_names
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

#' Internal helper to filter, subsample, and format stanfit data
#' @noRd
.get_stanfit_data_impl <- function(
  x,
  space = c("original", "transformed"),
  groups = NULL,
  categories = NULL,
  response_categories = NULL,
  n_samples = NULL,
  original_names = FALSE,
  phase = NULL
) {
  space <- match.arg(space)
  .assert_logical_scalar(original_names)
  if (!is.null(groups)) .assert_character(groups)
  if (!is.null(categories)) .assert_character(categories)
  if (!is.null(response_categories)) .assert_character(response_categories)
  if (!is.null(n_samples)) {
    .assert_numeric_scalar(n_samples)
    .assert_true(n_samples > 0L, msg = "n_samples must be a positive integer.")
    n_samples <- as.integer(n_samples)
  }

  data <- if (is.data.frame(x@data)) {
    x@data
  } else if (identical(space, "transformed")) {
    x@data$transformed
  } else {
    x@data$original
  }
  if (is.null(data) || !is.data.frame(data) || nrow(data) == 0L) {
    if (!is.null(groups) || !is.null(categories) || !is.null(response_categories)) {
      .stop("Object does not contain any data to filter.")
    }
    return(tibble::tibble())
  }

  # Renaming if original_names = TRUE
  if (isTRUE(original_names)) {
    orig <- get_original_variable_names(x)
    rename_map <- character(0)
    if (!is.null(orig$group) && !identical(orig$group, "group") && "group" %in% names(data)) {
      rename_map[orig$group] <- "group"
    }
    if (!is.null(orig$group_unique) && !identical(orig$group_unique, "group_unique") && !identical(orig$group_unique, orig$group) && "group_unique" %in% names(data)) {
      rename_map[orig$group_unique] <- "group_unique"
    }
    if (!is.null(orig$category) && !identical(orig$category, "category") && "category" %in% names(data)) {
      rename_map[orig$category] <- "category"
    }
    if (!is.null(orig$response_category) && !identical(orig$response_category, "response_category") && "response_category" %in% names(data)) {
      rename_map[orig$response_category] <- "response_category"
    }
    rename_map <- rename_map[!duplicated(names(rename_map))]
    if (length(rename_map) > 0L) {
      match_idx <- match(rename_map, names(data))
      valid <- !is.na(match_idx)
      names(data)[match_idx[valid]] <- names(rename_map)[valid]
    }
    grp_col <- if (!is.null(orig$group_unique) && orig$group_unique %in% names(data)) orig$group_unique else orig$group
    cat_col <- if (!is.null(orig$category) && orig$category %in% names(data)) orig$category else "category"
    resp_col <- if (!is.null(orig$response_category) && orig$response_category %in% names(data)) orig$response_category else "response_category"
  } else {
    grp_col <- if ("group_unique" %in% names(data)) "group_unique" else "group"
    cat_col <- "category"
    resp_col <- "response_category"
  }

  context <- if (!is.null(phase)) paste0(phase, " ") else ""

  # 1. Filter phase immediately if requested (message if phase contains zero data)
  if (!is.null(phase)) {
    .assert_character_scalar(phase)
    phase_mask <- data$Phase == phase
    if (!any(phase_mask)) {
      .message(sprintf("Object contains no %s data.", phase))
      return(tibble::tibble())
    }
    data <- data[phase_mask, , drop = FALSE]
  }

  # 2. Filter groups if requested (error if non-existing in object, message if zero data in phase)
  if (!is.null(groups)) {
    groups <- .validate_requested_labels(
      groups,
      get_group_labels(x, include_prior = FALSE),
      label_type = "group"
    )
    avail_grps <- unique(stats::na.omit(as.character(data[[grp_col]])))
    zero_data_grps <- setdiff(as.character(groups), avail_grps)
    if (length(zero_data_grps) > 0L) {
      .message(
        sprintf(
          "No %sdata found for group(s): %s",
          context,
          paste(zero_data_grps, collapse = ", ")
        )
      )
    }
    data <- data[as.character(data[[grp_col]]) %in% as.character(groups), , drop = FALSE]
  }

  # 3. Filter categories if requested (error if non-existing in object, message if zero data)
  if (!is.null(categories)) {
    categories <- .validate_requested_labels(
      categories,
      get_category_labels(x),
      label_type = "category"
    )
    if (cat_col %in% names(data)) {
      if (identical(phase, "exposure")) {
        avail_cats <- unique(stats::na.omit(as.character(data[[cat_col]])))
        zero_data_cats <- setdiff(as.character(categories), avail_cats)
        if (length(zero_data_cats) > 0L) {
          .message(
            sprintf(
              "No exposure data found for category(ies): %s",
              paste(zero_data_cats, collapse = ", ")
            )
          )
        }
        data <- data[as.character(data[[cat_col]]) %in% as.character(categories), , drop = FALSE]
      } else {
        # In get_data across phases, validate against exposure rows
        exp_mask <- data$Phase == "exposure"
        avail_cats <- unique(stats::na.omit(as.character(data[[cat_col]][exp_mask])))
        zero_data_cats <- setdiff(as.character(categories), avail_cats)
        if (length(zero_data_cats) > 0L && any(exp_mask)) {
          .message(
            sprintf(
              "No exposure data found for category(ies): %s",
              paste(zero_data_cats, collapse = ", ")
            )
          )
        }
        keep_exp <- exp_mask & (as.character(data[[cat_col]]) %in% as.character(categories))
        data <- data[!exp_mask | keep_exp, , drop = FALSE]
      }
    }
  }

  # 4. Filter response_categories if requested (error if non-existing, message if zero data)
  if (!is.null(response_categories)) {
    response_categories <- .validate_requested_labels(
      response_categories,
      get_response_category_labels(x),
      label_type = "response_category"
    )
    if (resp_col %in% names(data)) {
      if (identical(phase, "test")) {
        avail_rcats <- unique(stats::na.omit(as.character(data[[resp_col]])))
        zero_data_rcats <- setdiff(as.character(response_categories), avail_rcats)
        if (length(zero_data_rcats) > 0L) {
          .message(
            sprintf(
              "No test data found for response_category(ies): %s",
              paste(zero_data_rcats, collapse = ", ")
            )
          )
        }
        data <- data[as.character(data[[resp_col]]) %in% as.character(response_categories), , drop = FALSE]
      } else {
        # In get_data across phases, validate against test rows
        test_mask <- data$Phase == "test"
        avail_rcats <- unique(stats::na.omit(as.character(data[[resp_col]][test_mask])))
        zero_data_rcats <- setdiff(as.character(response_categories), avail_rcats)
        if (length(zero_data_rcats) > 0L && any(test_mask)) {
          .message(
            sprintf(
              "No test data found for response_category(ies): %s",
              paste(zero_data_rcats, collapse = ", ")
            )
          )
        }
        keep_test <- test_mask & (as.character(data[[resp_col]]) %in% as.character(response_categories))
        data <- data[!test_mask | keep_test, , drop = FALSE]
      }
    }
  }

  # 5. Subsampling
  if (!is.null(n_samples) && nrow(data) > n_samples) {
    data <- data[sample.int(nrow(data), n_samples), , drop = FALSE]
  }

  tibble::as_tibble(data)
}

#' @rdname get_data
#' @export
S7::method(get_data, IdealAdaptorStanfitInput) <- function(
  x,
  space = c("original", "transformed"),
  groups = NULL,
  categories = NULL,
  response_categories = NULL,
  n_samples = NULL,
  original_names = FALSE,
  ...
) {
  space <- match.arg(space)
  if (is.null(groups)) {
    lbls <- get_group_labels(x, include_prior = FALSE)
    if (length(lbls) > 0L) {
      groups <- lbls
    }
  }
  .get_stanfit_data_impl(
    x = x,
    space = space,
    groups = groups,
    categories = categories,
    response_categories = response_categories,
    n_samples = n_samples,
    original_names = original_names,
    phase = NULL
  )
}

#' @rdname get_data
#' @export
S7::method(get_exposure_data, IdealAdaptorStanfitInput) <- function(
  x,
  space = c("original", "transformed"),
  groups = NULL,
  categories = NULL,
  n_samples = NULL,
  original_names = FALSE,
  ...
) {
  space <- match.arg(space)
  .get_stanfit_data_impl(
    x = x,
    space = space,
    groups = groups,
    categories = categories,
    response_categories = NULL,
    n_samples = n_samples,
    original_names = original_names,
    phase = "exposure"
  )
}

#' @rdname get_data
#' @export
S7::method(get_test_data, IdealAdaptorStanfitInput) <- function(
  x,
  space = c("original", "transformed"),
  groups = NULL,
  response_categories = NULL,
  n_samples = NULL,
  original_names = FALSE,
  ...
) {
  space <- match.arg(space)
  .get_stanfit_data_impl(
    x = x,
    space = space,
    groups = groups,
    categories = NULL,
    response_categories = response_categories,
    n_samples = n_samples,
    original_names = original_names,
    phase = "test"
  )
}

S7::method(print, MVBU_Staninput) <- function(x, ...) {
  cls <- class(x)[1]
  cat("<", cls, ">\n", sep = "")
  cats <- get_category_labels(x)
  cues <- get_cue_labels(x)
  grps <- get_group_labels(x, include_prior = FALSE)
  if (length(cats) > 0) cat("  Categories (", length(cats), "): ", paste(cats, collapse = ", "), "\n", sep = "")
  if (length(cues) > 0) cat("  Cues (", length(cues), "): ", paste(cues, collapse = ", "), "\n", sep = "")
  if (length(grps) > 0) cat("  Groups (", length(grps), "): ", paste(grps, collapse = ", "), "\n", sep = "")
  invisible(x)
}

S7::method(print, IdealAdaptorStanfitInput) <- function(x, ...) {
  cls <- class(x)[1]
  cat("<", cls, ">\n", sep = "")
  cats <- get_category_labels(x)
  cues <- get_cue_labels(x)
  grps <- get_group_labels(x, include_prior = FALSE)
  if (length(cats) > 0) cat("  Categories (", length(cats), "): ", paste(cats, collapse = ", "), "\n", sep = "")
  if (length(cues) > 0) cat("  Cues (", length(cues), "): ", paste(cues, collapse = ", "), "\n", sep = "")
  if (length(grps) > 0) cat("  Groups (", length(grps), "): ", paste(grps, collapse = ", "), "\n", sep = "")
  orig_data <- get_data(x, space = "original")
  if (!is.null(orig_data) && is.data.frame(orig_data)) {
    cat("  Data: ", nrow(orig_data), " rows\n", sep = "")
  }
  invisible(x)
}

#' An S7 class for IdealAdaptorStanfitInput summaries
#'
#' @keywords internal
Summary_IdealAdaptorStanfitInput <- S7::new_class(
  "Summary_IdealAdaptorStanfitInput",
  parent = MVBU_Object,
  properties = list(
    object_class = S7::class_character,
    categories = S7::class_character,
    cues = S7::class_character,
    groups = S7::class_character,
    exposure_statistics = S7::class_any,
    test_summary = S7::class_list
  )
)

#' Summarize an IdealAdaptorStanfitInput object
#'
#' Provides a detailed overview of the exposure sufficient statistics across
#' group and category combinations as well as key properties of the test dataset.
#'
#' @param object An \code{\link{IdealAdaptorStanfitInput}} object.
#' @param categories Optional character vector of category labels to summarize.
#'   Defaults to all categories in \code{object}.
#' @param groups Optional character vector of group labels to summarize.
#'   Defaults to all groups in \code{object}.
#' @param ... Additional arguments.
#'
#' @rdname summary-methods
#' @exportS3Method base::summary
summary.IdealAdaptorStanfitInput <- function(
  object,
  categories = NULL,
  groups = NULL,
  ...
) {
  cats <- .validate_requested_labels(
    categories,
    get_category_labels(object),
    label_type = "category"
  )
  grps <- .validate_requested_labels(
    groups,
    get_group_labels(object, include_prior = FALSE),
    label_type = "group"
  )

  cues <- get_cue_labels(object)

  # Exposure sufficient statistics
  exp_stats <- tryCatch(
    get_exposure_category_statistics(
      object,
      categories = cats,
      groups = grps,
      stats = c("n", "mean", "cov"),
      simplify = FALSE
    ),
    error = function(e) tibble::tibble()
  )

  if (nrow(exp_stats) > 0) {
    if ("n" %in% names(exp_stats)) {
      exp_stats <- exp_stats[exp_stats$n > 0L, , drop = FALSE]
    }
    for (col in names(exp_stats)) {
      if (is.list(exp_stats[[col]])) {
        all_scalars <- all(vapply(exp_stats[[col]], function(x) length(x) == 1L && is.numeric(x), logical(1)))
        if (all_scalars) {
          exp_stats[[col]] <- vapply(exp_stats[[col]], as.numeric, numeric(1))
        } else {
          exp_stats[[col]] <- vapply(exp_stats[[col]], function(x) {
            if (is.matrix(x)) {
              paste0("matrix(", nrow(x), ", ", ncol(x), ")")
            } else if (is.numeric(x)) {
              paste0("(", paste(format(x, digits = 3), collapse = ", "), ")")
            } else {
              as.character(x)
            }
          }, character(1))
        }
      }
    }
  }

  # Test data summary
  test_data <- tryCatch(
    get_test_data(object, groups = grps),
    error = function(e) tibble::tibble()
  )

  test_by_group <- NULL
  if (nrow(test_data) > 0 && "group" %in% names(test_data)) {
    test_by_group <- table(test_data$group)
    test_by_group <- test_by_group[names(test_by_group) %in% grps]
    n_test_obs <- sum(test_by_group)
  } else {
    n_test_obs <- 0L
  }

  cue_ranges <- list()
  if (nrow(test_data) > 0 && "group" %in% names(test_data)) {
    for (cue in cues) {
      if (cue %in% names(test_data) && is.numeric(test_data[[cue]])) {
        cue_ranges[[cue]] <- list()
        for (grp in grps) {
          sub_df <- test_data[test_data$group == grp, , drop = FALSE]
          if (nrow(sub_df) > 0) {
            cue_ranges[[cue]][[grp]] <- range(sub_df[[cue]], na.rm = TRUE)
          }
        }
      }
    }
  }

  Summary_IdealAdaptorStanfitInput(
    object_class = class(object)[1],
    categories = cats,
    cues = cues,
    groups = grps,
    exposure_statistics = exp_stats,
    test_summary = list(
      n_observations = n_test_obs,
      by_group = test_by_group,
      cue_ranges = cue_ranges
    )
  )
}
S7::method(
  summary,
  IdealAdaptorStanfitInput
) <- summary.IdealAdaptorStanfitInput

S7::method(print, Summary_IdealAdaptorStanfitInput) <- function(x, ...) {
  cat("<", x@object_class, " summary>\n", sep = "")
  if (length(x@categories) > 0) cat("  Categories (", length(x@categories), "): ", paste(x@categories, collapse = ", "), "\n", sep = "")
  if (length(x@cues) > 0) cat("  Cues (", length(x@cues), "): ", paste(x@cues, collapse = ", "), "\n", sep = "")
  if (length(x@groups) > 0) cat("  Groups (", length(x@groups), "): ", paste(x@groups, collapse = ", "), "\n", sep = "")

  cat("\nExposure sufficient statistics:\n")
  if (nrow(x@exposure_statistics) > 0) {
    print(x@exposure_statistics, ...)
  } else {
    cat("  (No exposure statistics found)\n")
  }

  cat("\nTest data summary:\n")
  cat("  Total observations: ", x@test_summary$n_observations, "\n", sep = "")
  if (!is.null(x@test_summary$by_group) && length(x@test_summary$by_group) > 0) {
    grp_str <- paste(paste(names(x@test_summary$by_group), x@test_summary$by_group, sep = ": "), collapse = ", ")
    cat("  Observations by group: ", grp_str, "\n", sep = "")
  }
  if (length(x@test_summary$cue_ranges) > 0) {
    cat("  Cue ranges by group:\n")
    for (cue in names(x@test_summary$cue_ranges)) {
      grp_rngs <- x@test_summary$cue_ranges[[cue]]
      rng_strs <- vapply(names(grp_rngs), function(grp) {
        rng <- grp_rngs[[grp]]
        sprintf("%s: [%s, %s]", grp, format(rng[1], digits = 3), format(rng[2], digits = 3))
      }, character(1))
      cat("    ", cue, ": ", paste(rng_strs, collapse = ", "), "\n", sep = "")
    }
  }
  invisible(x)
}

#' @rdname get_family
#' @export
S7::method(get_family_name, IdealAdaptorStaninput) <- function(x, ...) {
  x@family
}

#' @rdname get_family
#' @export
S7::method(get_family, IdealAdaptorStaninput) <- function(x, ...) {
  .get_model_family(x@family)
}
