#' @include S7-core-classes.R
#' @include S7-generics.R
#' @include S7-transform-information.R
#' @include S7-staninput.R
#' @include S7-stanfit-input.R
#' @include S7-stanfit.R
#' @include S7-stanfit-methods.R
NULL

S7::method(get_staninput, IdealAdaptorStanfitInput) <- function(x) {
  x@staninput
}

S7::method(
  set_staninput,
  list(S7::class_any, S7::class_any)
) <- function(x, staninput) {
  .stop("x must be an IdealAdaptorStanfit or IdealAdaptorStanfitInput object.")
}

S7::method(
  set_staninput,
  list(MVBU_Stanfit, S7::class_any)
) <- function(x, staninput) {
  if (!S7::S7_inherits(staninput, MVBU_Staninput)) {
    .stop("staninput must be an MVBU_Staninput object.")
  }

  x@staninput <- staninput
  x
}

S7::method(
  set_staninput,
  list(IdealAdaptorStanfitInput, S7::class_any)
) <- function(x, staninput) {
  if (!S7::S7_inherits(staninput, MVBU_Staninput)) {
    .stop("staninput must be an MVBU_Staninput object.")
  }

  x@staninput <- staninput
  x
}

S7::method(get_transform_information, IdealAdaptorStanfitInput) <- function(x) {
  x@transform_information
}

S7::method(get_cue_labels, IdealAdaptorStanfitInput) <- function(
  x,
  indices = NULL,
  ...
) {
  cues <- if (!is.null(x@metadata$label_information$cue)) {
    x@metadata$label_information$cue
  } else {
    character(0)
  }
  if (!is.null(indices)) cues[indices] else cues
}

S7::method(get_category_labels, IdealAdaptorStanfitInput) <- function(
  x,
  indices = NULL,
  ...
) {
  cats <- if (!is.null(x@metadata$label_information$category)) {
    x@metadata$label_information$category
  } else {
    character(0)
  }
  if (!is.null(indices)) cats[indices] else cats
}

S7::method(get_group_labels, IdealAdaptorStanfitInput) <- function(
  x,
  indices = NULL,
  include_prior = FALSE,
  ...
) {
  label_info <- if (!is.null(x@metadata$label_information)) {
    x@metadata$label_information
  } else {
    list()
  }
  grp_unique_attr <- attr(x@data, "group.unique")
  grp_col <- if (!is.null(grp_unique_attr) && grp_unique_attr %in% names(x@data)) {
    grp_unique_attr
  } else {
    attr(x@data, "group")
  }
  expected_levels <- if (!is.null(grp_col) && is.factor(x@data[[grp_col]])) {
    levels(x@data[[grp_col]])
  } else {
    character(0)
  }
  grps <- if (length(label_info$group) > 0 &&
    (length(expected_levels) == 0 || length(label_info$group) == length(expected_levels))) {
    label_info$group
  } else if (length(expected_levels) > 0) {
    expected_levels
  } else if (length(label_info$group) > 0) {
    label_info$group
  } else {
    character(0)
  }
  if (include_prior) grps <- append("prior", grps)
  if (!is.null(indices)) grps[indices] else grps
}

S7::method(get_labels, IdealAdaptorStanfitInput) <- function(x, ...) {
  list(
    cue = get_cue_labels(x, ...),
    category = get_category_labels(x, ...),
    group = get_group_labels(x, ...)
  )
}

S7::method(get_model_type, IdealAdaptorStaninput) <- function(x) {
  if (S7::S7_inherits(x, NIX_IdealAdaptorStaninput)) {
    "NIX_ideal_adaptor"
  } else if (S7::S7_inherits(x, MNIX_IdealAdaptorStaninput)) {
    "MNIX_ideal_adaptor"
  } else if (S7::S7_inherits(x, NIW_IdealAdaptorStaninput)) {
    "NIW_ideal_adaptor"
  } else {
    "ideal_adaptor"
  }
}

S7::method(get_model_type, IdealAdaptorStanfitInput) <- function(x) {
  if (!is.null(x@staninput)) {
    get_model_type(x@staninput)
  } else {
    "ideal_adaptor"
  }
}

S7::method(get_data, IdealAdaptorStanfitInput) <- function(
  x,
  groups = get_group_labels(x, include_prior = FALSE),
  .rename_to_MVB_default = FALSE,
  ...
) {
  data <- x@data
  group.unique <- attr(data, "group.unique")
  if (.rename_to_MVB_default) {
    group <- attr(data, "group")
    category <- attr(data, "category")
    cues <- attr(data, "cues")
    response <- attr(data, "response")

    data <- data %>%
      dplyr::rename(
        group.unique = !!rlang::sym(group.unique),
        group = !!rlang::sym(group),
        category = !!rlang::sym(category),
        response = !!rlang::sym(response)
      )

    for (c in seq_along(cues)) {
      data <- data %>%
        dplyr::rename(!!rlang::sym(paste0("cue", c)) := !!rlang::sym(cues[c]))
    }
    group.unique <- "group.unique"
  }

  if (!is.null(group.unique) && group.unique %in% names(data)) {
    data %>% dplyr::filter(!!rlang::sym(group.unique) %in% groups)
  } else if ("group" %in% names(data)) {
    data %>% dplyr::filter(.data$group %in% groups)
  } else {
    data
  }
}

S7::method(get_exposure_data, IdealAdaptorStanfitInput) <- function(
  x,
  groups = get_group_labels(x, include_prior = FALSE),
  ...
) {
  get_data(x, groups = groups, ...) %>%
    dplyr::filter(.data$Phase == "exposure")
}

S7::method(get_test_data, IdealAdaptorStanfitInput) <- function(
  x,
  groups = get_group_labels(x, include_prior = FALSE),
  .recover_from_staninput = FALSE,
  ...
) {
  has_data <- !is.null(x@data) && is.data.frame(x@data) && nrow(x@data) > 0 && "Phase" %in% names(x@data)
  if (has_data) {
    df <- get_data(x, groups = groups, ...) %>%
      dplyr::filter(.data$Phase == "test")
    if (!"group" %in% names(df)) {
      grp_col <- attr(x@data, "group.unique") %||% attr(x@data, "group")
      if (!is.null(grp_col) && grp_col %in% names(df)) {
        df$group <- df[[grp_col]]
      }
    }
    df
  } else if (.recover_from_staninput && !is.null(get_staninput(x))) {
    stanvals <- get_staninput(x)@values
    df <- tibble::as_tibble(
      cbind(stanvals$x_test, stanvals$z_test_counts),
      .name_repair = "minimal"
    )
    df <- df %>%
      dplyr::mutate(
        group.id = stanvals$y_test,
        group = factor(
          attr(stanvals$y_test, "levels")[.data$group.id],
          levels = attr(stanvals$y_test, "levels")
        )
      ) %>%
      dplyr::filter(.data$group %in% groups)
    df
  } else {
    tibble::tibble()
  }
}

S7::method(
  get_exposure_category_statistic,
  IdealAdaptorStanfitInput
) <- function(
  x,
  categories = get_category_labels(x),
  groups = get_group_labels(x, include_prior = FALSE),
  statistic = c("n", "mean", "css", "uss", "cov"),
  untransform_cues = FALSE,
  ...
) {
  .compute_exposure_category_statistic(
    x,
    categories = categories,
    groups = groups,
    statistic = statistic,
    untransform_cues = untransform_cues,
    ...
  )
}

S7::method(
  get_exposure_category_mean,
  IdealAdaptorStanfitInput
) <- function(x, ...) {
  get_exposure_category_statistic(x, ..., statistic = "mean")
}

S7::method(
  get_exposure_category_css,
  IdealAdaptorStanfitInput
) <- function(x, ...) {
  get_exposure_category_statistic(x, ..., statistic = "css")
}

S7::method(
  get_exposure_category_uss,
  IdealAdaptorStanfitInput
) <- function(x, ...) {
  get_exposure_category_statistic(x, ..., statistic = "uss")
}

S7::method(
  get_exposure_category_cov,
  IdealAdaptorStanfitInput
) <- function(x, ...) {
  get_exposure_category_statistic(x, ..., statistic = "cov")
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
  if (!is.null(x@data) && is.data.frame(x@data)) {
    cat("  Data: ", nrow(x@data), " rows\n", sep = "")
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
#' @param ... Additional arguments.
#'
#' @return An object of class \code{Summary_IdealAdaptorStanfitInput} containing
#'   exposure statistics and test data summaries.
S7::method(summary, IdealAdaptorStanfitInput) <- function(object, ...) {
  cats <- get_category_labels(object)
  cues <- get_cue_labels(object)
  grps <- get_group_labels(object, include_prior = FALSE)

  # Exposure sufficient statistics
  exp_stats <- tryCatch(
    get_exposure_category_statistic(object, statistic = c("n", "mean", "cov")),
    error = function(e) tibble::tibble()
  )

  if (nrow(exp_stats) > 0) {
    for (col in names(exp_stats)) {
      if (is.list(exp_stats[[col]])) {
        all_scalars <- all(vapply(exp_stats[[col]], function(x) length(x) == 1L && is.numeric(x), logical(1)))
        if (all_scalars) {
          exp_stats[[col]] <- vapply(exp_stats[[col]], as.numeric, numeric(1))
        } else {
          exp_stats[[col]] <- vapply(exp_stats[[col]], function(x) {
            if (is.matrix(x)) paste0("matrix(", nrow(x), ", ", ncol(x), ")")
            else if (is.numeric(x)) paste0("(", paste(format(x, digits = 3), collapse = ", "), ")")
            else as.character(x)
          }, character(1))
        }
      }
    }
  }

  # Test data summary
  test_data <- tryCatch(
    get_test_data(object, .recover_from_staninput = TRUE),
    error = function(e) tibble::tibble()
  )

  vals <- if (!is.null(object@staninput)) object@staninput@values else NULL
  n_test_obs <- if (!is.null(vals$z_test_counts)) {
    sum(vals$z_test_counts)
  } else if (!is.null(vals$N_obs_test)) {
    vals$N_obs_test
  } else if (nrow(test_data) > 0) {
    nrow(test_data)
  } else {
    0L
  }

  test_by_group <- NULL
  if (!is.null(vals$y_test) && length(vals$y_test) > 0) {
    grp_vector <- if (is.numeric(vals$y_test)) grps[vals$y_test] else vals$y_test
    if (!is.null(vals$z_test_counts)) {
      row_counts <- rowSums(vals$z_test_counts)
      test_by_group <- tapply(row_counts, grp_vector, sum)
    } else {
      test_by_group <- table(grp_vector)
    }
  } else if (nrow(test_data) > 0 && "group" %in% names(test_data)) {
    test_by_group <- table(test_data$group)
  }

  cue_ranges <- list()
  if (!is.null(vals$x_test) && !is.null(vals$y_test)) {
    x_mat <- as.matrix(vals$x_test)
    grp_vector <- if (is.numeric(vals$y_test)) grps[vals$y_test] else vals$y_test
    for (i in seq_along(cues)) {
      cue_name <- cues[i]
      if (i <= ncol(x_mat)) {
        cue_ranges[[cue_name]] <- list()
        for (grp in grps) {
          idx <- which(grp_vector == grp)
          if (length(idx) > 0) {
            cue_ranges[[cue_name]][[grp]] <- range(x_mat[idx, i], na.rm = TRUE)
          }
        }
      }
    }
  } else if (nrow(test_data) > 0 && "group" %in% names(test_data)) {
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
