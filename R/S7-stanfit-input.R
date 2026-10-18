#' @include asserts.R
#' @include S7-class.R
#' @include S7-transform-information.R
#' @include S7-staninput.R
#' @include S7-prior.R
NULL

#' Specify control parameters for new_ideal_adaptor_stanfit_input()
#'
#' This function is used to specify control parameters for the
#' `new_ideal_adaptor_stanfit_input()` function, and to provide reasonable
#' defaults for any of the unspecified parameters.
#'
#' @param tau_scale `r lifecycle::badge("deprecated")` Deprecated. Use
#'   `priors = prior(..., class = "tau")` in \code{\link{fit_ideal_adaptor}}
#'   instead.
#' @param L_omega_eta `r lifecycle::badge("deprecated")` Deprecated. Use
#'   `priors = prior(..., class = "cor")` in \code{\link{fit_ideal_adaptor}}
#'   instead.
#' @param split_loglik_per_observation Split log likelihood per observation.
#' @param transform_type Affine transformation type applied to the data.
#'
#' @return A list of control parameters.
#' @rdname MVBU_StanfitInput
#' @export
control_staninput <- function(
  tau_scale = NULL,
  L_omega_eta = NULL,
  split_loglik_per_observation = 0,
  transform_type = c(
    "identity", "center", "standardize", "PCA whiten", "ZCA whiten"
  )[3]
) {
  if (!is.null(tau_scale)) {
    lifecycle::deprecate_warn(
      "0.9.0",
      "control_staninput(tau_scale)",
      details = paste0(
        "Use priors = prior(..., class = 'tau') in fit_ideal_adaptor() ",
        "instead."
      )
    )
  }
  if (!is.null(L_omega_eta)) {
    lifecycle::deprecate_warn(
      "0.9.0",
      "control_staninput(L_omega_eta)",
      details = paste0(
        "Use priors = prior(..., class = 'cor') in fit_ideal_adaptor() ",
        "instead."
      )
    )
  }
  list(
    tau_scale = tau_scale,
    L_omega_eta = L_omega_eta,
    split_loglik_per_observation = split_loglik_per_observation,
    transform_type = transform_type
  )
}

#' MVBeliefUpdatr Stanfit Input Classes and Constructors
#'
#' An S7 Class for Ideal-Adaptor Fit-Input Objects
#'
#' @param data A named \code{list} with \code{$original} and
#'   \code{$transformed} data frames, or a single \code{data.frame} containing
#'   input data. Default: \code{list()}.
#' @param transform_information An \code{\link{MVBU_TransformInformation}}
#'   object. Default: \code{NULL}.
#' @param metadata A named \code{list} containing \code{label_information}.
#'   Default: \code{list()}.
#'
#' @name MVBU_StanfitInput
#' @rdname MVBU_StanfitInput
#' @docType class
#' @export
IdealAdaptorStanfitInput <- S7::new_class(
  "IdealAdaptorStanfitInput",
  package = NULL,
  parent = MVBU_Object,
  properties = list(
    data = S7::class_list,
    transform_information = MVBU_TransformInformation,
    metadata = S7::class_list
  ),
  constructor = function(
    data = list(),
    transform_information = NULL,
    metadata = list()
  ) {
    if (is.null(transform_information)) {
      transform_information <- MVBU_TransformInformation()
    }

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
    if (!is.null(metadata$original_variable_names) &&
      is.list(metadata$original_variable_names)) {
      orig <- metadata$original_variable_names
      metadata$original_variable_names <- list(
        group = if (!is.null(orig$group)) {
          as.character(orig$group)
        } else {
          "group"
        },
        group_unique = if (!is.null(orig$group_unique)) {
          as.character(orig$group_unique)
        } else if (!is.null(orig$group.unique)) {
          as.character(orig$group.unique)
        } else {
          "group"
        },
        category = if (!is.null(orig$category)) {
          as.character(orig$category)
        } else {
          "category"
        },
        response_category = if (!is.null(orig$response_category)) {
          as.character(orig$response_category)
        } else {
          "response_category"
        },
        cues = if (!is.null(orig$cues)) {
          as.character(orig$cues)
        } else {
          character(0)
        }
      )
    }

    S7::new_object(
      MVBU_Object(),
      data = data,
      transform_information = transform_information,
      metadata = metadata
    )
  },
  validator = function(self) {
    if (!is.list(self@data) || is.data.frame(self@data)) {
      return("`data` must be a list containing $original and $transformed.")
    }
    if (!inherits(self@data$original, "data.frame")) {
      return("`data$original` must inherit from data.frame.")
    }
    if (!inherits(self@data$transformed, "data.frame")) {
      return("`data$transformed` must inherit from data.frame.")
    }
    if (!is.null(self@transform_information) &&
      !S7::S7_inherits(self@transform_information, MVBU_TransformInformation)) {
      return(
        "`transform_information` must inherit from MVBU_TransformInformation"
      )
    }
    if (!is.list(self@metadata)) {
      return("`metadata` must be a list")
    }
    NULL
  }
)

#' Construct a streamlined ideal-adaptor Stanfit input object
#'
#' This function composes the exposure and test data in a compact form that can
#' be used as input to the ideal-adaptor Stan models. It accepts empty exposure
#' data as long as test data are present, and returns an S7 fit-input object
#' with the prepared data (in both original and transformed spaces) and
#' transformation metadata.
#'
#' Exposure and test data are checked for the required cue, group, category,
#' and response columns. The exposure and test data are transformed to
#' facilitate model-fitting. This transformation is undone in the stanfit object
#' returned by the Stan model, so that the model parameters are returned in the
#' original cue space.
#'
#' It is important to use `group` to identify individuals that had a specific
#' exposure (or no exposure at all) and specific test trials. You should not
#' use `group` to identify exposure conditions. Setting `group` to an exposure
#' condition results in exposure that concatenates the exposure observations
#' from all subjects in that condition. Typically, this is not what users
#' intend, as it models exposure to the combination of exposure tokens across
#' all subjects, rather than exposure to one set of those exposure tokens. To
#' achieve this intended outcome, use `group_unique` to identify groups with
#' identical exposure. This will correctly use only one unique instance of the
#' observations that any level of `group` receives during exposure.
#'
#' @param exposure A data frame or tibble with exposure data. Each row is
#'   assumed to contain one observation.
#' @param test A data frame or tibble with test data. Each row is assumed to
#'   contain one observation.
#' @param cues Names of columns with cue values. These columns must be present
#'   in both `exposure` and `test`.
#' @param category Name of the column in `exposure` that stores the category
#'   information. (default: "category")
#' @param response Name of the column in `test` that stores the response
#'   information. (default: "response_category")
#' @param group Name of the column in `exposure` and `test` that stores the
#'   grouping information. (default: "group")
#' @param group_unique Optional grouping column used to collapse repeated
#'   exposure conditions.
#' @param control A list of control parameters for the constructor, such as
#'   `list(transform_type = "standardize")` or produced by
#'   `control_staninput()`.
#' @param ... Additional arguments passed to methods or deprecated arguments.
#' @return An object of class `IdealAdaptorStanfitInput` containing the
#'   prepared data (original and transformed) and transform metadata.
#' @export
new_ideal_adaptor_stanfit_input <- function(
  exposure,
  test,
  cues,
  category = "category",
  response = "response_category",
  group = "group",
  group_unique = NULL,
  control = list(transform_type = "standardize"),
  ...
) {
  dots <- list(...)
  if ("group.unique" %in% names(dots)) {
    warning(
      paste0(
        "Argument `group.unique` is deprecated; ",
        "please use `group_unique` instead."
      ),
      call. = FALSE
    )
    if (is.null(group_unique)) {
      group_unique <- dots$group.unique
    }
  }

  if ("stanmodel" %in% names(dots)) {
    lifecycle::deprecate_warn(
      "0.9.0",
      "new_ideal_adaptor_stanfit_input(stanmodel)",
      "fit_ideal_adaptor(stanmodel)"
    )
  }
  if ("fixed_parameters" %in% names(dots)) {
    lifecycle::deprecate_warn(
      "0.9.0",
      "new_ideal_adaptor_stanfit_input(fixed_parameters)",
      "fit_ideal_adaptor(fixed_parameters)"
    )
  }

  .assert_list(control)
  transform_type <- if ("transform_type" %in% names(control)) {
    control$transform_type
  } else {
    "standardize"
  }
  .assert_all(
    is.character(transform_type),
    length(transform_type) == 1L,
    msg = "transform_type must be a single character value."
  )
  .assert_true(
    transform_type %in% c(
      "identity", "center", "standardize", "PCA whiten", "ZCA whiten"
    ),
    msg = paste0(
      "transform_type must be one of identity, center, standardize, ",
      "PCA whiten, or ZCA whiten."
    )
  )

  .assert_non_NA_character(cues)
  cues <- unique(cues)

  exposure <- .prepare_staninput_frame(
    exposure,
    cues = cues,
    category = category,
    response = NULL,
    group = group,
    group_unique = group_unique
  )
  test <- .prepare_staninput_frame(
    test,
    cues = cues,
    category = NULL,
    response = response,
    group = group,
    group_unique = group_unique
  )
  .assert_true(
    nrow(test) >= 1L,
    msg = "new_ideal_adaptor_stanfit_input requires non-empty test data."
  )

  # Handle grouping structure if a unique group column is provided.
  has_unique_groups <- !is.null(group_unique) &&
    group_unique %in% names(exposure)
  if (has_unique_groups) {
    .validate_unique_group_mapping(exposure, group, group_unique)
    .check_unique_group_identity(
      exposure = exposure,
      cues = cues,
      category = category,
      group = group,
      group_unique = group_unique
    )
  }

  if (!is.null(category) && !is.null(response)) {
    shared_category_labels <- union(
      levels(exposure[[category]]),
      levels(test[[response]])
    )
    exposure[[category]] <- factor(
      exposure[[category]],
      levels = shared_category_labels
    )
    test[[response]] <- factor(
      test[[response]],
      levels = shared_category_labels
    )
  }
  shared_group_labels <- union(
    levels(exposure[[group]]),
    levels(test[[group]])
  )
  exposure[[group]] <- factor(
    exposure[[group]],
    levels = shared_group_labels
  )
  test[[group]] <- factor(test[[group]], levels = shared_group_labels)

  if (has_unique_groups) {
    shared_unique_labels <- union(
      levels(factor(exposure[[group_unique]])),
      if (group_unique %in% names(test)) {
        levels(factor(test[[group_unique]]))
      } else {
        character(0)
      }
    )
    exposure[[group_unique]] <- factor(
      exposure[[group_unique]],
      levels = shared_unique_labels
    )
    if (group_unique %in% names(test)) {
      test[[group_unique]] <- factor(
        test[[group_unique]],
        levels = shared_unique_labels
      )
    }
    group_labels <- shared_unique_labels
  } else {
    group_labels <- shared_group_labels
  }

  category_labels <- if (!is.null(category) && !is.null(response)) {
    shared_category_labels
  } else if (!is.null(category) && category %in% names(exposure)) {
    levels(exposure[[category]])
  } else {
    character(0)
  }
  response_labels <- if (!is.null(test) && !is.null(response) &&
    response %in% names(test)) {
    levels(test[[response]])
  } else {
    category_labels
  }
  if (length(category_labels) == 0L) {
    category_labels <- response_labels
  }
  if (length(response_labels) == 0L) {
    response_labels <- category_labels
  }

  transform <- get_affine_transform(exposure, cues, transform_type)
  exposure_transformed <- transform@transform.function(
    exposure,
    return_type = "replace"
  )
  test_transformed <- transform@transform.function(
    test,
    return_type = "replace"
  )

  data <- .build_stanfit_input_data(
    exposure = exposure,
    test = test,
    category = category,
    response = response,
    group = group,
    group_unique = group_unique,
    cues = cues
  )
  data_transformed <- .build_stanfit_input_data(
    exposure = exposure_transformed,
    test = test_transformed,
    category = category,
    response = response,
    group = group,
    group_unique = group_unique,
    cues = cues
  )

  orig_var_names <- list(
    group = group,
    group_unique = if (!is.null(group_unique)) group_unique else group,
    category = category,
    response_category = response,
    cues = cues
  )

  suff_stats <- tryCatch(
    get_sufficient_category_statistics(
      exposure,
      cues = cues,
      category = category,
      group = group,
      simplify = FALSE
    ),
    error = function(e) NULL
  )

  IdealAdaptorStanfitInput(
    data = list(
      original = data,
      transformed = data_transformed
    ),
    transform_information = transform,
    metadata = list(
      label_information = list(
        cue = cues,
        category = category_labels,
        response_category = response_labels,
        group = group_labels
      ),
      original_variable_names = orig_var_names,
      exposure_category_statistics = suff_stats
    )
  )
}

#' Prepare a data frame for Stan input construction
#'
#' Internal helper that retains the required columns, coerces cue values to
#' numeric data, removes rows with missing values, and converts the required
#' columns to factors before downstream processing.
#'
#' @param data Input data frame.
#' @param cues Character vector of cue columns.
#' @param category Optional category column.
#' @param response Optional response column.
#' @param group Grouping column.
#' @param group_unique Optional grouping column used to collapse repeated
#'   exposure conditions.
#' @keywords internal
#' @noRd
.prepare_staninput_frame <- function(
  data, cues, category, response, group, group_unique = NULL
) {
  arg_name <- deparse(substitute(data))
  data <- as.data.frame(data)
  .assert_non_NA_character(cues)
  required_cols <- c(group)
  if (!is.null(category)) required_cols <- c(required_cols, category)
  if (!is.null(response)) required_cols <- c(required_cols, response)
  if (!is.null(group_unique)) required_cols <- c(required_cols, group_unique)
  required_cols <- unique(required_cols)
  missing_cols <- setdiff(required_cols, names(data))
  .assert_true(
    length(missing_cols) == 0L,
    msg = sprintf(
      "Missing columns in %s data: %s",
      arg_name, paste(missing_cols, collapse = ", ")
    )
  )
  missing_cues <- setdiff(cues, names(data))
  .assert_true(
    length(missing_cues) == 0L,
    msg = sprintf(
      "Missing cue columns in %s data: %s",
      arg_name, paste(missing_cues, collapse = ", ")
    )
  )

  keep_cols <- unique(c(group, cues, required_cols))
  data <- data[, keep_cols, drop = FALSE]

  factor_cols <- required_cols
  preserved_levels <- lapply(data[, factor_cols, drop = FALSE], function(x) {
    if (is.factor(x)) levels(x) else NULL
  })
  data[factor_cols] <- lapply(data[factor_cols], function(x) {
    values <- as.character(x)
    values[is.na(values)] <- NA
    values
  })
  data[, cues] <- lapply(data[, cues, drop = FALSE], function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  })
  complete_rows <- stats::complete.cases(
    data[, c(cues, required_cols), drop = FALSE]
  )
  data <- data[complete_rows, , drop = FALSE]

  for (col_name in factor_cols) {
    current_levels <- preserved_levels[[col_name]]
    if (!is.null(current_levels) && length(current_levels) > 0) {
      data[[col_name]] <- factor(data[[col_name]], levels = current_levels)
    } else {
      data[[col_name]] <- factor(data[[col_name]])
    }
  }

  data
}

#' Align factor levels across exposure and test data
#'
#' Internal helper that ensures the grouping and category factors are aligned
#' with the levels present in the reference data.
#'
#' @param data Input data frame to align.
#' @param reference Reference data frame providing the target factor levels.
#' @param category Optional category column.
#' @param response Optional response column.
#' @param group Grouping column.
#' @keywords internal
#' @noRd
.align_staninput_factors <- function(
  data, reference, category, response, group
) {
  if (!is.null(category) && !is.null(response)) {
    category_labels <- levels(reference[[category]])
    if (length(category_labels) == 0) {
      category_labels <- levels(reference[[response]])
    }
    if (length(category_labels) > 0) {
      data[[category]] <- factor(data[[category]], levels = category_labels)
      data[[response]] <- factor(data[[response]], levels = category_labels)
    }
  }
  group_labels <- levels(reference[[group]])
  if (length(group_labels) == 0) {
    group_labels <- levels(data[[group]])
  }
  if (length(group_labels) > 0) {
    data[[group]] <- factor(data[[group]], levels = group_labels)
  }
  data
}

#' Build the combined exposure and test data frame
#'
#' Internal helper that stacks the processed exposure and test rows into a
#' single data frame annotated with a Phase column.
#'
#' @keywords internal
#' @noRd
.build_stanfit_input_data <- function(
  exposure, test, category, response, group, group_unique, cues
) {
  if (nrow(exposure) > 0L) {
    exp_grp_unique <- if (!is.null(group_unique) &&
      group_unique %in% names(exposure)) {
      exposure[[group_unique]]
    } else {
      exposure[[group]]
    }
    exp_cat <- if (category %in% names(exposure)) {
      exposure[[category]]
    } else {
      rep(NA, nrow(exposure))
    }
    exp_resp <- if (response %in% names(exposure)) {
      exposure[[response]]
    } else {
      rep(NA, nrow(exposure))
    }
    exposure_data <- data.frame(
      group = exposure[[group]],
      group_unique = exp_grp_unique,
      category = exp_cat,
      response_category = exp_resp,
      stringsAsFactors = FALSE
    )
    for (c in cues) {
      exposure_data[[c]] <- exposure[[c]]
    }
    exposure_data$Phase <- "exposure"
  } else {
    exposure_data <- data.frame(
      group = character(0),
      group_unique = character(0),
      category = character(0),
      response_category = character(0),
      stringsAsFactors = FALSE
    )
    for (c in cues) {
      exposure_data[[c]] <- numeric(0)
    }
    exposure_data$Phase <- character(0)
  }

  if (nrow(test) > 0L) {
    test_grp_unique <- if (!is.null(group_unique) &&
      group_unique %in% names(test)) {
      test[[group_unique]]
    } else {
      test[[group]]
    }
    test_cat <- if (category %in% names(test)) {
      test[[category]]
    } else {
      rep(NA, nrow(test))
    }
    test_resp <- if (response %in% names(test)) {
      test[[response]]
    } else {
      rep(NA, nrow(test))
    }
    test_data <- data.frame(
      group = test[[group]],
      group_unique = test_grp_unique,
      category = test_cat,
      response_category = test_resp,
      stringsAsFactors = FALSE
    )
    for (c in cues) {
      test_data[[c]] <- test[[c]]
    }
    test_data$Phase <- "test"
  } else {
    test_data <- data.frame(
      group = character(0),
      group_unique = character(0),
      category = character(0),
      response_category = character(0),
      stringsAsFactors = FALSE
    )
    for (c in cues) {
      test_data[[c]] <- numeric(0)
    }
    test_data$Phase <- character(0)
  }

  data <- rbind(exposure_data, test_data)
  data$Phase <- factor(data$Phase, levels = c("exposure", "test"))
  data
}

#' Validate unique group mapping consistency
#'
#' @param exposure Exposure data frame.
#' @param group Grouping column name.
#' @param group_unique Unique grouping column name.
#' @keywords internal
#' @noRd
.validate_unique_group_mapping <- function(exposure, group, group_unique) {
  .assert_data_contains_cols(exposure, group_unique)
  group_to_unique <- tapply(
    as.character(exposure[[group_unique]]),
    exposure[[group]],
    function(u) unique(stats::na.omit(u))
  )
  for (g in names(group_to_unique)) {
    u <- group_to_unique[[g]]
    if (length(u) > 1) {
      .stop(
        sprintf(
          "Group '%s' is mapped to multiple unique groups: %s",
          g,
          paste(u, collapse = ", ")
        )
      )
    }
  }
}

#' Check exposure statistics across unique groups
#'
#' @param exposure Exposure data frame.
#' @param cues Character vector of cue names.
#' @param category Category column name.
#' @param group Grouping column name.
#' @param group_unique Unique grouping column name.
#' @param tol Tolerance for floating point comparisons.
#' @keywords internal
#' @noRd
.check_unique_group_identity <- function(
  exposure,
  cues,
  category,
  group,
  group_unique,
  tol = 1e-6
) {
  unique_labels <- levels(factor(exposure[[group_unique]]))
  category_labels <- levels(factor(exposure[[category]]))
  group_to_unique <- tapply(
    as.character(exposure[[group_unique]]),
    exposure[[group]],
    function(u) unique(stats::na.omit(u))
  )

  mismatched_uniques <- character(0)
  mismatch_details <- list()

  for (u in unique_labels) {
    groups_in_u <- names(group_to_unique)[
      vapply(
        group_to_unique,
        function(x) length(x) == 1L && identical(x, u),
        logical(1)
      )
    ]
    if (length(groups_in_u) > 1L) {
      extract_grp_stats <- function(grp_name) {
        sub_grp <- exposure[exposure[[group]] == grp_name, , drop = FALSE]
        s <- list()
        for (cat_name in category_labels) {
          sub_cat <- sub_grp[sub_grp[[category]] == cat_name, , drop = FALSE]
          n_obs <- nrow(sub_cat)
          if (n_obs == 0L) {
            m_obs <- rep(0, length(cues))
            sp_obs <- matrix(0, nrow = length(cues), ncol = length(cues))
          } else {
            vals_mat <- as.matrix(sub_cat[, cues, drop = FALSE])
            m_obs <- .colMeans(vals_mat)
            if (n_obs > 1L) {
              centered <- sweep(vals_mat, 2, m_obs, "-")
              sp_obs <- crossprod(centered)
            } else {
              sp_obs <- matrix(0, nrow = length(cues), ncol = length(cues))
            }
          }
          s[[cat_name]] <- list(n = n_obs, mean = m_obs, spread = sp_obs)
        }
        s
      }

      are_stats_equal <- function(s1, s2) {
        for (cat_name in names(s1)) {
          if (s1[[cat_name]]$n != s2[[cat_name]]$n) {
            return(FALSE)
          }
          if (s1[[cat_name]]$n > 0L) {
            if (max(abs(s1[[cat_name]]$mean - s2[[cat_name]]$mean)) > tol) {
              return(FALSE)
            }
            if (s1[[cat_name]]$n > 1L) {
              if (
                max(abs(s1[[cat_name]]$spread - s2[[cat_name]]$spread)) > tol
              ) {
                return(FALSE)
              }
            }
          }
        }
        TRUE
      }

      eq_sets <- list()
      for (k in seq_along(groups_in_u)) {
        g_name <- groups_in_u[k]
        g_stat <- extract_grp_stats(g_name)
        matched <- FALSE
        for (set_idx in seq_along(eq_sets)) {
          if (are_stats_equal(eq_sets[[set_idx]]$stats, g_stat)) {
            eq_sets[[set_idx]]$groups <- c(eq_sets[[set_idx]]$groups, g_name)
            matched <- TRUE
            break
          }
        }
        if (!matched) {
          eq_sets[[length(eq_sets) + 1L]] <- list(
            groups = g_name,
            stats = g_stat
          )
        }
      }

      if (length(eq_sets) > 1L) {
        mismatched_uniques <- c(mismatched_uniques, u)
        mismatch_details[[u]] <- eq_sets
      }
    }
  }

  if (length(mismatched_uniques) > 0L) {
    first_u <- mismatched_uniques[1]
    first_sets <- mismatch_details[[first_u]]

    set1_groups <- paste(
      sprintf('"%s"', first_sets[[1]]$groups),
      collapse = ", "
    )
    set2_groups <- paste(
      sprintf('"%s"', first_sets[[2]]$groups),
      collapse = ", "
    )

    format_stats <- function(s) {
      vapply(names(s), function(cat) {
        st <- s[[cat]]
        mean_str <- paste(round(st$mean, 3), collapse = ", ")
        sprintf("category '%s' (N=%d, mean=[%s])", cat, st$n, mean_str)
      }, character(1))
    }
    s1_str <- paste(format_stats(first_sets[[1]]$stats), collapse = "; ")
    s2_str <- paste(format_stats(first_sets[[2]]$stats), collapse = "; ")

    msg <- sprintf(
      paste0(
        "Non-identical exposure sufficient statistics found within ",
        "unique group(s): %s.\n",
        "Details for unique group '%s': contains %d distinct sets of ",
        "sufficient statistics across groups.\n",
        "Example difference:\n",
        "  - In group(s) %s sufficient statistics are: %s\n",
        "  - In group(s) %s sufficient statistics are: %s"
      ),
      paste(mismatched_uniques, collapse = ", "),
      first_u,
      length(first_sets),
      set1_groups,
      s1_str,
      set2_groups,
      s2_str
    )
    .stop(msg)
  }

  invisible(TRUE)
}
