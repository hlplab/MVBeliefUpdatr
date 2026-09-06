#' @include asserts.R
#' @include S7-core-classes.R
#' @include S7-generics.R
#' @include S7-transform-information.R
#' @include S7-staninput.R
#' @include S7-stanfit-input.R
#' @include S7-stanfit.R
NULL

get_ideal_adaptor_stanfit_constructor <- function(staninput = NULL) {
  if (!is.null(staninput)) {
    if (S7::S7_inherits(staninput, NIX_IdealAdaptorStaninput)) {
      NIX_IdealAdaptorStanfit
    } else if (S7::S7_inherits(staninput, MNIX_IdealAdaptorStaninput)) {
      MNIX_IdealAdaptorStanfit
    } else if (S7::S7_inherits(staninput, NIW_IdealAdaptorStaninput)) {
      NIW_IdealAdaptorStanfit
    } else {
      IdealAdaptorStanfit
    }
  } else {
    IdealAdaptorStanfit
  }
}

#' @rdname get_stanfit
#' @export
S7::method(get_stanfit, S7::class_any) <- function(x) {
  .stop("x must be an IdealAdaptorStanfit object.")
}

#' @rdname get_stanfit
#' @export
S7::method(get_stanfit, MVBU_Stanfit) <- function(x) {
  x@stanfit
}

#' @rdname get_stanfit
#' @export
S7::method(set_stanfit, list(S7::class_any, S7::class_any)) <- function(x, stanfit) {
  .stop("x must be an IdealAdaptorStanfit object.")
}

#' @rdname get_stanfit
#' @export
S7::method(set_stanfit, list(MVBU_Stanfit, S7::class_any)) <- function(x, stanfit) {
  # no assertions for stanfit here since the @<- assignment operator applied to S7 objects
  # will automatically call the validator for the class, which already check that the stanfit
  # is valid.
  x@stanfit <- stanfit
  x
}

#' @rdname get_parameters
#' @export
S7::method(get_parameter_names, MVBU_Stanfit) <- function(x, original_pars = FALSE, ...) {
  stanfit <- get_stanfit(x)
  if (is.null(stanfit)) {
    return(character(0))
  }
  if (original_pars) stanfit@model_pars else names(stanfit)
}

#' @rdname get_parameters
#' @export
S7::method(get_parameter_names, S7::class_any) <- function(x, original_pars = FALSE, ...) {
  if (inherits(x, "stanfit")) {
    if (original_pars) x@model_pars else names(x)
  } else {
    .stop("get_parameter_names is not implemented for objects of class ", class(x)[1])
  }
}


#' @rdname get_staninput
#' @export
S7::method(get_staninput, S7::class_any) <- function(x) {
  .stop("x must be an IdealAdaptorStanfit or IdealAdaptorStanfitInput object.")
}

#' @rdname get_staninput
#' @export
S7::method(get_staninput, MVBU_Stanfit) <- function(x) {
  x@staninput
}

#' @rdname get_transform_information
#' @export
S7::method(get_transform_information, S7::class_any) <- function(x) {
  .stop("x must be an IdealAdaptorStanfit or IdealAdaptorStanfitInput object.")
}

#' @rdname get_transform_information
#' @export
S7::method(get_transform_information, MVBU_Stanfit) <- function(x) {
  x@transform_information
}

#' @rdname get_cue_labels
#' @export
S7::method(get_cue_labels, MVBU_Stanfit) <- function(x, indices = NULL, ...) {
  label_info <- if (!is.null(x@metadata$label_information)) {
    x@metadata$label_information
  } else {
    list()
  }
  cues <- if (length(label_info$cue) > 0) {
    label_info$cue
  } else if (!is.null(attr(x@data, "cues"))) {
    as.character(attr(x@data, "cues"))
  } else {
    character(0)
  }
  if (!is.null(indices)) cues[indices] else cues
}

#' @rdname get_category_labels
#' @export
S7::method(get_category_labels, MVBU_Stanfit) <- function(x, indices = NULL, ...) {
  label_info <- if (!is.null(x@metadata$label_information)) {
    x@metadata$label_information
  } else {
    list()
  }
  cat_attr <- attr(x@data, "category")
  cats <- if (length(label_info$category) > 0) {
    label_info$category
  } else if (!is.null(cat_attr) && is.factor(x@data[[cat_attr]])) {
    levels(x@data[[cat_attr]])
  } else {
    character(0)
  }
  if (!is.null(indices)) cats[indices] else cats
}

#' @rdname get_group_labels
#' @export
S7::method(
  get_group_labels,
  MVBU_Stanfit
) <- function(x, indices = NULL, include_prior = FALSE, ...) {
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

#' @rdname get_labels
#' @export
S7::method(get_labels, MVBU_Stanfit) <- function(x, ...) {
  list(
    cue = get_cue_labels(x, ...),
    category = get_category_labels(x, ...),
    group = get_group_labels(x, ...)
  )
}

#' @rdname get_model_type
#' @export
S7::method(get_model_type, MVBU_Stanfit) <- function(x) {
  if (!is.null(x@stanfit) && length(x@stanfit@model_name) > 0) {
    x@stanfit@model_name
  } else {
    "MVBU_Stanfit type not available."
  }
}

#' @rdname get_data
#' @export
S7::method(get_data, MVBU_Stanfit) <- function(
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

#' @rdname get_exposure_data
#' @export
S7::method(get_exposure_data, MVBU_Stanfit) <- function(
  x,
  groups = get_group_labels(x, include_prior = FALSE),
  ...
) {
  get_data(x, groups = groups, ...) %>%
    dplyr::filter(.data$Phase == "exposure")
}

#' @rdname get_test_data
#' @export
S7::method(get_test_data, MVBU_Stanfit) <- function(
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

.compute_exposure_category_statistic <- function(
  x,
  categories = get_category_labels(x),
  groups = get_group_labels(x, include_prior = FALSE),
  statistic = c("n", "mean", "css", "uss", "cov"),
  untransform_cues = FALSE,
  ...
) {
  .assert_that(
    all(statistic %in% c("n", "mean", "css", "uss", "cov")),
    msg = "statistic must be one of 'n', 'mean', 'css', 'uss', or 'cov'."
  )
  .assert_that(
    any(is.factor(categories), is.character(categories), is.numeric(categories))
  )
  .assert_that(
    any(is.factor(groups), is.character(groups), is.numeric(groups))
  )
  avail_cats <- get_category_labels(x)
  .assert_that(
    all(categories %in% avail_cats),
    msg = paste(
      "Some categories not found in the exposure data:",
      paste(setdiff(categories, avail_cats), collapse = ", ")
    )
  )
  avail_grps <- get_group_labels(x, include_prior = FALSE)
  .assert_that(
    all(groups %in% avail_grps),
    msg = paste(
      "Some groups not found in the exposure data:",
      paste(setdiff(groups, avail_grps), collapse = ", ")
    )
  )

  staninput <- get_staninput(x)@values
  category_names <- get_category_labels(x)
  group_names <- get_group_labels(x, include_prior = FALSE)
  cue_names <- get_cue_labels(x)
  stanmodelname <- get_model_type(x)

  df <- NULL

  # Get counts n
  if (any(untransform_cues, c("n", "css", "cov") %in% statistic)) {
    n <- staninput$N_exposure
    d <- dim(n)
    if (!length(category_names)) {
      category_names <- paste0("category_", seq_len(d[1]))
    }
    if (!length(group_names)) {
      group_names <- paste0("group_", seq_len(d[2]))
    }
    dn <- list(category = category_names, group = group_names)

    df.n <- tibble::tibble()
    for (c in 1:d[1]) {
      for (g in 1:d[2]) {
        df.n <- dplyr::bind_rows(
          df.n,
          tibble::tibble(
            group = dn[[2]][g],
            category = dn[[1]][c],
            n = n[c, g]
          )
        )
      }
    }

    df <- if (!is.null(df)) {
      dplyr::left_join(df, df.n, by = c("group", "category"))
    } else {
      df.n
    }
  }

  # Get central tendencies m
  if (any(untransform_cues, c("mean", "css", "cov") %in% statistic)) {
    m <- staninput$x_mean_exposure
    d <- dim(m)
    if (!length(category_names)) {
      category_names <- paste0("category_", seq_len(d[1]))
    }
    if (!length(group_names)) {
      group_names <- paste0("group_", seq_len(d[2]))
    }
    if (grepl("^NIX", stanmodelname)) {
      dn <- list(category = category_names, group = group_names)
    } else {
      if (!length(cue_names)) {
        cue_names <- paste0("cues", seq_len(d[3]))
      }
      dn <- list(
        category = category_names,
        group = group_names,
        cue = cue_names
      )
    }

    df.m <- tibble::tibble()
    for (c in 1:d[1]) {
      for (g in 1:d[2]) {
        if (grepl("^NIX", stanmodelname)) {
          df.m <- dplyr::bind_rows(
            df.m,
            tibble::tibble(
              group = dn[[2]][g],
              category = dn[[1]][c],
              value = m[c, g]
            )
          )
        } else {
          for (f in 1:d[3]) {
            df.m <- dplyr::bind_rows(
              df.m,
              tibble::tibble(
                group = dn[[2]][g],
                category = dn[[1]][c],
                cue = dn[[3]][f],
                value = m[c, g, f]
              )
            )
          }
        }
      }
    }

    if (grepl("^NIX", stanmodelname)) {
      df.m <- dplyr::mutate(df.m, mean = .data$value) %>%
        dplyr::select(-.data$value)
    } else {
      df.m <- df.m %>%
        tidyr::pivot_wider(names_from = "cue", values_from = "value") %>%
        make_vector_column(
          cols = dn[[3]],
          vector_col = "mean",
          .keep = "unused"
        )
    }

    df <- if (!is.null(df)) {
      dplyr::left_join(df, df.m, by = c("group", "category"))
    } else {
      df.m
    }
  }

  # Get scatter or covariance matrices s
  if (any(untransform_cues, c("uss", "css", "cov") %in% statistic)) {
    # TO DO: Stan programs should be updated to return either css or uss
    # for all types of models (or all stats), rather than storing different stats for each model.
    # The stancode could then transform the input data to the correct quantities. This would make
    # the handling here a lot easier.
    if (grepl("^NIX", stanmodelname)) {
      .stop(
        "Extraction of uss, css, or cov not yet implemented for NIX models."
      )
    } else if (grepl("^NIW", stanmodelname)) {
      s <- staninput$x_ss_exposure
    } else if (grepl("^MNIX", stanmodelname)) {
      .stop(
        "Extraction of uss, css, or cov not yet implemented for MNIX models."
      )
    } else {
      .stop(
        "Unrecognized model. No method available to extract category variance."
      )
    }

    d <- dim(s)
    if (!length(category_names)) {
      category_names <- paste0("category_", seq_len(d[1]))
    }
    if (!length(group_names)) {
      group_names <- paste0("group_", seq_len(d[2]))
    }
    if (!length(cue_names)) {
      cue_names <- paste0("cues", seq_len(d[3]))
    }
    dn <- list(
      category = category_names,
      group = group_names,
      cue = cue_names,
      cue2 = cue_names
    )

    df.s <- tibble::tibble()
    for (c in 1:d[1]) {
      for (g in 1:d[2]) {
        for (f1 in 1:d[3]) {
          for (f2 in 1:d[4]) {
            df.s <- dplyr::bind_rows(
              df.s,
              tibble::tibble(
                group = dn[[2]][g],
                category = dn[[1]][c],
                cue = dn[[3]][f1],
                cue2 = dn[[4]][f2],
                value = s[c, g, f1, f2]
              )
            )
          }
        }
      }
    }

    df.s <- df.s %>%
      dplyr::group_by(.data$category, .data$group) %>%
      dplyr::summarise(
        uss = list(matrix(.data$value, nrow = sqrt(length(.data$value)))),
        .groups = "drop"
      )

    df <- if (!is.null(df)) {
      dplyr::left_join(df, df.s, by = c("group", "category"))
    } else {
      df.s
    }
  }

  if (any(untransform_cues, c("css", "cov") %in% statistic)) {
    df <- dplyr::mutate(
      df,
      css = purrr::pmap(list(.data$uss, .data$n, .data$mean), uss2css)
    )
  }

  if (any(untransform_cues, c("cov") %in% statistic)) {
    df <- dplyr::mutate(
      df,
      cov = purrr::map2(.data$css, .data$n, css2cov)
    )
  }

  if (untransform_cues) {
    trans_info <- get_transform_information(x)
    if ("cov" %in% statistic) {
      df <- dplyr::mutate(
        df,
        cov = purrr::map(.data$cov, ~ untransform_category_cov(.x, trans_info))
      )
    }
    if (any(c("css", "uss") %in% statistic)) {
      df <- dplyr::mutate(
        df,
        css = purrr::map2(.data$cov, .data$n, cov2css)
      )
    }
    if (any(c("uss", "mean") %in% statistic)) {
      df <- dplyr::mutate(
        df,
        mean = purrr::map(
          .data$mean,
          ~ untransform_category_mean(.x, trans_info)
        )
      )
    }
    if ("uss" %in% statistic) {
      df <- dplyr::mutate(
        df,
        uss = purrr::pmap(list(.data$cov, .data$n, .data$mean), css2uss)
      )
    }
  }

  df <- df %>%
    dplyr::select(dplyr::all_of(c("group", "category", statistic))) %>%
    dplyr::filter(
      .data$group %in% groups,
      .data$category %in% categories
    ) %>%
    dplyr::mutate(
      category = factor(.data$category, levels = categories),
      group = factor(.data$group, levels = groups)
    )

  if (nrow(df) == 1 && length(statistic) == 1) {
    df <- df[[statistic]][[1]]
  }

  df
}

#' @rdname get_exposure_category_statistic
#' @export
S7::method(get_exposure_category_statistic, MVBU_Stanfit) <- function(
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

#' @rdname get_exposure_category_statistic
#' @export
S7::method(get_exposure_category_mean, MVBU_Stanfit) <- function(x, ...) {
  get_exposure_category_statistic(x, ..., statistic = "mean")
}

#' @rdname get_exposure_category_statistic
#' @export
S7::method(get_exposure_category_css, MVBU_Stanfit) <- function(x, ...) {
  get_exposure_category_statistic(x, ..., statistic = "css")
}

#' @rdname get_exposure_category_statistic
#' @export
S7::method(get_exposure_category_uss, MVBU_Stanfit) <- function(x, ...) {
  get_exposure_category_statistic(x, ..., statistic = "uss")
}

#' @rdname get_exposure_category_statistic
#' @export
S7::method(get_exposure_category_cov, MVBU_Stanfit) <- function(x, ...) {
  get_exposure_category_statistic(x, ..., statistic = "cov")
}

.nest_draw_cues <- function(d.pars) {
  if (!all(c("cue", "cue2") %in% names(d.pars))) {
    return(d.pars)
  }

  group_cols <- setdiff(names(d.pars), c("cue", "cue2", "m", "S"))
  cues_order <- unique(d.pars$cue)
  d_dim <- length(cues_order)

  d.pars %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::arrange(.data$cue, .data$cue2, .by_group = TRUE) %>%
    dplyr::summarise(
      m = list({
        m_vec <- unique(.data$m)
        names(m_vec) <- unique(.data$cue)
        m_vec
      }),
      S = list({
        matrix(.data$S, nrow = d_dim, ncol = d_dim, dimnames = list(cues_order, cues_order))
      }),
      .groups = "drop"
    ) %>%
    dplyr::relocate(dplyr::starts_with(c("lapse_", "prior")), .after = "S")
}

.unnest_draw_cues <- function(d.pars) {
  if (all(c("cue", "cue2") %in% names(d.pars))) {
    return(d.pars)
  }

  cue_labels <- if ("m" %in% names(d.pars) && length(d.pars$m) > 0 && !is.null(names(d.pars$m[[1]]))) {
    names(d.pars$m[[1]])
  } else if ("S" %in% names(d.pars) && length(d.pars$S) > 0 && !is.null(colnames(d.pars$S[[1]]))) {
    colnames(d.pars$S[[1]])
  } else {
    NULL
  }

  d.pars <- d.pars %>%
    tidyr::unnest(c("m", "S")) %>%
    dplyr::group_by(dplyr::across(-dplyr::any_of(c("m", "S")))) %>%
    dplyr::mutate(cue = cue_labels)

  for (i in seq_along(cue_labels)) {
    d.pars <- d.pars %>% dplyr::mutate(!!rlang::sym(cue_labels[i]) := (.data$S)[, i])
  }

  d.pars %>%
    dplyr::select(-"S") %>%
    tidyr::pivot_longer(cols = dplyr::all_of(cue_labels), values_to = "S", names_to = "cue2") %>%
    dplyr::ungroup() %>%
    dplyr::relocate(.data$cue, .data$cue2, .after = dplyr::any_of("nu"))
}

#' @rdname get_draws
#' @export
S7::method(get_draws, MVBU_Stanfit) <- function(
  fit,
  categories = get_category_labels(fit),
  groups = get_group_labels(fit, include_prior = TRUE),
  which = if ("prior" %in% groups) {
    if (length(groups) > 1) "both" else "prior"
  } else {
    "posterior"
  },
  ndraws = NULL,
  untransform_cues = FALSE,
  summarize = FALSE,
  nest = TRUE,
  seed = if (!is.null(ndraws)) runif(1, -1e6, 1e6) else NULL,
  ...
) {
  dots <- list(...)
  if ("wide" %in% names(dots)) {
    lifecycle::deprecate_warn("0.0.9", "get_draws(wide = )")
  }
  .assert_contains_draws(fit)
  .assert_that(
    any(is.factor(categories), is.character(categories), is.numeric(categories))
  )
  .assert_that(
    any(is.factor(groups), is.character(groups), is.numeric(groups))
  )
  avail_cats <- get_category_labels(fit)
  .assert_that(
    all(categories %in% avail_cats),
    msg = paste(
      "Some categories not found in model:",
      paste(setdiff(categories, avail_cats), collapse = ", ")
    )
  )
  avail_grps <- get_group_labels(fit, include_prior = TRUE)
  .assert_that(
    all(groups %in% avail_grps),
    msg = paste(
      "Some groups not found in model:",
      paste(setdiff(groups, avail_grps), collapse = ", ")
    )
  )
  .assert_that(
    which %in% c("prior", "posterior", "both"),
    msg = "which must be one of 'prior', 'posterior', or 'both'."
  )
  .assert_that(
    any(is.null(ndraws), .is_scalar_count(ndraws)),
    msg = "If not NULL, ndraw must be a count."
  )
  .assert_that(
    any(is.null(ndraws), !is.null(seed)),
    msg = "If ndraws is not NULL, seed must be specified."
  )
  .assert_that(.is_scalar_logical(summarize))

  if ("prior" %in% groups && length(groups) > 1) {
    d.prior <- get_draws(
      fit = fit,
      categories = categories,
      groups = "prior",
      ndraws = ndraws,
      untransform_cues = untransform_cues,
      summarize = summarize,
      nest = nest,
      seed = seed,
      ...
    )
    d.posterior <- get_draws(
      fit = fit,
      categories = categories,
      groups = setdiff(groups, "prior"),
      ndraws = ndraws,
      untransform_cues = untransform_cues,
      summarize = summarize,
      nest = nest,
      seed = seed,
      ...
    )
    d.pars <- rbind(d.prior, d.posterior) %>%
      dplyr::mutate(
        group = factor(
          .data$group,
          levels = c(levels(d.prior$group), levels(d.posterior$group))
        )
      )
    return(d.pars)
  }

  postfix <- if ("prior" %in% groups) "_0" else "_n"
  kappa <- paste0("kappa", postfix)
  nu <- paste0("nu", postfix)
  m <- paste0("m", postfix)
  S <- paste0("S", postfix)

  pars.index <- if ("prior" %in% groups) "category" else c("category", "group")

  stanfit <- get_stanfit(fit)
  var_names <- names(stanfit)
  is_nix_1d <- "m_0[1]" %in% var_names || "m_n[1,1]" %in% var_names

  if (is_nix_1d) {
    cue_lab <- get_cue_labels(fit)
    if (length(cue_lab) == 0) cue_lab <- "cue"
    cue_lab <- cue_lab[1]

    if ("prior" %in% groups) {
      d.pars <- stanfit %>%
        tidybayes::spread_draws(
          !!rlang::sym(kappa),
          !!rlang::sym(nu),
          (!!rlang::sym(m))[!!!rlang::syms(pars.index)],
          (!!rlang::sym(S))[!!!rlang::syms(pars.index)],
          lapse_rate,
          ndraws = ndraws,
          seed = seed
        ) %>%
        dplyr::mutate(cue = cue_lab, cue2 = cue_lab)
    } else {
      d.pars <- stanfit %>%
        tidybayes::spread_draws(
          (!!rlang::sym(kappa))[!!!rlang::syms(pars.index)],
          (!!rlang::sym(nu))[!!!rlang::syms(pars.index)],
          (!!rlang::sym(m))[!!!rlang::syms(pars.index)],
          (!!rlang::sym(S))[!!!rlang::syms(pars.index)],
          lapse_rate,
          ndraws = ndraws,
          seed = seed
        ) %>%
        dplyr::mutate(cue = cue_lab, cue2 = cue_lab)
    }
  } else {
    if ("prior" %in% groups) {
      d.pars <- stanfit %>%
        tidybayes::spread_draws(
          !!rlang::sym(kappa),
          !!rlang::sym(nu),
          (!!rlang::sym(m))[!!!rlang::syms(pars.index), cue],
          (!!rlang::sym(S))[!!!rlang::syms(pars.index), cue, cue2],
          lapse_rate,
          ndraws = ndraws,
          seed = seed
        )
    } else {
      d.pars <- stanfit %>%
        tidybayes::spread_draws(
          (!!rlang::sym(kappa))[!!!rlang::syms(pars.index)],
          (!!rlang::sym(nu))[!!!rlang::syms(pars.index)],
          (!!rlang::sym(m))[!!!rlang::syms(pars.index), cue],
          (!!rlang::sym(S))[!!!rlang::syms(pars.index), cue, cue2],
          lapse_rate,
          ndraws = ndraws,
          seed = seed
        )
    }
  }

  d.pars <- d.pars %>%
    dplyr::rename_with(~ sub(postfix, "", .), dplyr::ends_with(postfix))

  if (summarize) {
    d.pars <- d.pars %>%
      dplyr::group_by(!!!rlang::syms(pars.index), .data$cue, .data$cue2) %>%
      dplyr::summarise(
        dplyr::across(
          c("kappa", "nu", "m", "S", "lapse_rate"),
          mean
        ),
        .groups = "drop"
      ) %>%
      dplyr::mutate(.chain = "all", .iteration = "all", .draw = "all")
  }

  if ("prior" %in% groups) {
    d.pars <- dplyr::mutate(d.pars, group = "prior")
  }

  d.pars <- d.pars %>%
    dplyr::filter(
      .data$group %in% groups,
      .data$category %in% categories
    ) %>%
    dplyr::relocate(
      dplyr::any_of(c(
        ".chain",
        ".iteration",
        ".draw",
        "group",
        "category",
        "kappa",
        "nu",
        "cue",
        "cue2",
        "m",
        "S",
        "lapse_rate"
      ))
    )

  check_ok <- d.pars %>%
    dplyr::ungroup() %>%
    dplyr::summarise(
      dplyr::across(
        c("kappa", "nu", "m", "S", "lapse_rate"),
        ~ all(!is.na(.x) & !is.nan(.x) & !is.infinite(.x))
      )
    ) %>%
    unlist() %>%
    all()

  if (!check_ok) {
    .stop(
      paste0(
        "Some draws are not well-formed (NA, NaN, infinite values). ",
        "This likely means that there was an issue during fitting."
      )
    )
  }

  if (untransform_cues) {
    transform_info <- get_transform_information(fit)
    if (!is.null(transform_info)) {
      d.pars <- .nest_draw_cues(d.pars)
      d.pars$m <- lapply(d.pars$m, untransform_category_mean, transform = transform_info)
      d.pars$S <- lapply(d.pars$S, untransform_category_cov, transform = transform_info)
      if (!nest) {
        d.pars <- .unnest_draw_cues(d.pars)
      }
    }
  } else if (nest) {
    d.pars <- .nest_draw_cues(d.pars)
  }

  d.pars <- d.pars %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      category = factor(.data$category, levels = categories),
      group = factor(.data$group, levels = groups)
    )

  d.pars
}

#' @rdname get_parameters
#' @export
S7::method(get_expected_category_statistic, MVBU_Stanfit) <- function(
  x,
  categories = get_category_labels(x),
  groups = get_group_labels(x, include_prior = TRUE),
  statistic = c("mu", "Sigma"),
  ...
) {
  .assert_that(all(statistic %in% c("mu", "Sigma")))
  .assert_that(
    any(is.factor(categories), is.character(categories), is.numeric(categories))
  )
  .assert_that(
    any(is.factor(groups), is.character(groups), is.numeric(groups))
  )
  avail_cats <- get_category_labels(x)
  .assert_that(
    all(categories %in% avail_cats),
    msg = paste(
      "Some categories not found in model:",
      paste(setdiff(categories, avail_cats), collapse = ", ")
    )
  )
  avail_grps <- get_group_labels(x, include_prior = TRUE)
  .assert_that(
    all(groups %in% avail_grps),
    msg = paste(
      "Some groups not found in model:",
      paste(setdiff(groups, avail_grps), collapse = ", ")
    )
  )

  draws <- get_draws(
    x,
    categories = categories,
    groups = groups,
    nest = TRUE,
    summarize = FALSE,
    ...
  ) %>%
    dplyr::mutate(Sigma = get_expected_Sigma_from_S(.data$S, .data$nu)) %>%
    dplyr::group_by(.data$group, .data$category) %>%
    dplyr::summarise(
      mu.mean = list(purrr::reduce(.data$m, `+`) / length(.data$m)),
      Sigma.mean = list(purrr::reduce(.data$Sigma, `+`) / length(.data$Sigma)),
      .groups = "drop"
    ) %>%
    dplyr::select(
      dplyr::all_of(c("group", "category", paste0(statistic, ".mean")))
    ) %>%
    dplyr::mutate(
      category = factor(.data$category, levels = categories),
      group = factor(.data$group, levels = groups)
    )

  if (nrow(draws) == 1 && length(statistic) == 1) {
    draws <- draws[[paste0(statistic, ".mean")]][[1]]
  }

  draws
}

#' @rdname get_parameters
#' @export
S7::method(get_expected_mu, MVBU_Stanfit) <- function(x, ...) {
  get_expected_category_statistic(x, statistic = "mu", ...)
}

#' @rdname get_parameters
#' @export
S7::method(get_expected_sigma, MVBU_Stanfit) <- function(x, ...) {
  get_expected_category_statistic(x, statistic = "Sigma", ...)
}

# -----------------------------------------------------------------------------
# Transform and Draw Accessors
# -----------------------------------------------------------------------------

#' @rdname get_transform_function
#' @export
S7::method(get_transform_function, S7::class_any) <- function(x, ...) {
  info <- get_transform_information(x)
  if (S7::S7_inherits(info, MVBU_TransformInformation)) {
    return(info@`transform.function`)
  }
  identity
}

#' @rdname get_transform_function
#' @export
S7::method(get_untransform_function, S7::class_any) <- function(x, ...) {
  info <- get_transform_information(x)
  if (S7::S7_inherits(info, MVBU_TransformInformation)) {
    return(info@`untransform.function`)
  }
  identity
}

#' @rdname get_number_of_draws
#' @export
S7::method(get_number_of_draws, S7::class_any) <- function(fit) {
  stanfit <- get_stanfit(fit)
  if (is.null(stanfit)) {
    return(0L)
  }
  if (isS4(stanfit) && .hasSlot(stanfit, "sim") && length(stanfit@sim$samples) > 0) {
    return(length(stanfit@sim$samples[[1]][[1]]))
  }
  posterior::ndraws(posterior::as_draws(stanfit))
}

#' @rdname get_number_of_draws
#' @export
S7::method(get_random_draw_indices, S7::class_any) <- function(fit, ndraws) {
  n.all.draws <- get_number_of_draws(fit)
  .assert_that(
    ndraws <= n.all.draws,
    msg = paste0(
      "Cannot return ", ndraws, " draws because there are only ",
      n.all.draws, " in the object."
    )
  )
  sample(seq_len(n.all.draws), size = ndraws)
}


#' Get categorization function from a model object
#'
#' Returns a categorization function that can be used to categorize a set of
#' cues into categories.
#'
#' @param x A model object.
#' @param lapse_treatment Should the consequences of attentional lapses be
#'   included in the categorization function ("marginalize") or not
#'   ("no_lapses")? (default: "marginalize")
#' @param groups Character vector of groups to include.
#' @param ... Optionally, additional arguments handed to \code{\link{get_draws}}.
#'
#' @return A tibble with categorization functions per group and draw.
#' @export
get_categorization_function <- function(
  x,
  lapse_treatment = c("no_lapses", "sample", "marginalize")[3],
  groups = get_group_labels(x, include_prior = TRUE),
  ...
) {
  d.pars <- get_draws(
    x,
    groups = groups,
    summarize = FALSE,
    ...
  )

  d.pars %>%
    dplyr::group_by(.data$group, .data$.draw) %>%
    dplyr::group_modify(~ {
      tibble::tibble(
        f = list(
          get_categorization_function_from_stanfit_draws(
            .x,
            noise_treatment = "no_noise",
            lapse_treatment = lapse_treatment
          )
        )
      )
    })
}

#' @keywords internal
#' @noRd
get_categorization_function_from_stanfit_draws <- function(x, ...) {
  get_NIW_categorization_function(
    ms = x$m,
    Ss = x$S,
    kappas = x$kappa,
    nus = x$nu,
    lapse_rate = unlist(x$lapse_rate)[1],
    ...
  )
}

S7::method(print, MVBU_Stanfit) <- function(x, ...) {
  cls <- class(x)[1]
  cat("<", cls, ">\n", sep = "")
  cat("  Model type: ", get_model_type(x), "\n", sep = "")
  cats <- get_category_labels(x)
  cues <- get_cue_labels(x)
  grps <- get_group_labels(x, include_prior = FALSE)
  cat("  Categories (", length(cats), "): ", paste(cats, collapse = ", "), "\n", sep = "")
  cat("  Cues (", length(cues), "): ", paste(cues, collapse = ", "), "\n", sep = "")
  if (length(grps) > 0) cat("  Groups (", length(grps), "): ", paste(grps, collapse = ", "), "\n", sep = "")
  n_draws <- get_number_of_draws(x)
  cat("  Draws: ", n_draws, "\n", sep = "")
  invisible(x)
}

#' An S7 class for MVBU_Stanfit summaries
#'
#' @keywords internal
Summary_MVBU_Stanfit <- S7::new_class(
  "Summary_MVBU_Stanfit",
  parent = MVBU_Object,
  properties = list(
    fitted = S7::class_any,
    fixed = S7::class_any
  )
)

S7::method(print, Summary_MVBU_Stanfit) <- function(x, ...) {
  if (!is.null(x@fixed) && nrow(x@fixed) > 0) {
    cat("Fixed parameters:\n")
    print(as.data.frame(x@fixed), row.names = FALSE, max = nrow(x@fixed) * 100, ...)
    if (!is.null(x@fitted) && nrow(x@fitted) > 0) {
      cat("\n")
    }
  }

  if (!is.null(x@fitted) && nrow(x@fitted) > 0) {
    cat("Fitted parameters:\n")
    print(as.data.frame(x@fitted), row.names = FALSE, max = nrow(x@fitted) * 100, ...)
  } else {
    cat("No fitted parameters.\n")
  }
  invisible(x)
}

S7::method(as.data.frame, Summary_MVBU_Stanfit) <- function(x, ...) {
  as.data.frame(x@fitted, ...)
}

#' @export
head.Summary_MVBU_Stanfit <- function(x, n = 6L, ...) {
  utils::head(x@fitted, n = n, ...)
}

.extract_fixed_parameters <- function(object) {
  staninput <- get_staninput(object)
  stanvals <- if (!is.null(staninput)) staninput@values else list()
  category_levels <- get_category_labels(object)
  cue_levels <- get_cue_labels(object)

  rows <- list()

  # 1. lapse_rate
  if (isTRUE(stanvals$lapse_rate_known == 1) || isTRUE(stanvals$lapse_rate_known == 1L)) {
    rows[[length(rows) + 1L]] <- tibble::tibble(
      Parameter = "lapse_rate",
      `Dist.` = "",
      Group = "",
      Category = "",
      Cue1 = "",
      Cue2 = "",
      Value = as.numeric(stanvals$lapse_rate_data)
    )
  }

  # 2. mu_0 (indirectly fixing m_0)
  if (isTRUE(stanvals$mu_0_known == 1) || isTRUE(stanvals$mu_0_known == 1L)) {
    mu_mat <- as.matrix(stanvals$mu_0_data)
    K <- nrow(mu_mat)
    M <- ncol(mu_mat)
    for (k in seq_len(K)) {
      for (m in seq_len(M)) {
        rows[[length(rows) + 1L]] <- tibble::tibble(
          Parameter = "m",
          `Dist.` = "prior",
          Group = "",
          Category = if (k <= length(category_levels)) category_levels[k] else as.character(k),
          Cue1 = if (m <= length(cue_levels)) cue_levels[m] else as.character(m),
          Cue2 = "",
          Value = mu_mat[k, m]
        )
      }
    }
  }

  # 3. Sigma_0 (indirectly fixing S_0)
  if (isTRUE(stanvals$Sigma_0_known == 1) || isTRUE(stanvals$Sigma_0_known == 1L)) {
    sig_data <- stanvals$Sigma_0_data
    if (is.array(sig_data) && length(dim(sig_data)) == 3) {
      for (k in seq_len(dim(sig_data)[1])) {
        for (m1 in seq_len(dim(sig_data)[2])) {
          for (m2 in seq_len(dim(sig_data)[3])) {
            rows[[length(rows) + 1L]] <- tibble::tibble(
              Parameter = "S",
              `Dist.` = "prior",
              Group = "",
              Category = if (k <= length(category_levels)) category_levels[k] else as.character(k),
              Cue1 = if (m1 <= length(cue_levels)) cue_levels[m1] else as.character(m1),
              Cue2 = if (m2 <= length(cue_levels)) cue_levels[m2] else as.character(m2),
              Value = sig_data[k, m1, m2]
            )
          }
        }
      }
    } else if (is.matrix(sig_data)) {
      for (k in seq_len(nrow(sig_data))) {
        for (m in seq_len(ncol(sig_data))) {
          rows[[length(rows) + 1L]] <- tibble::tibble(
            Parameter = "S",
            `Dist.` = "prior",
            Group = "",
            Category = if (k <= length(category_levels)) category_levels[k] else as.character(k),
            Cue1 = if (m <= length(cue_levels)) cue_levels[m] else as.character(m),
            Cue2 = if (m <= length(cue_levels)) cue_levels[m] else as.character(m),
            Value = sig_data[k, m]
          )
        }
      }
    }
  }

  if (length(rows) > 0) {
    dplyr::bind_rows(rows)
  } else {
    tibble::tibble(
      Parameter = character(0),
      `Dist.` = character(0),
      Group = character(0),
      Category = character(0),
      Cue1 = character(0),
      Cue2 = character(0),
      Value = numeric(0)
    )
  }
}

#' Summarize an MVBU Stanfit object
#'
#' Specifies reasonable defaults for the parameters to be summarized for the stanfit object.
#'
#' @param object An \code{\link{MVBU_Stanfit}} object.
#' @param pars A character vector of parameter names to be summarized. If `NULL`, all model
#'   parameters are summarized. (default: `NULL`)
#' @param sufficient_pars_only Should only the sufficient parameters be summarized? (default: `TRUE`)
#' @param indices_as_names Should the indices of the parameters be translated into level names in the summary?
#'   This output is substantially more readable, with names shown in separate columns (distribution, group,
#'   category, cue). (default: `TRUE`)
#' @param include_transformed_pars Should transformed parameters be included in the summary? (default: `FALSE`)
#' @param ... Additional arguments passed to \code{rstan::summary}.
#'
#' @return An object of class \code{Summary_MVBU_Stanfit} containing parameter summaries and diagnostic statistics.
S7::method(summary, MVBU_Stanfit) <- function(
  object,
  pars = NULL,
  sufficient_pars_only = TRUE,
  indices_as_names = TRUE,
  include_transformed_pars = FALSE,
  ...
) {
  stanfit <- get_stanfit(object)
  if (is.null(stanfit)) {
    message("No Stanfit object found in MVBU_Stanfit object. Printing object instead.")
    print(object, ...)
    return(invisible(object))
  }
  .assert_contains_draws(stanfit)
  if (is.null(pars)) {
    pars <- names(stanfit)
    pars <- grep("^((kappa|nu|m|S|cue_weight)_|lapse_rate|p_category|Sigma_noise)", pars, value = TRUE)
    pars <- grep("^m_0_(tau|L_omega|cov)", pars, value = TRUE, invert = TRUE)
    pars <- grep("^((m|S)_0|lapse_rate)_param", pars, value = TRUE, invert = TRUE)
    if (!include_transformed_pars) pars <- grep("_transformed", pars, value = TRUE, invert = TRUE)
  }

  raw_sum <- rstan::summary(stanfit, pars = pars, ...)$summary
  if (is.null(raw_sum) || nrow(raw_sum) == 0) {
    return(Summary_MVBU_Stanfit(fitted = as.data.frame(raw_sum), fixed = .extract_fixed_parameters(object)))
  }

  # Sort and filter output
  full_summary <-
    raw_sum %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    dplyr::mutate(
      name = factor(
        gsub("^(kappa|nu|m|S|lapse_rate|p_category|cue_weight|Sigma_noise).*$", "\\1", .data$Parameter),
        levels = c("kappa", "nu", "m", "S", "cue_weight", "lapse_rate", "p_category", "Sigma_noise")
      ),
      distribution = gsub("^.*_(0|n).*$", "\\1", .data$Parameter),
      index = gsub("^.*_(0|n)?\\[(.*)\\]$", "\\2", .data$Parameter),
      index = ifelse(.data$index == .data$Parameter, 1, .data$index)
    ) %>%
    {
      if (sufficient_pars_only) dplyr::filter(., .data$distribution == "0" | .data$name %in% c("cue_weight", "lapse_rate", "p_category", "Sigma_noise")) else .
    } %>%
    tidyr::separate(.data$index, into = c("i1", "i2", "i3", "i4"), sep = ",", fill = "right") %>%
    dplyr::mutate(dplyr::across(c("i1", "i2", "i3", "i4"), as.integer)) %>%
    dplyr::arrange(.data$distribution, .data$name, .data$i1, .data$i2, .data$i3, .data$i4)

  category_levels <- get_category_labels(object)
  group_levels <- get_group_labels(object, include_prior = FALSE)
  cue_levels <- get_cue_labels(object)
  if (indices_as_names) {
    full_summary <-
      full_summary %>%
      dplyr::mutate(
        Parameter = as.character(.data$name),
        `Dist.` = dplyr::case_when(
          .data$distribution == "0" ~ "prior",
          .data$distribution == "n" ~ "posterior",
          TRUE ~ ""
        ),
        Group = dplyr::case_when(
          .data$name %in% c("kappa", "nu", "m", "S") & .data$`Dist.` == "posterior" ~ group_levels[.data$i2],
          .data$name == "cue_weight" ~ group_levels[.data$i1],
          TRUE ~ ""
        ),
        Category = dplyr::case_when(
          .data$name %in% c("m", "S") & .data$`Dist.` == "prior" ~ category_levels[.data$i1],
          .data$name %in% c("kappa", "nu", "m", "S") & .data$`Dist.` == "posterior" ~ category_levels[.data$i1],
          .data$name == "p_category" ~ category_levels[.data$i1],
          TRUE ~ ""
        ),
        Cue1 = dplyr::case_when(
          .data$name %in% c("m", "S") & .data$`Dist.` == "prior" ~ cue_levels[.data$i2],
          .data$name %in% c("m", "S") & .data$`Dist.` == "posterior" ~ cue_levels[.data$i3],
          .data$name == "cue_weight" ~ cue_levels[.data$i2],
          TRUE ~ ""
        ),
        Cue2 = dplyr::case_when(
          .data$name %in% c("S") & .data$`Dist.` == "prior" ~ cue_levels[.data$i3],
          .data$name %in% c("S") & .data$`Dist.` == "posterior" ~ cue_levels[.data$i4],
          TRUE ~ ""
        )
      ) %>%
      dplyr::relocate(
        tidyselect::all_of(c("Parameter", "Dist.", "Group", "Category", "Cue1", "Cue2")),
        tidyselect::everything()
      )
  } else {
    full_summary <- full_summary %>%
      dplyr::relocate(tidyselect::all_of("Parameter"), tidyselect::everything())
  }

  full_summary <-
    full_summary %>%
    dplyr::select(-dplyr::any_of(c("name", "distribution", "i1", "i2", "i3", "i4")))

  Rhats <- full_summary[["Rhat"]]
  if (!is.null(Rhats) && any(Rhats > 1.05, na.rm = TRUE)) {
    .warning(
      "Parts of the model have not converged (some Rhats are > 1.05). ",
      "Be careful when analysing the results! We recommend running ",
      "more iterations and/or setting stronger priors."
    )
  }
  div_trans <- tryCatch(sum(nuts_params(object, pars = "divergent__")$Value), error = function(e) 0)
  adapt_delta <- tryCatch(control_params(object)$adapt_delta, error = function(e) NULL)
  if (div_trans > 0) {
    .warning(
      "There were ", div_trans, " divergent transitions after warmup. ",
      if (!is.null(adapt_delta)) paste0("Increasing adapt_delta above ", adapt_delta, " may help. ") else "",
      "See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup"
    )
  }

  fixed_summary <- .extract_fixed_parameters(object)

  if (nrow(fixed_summary) > 0 && nrow(full_summary) > 0) {
    drop_idx <- integer(0)
    for (i in seq_len(nrow(fixed_summary))) {
      fp <- fixed_summary[i, ]
      matches <- which(
        full_summary$Parameter == fp$Parameter &
          (fp$`Dist.` == "" | full_summary$`Dist.` == fp$`Dist.`) &
          (fp$Category == "" | full_summary$Category == fp$Category) &
          (fp$Cue1 == "" | full_summary$Cue1 == fp$Cue1) &
          (fp$Cue2 == "" | full_summary$Cue2 == fp$Cue2)
      )
      drop_idx <- c(drop_idx, matches)
    }
    if (length(drop_idx) > 0) {
      full_summary <- full_summary[-unique(drop_idx), , drop = FALSE]
    }
  }

  return(Summary_MVBU_Stanfit(fitted = full_summary, fixed = fixed_summary))
}

#' loo for MVBU stanfit objects
#'
#' \code{loo} method for \code{\link{MVBU_Stanfit}} objects. For details,
#' see [loo::loo].
#'
#' @param x An \code{\link{MVBU_Stanfit}} object.
#' @param pars Parameter name containing the pointwise log likelihood. (default: \code{"log_lik"})
#' @param ... Additional arguments passed to \code{loo::loo.array}.
#' @param save_psis Logical scalar passed to \code{loo::loo.array}. (default: \code{FALSE})
#' @param cores Number of cores for parallel execution. (default: \code{getOption("mc.cores", 1)})
#'
#' @return A \code{loo} object.
#'
#' @method loo MVBU_Stanfit
#'
#' @importFrom loo loo extract_log_lik relative_eff loo.array
#' @export
loo.MVBU_Stanfit <- function(
  x,
  pars = "log_lik",
  ...,
  save_psis = FALSE,
  cores = getOption("mc.cores", 1)
) {
  .assert_true(length(pars) == 1L, msg = "pars must contain exactly one element.")
  stanfit <- get_stanfit(x)
  LLarray <- loo::extract_log_lik(
    stanfit = stanfit,
    parameter_name = pars,
    merge_chains = FALSE
  )
  r_eff <- loo::relative_eff(x = exp(LLarray), cores = cores)
  loo::loo.array(
    LLarray,
    r_eff = r_eff,
    cores = cores,
    save_psis = save_psis
  )
}

#' @rdname evaluate_model
#' @export
S7::method(evaluate_model, MVBU_Stanfit) <- function(
  model,
  x = NULL,
  response_category = NULL,
  method = "log_lik",
  decision_rule = if (identical(method, "accuracy")) "criterion" else "proportional",
  return_by_x = FALSE,
  ...
) {
  pars_sum <- get_draws(model, summarize = TRUE, nest = TRUE, untransform_cues = TRUE)
  avail_groups <- if (is.factor(pars_sum$group)) levels(pars_sum$group) else unique(pars_sum$group)
  post_groups <- setdiff(avail_groups, "prior")
  target_group <- if (length(post_groups) > 0) post_groups[1] else "prior"
  pars_group <- pars_sum[pars_sum$group == target_group, ]

  cues <- get_cue_labels(model)
  cats <- get_category_labels(model)
  is_nix <- S7::S7_inherits(model, NIX_IdealAdaptorStanfit)

  cat_reps <- lapply(cats, function(cat_name) {
    row_match <- which(pars_group$category == cat_name)
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

  lapse <- if ("lapse_rate" %in% names(pars_group)) pars_group$lapse_rate[1] else 0

  cog_model <- if (is_nix) {
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

  if (is.null(x) || is.null(response_category)) {
    test_df <- get_test_data(model, .recover_from_staninput = TRUE)
    if (is.null(test_df) || nrow(test_df) == 0) {
      .stop("No test data found in the stanfit object. Please supply x and response_category.")
    }
    resp_col <- if ("response" %in% names(test_df)) {
      "response"
    } else if ("category" %in% names(test_df)) {
      "category"
    } else {
      names(test_df)[1]
    }
    x_mat <- as.matrix(test_df[, cues, drop = FALSE])
    response_category <- test_df[[resp_col]]
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

#' @rdname sample_observations
#' @export
S7::method(sample_observations, MVBU_Stanfit) <- function(x, n = 1L, with_replacement = TRUE, randomize_order = TRUE, ...) {
  .assert_true(S7::S7_inherits(x, MVBU_Stanfit), msg = "x must be an MVBU_Stanfit object.")
  data_df <- x@data
  .assert_true(!is.null(data_df) && nrow(data_df) > 0L, msg = "Stanfit model object contains no data.")
  
  idx <- sample(seq_len(nrow(data_df)), size = as.integer(n), replace = isTRUE(with_replacement))
  if (isFALSE(randomize_order)) idx <- sort(idx)
  sampled_df <- data_df[idx, , drop = FALSE]
  rownames(sampled_df) <- NULL
  sampled_df
}


