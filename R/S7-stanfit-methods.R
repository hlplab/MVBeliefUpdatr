#' @include asserts.R
#' @include S7-family.R
#' @include S7-class.R
#' @include S7-generics.R
#' @include S7-transform-information.R
#' @include S7-staninput.R
#' @include S7-stanfit-input.R
#' @include S7-stanfit.R
#' @include S7-stanfit-posterior.R
#' @include S7-class-nix.R
#' @include S7-class-mnix.R
#' @include S7-class-niw.R
NULL

.get_ideal_adaptor_stanfit_constructor <- function(staninput = NULL,
                                                   stanmodel = NULL) {
  if (is.null(stanmodel) && !is.null(staninput)) {
    fam <- tryCatch(
      get_family_name(staninput),
      error = function(e) character(0)
    )
    if (length(fam) > 0 && nzchar(fam)) {
      stanmodel <- .to_stanmodel(fam)
    }
  }

  if (!is.null(stanmodel) && nzchar(stanmodel)) {
    ctor <- .get_stanfit_constructor(stanmodel)
    if (!is.null(ctor)) {
      return(ctor)
    }
    .stop(
      "Unrecognized stanmodel '", stanmodel,
      "'. Cannot determine concrete IdealAdaptorStanfit constructor."
    )
  }

  IdealAdaptorStanfit
}

#' @rdname get_stanfit
#' @export
S7::method(get_stanfit, S7::class_any) <- function(x, ...) {
  .stop("x must be an IdealAdaptorStanfit object.")
}

#' @rdname get_stanfit
#' @export
S7::method(get_stanfit, MVBU_Stanfit) <- function(x, ...) {
  x@stanfit
}

#' @rdname get_model_family
S7::method(
  get_model_family,
  S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior)
) <- function(x, ...) {
  if (!is.null(x@metadata$model_family)) {
    return(x@metadata$model_family)
  }
  if (!is.null(x@metadata$model_name)) {
    return(x@metadata$model_name)
  }
  if (!is.null(x@metadata$stanmodel)) {
    return(x@metadata$stanmodel)
  }
  class(x)[1]
}

#' @rdname get_family_name
S7::method(
  get_family_name,
  S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior)
) <- function(x, ...) {
  fam_raw <- get_model_family(x, ...)
  .get_family_from_class_name(fam_raw)
}

#' @rdname get_family
S7::method(
  get_family,
  S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior)
) <- function(x, ...) {
  fam_name <- get_family_name(x, ...)
  .get_model_family(fam_name)
}


#' @rdname get_stanfit
#' @export
S7::method(set_stanfit, list(S7::class_any, S7::class_any)) <- function(
  x,
  stanfit,
  ...
) {
  .stop("x must be an IdealAdaptorStanfit object.")
}

#' @rdname get_stanfit
#' @export
S7::method(set_stanfit, list(MVBU_Stanfit, S7::class_any)) <- function(
  x,
  stanfit,
  ...
) {
  # no assertions for stanfit here since the @<- assignment operator applied to
  # S7 objects will automatically call the validator for the class, which
  # already checks that the stanfit is valid.
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
S7::method(get_staninput, S7::class_any) <- function(x, ...) {
  .stop("x must be an IdealAdaptorStanfit object.")
}

#' @rdname get_staninput
S7::method(get_staninput, MVBU_Stanfit) <- function(x, ...) {
  x@staninput
}

#' @rdname get_transform_information
#' @export
S7::method(get_transform_information, S7::class_any) <- function(x, ...) {
  .stop(paste0(
    "x must be an MVBU_Stanfit, MVBU_StanfitPosterior, or ",
    "IdealAdaptorStanfitInput object."
  ))
}

#' @rdname get_transform_information
S7::method(
  get_transform_information,
  S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior, IdealAdaptorStanfitInput)
) <- function(x, ...) {
  x@transform_information
}

#' @rdname get_cue_labels
S7::method(
  get_cue_labels,
  S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior)
) <- function(x, indices = NULL, ...) {
  cues <- if (!is.null(x@metadata$label_information$cue)) {
    x@metadata$label_information$cue
  } else {
    character(0)
  }
  if (!is.null(indices)) cues[indices] else cues
}

#' @rdname get_category_labels
S7::method(
  get_category_labels,
  S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior)
) <- function(x, indices = NULL, ...) {
  cats <- if (!is.null(x@metadata$label_information$category)) {
    x@metadata$label_information$category
  } else {
    character(0)
  }
  if (!is.null(indices)) cats[indices] else cats
}

#' @rdname get_category_labels
S7::method(
  get_response_category_labels,
  S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior)
) <- function(
  x,
  indices = NULL,
  ...
) {
  rcats <- if (!is.null(x@metadata$label_information$response_category)) {
    x@metadata$label_information$response_category
  } else {
    character(0)
  }
  if (!is.null(indices)) rcats[indices] else rcats
}

#' @rdname get_group_labels
S7::method(
  get_group_labels,
  S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior)
) <- function(x, indices = NULL, include_prior = FALSE, ...) {
  grps <- if (!is.null(x@metadata$label_information$group)) {
    x@metadata$label_information$group
  } else {
    character(0)
  }
  if (include_prior) grps <- append("prior", grps)
  if (!is.null(indices)) grps[indices] else grps
}

#' @rdname get_labels
S7::method(
  get_labels,
  S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior)
) <- function(x, ...) {
  list(
    cue = get_cue_labels(x, ...),
    category = get_category_labels(x, ...),
    response_category = get_response_category_labels(x, ...),
    group = get_group_labels(x, ...)
  )
}

#' @rdname get_original_variable_names
S7::method(
  get_original_variable_names,
  S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior)
) <- function(
  x,
  variable = c(
    "group", "group_unique", "category", "response_category", "cues"
  ),
  ...
) {
  orig <- x@metadata$original_variable_names
  if (
    is.null(orig) &&
      S7::S7_inherits(x, MVBU_Stanfit) &&
      !is.null(get_staninput(x))
  ) {
    return(
      get_original_variable_names(get_staninput(x), variable = variable, ...)
    )
  }
  if (is.null(orig)) {
    .stop(
      "Original variable names are not available in this object's metadata."
    )
  }
  if (missing(variable)) {
    return(orig)
  }
  valid_vars <- c(
    "group", "group_unique", "category", "response_category", "cues"
  )
  selected <- match.arg(variable, valid_vars, several.ok = TRUE)
  if (length(selected) == 1L) {
    orig[[selected]]
  } else {
    orig[selected]
  }
}

#' @rdname get_data
S7::method(
  get_data,
  S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior)
) <- function(
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
S7::method(
  get_exposure_data,
  S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior)
) <- function(
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
S7::method(
  get_test_data,
  S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior)
) <- function(
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


#' @rdname get_exposure_category_statistics
S7::method(
  get_exposure_category_statistics,
  S7::new_union(IdealAdaptorStanfitInput, MVBU_Stanfit, MVBU_StanfitPosterior)
) <- function(
  x,
  categories = NULL,
  groups = NULL,
  stats = c("n", "mean", "cov", "css", "uss"),
  simplify = TRUE,
  ...
) {
  valid_stats <- c("n", "mean", "cov", "css", "uss")
  stats <- match.arg(stats, valid_stats, several.ok = TRUE)

  avail_cats <- get_category_labels(x)
  if (is.null(categories)) {
    categories <- avail_cats
  } else {
    categories <- .validate_requested_labels(
      categories,
      avail_cats,
      label_type = "category"
    )
  }

  avail_grps <- get_group_labels(x, include_prior = FALSE)
  if (is.null(groups)) {
    groups <- avail_grps
  } else {
    groups <- .validate_requested_labels(
      groups,
      avail_grps,
      label_type = "group"
    )
  }

  meta <- get_metadata(x)
  target_stats <- meta$exposure_category_statistics
  .assert_that(
    !is.null(target_stats) && is.data.frame(target_stats),
    msg = "Exposure category statistics are not available in this object."
  )

  res <- as.data.frame(target_stats)
  if ("category" %in% names(res) && !is.null(categories)) {
    res <- res[res$category %in% categories, , drop = FALSE]
  }
  if ("group" %in% names(res) && !is.null(groups)) {
    res <- res[res$group %in% groups, , drop = FALSE]
  }
  if ("cov" %in% stats && !"cov" %in% names(res) &&
    all(c("css", "n") %in% names(res))) {
    res$cov <- mapply(function(css, n) {
      if (is.matrix(css)) {
        css2cov(css, n)
      } else {
        k <- length(css)
        c_mat <- if (n > 1L) {
          diag(css / (n - 1L), k, k)
        } else {
          matrix(NA_real_, k, k)
        }
        dimnames(c_mat) <- list(names(css), names(css))
        c_mat
      }
    }, res$css, res$n, SIMPLIFY = FALSE)
  }
  if ("uss" %in% stats && !"uss" %in% names(res) &&
    all(c("css", "n", "mean") %in% names(res))) {
    res$uss <- mapply(function(css, n, m) {
      if (is.matrix(css)) {
        css2uss(css, n, m)
      } else {
        u_vec <- css + n * (m^2)
        names(u_vec) <- names(css)
        u_vec
      }
    }, res$css, res$n, res$mean, SIMPLIFY = FALSE)
  }
  if ("css" %in% stats && !"css" %in% names(res) &&
    all(c("cov", "n") %in% names(res))) {
    res$css <- mapply(function(cov, n) {
      if (is.matrix(cov)) {
        cov2css(cov, n)
      } else {
        cov * (n - 1L)
      }
    }, res$cov, res$n, SIMPLIFY = FALSE)
  }
  .assert_that(
    all(stats %in% names(res)),
    msg = "Exposure category statistics are not available in this object."
  )

  cols <- c(intersect(c("group", "category"), names(res)), stats)
  df <- tibble::as_tibble(res[, cols, drop = FALSE])
  if (isTRUE(simplify) && nrow(df) == 1L && length(stats) == 1L) {
    return(df[[stats]][[1L]])
  }
  df
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
        if (length(.data$S) == d_dim && d_dim > 1) {
          mat <- diag(.data$S, nrow = d_dim, ncol = d_dim)
          dimnames(mat) <- list(cues_order, cues_order)
          mat
        } else {
          matrix(
            .data$S,
            nrow = d_dim,
            ncol = d_dim,
            dimnames = list(cues_order, cues_order)
          )
        }
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
  if (is.null(groups)) groups <- get_group_labels(fit, include_prior = TRUE)
  if (missing(which)) {
    which <- if ("prior" %in% groups) {
      if (length(groups) > 1) "both" else "prior"
    } else {
      "posterior"
    }
  }
  if (is.null(seed) && !is.null(ndraws)) {
    seed <- stats::runif(1, -1e6, 1e6)
  }
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
  categories <- .validate_requested_labels(
    categories,
    get_category_labels(fit),
    label_type = "category"
  )
  groups <- .validate_requested_labels(
    groups,
    get_group_labels(fit, include_prior = TRUE),
    label_type = "group"
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

  if (!is.null(ndraws)) {
    n_avail <- get_number_of_draws(fit)
    if (n_avail > 0L && ndraws > n_avail) {
      ndraws <- n_avail
    }
  }

  if ("prior" %in% groups && length(groups) > 1) {
    d.prior <- get_draws(
      fit = fit,
      categories = categories,
      groups = "prior",
      ndraws = ndraws,
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
  fam <- get_model_family(fit)
  is_nix <- grepl("NIX", fam, ignore.case = TRUE) &&
    !grepl("MNIX", fam, ignore.case = TRUE)
  is_mnix <- grepl("MNIX", fam, ignore.case = TRUE)
  is_niw <- grepl("NIW", fam, ignore.case = TRUE)

  stanfit_pars <- dimnames(stanfit)$parameters
  if (is.null(stanfit_pars) && inherits(stanfit, "stanfit")) {
    stanfit_pars <- stanfit@model_pars
  }
  has_indexed_kappa <- any(startsWith(stanfit_pars, paste0(kappa, "[")))
  has_indexed_nu <- any(startsWith(stanfit_pars, paste0(nu, "[")))

  kappa_expr <- if (has_indexed_kappa) {
    rlang::expr((!!rlang::sym(kappa))[!!!rlang::syms(pars.index)])
  } else {
    rlang::expr(!!rlang::sym(kappa))
  }
  nu_expr <- if (has_indexed_nu) {
    rlang::expr((!!rlang::sym(nu))[!!!rlang::syms(pars.index)])
  } else {
    rlang::expr(!!rlang::sym(nu))
  }

  if (is_nix) {
    cue_lab <- get_cue_labels(fit)
    if (length(cue_lab) == 0) cue_lab <- "cue"
    cue_lab <- cue_lab[1]

    d.pars <- stanfit %>%
      tidybayes::spread_draws(
        !!kappa_expr,
        !!nu_expr,
        (!!rlang::sym(m))[!!!rlang::syms(pars.index)],
        (!!rlang::sym(S))[!!!rlang::syms(pars.index)],
        lapse_rate,
        ndraws = ndraws,
        seed = seed
      ) %>%
      dplyr::mutate(cue = cue_lab, cue2 = cue_lab)
  } else if (is_mnix) {
    d.pars <- stanfit %>%
      tidybayes::spread_draws(
        !!kappa_expr,
        !!nu_expr,
        (!!rlang::sym(m))[!!!rlang::syms(pars.index), cue],
        (!!rlang::sym(S))[!!!rlang::syms(pars.index), cue],
        lapse_rate,
        ndraws = ndraws,
        seed = seed
      ) %>%
      dplyr::mutate(cue2 = .data$cue)
  } else if (is_niw) {
    d.pars <- stanfit %>%
      tidybayes::spread_draws(
        !!kappa_expr,
        !!nu_expr,
        (!!rlang::sym(m))[!!!rlang::syms(pars.index), cue],
        (!!rlang::sym(S))[!!!rlang::syms(pars.index), cue, cue2],
        lapse_rate,
        ndraws = ndraws,
        seed = seed
      )
  } else {
    stop("Unknown model family.")
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

  if (nest) {
    d.pars <- .nest_draw_cues(d.pars)
  }

  if (nest && "S" %in% names(d.pars) && "nu" %in% names(d.pars)) {
    d.pars$Sigma_exp <- get_expected_Sigma_from_S(d.pars$S, d.pars$nu)
    if ("kappa" %in% names(d.pars)) {
      d.pars$Sigma_marg <- get_marginal_Sigma_from_S(d.pars$S, d.pars$nu, d.pars$kappa)
    }
  } else if (!nest && "S" %in% names(d.pars) && "nu" %in% names(d.pars)) {
    D <- length(unique(d.pars$cue))
    d.pars$Sigma_exp <- d.pars$S / (d.pars$nu - D - 1)
    if ("kappa" %in% names(d.pars)) {
      d.pars$Sigma_marg <- ((d.pars$kappa + 1) / d.pars$kappa) * d.pars$Sigma_exp
    }
  }

  d.pars <- d.pars %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      category = factor(.data$category, levels = categories),
      group = factor(.data$group, levels = groups)
    )

  d.pars
}

#' @rdname get_expected_category_statistic
#' @export
S7::method(get_expected_category_statistic, MVBU_Stanfit) <- function(
  x,
  categories = NULL,
  groups = NULL,
  statistic = c("mu", "Sigma"),
  ...
) {
  if (is.null(categories)) categories <- get_category_labels(x)
  if (is.null(groups)) groups <- get_group_labels(x, include_prior = TRUE)
  .assert_that(all(statistic %in% c("mu", "Sigma")))
  .assert_that(
    any(is.factor(categories), is.character(categories), is.numeric(categories))
  )
  .assert_that(
    any(is.factor(groups), is.character(groups), is.numeric(groups))
  )
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
      dplyr::filter(.data$category %in% categories, .data$group %in% groups) %>%
      dplyr::select(dplyr::all_of(c("group", "category", paste0(statistic, ".mean")))) %>%
      dplyr::mutate(
        category = factor(.data$category, levels = categories),
        group = factor(.data$group, levels = groups)
      )
  } else {
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
      dplyr::mutate(
        category = factor(.data$category, levels = categories),
        group = factor(.data$group, levels = groups)
      )

    # Populate cache if computing on full default domain
    if (setequal(categories, avail_cats) && setequal(groups, avail_grps)) {
      if (is.null(x@cache$summary)) x@cache$summary <- list()
      x@cache$summary$expected_moments <- draws
    }

    draws <- draws %>%
      dplyr::select(
        dplyr::all_of(c("group", "category", paste0(statistic, ".mean")))
      )
  }

  if (nrow(draws) == 1 && length(statistic) == 1) {
    draws <- draws[[paste0(statistic, ".mean")]][[1]]
  }

  draws
}

#' @rdname get_marginal_category_statistic
#' @export
S7::method(get_marginal_category_statistic, MVBU_Stanfit) <- function(
  x,
  categories = NULL,
  groups = NULL,
  statistic = c("mu", "Sigma"),
  ...
) {
  if (is.null(categories)) categories <- get_category_labels(x)
  if (is.null(groups)) groups <- get_group_labels(x, include_prior = TRUE)
  .assert_that(all(statistic %in% c("mu", "Sigma")))
  .assert_that(
    any(is.factor(categories), is.character(categories), is.numeric(categories))
  )
  .assert_that(
    any(is.factor(groups), is.character(groups), is.numeric(groups))
  )
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
      dplyr::filter(.data$category %in% categories, .data$group %in% groups) %>%
      dplyr::select(dplyr::all_of(c("group", "category", paste0(statistic, ".mean")))) %>%
      dplyr::mutate(
        category = factor(.data$category, levels = categories),
        group = factor(.data$group, levels = groups)
      )
  } else {
    d_raw <- get_draws(
      x,
      categories = categories,
      groups = groups,
      nest = TRUE,
      summarize = FALSE,
      ...
    )
    if (!"Sigma_marg" %in% names(d_raw) && "S" %in% names(d_raw) && "nu" %in% names(d_raw)) {
      d_raw$Sigma_marg <- if ("kappa" %in% names(d_raw)) {
        get_marginal_Sigma_from_S(d_raw$S, d_raw$nu, d_raw$kappa)
      } else {
        get_expected_Sigma_from_S(d_raw$S, d_raw$nu)
      }
    }
    draws <- d_raw %>%
      dplyr::group_by(.data$group, .data$category) %>%
      dplyr::summarise(
        mu.mean = list(purrr::reduce(.data$m, `+`) / length(.data$m)),
        Sigma.mean = list(purrr::reduce(.data$Sigma_marg, `+`) / length(.data$Sigma_marg)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        category = factor(.data$category, levels = categories),
        group = factor(.data$group, levels = groups)
      )

    # Populate cache if computing on full default domain
    if (setequal(categories, avail_cats) && setequal(groups, avail_grps)) {
      if (is.null(x@cache$summary)) x@cache$summary <- list()
      x@cache$summary$marginal_moments <- draws
    }

    draws <- draws %>%
      dplyr::select(
        dplyr::all_of(c("group", "category", paste0(statistic, ".mean")))
      )
  }

  if (nrow(draws) == 1 && length(statistic) == 1) {
    draws <- draws[[paste0(statistic, ".mean")]][[1]]
  }

  draws
}

#' @rdname get_parameters
S7::method(
  get_expected_mu,
  S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior)
) <- function(x, ...) {
  get_expected_category_statistic(x, statistic = "mu", ...)
}

#' @rdname get_parameters
S7::method(
  get_expected_sigma,
  S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior)
) <- function(x, ...) {
  get_expected_category_statistic(x, statistic = "Sigma", ...)
}

#' @rdname get_parameters
S7::method(
  get_marginal_mu,
  S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior)
) <- function(x, ...) {
  get_marginal_category_statistic(x, statistic = "mu", ...)
}

#' @rdname get_parameters
S7::method(
  get_marginal_sigma,
  S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior)
) <- function(x, ...) {
  get_marginal_category_statistic(x, statistic = "Sigma", ...)
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
    return(invert_transform(info)@transform.function)
  }
  identity
}

#' @rdname get_number_of_draws
#' @export
S7::method(get_number_of_draws, S7::class_any) <- function(fit, ...) {
  stanfit <- get_stanfit(fit)
  if (is.null(stanfit)) {
    return(0L)
  }
  posterior::ndraws(posterior::as_draws(stanfit))
}

#' @rdname get_number_of_draws
#' @export
S7::method(get_random_draw_indices, S7::class_any) <- function(
  fit,
  ndraws = NULL,
  ...
) {
  n.all.draws <- get_number_of_draws(fit)
  if (is.null(ndraws)) {
    return(seq_len(n.all.draws))
  }
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
          .get_categorization_function_from_stanfit_draws(
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
.get_categorization_function_from_stanfit_draws <- function(x, ...) {
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
  cat("  Model family: ", get_family_name(x), "\n", sep = "")
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
    header = S7::class_any,
    fitted = S7::class_any,
    fixed = S7::class_any,
    expected_moments = S7::class_any,
    marginal_moments = S7::class_any,
    high_rhats = S7::class_any,
    warnings = S7::class_any
  ),
  constructor = function(
    fitted = NULL,
    fixed = NULL,
    expected_moments = NULL,
    marginal_moments = NULL,
    high_rhats = NULL,
    header = NULL,
    warnings = NULL
  ) {
    S7::new_object(
      MVBU_Object(),
      header = header,
      fitted = fitted,
      fixed = fixed,
      expected_moments = expected_moments,
      marginal_moments = marginal_moments,
      high_rhats = high_rhats,
      warnings = warnings
    )
  }
)

S7::method(print, Summary_MVBU_Stanfit) <- function(x, ...) {
  if (!is.null(x@header)) {
    cat(" // MVBeliefUpdatr Stanfit Model Summary\n")
    if (!is.null(x@header$family)) cat(" Model Family: ", x@header$family, "\n", sep = "")
    if (!is.null(x@header$cues)) cat(" Cues:         ", x@header$cues, "\n", sep = "")
    if (!is.null(x@header$categories)) cat(" Categories:   ", x@header$categories, "\n", sep = "")
    if (!is.null(x@header$observations)) cat(" Observations: ", x@header$observations, "\n", sep = "")
    if (!is.null(x@header$draws)) {
      cat(" Draws:        ", x@header$draws$line1, "\n", sep = "")
      if (!is.null(x@header$draws$line2)) cat("               ", x@header$draws$line2, "\n", sep = "")
    }
    cat("\n")
  }

  if (!is.null(x@fixed) && nrow(x@fixed) > 0) {
    cat("Fixed parameters:\n")
    print(
      as.data.frame(x@fixed),
      row.names = FALSE,
      max = nrow(x@fixed) * 100,
      ...
    )
    if (!is.null(x@fitted) && nrow(x@fitted) > 0) {
      cat("\n")
    }
  }

  if (!is.null(x@fitted) && nrow(x@fitted) > 0) {
    cat("Fitted parameters:\n")
    print(
      as.data.frame(x@fitted),
      row.names = FALSE,
      max = nrow(x@fitted) * 100,
      ...
    )
  } else {
    cat("No fitted parameters.\n")
  }

  moments <- if (is.null(x@expected_moments)) {
    x@marginal_moments
  } else if (is.null(x@marginal_moments)) {
    x@expected_moments
  } else {
    dplyr::bind_rows(
      x@expected_moments,
      x@marginal_moments[x@marginal_moments$Parameter != "mu", ]
    )
  }

  if (!is.null(moments) && nrow(moments) > 0) {
    param_levels <- c("mu", "Sigma_exp", "Sigma_marg")
    order_cols <- if ("Cue" %in% names(moments)) {
      list(
        moments$Group,
        moments$Category,
        match(moments$Parameter, param_levels),
        moments$Cue
      )
    } else {
      list(
        moments$Group,
        moments$Category,
        match(moments$Parameter, param_levels),
        moments$Cue1,
        moments$Cue2
      )
    }
    moments <- moments[do.call(order, order_cols), , drop = FALSE]
    cat("\nCategory moments (distributions across draws):\n")
    print(
      as.data.frame(moments),
      row.names = FALSE,
      max = nrow(moments) * 100,
      ...
    )
  }

  div_trans <- if (!is.null(x@warnings$divergent)) x@warnings$divergent else 0
  has_high_rhats <- !is.null(x@high_rhats) && nrow(x@high_rhats) > 0
  if (div_trans > 0 || has_high_rhats) {
    cat("\nWarnings:\n")
    if (div_trans > 0) {
      cat("  * Warning: There were ", div_trans, " divergent transitions after warmup.\n", sep = "")
    }
    if (has_high_rhats) {
      max_rhat <- round(max(x@high_rhats$Rhat, na.rm = TRUE), 2)
      cat("  * Warning: The largest R-hat is ", sprintf("%.2f", max_rhat), ", indicating potential convergence issues.\n", sep = "")
      cat("\nParameters with Rhats > 1.05:\n")
      print(
        as.data.frame(x@high_rhats),
        row.names = FALSE,
        max = nrow(x@high_rhats) * 100,
        ...
      )
    }
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

# Summarize parameter posterior draws into mean, standard deviation, and credible quantiles.
# Handles mean vector (mu) and covariance matrices (expected Sigma_exp and marginal Sigma_marg)
# across all categories and experimental groups.
.summarize_moment_draws <- function(
  d_draws,
  cues,
  type = c("expected", "marginal"),
  probs = c(0.025, 0.5, 0.975),
  model_family = NULL
) {
  if (is.null(d_draws) || nrow(d_draws) == 0L) {
    return(NULL)
  }
  type <- match.arg(type, c("expected", "marginal"), several.ok = TRUE)

  has_exp <- "Sigma_exp" %in% names(d_draws) || "Sigma" %in% names(d_draws)
  has_marg <- "Sigma_marg" %in% names(d_draws)
  exp_col <- if ("Sigma_exp" %in% names(d_draws)) "Sigma_exp" else "Sigma"

  include_exp <- "expected" %in% type && has_exp
  include_marg <- "marginal" %in% type && has_marg

  if (!"m" %in% names(d_draws)) {
    return(NULL)
  }

  fam_clean <- if (!is.null(model_family)) sub("_.*$", "", model_family) else ""
  is_nix_or_mnix <- isTRUE(fam_clean %in% c("NIX", "MNIX"))

  grp_col <- if ("group" %in% names(d_draws)) {
    as.character(d_draws$group)
  } else {
    rep("all", nrow(d_draws))
  }
  cat_col <- as.character(d_draws$category)
  combos <- unique(
    data.frame(group = grp_col, category = cat_col, stringsAsFactors = FALSE)
  )

  rows <- list()
  n_cues <- length(cues)
  prob_names <- if (!is.null(names(probs))) {
    names(probs)
  } else {
    paste0(probs * 100, "%")
  }

  for (k in seq_len(nrow(combos))) {
    grp_k <- combos$group[k]
    cat_k <- combos$category[k]
    idx <- which(grp_col == grp_k & cat_col == cat_k)
    if (length(idx) == 0L) next

    sub_m <- d_draws$m[idx]

    # Mean vector draws
    for (i in seq_len(n_cues)) {
      vals <- vapply(sub_m, function(vec) {
        if (is.numeric(vec) && length(vec) >= i) vec[i] else NA_real_
      }, numeric(1L))
      q <- stats::quantile(vals, probs = probs, na.rm = TRUE)
      row_df <- if (is_nix_or_mnix) {
        tibble::tibble(
          Parameter = "mu",
          Group = grp_k,
          Category = cat_k,
          Cue = cues[i],
          mean = mean(vals, na.rm = TRUE),
          sd = stats::sd(vals, na.rm = TRUE)
        )
      } else {
        tibble::tibble(
          Parameter = "mu",
          Group = grp_k,
          Category = cat_k,
          Cue1 = cues[i],
          Cue2 = "",
          mean = mean(vals, na.rm = TRUE),
          sd = stats::sd(vals, na.rm = TRUE)
        )
      }
      for (p_idx in seq_along(probs)) {
        row_df[[prob_names[p_idx]]] <- q[p_idx]
      }
      rows[[length(rows) + 1L]] <- row_df
    }

    # Expected covariance matrix draws
    if (include_exp) {
      sub_sig_exp <- d_draws[[exp_col]][idx]
      for (i in seq_len(n_cues)) {
        j_range <- if (is_nix_or_mnix) i else i:n_cues
        for (j in j_range) {
          vals <- vapply(sub_sig_exp, function(mat) {
            m <- as.matrix(mat)
            if (nrow(m) >= i && ncol(m) >= j) m[i, j] else NA_real_
          }, numeric(1L))
          q <- stats::quantile(vals, probs = probs, na.rm = TRUE)
          row_df <- if (is_nix_or_mnix) {
            tibble::tibble(
              Parameter = "Sigma_exp",
              Group = grp_k,
              Category = cat_k,
              Cue = cues[i],
              mean = mean(vals, na.rm = TRUE),
              sd = stats::sd(vals, na.rm = TRUE)
            )
          } else {
            tibble::tibble(
              Parameter = "Sigma_exp",
              Group = grp_k,
              Category = cat_k,
              Cue1 = cues[i],
              Cue2 = cues[j],
              mean = mean(vals, na.rm = TRUE),
              sd = stats::sd(vals, na.rm = TRUE)
            )
          }
          for (p_idx in seq_along(probs)) {
            row_df[[prob_names[p_idx]]] <- q[p_idx]
          }
          rows[[length(rows) + 1L]] <- row_df
        }
      }
    }

    # Marginal covariance matrix draws
    if (include_marg) {
      sub_sig_marg <- d_draws[["Sigma_marg"]][idx]
      for (i in seq_len(n_cues)) {
        j_range <- if (is_nix_or_mnix) i else i:n_cues
        for (j in j_range) {
          vals <- vapply(sub_sig_marg, function(mat) {
            m <- as.matrix(mat)
            if (nrow(m) >= i && ncol(m) >= j) m[i, j] else NA_real_
          }, numeric(1L))
          q <- stats::quantile(vals, probs = probs, na.rm = TRUE)
          row_df <- if (is_nix_or_mnix) {
            tibble::tibble(
              Parameter = "Sigma_marg",
              Group = grp_k,
              Category = cat_k,
              Cue = cues[i],
              mean = mean(vals, na.rm = TRUE),
              sd = stats::sd(vals, na.rm = TRUE)
            )
          } else {
            tibble::tibble(
              Parameter = "Sigma_marg",
              Group = grp_k,
              Category = cat_k,
              Cue1 = cues[i],
              Cue2 = cues[j],
              mean = mean(vals, na.rm = TRUE),
              sd = stats::sd(vals, na.rm = TRUE)
            )
          }
          for (p_idx in seq_along(probs)) {
            row_df[[prob_names[p_idx]]] <- q[p_idx]
          }
          rows[[length(rows) + 1L]] <- row_df
        }
      }
    }
  }

  if (length(rows) == 0L) {
    return(NULL)
  }
  dplyr::bind_rows(rows)
}

#' Format fixed parameters as a summary table
#'
#' @param fixed_params An \code{\link{MVBU_FixedParameter}} object or NULL.
#' @param object Optional \code{MVBU_Stanfit} or model object providing label
#'   context (unused, kept for signature consistency).
#' @return A \code{\link[tibble]{tibble}} with fixed parameters formatted
#'   for display.
#' @keywords internal
#' @noRd
.format_fixed_parameters_table <- function(fixed_params, object = NULL) {
  if (is.null(fixed_params) || length(fixed_params) == 0L) {
    return(tibble::tibble())
  }
  if (S7::S7_inherits(fixed_params, MVBU_FixedParameter)) {
    return(tibble::as_tibble(fixed_params))
  }
  tibble::tibble()
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
#' @param probs Numeric vector of probabilities for credible interval quantiles.
#'   (default: \code{c(0.025, 0.5, 0.975)})
#' @param ... Additional arguments passed to \code{rstan::summary}.
#'
#' @rdname summary-methods
#' @exportS3Method base::summary
summary.MVBU_Stanfit <- function(
  object,
  pars = NULL,
  sufficient_pars_only = TRUE,
  indices_as_names = TRUE,
  include_transformed_pars = FALSE,
  probs = c(0.025, 0.5, 0.975),
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
    pars <- grep(
      "^((kappa|nu|m|S|weights)_|weights\\[|lapse_rate|lapse_bias|prior_p|Sigma_noise)",
      pars,
      value = TRUE
    )
    pars <- grep("^m_0_(tau|L_omega|cov)", pars, value = TRUE, invert = TRUE)
    pars <- grep(
      "^((m|S)_0|lapse_rate|weights)_param",
      pars,
      value = TRUE,
      invert = TRUE
    )
    if (!include_transformed_pars) {
      pars <- grep("_transformed", pars, value = TRUE, invert = TRUE)
    }
    fixed_names <- names(get_parameters(object, space = "original")$fixed)
    if (length(fixed_names) > 0) {
      fixed_pat <- paste0(
        "^(", paste(fixed_names, collapse = "|"), ")(\\[|_|$)"
      )
      pars <- grep(fixed_pat, pars, value = TRUE, invert = TRUE)
    }
  }

  raw_sum <- rstan::summary(stanfit, pars = pars, probs = probs, ...)$summary
  fixed_params <- get_parameters(object, space = "original")$fixed
  fixed_summary <- .format_fixed_parameters_table(fixed_params, object = object)

  if (is.null(raw_sum) || nrow(raw_sum) == 0) {
    return(Summary_MVBU_Stanfit(
      fitted = as.data.frame(raw_sum),
      fixed = fixed_summary
    ))
  }

  # Sort and filter output
  raw_df <-
    raw_sum %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    dplyr::mutate(
      name = factor(
        gsub(
          "^(kappa|nu|m|S|weights|lapse_rate|lapse_bias|prior_p|Sigma_noise).*$",
          "\\1",
          .data$Parameter
        ),
        levels = c(
          "weights", "prior_p", "lapse_rate", "lapse_bias",
          "kappa", "nu", "m", "S", "Sigma_noise"
        )
      ),
      distribution = gsub("^.*_(0|n).*$", "\\1", .data$Parameter),
      index = gsub("^.*_(0|n)?\\[(.*)\\]$", "\\2", .data$Parameter),
      index = gsub("^(weights|lapse_bias|prior_p)\\[(.*)\\]$", "\\2", .data$index),
      index = ifelse(.data$index == .data$Parameter, 1, .data$index)
    ) %>%
    tidyr::separate(
      .data$index,
      into = c("i1", "i2", "i3", "i4"),
      sep = ",",
      fill = "right"
    ) %>%
    dplyr::mutate(dplyr::across(c("i1", "i2", "i3", "i4"), as.integer)) %>%
    dplyr::arrange(
      .data$distribution,
      .data$name,
      .data$i1,
      .data$i2,
      .data$i3,
      .data$i4
    )

  category_labels <- get_category_labels(object)
  group_labels <- get_group_labels(object, include_prior = FALSE)
  cue_labels <- get_cue_labels(object)
  fam <- sub("_.*$", "", get_model_family(object))
  is_nix_or_mnix <- fam %in% c("NIX", "MNIX")

  if (indices_as_names) {
    full_summary <-
      raw_df %>%
      dplyr::mutate(
        Parameter = as.character(.data$name),
        `Dist.` = dplyr::case_when(
          .data$distribution == "0" ~ "prior",
          .data$distribution == "n" ~ "posterior",
          TRUE ~ ""
        ),
        Group = dplyr::case_when(
          .data$name %in% c("kappa", "nu", "m", "S") &
            .data$`Dist.` == "posterior" ~ group_labels[.data$i2],
          .data$name == "weights" & !is.na(.data$i2) ~ group_labels[.data$i1],
          TRUE ~ ""
        ),
        Category = dplyr::case_when(
          .data$name %in% c("m", "S") &
            .data$`Dist.` == "prior" ~ category_labels[.data$i1],
          .data$name %in% c("kappa", "nu", "m", "S") &
            .data$`Dist.` == "posterior" ~ category_labels[.data$i1],
          .data$name %in% c("prior_p", "lapse_bias") ~
            category_labels[.data$i1],
          TRUE ~ ""
        )
      )

    if (fam == "NIX") {
      full_summary <- full_summary %>%
        dplyr::mutate(
          Cue = dplyr::case_when(
            .data$name %in% c("m", "S") ~ cue_labels[1],
            .data$name == "weights" ~ cue_labels[1],
            TRUE ~ ""
          )
        ) %>%
        dplyr::relocate(
          tidyselect::all_of(
            c("Parameter", "Dist.", "Group", "Category", "Cue")
          ),
          tidyselect::everything()
        )
    } else if (fam == "MNIX") {
      full_summary <- full_summary %>%
        dplyr::mutate(
          Cue = dplyr::case_when(
            .data$name %in% c("m", "S") &
              .data$`Dist.` == "prior" ~ cue_labels[.data$i2],
            .data$name %in% c("m", "S") &
              .data$`Dist.` == "posterior" ~ cue_labels[.data$i3],
            .data$name == "weights" & is.na(.data$i2) ~
              cue_labels[.data$i1],
            .data$name == "weights" &
              !is.na(.data$i2) &
              is.na(.data$i3) ~ cue_labels[.data$i2],
            .data$name == "weights" & !is.na(.data$i3) ~
              cue_labels[.data$i3],
            TRUE ~ ""
          )
        ) %>%
        dplyr::relocate(
          tidyselect::all_of(
            c("Parameter", "Dist.", "Group", "Category", "Cue")
          ),
          tidyselect::everything()
        )
    } else {
      full_summary <- full_summary %>%
        dplyr::mutate(
          Cue1 = dplyr::case_when(
            .data$name %in% c("m", "S") &
              .data$`Dist.` == "prior" ~ cue_labels[.data$i2],
            .data$name %in% c("m", "S") &
              .data$`Dist.` == "posterior" ~ cue_labels[.data$i3],
            .data$name == "weights" ~ cue_labels[.data$i1],
            TRUE ~ ""
          ),
          Cue2 = dplyr::case_when(
            .data$name %in% c("S") &
              .data$`Dist.` == "prior" ~ cue_labels[.data$i3],
            .data$name %in% c("S") &
              .data$`Dist.` == "posterior" ~ cue_labels[.data$i4],
            TRUE ~ ""
          )
        ) %>%
        dplyr::relocate(
          tidyselect::all_of(
            c("Parameter", "Dist.", "Group", "Category", "Cue1", "Cue2")
          ),
          tidyselect::everything()
        )
    }
  } else {
    full_summary <- raw_df %>%
      dplyr::relocate(tidyselect::all_of("Parameter"), tidyselect::everything())
  }

  # Build high_rhats table across all parameters before sufficient_pars_only
  Rhats <- full_summary[["Rhat"]]
  high_rhat_df <- NULL
  if (!is.null(Rhats) && any(Rhats > 1.05, na.rm = TRUE)) {
    .warning(
      "Parts of the model have not converged (some Rhats are > 1.05). ",
      "Be careful when analysing the results! We recommend running ",
      "more iterations and/or setting stronger priors."
    )
    cue_cols <- if (is_nix_or_mnix) "Cue" else c("Cue1", "Cue2")
    id_cols <- if (indices_as_names) {
      c("Parameter", "Dist.", "Group", "Category", cue_cols)
    } else {
      "Parameter"
    }
    high_rhat_df <- full_summary[
      !is.na(full_summary$Rhat) & full_summary$Rhat > 1.05,
      c(id_cols, "Rhat"),
      drop = FALSE
    ]
  }

  if (sufficient_pars_only) {
    full_summary <- full_summary %>%
      dplyr::filter(
        .data$distribution == "0" |
          .data$name %in% c(
            "weights", "lapse_rate", "lapse_bias", "prior_p", "Sigma_noise"
          )
      )
  }

  full_summary <-
    full_summary %>%
    dplyr::select(
      -dplyr::any_of(c("name", "distribution", "i1", "i2", "i3", "i4"))
    )

  div_trans <- tryCatch(
    sum(nuts_params(object, pars = "divergent__")$Value),
    error = function(e) 0
  )
  adapt_delta <- tryCatch(
    control_params(object)$adapt_delta,
    error = function(e) NULL
  )
  if (div_trans > 0) {
    .warning(
      "There were ", div_trans, " divergent transitions after warmup. ",
      if (!is.null(adapt_delta)) {
        paste0("Increasing adapt_delta above ", adapt_delta, " may help. ")
      } else {
        ""
      },
      "See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup"
    )
  }

  if (nrow(fixed_summary) > 0 && nrow(full_summary) > 0) {
    drop_idx <- integer(0)
    for (i in seq_len(nrow(fixed_summary))) {
      fp <- fixed_summary[i, ]
      cue_match <- if ("Cue" %in% names(fp)) {
        fp$Cue == "" | full_summary$Cue == fp$Cue
      } else {
        (fp$Cue1 == "" | full_summary$Cue1 == fp$Cue1) &
          (fp$Cue2 == "" | full_summary$Cue2 == fp$Cue2)
      }
      matches <- which(
        full_summary$Parameter == fp$Parameter &
          (fp$`Dist.` == "" | full_summary$`Dist.` == fp$`Dist.`) &
          (fp$Category == "" | full_summary$Category == fp$Category) &
          cue_match
      )
      drop_idx <- c(drop_idx, matches)
    }
    if (length(drop_idx) > 0) {
      full_summary <- full_summary[-unique(drop_idx), , drop = FALSE]
    }
  }

  # If cached summary stats are present with matching probs, use them;
  # otherwise compute dynamically
  cache <- .get_cache(object)
  exp_moments <- NULL
  marg_moments <- NULL

  prob_names <- if (!is.null(names(probs))) {
    names(probs)
  } else {
    paste0(probs * 100, "%")
  }
  if (
    !is.null(cache$summary$expected_moments) &&
      all(prob_names %in% names(cache$summary$expected_moments))
  ) {
    exp_moments <- cache$summary$expected_moments
  }
  if (
    !is.null(cache$summary$marginal_moments) &&
      all(prob_names %in% names(cache$summary$marginal_moments))
  ) {
    marg_moments <- cache$summary$marginal_moments
  }

  if (is.null(exp_moments) || is.null(marg_moments)) {
    d_draws <- tryCatch(
      get_draws(object, nest = TRUE, summarize = FALSE),
      error = function(e) NULL
    )
    if (!is.null(d_draws) && nrow(d_draws) > 0L) {
      moments <- .summarize_moment_draws(
        d_draws,
        cue_labels,
        type = c("expected", "marginal"),
        probs = probs,
        model_family = fam
      )
      if (!is.null(moments)) {
        if (is.null(exp_moments)) {
          exp_moments <- moments[
            moments$Parameter %in% c("mu", "Sigma_exp"),
          ]
        }
        if (is.null(marg_moments)) {
          marg_moments <- moments[
            moments$Parameter %in% c("mu", "Sigma_marg"),
          ]
        }
      }
    }
  }

  # Build header metadata
  model_cls <- class(object)[1]
  fam_clean <- gsub("Stanfit$", "", model_cls)
  fam_clean <- gsub("_", " ", fam_clean)
  cues <- get_cue_labels(object)
  cats <- get_category_labels(object)
  grps <- get_group_labels(object, include_prior = FALSE)

  orig_data <- get_data(object, space = "original")
  n_obs <- if (!is.null(object@staninput) && !is.null(object@staninput@values$n_test)) {
    object@staninput@values$n_test
  } else if (!is.null(orig_data) && nrow(orig_data) > 0) {
    nrow(orig_data)
  } else {
    0
  }
  L_groups <- if (!is.null(object@staninput) && !is.null(object@staninput@values$L)) {
    object@staninput@values$L
  } else {
    max(1L, length(grps))
  }
  obs_str <- if (L_groups > 1) {
    paste0("Total test observations = ", n_obs, " across L = ", L_groups, " exposure groups")
  } else {
    paste0("Total test observations = ", n_obs)
  }

  chains <- tryCatch(stanfit@sim$chains, error = function(e) 1)
  iter <- tryCatch(stanfit@sim$iter, error = function(e) 0)
  warmup <- tryCatch(stanfit@sim$warmup, error = function(e) 0)
  thin <- tryCatch(stanfit@sim$thin, error = function(e) 1)
  total_draws <- get_number_of_draws(object)

  draws_line1 <- paste0(chains, " chains, each with iter = ", iter, "; warmup = ", warmup, "; thin = ", thin)
  draws_line2 <- paste0("Total post-warmup draws = ", total_draws)

  header_info <- list(
    family = fam_clean,
    cues = paste0(paste(cues, collapse = ", "), " (K = ", length(cues), ")"),
    categories = paste0(paste(cats, collapse = ", "), " (M = ", length(cats), ")"),
    observations = obs_str,
    draws = list(line1 = draws_line1, line2 = draws_line2)
  )

  return(Summary_MVBU_Stanfit(
    fitted = full_summary,
    fixed = fixed_summary,
    expected_moments = exp_moments,
    marginal_moments = marg_moments,
    high_rhats = high_rhat_df,
    header = header_info,
    warnings = list(divergent = div_trans)
  ))
}
S7::method(summary, MVBU_Stanfit) <- summary.MVBU_Stanfit

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

.combine_multigroup_evaluation <- function(res_list, method, return_by_x) {
  if (isTRUE(return_by_x)) {
    if (length(method) == 1L) {
      dfs <- lapply(res_list, function(it) {
        df <- it$eval
        df$group <- it$group
        df[, c("group", setdiff(names(df), "group")), drop = FALSE]
      })
      return(dplyr::bind_rows(dfs))
    }
    out <- list()
    for (m in method) {
      dfs <- lapply(res_list, function(it) {
        df <- it$eval[[m]]
        df$group <- it$group
        df[, c("group", setdiff(names(df), "group")), drop = FALSE]
      })
      out[[m]] <- dplyr::bind_rows(dfs)
    }
    return(out)
  }

  total_n <- sum(vapply(res_list, function(it) it$n, integer(1)))
  calc_method_val <- function(m) {
    if (identical(m, "accuracy")) {
      sum(vapply(res_list, function(it) {
        val <- if (length(method) == 1L) it$eval else it$eval[[m]]
        val * it$n
      }, numeric(1))) / total_n
    } else {
      sum(vapply(res_list, function(it) {
        if (length(method) == 1L) it$eval else it$eval[[m]]
      }, numeric(1)))
    }
  }

  if (length(method) == 1L) {
    return(calc_method_val(method))
  }

  out <- list()
  for (m in method) {
    out[[m]] <- calc_method_val(m)
  }
  out
}

#' @rdname evaluate_model
S7::method(
  evaluate_model,
  S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior)
) <- function(
  model,
  x = NULL,
  response_category = NULL,
  method = "log_lik",
  decision_rule = if (identical(method, "accuracy")) "criterion" else "proportional",
  return_by_x = FALSE,
  ...
) {
  post_obj <- as_MVBU_stanfit_posterior(model)
  dots <- list(...)
  group_arg <- dots$group
  cues <- get_cue_labels(post_obj)

  if (is.null(x) || is.null(response_category)) {
    test_df <- get_test_data(post_obj)
    if (is.null(test_df) || nrow(test_df) == 0L) {
      .stop(
        "No test data found in the model object. ",
        "Please supply x and response_category."
      )
    }
    if (!"response_category" %in% names(test_df)) {
      .stop("No 'response_category' column found in test data frame.")
    }
    if (!is.null(group_arg)) {
      test_df <- test_df[test_df$group %in% group_arg, , drop = FALSE]
      if (nrow(test_df) == 0L) {
        .stop(
          "No test data found for specified group: ",
          paste(group_arg, collapse = ", ")
        )
      }
    }
    has_groups <- "group" %in% names(test_df) && length(unique(test_df$group)) > 1L
    if (has_groups) {
      unique_grps <- unique(test_df$group)
      res_list <- lapply(unique_grps, function(g) {
        sub_df <- test_df[test_df$group == g, , drop = FALSE]
        x_mat <- as.matrix(sub_df[, cues, drop = FALSE])
        resp <- sub_df[["response_category"]]
        valid_idx <- !is.na(resp)
        if (any(!valid_idx)) {
          x_mat <- x_mat[valid_idx, , drop = FALSE]
          resp <- resp[valid_idx]
        }
        cog_model_g <- .as_cognitive_model_from_posterior(post_obj, group = g)
        ev <- evaluate_model(
          cog_model_g,
          x = x_mat,
          response_category = resp,
          method = method,
          decision_rule = decision_rule,
          return_by_x = return_by_x,
          ...
        )
        list(group = g, n = length(resp), eval = ev)
      })
      return(.combine_multigroup_evaluation(
        res_list,
        method = method,
        return_by_x = return_by_x
      ))
    }

    grp <- if ("group" %in% names(test_df)) unique(test_df$group)[1] else group_arg
    cog_model <- .as_cognitive_model_from_posterior(post_obj, group = grp)
    x_mat <- as.matrix(test_df[, cues, drop = FALSE])
    response_category <- test_df[["response_category"]]
    valid_idx <- !is.na(response_category)
    if (any(!valid_idx)) {
      x_mat <- x_mat[valid_idx, , drop = FALSE]
      response_category <- response_category[valid_idx]
    }
    x <- x_mat
    return(evaluate_model(
      cog_model,
      x = x,
      response_category = response_category,
      method = method,
      decision_rule = decision_rule,
      return_by_x = return_by_x,
      ...
    ))
  }

  if (is.null(group_arg)) {
    model_groups <- setdiff(get_group_labels(post_obj), "prior")
    if (length(model_groups) == 0L) model_groups <- get_group_labels(post_obj)

    if (is.data.frame(x) &&
      "group" %in% names(x) &&
      length(unique(x$group)) > 1L) {
      unique_grps <- unique(x$group)
      if (is.data.frame(response_category) &&
        "response_category" %in% names(response_category)) {
        response_category <- response_category$response_category
      }
      res_list <- lapply(unique_grps, function(g) {
        g_idx <- which(x$group == g)
        sub_x <- as.matrix(x[g_idx, cues, drop = FALSE])
        sub_resp <- response_category[g_idx]
        cog_model_g <- .as_cognitive_model_from_posterior(post_obj, group = g)
        ev <- evaluate_model(
          cog_model_g,
          x = sub_x,
          response_category = sub_resp,
          method = method,
          decision_rule = decision_rule,
          return_by_x = return_by_x,
          ...
        )
        list(group = g, n = length(sub_resp), eval = ev)
      })
      return(.combine_multigroup_evaluation(
        res_list,
        method = method,
        return_by_x = return_by_x
      ))
    }

    if (length(model_groups) > 1L &&
      (!is.data.frame(x) || !"group" %in% names(x))) {
      if (is.data.frame(response_category) &&
        "response_category" %in% names(response_category)) {
        response_category <- response_category$response_category
      }
      res_list <- lapply(model_groups, function(g) {
        cog_model_g <- .as_cognitive_model_from_posterior(post_obj, group = g)
        ev <- evaluate_model(
          cog_model_g,
          x = x,
          response_category = response_category,
          method = method,
          decision_rule = decision_rule,
          return_by_x = return_by_x,
          ...
        )
        list(group = g, n = length(response_category), eval = ev)
      })
      return(.combine_multigroup_evaluation(
        res_list,
        method = method,
        return_by_x = return_by_x
      ))
    }
  }

  cog_model <- .as_cognitive_model_from_posterior(post_obj, group = group_arg)
  evaluate_model(
    cog_model,
    x = x,
    response_category = response_category,
    method = method,
    decision_rule = decision_rule,
    return_by_x = return_by_x,
    ...
  )
}

#' @rdname sample_observations
#' @export
S7::method(sample_observations, MVBU_Stanfit) <- function(
  x,
  n = 1L,
  with_replacement = TRUE,
  randomize_order = TRUE,
  ...
) {
  .assert_true(
    S7::S7_inherits(x, MVBU_Stanfit),
    msg = "x must be an MVBU_Stanfit object."
  )
  data_df <- get_data(x)
  .assert_true(
    !is.null(data_df) && nrow(data_df) > 0L,
    msg = "Stanfit model object contains no data."
  )

  idx <- sample(
    seq_len(nrow(data_df)),
    size = as.integer(n),
    replace = isTRUE(with_replacement)
  )
  if (isFALSE(randomize_order)) idx <- sort(idx)
  sampled_df <- data_df[idx, , drop = FALSE]
  rownames(sampled_df) <- NULL
  sampled_df
}

#' @rdname reconstruct_update_history
#' @param uncertainty_treatment Character string specifying treatment of posterior parameter
#'   uncertainty: \code{"marginalize"} (default) reconstructs trajectories across posterior
#'   draws and summarizes them, while \code{"discard"} reconstructs a single trajectory
#'   from the prior point estimate.
#' @param ndraws Number of posterior draws to use when \code{uncertainty_treatment = "marginalize"}.
#'   Defaults to \code{20L}.
#' @param step_size Observation step size between checkpoints.
#' @param groups Groups to reconstruct.
#' @param categories Categories to reconstruct.
#' @param parallel Logical; whether to parallelize draw reconstruction.
#' @param seed Random seed.
#'
# Internal helper: Create an ideal adaptor category representation dynamically
# Factory dispatcher that constructs a category representation (NIX, MNIX, NIW)
# without hardcoding if/else checks across calling functions.
.create_ideal_adaptor_representation <- function(
  model_family,
  category_labels,
  cue_labels,
  m,
  s,
  kappa,
  nu
) {
  fam <- toupper(as.character(model_family))

  if (grepl("NIX", fam) && !grepl("MNIX", fam)) {
    new_nix_category_representation(
      category_labels = category_labels,
      cue_labels = cue_labels,
      m = as.numeric(m)[1],
      sigma2 = as.numeric(s)[1] / as.numeric(nu),
      kappa = as.numeric(kappa),
      nu = as.numeric(nu)
    )
  } else if (grepl("MNIX", fam)) {
    k_comp <- length(as.vector(m))
    new_mnix_category_representation(
      category_labels = category_labels,
      cue_labels = cue_labels,
      m = as.vector(m),
      sigma2 = as.vector(s) / as.numeric(nu),
      kappa = rep_len(as.numeric(kappa), k_comp),
      nu = rep_len(as.numeric(nu), k_comp)
    )
  } else if (grepl("NIW", fam)) {
    new_niw_category_representation(
      category_labels = category_labels,
      cue_labels = cue_labels,
      m = as.vector(m),
      S = as.matrix(s),
      kappa = as.numeric(kappa),
      nu = as.numeric(nu)
    )
  } else {
    .stop(sprintf(
      "Unsupported ideal adaptor model family for representation: '%s'.",
      model_family
    ))
  }
}

#' @export
S7::method(reconstruct_update_history, MVBU_Stanfit) <- function(
  object,
  uncertainty_treatment = c("marginalize", "discard"),
  ndraws = 20L,
  step_size = 10L,
  groups = NULL,
  categories = NULL,
  parallel = FALSE,
  seed = 42L,
  ...
) {
  dots <- list(...)
  if (!missing(uncertainty_treatment)) {
    uncertainty_treatment <- match.arg(uncertainty_treatment, c("marginalize", "discard"))
  } else {
    uncertainty_treatment <- "marginalize"
  }

  cues <- get_cue_labels(object)
  cats <- if (!is.null(categories)) categories else get_category_labels(object)
  model_fam <- get_model_family(object)

  # Extract exposure data from Stanfit object
  exp_df <- get_exposure_data(object)
  if (!is.null(groups) && "group" %in% names(exp_df)) {
    exp_df <- exp_df[exp_df$group %in% groups, , drop = FALSE]
  }
  if (!is.null(categories) && "category" %in% names(exp_df)) {
    exp_df <- exp_df[exp_df$category %in% categories, , drop = FALSE]
  }

  n_total <- if (is.null(exp_df)) 0L else nrow(exp_df)
  step_size <- as.integer(step_size)
  if (n_total == 0L) {
    checkpoints <- 0L
  } else {
    checkpoints <- unique(c(0L, seq(step_size, n_total, by = step_size), n_total))
  }
  lbls <- paste0("N_", checkpoints)

  if (identical(uncertainty_treatment, "discard")) {
    # Extract expected prior parameters & instantiate prior ideal adaptor
    pars_sum <- get_draws(object, groups = "prior", summarize = TRUE, nest = TRUE)
    if (nrow(pars_sum) == 0L) {
      pars_sum <- get_draws(object, summarize = TRUE, nest = TRUE)
      pars_sum <- pars_sum[1:length(get_category_labels(object)), ]
    }

    cat_reps <- lapply(cats, function(cat_name) {
      row_match <- which(pars_sum$category == cat_name)
      if (length(row_match) == 0L) row_match <- 1L
      m_val <- pars_sum$m[[row_match]]
      s_val <- pars_sum$S[[row_match]]
      kappa_val <- pars_sum$kappa[row_match]
      nu_val <- pars_sum$nu[row_match]
      .create_ideal_adaptor_representation(
        model_family = model_fam,
        category_labels = cat_name,
        cue_labels = cues,
        m = m_val,
        s = s_val,
        kappa = kappa_val,
        nu = nu_val
      )
    })
    names(cat_reps) <- cats

    prior_model <- as_ideal_adaptor(
      new_category_representation_template(cat_reps)
    )

    if (n_total == 0L) {
      return(as_model_list(list(prior_model), model_labels = "Prior"))
    }

    model_seq <- vector("list", length(checkpoints))
    model_seq[[1]] <- prior_model

    curr <- prior_model
    for (idx in seq_along(checkpoints)[-1]) {
      start_i <- checkpoints[idx - 1] + 1L
      end_i <- checkpoints[idx]
      sub_exp <- exp_df[start_i:end_i, , drop = FALSE]
      curr <- update_template(curr, observations = sub_exp)
      model_seq[[idx]] <- curr
    }

    return(as_model_list(model_seq, model_labels = lbls, metadata = get_metadata(object)))
  } else {
    # Full posterior MCMC draws representation
    if (is.null(ndraws)) ndraws <- 20L
    n_avail <- get_number_of_draws(object)
    if (n_avail > 0L && ndraws > n_avail) {
      ndraws <- n_avail
    }

    d_prior <- get_draws(
      object,
      groups = "prior",
      ndraws = ndraws,
      nest = TRUE,
      seed = seed
    )
    if (nrow(d_prior) == 0L) {
      d_prior <- get_draws(
        object,
        ndraws = ndraws,
        nest = TRUE,
        seed = seed
      )
    }
    if (!is.null(categories)) {
      d_prior <- d_prior[d_prior$category %in% categories, , drop = FALSE]
    }
    draw_ids <- unique(d_prior$.draw)

    initial_models <- lapply(draw_ids, function(d_id) {
      sub_draw <- d_prior[d_prior$.draw == d_id, ]
      cat_reps <- lapply(cats, function(cat_name) {
        row_match <- which(sub_draw$category == cat_name)
        if (length(row_match) == 0L) row_match <- 1L
        m_val <- sub_draw$m[[row_match]]
        s_val <- sub_draw$S[[row_match]]
        kappa_val <- sub_draw$kappa[row_match]
        nu_val <- sub_draw$nu[row_match]
        .create_ideal_adaptor_representation(
          model_family = model_fam,
          category_labels = cat_name,
          cue_labels = cues,
          m = m_val,
          s = s_val,
          kappa = kappa_val,
          nu = nu_val
        )
      })
      names(cat_reps) <- cats
      as_ideal_adaptor(
        new_category_representation_template(cat_reps)
      )
    })

    update_one_draw <- function(model_i) {
      seq_i <- vector("list", length(checkpoints))
      seq_i[[1]] <- model_i
      curr_i <- model_i
      if (length(checkpoints) > 1L) {
        for (idx in seq_along(checkpoints)[-1]) {
          start_i <- checkpoints[idx - 1] + 1L
          end_i <- checkpoints[idx]
          sub_exp <- exp_df[start_i:end_i, , drop = FALSE]
          curr_i <- update_template(curr_i, observations = sub_exp)
          seq_i[[idx]] <- curr_i
        }
      }
      seq_i
    }

    if (isTRUE(parallel) && .Platform$OS.type != "windows") {
      draw_trajectories <- parallel::mclapply(initial_models, update_one_draw)
    } else {
      draw_trajectories <- lapply(initial_models, update_one_draw)
    }

    model_seq <- vector("list", length(checkpoints))
    target_grp <- if (!is.null(groups) && length(groups) > 0) groups[1] else "all"

    for (k in seq_along(checkpoints)) {
      models_at_k <- lapply(draw_trajectories, function(traj) traj[[k]])

      draw_rows <- list()
      for (d_idx in seq_along(models_at_k)) {
        mod <- models_at_k[[d_idx]]
        reps <- get_category_representations(mod)
        for (c_idx in seq_along(cats)) {
          rep_c <- reps[[c_idx]]
          draw_rows[[length(draw_rows) + 1L]] <- tibble::tibble(
            .draw = draw_ids[d_idx],
            group = target_grp,
            category = cats[c_idx],
            m = list(rep_c@m),
            S = list(rep_c@S),
            kappa = rep_c@kappa,
            nu = rep_c@nu,
            Sigma_exp = list(get_expected_sigma(rep_c)),
            Sigma = list(get_expected_sigma(rep_c)),
            Sigma_marg = list(get_marginal_sigma(rep_c))
          )
        }
      }
      draws_k <- dplyr::bind_rows(draw_rows)
      summary_k <- draws_k %>%
        dplyr::group_by(.data$group, .data$category) %>%
        dplyr::summarise(
          m = list(Reduce(`+`, .data$m) / length(.data$m)),
          S = list(Reduce(`+`, .data$S) / length(.data$S)),
          kappa = mean(.data$kappa),
          nu = mean(.data$nu),
          Sigma_exp = list(Reduce(`+`, .data$Sigma_exp) / length(.data$Sigma_exp)),
          Sigma = list(Reduce(`+`, .data$Sigma) / length(.data$Sigma)),
          Sigma_marg = list(Reduce(`+`, .data$Sigma_marg) / length(.data$Sigma_marg)),
          .groups = "drop"
        )

      meta_k <- get_metadata(object)
      meta_k$checkpoint <- checkpoints[k]
      meta_k$label_information <- list(
        cue = cues,
        category = cats,
        group = target_grp
      )

      model_seq[[k]] <- MVBU_StanfitPosterior(
        draws = list(raw = draws_k, summary = summary_k),
        metadata = meta_k
      )
    }

    return(as_model_list(model_seq, model_labels = lbls, metadata = get_metadata(object)))
  }
}

#' @rdname likelihood
S7::method(
  likelihood,
  list(
    S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior),
    S7::class_any,
    S7::class_any
  )
) <- function(
  x,
  new_data,
  categories,
  ...
) {
  post_obj <- as_MVBU_stanfit_posterior(x)
  dots <- list(...)
  group_arg <- dots$group

  if (is.data.frame(new_data) &&
    "group" %in% names(new_data) &&
    is.null(group_arg)) {
    unique_grps <- unique(new_data$group)
    if (length(unique_grps) > 1L) {
      if (missing(categories) || is.null(categories)) {
        categories <- get_category_labels(post_obj)
      }
      res <- matrix(
        NA_real_,
        nrow = nrow(new_data),
        ncol = length(categories),
        dimnames = list(NULL, categories)
      )
      cues <- get_cue_labels(post_obj)
      for (g in unique_grps) {
        idx <- which(new_data$group == g)
        cog_model_g <- .as_cognitive_model_from_posterior(post_obj, group = g)
        sub_data <- as.matrix(new_data[idx, cues, drop = FALSE])
        res[idx, ] <- likelihood(cog_model_g, sub_data, categories, ...)
      }
      return(res)
    }
    group_arg <- unique_grps[1]
  }

  cog_model <- .as_cognitive_model_from_posterior(post_obj, group = group_arg)
  if (missing(categories)) categories <- NULL
  likelihood(cog_model, new_data, categories, ...)
}

#' @rdname posterior
S7::method(
  posterior,
  list(
    S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior),
    S7::class_any,
    S7::class_any
  )
) <- function(
  x,
  new_data,
  categories,
  ...
) {
  post_obj <- as_MVBU_stanfit_posterior(x)
  dots <- list(...)
  group_arg <- dots$group

  if (is.data.frame(new_data) &&
    "group" %in% names(new_data) &&
    is.null(group_arg)) {
    unique_grps <- unique(new_data$group)
    if (length(unique_grps) > 1L) {
      if (missing(categories) || is.null(categories)) {
        categories <- get_category_labels(post_obj)
      }
      res <- matrix(
        NA_real_,
        nrow = nrow(new_data),
        ncol = length(categories),
        dimnames = list(NULL, categories)
      )
      for (g in unique_grps) {
        idx <- which(new_data$group == g)
        cog_model_g <- .as_cognitive_model_from_posterior(post_obj, group = g)
        res[idx, ] <- posterior(
          cog_model_g,
          new_data[idx, , drop = FALSE],
          categories,
          ...
        )
      }
      return(res)
    }
    group_arg <- unique_grps[1]
  }

  cog_model <- .as_cognitive_model_from_posterior(post_obj, group = group_arg)
  if (missing(categories)) categories <- NULL
  posterior(cog_model, new_data, categories, ...)
}

#' @rdname categorize
S7::method(
  categorize,
  list(
    S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior),
    S7::class_any,
    S7::class_any
  )
) <- function(
  x,
  new_data,
  decision_rule,
  simplify = NULL,
  ...
) {
  post_obj <- as_MVBU_stanfit_posterior(x)
  dots <- list(...)
  group_arg <- dots$group

  if (is.data.frame(new_data) &&
    "group" %in% names(new_data) &&
    is.null(group_arg)) {
    unique_grps <- unique(new_data$group)
    if (length(unique_grps) > 1L) {
      if (missing(decision_rule) || is.null(decision_rule)) {
        decision_rule <- "proportional"
      }
      if (is.null(simplify)) {
        simplify <- !identical(decision_rule, "proportional")
      }
      res_list <- vector("list", length(unique_grps))
      for (i in seq_along(unique_grps)) {
        g <- unique_grps[i]
        idx <- which(new_data$group == g)
        cog_model_g <- .as_cognitive_model_from_posterior(post_obj, group = g)
        val <- categorize(
          cog_model_g,
          new_data[idx, , drop = FALSE],
          decision_rule = decision_rule,
          simplify = simplify,
          ...
        )
        res_list[[i]] <- list(idx = idx, val = val)
      }
      if (isTRUE(simplify)) {
        out <- character(nrow(new_data))
        for (item in res_list) {
          out[item$idx] <- item$val
        }
        return(out)
      }
      out <- data.frame(
        response_category = character(nrow(new_data)),
        response_probability = numeric(nrow(new_data)),
        stringsAsFactors = FALSE
      )
      for (item in res_list) {
        out[item$idx, ] <- item$val
      }
      return(out)
    }
    group_arg <- unique_grps[1]
  }

  cog_model <- .as_cognitive_model_from_posterior(post_obj, group = group_arg)
  categorize(cog_model, new_data, decision_rule, simplify = simplify, ...)
}

#' @rdname get_category_likelihood_function
S7::method(
  get_category_likelihood_function,
  S7::new_union(MVBU_Stanfit, MVBU_StanfitPosterior)
) <- function(x, ...) {
  post_obj <- as_MVBU_stanfit_posterior(x)
  dots <- list(...)
  cog_model <- .as_cognitive_model_from_posterior(post_obj, group = dots$group)
  get_category_likelihood_function(cog_model, ...)
}


#' Summarise category parameter draws across MCMC samples
#'
#' Computes the posterior mean of location (mean vector \code{mu.mean})
#' and dispersion (covariance matrix \code{Sigma.mean}) per group and category
#' from a draws data frame (supporting NIW, NIX, and MNIX models).
#'
#' @param draws_df Data frame of draws containing at least \code{group}, \code{category},
#'   and mean draws \code{m}, plus \code{Sigma} or (\code{S} and \code{nu}) or \code{sigma2}.
#' @return A summarised data frame with columns \code{group}, \code{category},
#'   \code{mu.mean}, and \code{Sigma.mean}.
#' @noRd
#' @keywords internal
.summarise_category_parameter_draws <- function(draws_df) {
  d_raw <- draws_df

  # Convert conjugate scatter / df or variance to Sigma if needed
  if (!"Sigma" %in% names(d_raw)) {
    if ("S" %in% names(d_raw) && "nu" %in% names(d_raw)) {
      d_raw$Sigma <- get_expected_Sigma_from_S(d_raw$S, d_raw$nu)
    } else if ("sigma2" %in% names(d_raw)) {
      # For univariate NIX or MNIX draws
      d_raw$Sigma <- lapply(d_raw$sigma2, function(s2) {
        if (is.matrix(s2)) {
          s2
        } else if (length(s2) == 1L) {
          matrix(s2, 1L, 1L)
        } else {
          diag(as.numeric(s2), nrow = length(s2))
        }
      })
    }
  }

  # Split by group and category
  grp_col <- if ("group" %in% names(d_raw)) d_raw$group else rep("all", nrow(d_raw))
  cat_col <- d_raw$category
  split_factor <- interaction(grp_col, cat_col, drop = TRUE, lex.order = TRUE)
  splits <- split(seq_len(nrow(d_raw)), split_factor)

  res_list <- vector("list", length(splits))
  for (i in seq_along(splits)) {
    idx <- splits[[i]]
    m_list <- d_raw$m[idx]
    sig_list <- d_raw$Sigma[idx]
    n_draws <- length(idx)

    # Base R Reduce for summation without purrr
    m_sum <- Reduce(`+`, lapply(m_list, as.numeric))
    mu_mean <- m_sum / n_draws

    sig_sum <- Reduce(`+`, lapply(sig_list, as.matrix))
    sigma_mean <- sig_sum / n_draws

    res_list[[i]] <- data.frame(
      group = grp_col[idx[1L]],
      category = cat_col[idx[1L]],
      stringsAsFactors = FALSE
    )
    res_list[[i]]$mu.mean <- list(mu_mean)
    res_list[[i]]$Sigma.mean <- list(sigma_mean)
  }

  do.call(rbind, res_list)
}
