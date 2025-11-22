#' @import ggplot2 cowplot gganimate transformr av viridis
#' @importFrom ellipse ellipse
#' @importFrom mvtnorm dmvt
#' @importFrom scales trans_new
NULL

#' Checking that input to plotting function is compatible
check_compatibility_between_NIW_belief_and_data <- function(
    x,
    data.exposure,
    data.test,
    facet_rows_by, facet_cols_by, facet_wrap_by, animate_by
) {
  facet_rows_by = enquo(facet_rows_by)
  facet_cols_by = enquo(facet_cols_by)
  facet_wrap_by = enquo(facet_wrap_by)
  animate_by = enquo(animate_by)

  assert_that(is.NIW_belief(x))

  if (!quo_is_null(facet_rows_by)) {
    assert_that(quo_is_null(facet_wrap_by), msg = "Can only specify either facet_wrap_by or facet_rows_by/facet_cols_by.")
    assert_that(all(as_name(facet_rows_by) %in% names(x)),
                msg = paste(as_name(facet_rows_by), "not found in NIW_belief (x)."))
    assert_that(!all(!is.null(data.exposure), as_name(facet_rows_by) %nin% names(data.exposure)),
                msg = "When facet_rows_by is specified, it must be present in the exposure data.")
  }
  if (!quo_is_null(facet_cols_by)) {
    assert_that(all(as_name(facet_cols_by) %in% names(x)),
                msg = paste(as_name(facet_cols_by), "not found in NIW_belief (x)."))
    assert_that(!all(!is.null(data.exposure), as_name(facet_cols_by) %nin% names(data.exposure)),
                msg = "When facet_cols_by is specified, it must be present in the exposure data.")
  }
  if (!quo_is_null(facet_wrap_by)) {
    assert_that(all(as_name(facet_wrap_by) %in% names(x)),
                msg = paste(as_name(facet_wrap_by), "not found in NIW_belief (x)."))
    assert_that(!all(!is.null(data.exposure), as_name(facet_wrap_by) %nin% names(data.exposure)),
                msg = "When facet_wrap_by is specified, it must be present in the exposure data.")
  }
  if (!quo_is_null(animate_by)) {
    assert_that(all(as_name(animate_by) %in% names(x)),
                msg = paste(as_name(animate_by), "not found in NIW_belief (x)."))
    assert_that(!all(!is.null(data.exposure), as_name(animate_by) %nin% names(data.exposure)),
                msg = "When animate_by is specified, it must be present in the exposure data.")
  }

  cue.labels = get_cue_labels_from_model(x)
  assert_that(!all(!is.null(data.exposure), cue.labels %nin% names(data.exposure)),
              msg = "Can't plot exposure data: cue names in exposure data must match those in the NIW belief object.")
  assert_that(!all(!is.null(data.exposure), "category" %nin% names(data.exposure)),
              msg = "Can't plot exposure data: exposure data does not contain column category.")
  assert_that(!all(!is.null(data.test), is.null(data.test)),
              msg = "Can't plot test data: No test data provided.")
  assert_that(!all(!is.null(data.test), cue.labels %nin% names(data.test)),
              msg = "Can't plot test data: cue names in test data must match those in the NIW belief object.")

  return(TRUE)
}

#' Plot expected categories of ideal adaptor stanfit
#'
#' Plot bivariate Gaussian categories expected given the parameters inferred by incremental Bayesian belief-
#' updating (IBBU) for an ideal adaptor. Specifically, the categories are derived by marginalizing over the uncertainty represented
#' by the posterior MCMC samples. Two methods are available (specified by `type`), which differ in their
#' computational demands and speed.
#'
#' @aliases plot_expected_ibbu_stanfit_categories
#'
#' @param model An ideal adaptor stanfit object.
#' @param type Either `"contour"` or `"density"`, specifying the type of plot. Note that the contour plot is *much*
#' faster. It simply gets the expected values of \code{mu} (based on the NIW parameter \code{m}) and \code{Sigma}
#' (based on the NIW parameters \code{S} and \code{nu}) at each MCMC draw, and then averages over
#' all MCMC draws. The plotted categories represent those means of the expected \code{mu} and \code{Sigma}. The
#' density plot instead calculates the posterior predictive for each MCMC draw (i.e, the multivariate Student-T
#' density based on the NIW parameters \code{m, S, kappa, nu}), and then averages those densities. Since this is
#' done for *all* points defined by the data.grid this can be rather computationally expensive and slow.
#' @param categories,groups,cues Character vector of categories, groups, and cues to be plotted.
#' (default: all categories, groups, and cues in the model will be plotted)
#' @param plot.test,plot.exposure Should the test and/or exposure stimuli be plotted? (default: `TRUE` for `plot.test`,
#' `FALSE` for `plot.exposure`) The test items are plotted as black points. The exposure mean is plotted as point,
#' and the .95 interval of cue distributions during exposure are plotted as dashed ellipse in the same color as the
#' expected categories.
#' @param ndraws Number of draws from posterior to use for plot, or `NULL` if all draws are to be returned. (default: `NULL`)
#' @param annotate_inferred_category_means Character vector indicating whether the location and value of the mean be
#' indicated through data rugs (`"rug"`) and/or text labels (`"text"`)? Set to NULL to ignore. (default: `c("rug", "text")`)
#' @param untransform_cues DEPRECATED. Should m_0 and S_0 be transformed back into the original cue space? (default: `FALSE`)
#' @param levels Used only if `type` is `"contour"`. levels The cumulative probability levels that should be plotted (using
#' `geom_polygon()`) around the mean. By default the most transparent ellipse still drawn corresponds to .95.
#' @param xlim,ylim For density plots. Limits for the x- and y-axis.
#' @param resolution For density plots. How many steps along x and y should be calculated? Note that computational
#' complexity increases quadratically with resolution. (default: 25)
#' @param category.colors Vector of colors of same length as categories.
#' @param data.grid.xlim,data.grid.ylim,data.grid.resolution Used only if `type` is `"density"`. Limits for x- and y-axis as
#' well as resolution of the data.grid, defining the range over which the posterior predictive (multivariate Student-T density)
#' is calculated. Note that the number of densities to calculate is a *quadratic* function of `data.grid.resolution`. The default
#' for `data.grid.resolution` is 10, corresponding to 100 densities to be calculated for each MCMC draw.
#'
#'
#' @return ggplot object.
#'
#' @details
#' Typically, the categories, groups, and cues
#' are automatically added to the fit during the creation of the fit. If necessary, however, it is possible to use
#' [tidybayes::recover_types] on the stanfit object to add or change these levels later.
#'
#'
#' @seealso TBD
#' @keywords TBD
#'
#' @importFrom purrr map_dbl
#' @export
plot_expected_categories <- function(model, ...) {
  UseMethod("plot_expected_categories")
}

#' @rdname plot_expected_categories
#' @export
plot_expected_categories_contour <- function(model, ...) {
  UseMethod("plot_expected_categories_contour")
}

#' @rdname plot_expected_categories
#' @export
plot_expected_categories_density <- function(model, ...) {
  UseMethod("plot_expected_categories_density")
}

#' @aliases plot_expected_categories_contour_2D
#' @rdname plot_expected_categories
#' @export
plot_expected_categories_contour2D <- function(model, ...) {
  UseMethod("plot_expected_categories_contour2D")
}

#' @aliases plot_expected_categories_density_1D
#' @rdname plot_expected_categories
#' @export
plot_expected_categories_density1D <- function(model, ...) {
  UseMethod("plot_expected_categories_density1D")
}

#' @aliases plot_expected_categories_density_2D
#' @rdname plot_expected_categories
#' @export
plot_expected_categories_density2D <- function(model, ...) {
  UseMethod("plot_expected_categories_density2D")
}

#' @rdname plot_expected_categories
#' @export
plot_expected_categories.ideal_adaptor_stanfit <- function(
    model,
    type,
    ...
) {
  assert_that(all(type %in% c("contour", "density"), length(type) == 1))
  if (type == "contour")
    plot_expected_categories_contour(model = model, ...)
  else {
    plot_expected_categories_density(model = model, ...)
  }
}

#' @aliases plot_expected_ibbu_stanfit_categories_contour
#' @rdname plot_expected_categories
#' @export
plot_expected_categories_contour.ideal_adaptor_stanfit <- function(
    model,
    cues = get_cue_levels(model),
    ...
) {
  if (length(cues) == 1) {
    warning("Contour plots are only supported when at least 2 cues are selected for plotting.")
  } else if (length(cues) == 2) {
    plot_expected_categories_contour2D(model = model, cues = cues, ...)
  } else {
    warning("Contour plots for more than 2 are not yet supported.")
  }
}

#' @aliases plot_expected_ibbu_stanfit_categories_density
#' @rdname plot_expected_categories
#' @export
plot_expected_categories_density.ideal_adaptor_stanfit <- function(
    model,
    cues = get_cue_levels(model),
    ...
) {
  if (length(cues) == 1) {
    plot_expected_categories_density1D(model = model, cues = cues, ...)
  } else if (length(cues) == 2) {
    plot_expected_categories_density2D(model = model, cues = cues, ...)
  } else {
    warning("Density plots for more than 2 are not yet supported.")
  }
}

#' Plot expected bivariate (2D) category likelihoods
#'
#' Plot bivariate Gaussian categories expected given NIW belief(s). One NIW belief describes the uncertainty about the
#' category statistics of all categories. This includes the m (the mean of category means \eqn{\mu}), S (the scattermatrix),
#' kappa (the strength of the belief in m) and nu (the strength of the belief in S).
#'
#' It is possible to hand more than one NIW belief to this function, and to facet or animate by variables that uniquely
#' identify the different beliefs. For example, one can plot
#' different priors for different talkers (grouping by talker), or different posteriors for different exposure conditions
#' (grouping by exposure condition), the incremental updating of NIW beliefs (grouping by observations), or any combinations
#' of these.
#'
#' @param x An \code{\link{NIW_belief}} or \code{\link{NIW_beliefs}} object.
#' @param levels Levels of the confidence ellipses. (default: .5, .66, .8, .9., and .95)
#' @param data.exposure Optional \code{tibble} or \code{data.frame} that contains exposure data to be plotted. (default: `NULL`)
#' @param data.test Optional \code{tibble} or \code{data.frame} that contains test data to be plotted. (default: `NULL`)
#' @param facet_rows_by,facet_cols_by,facet_wrap_by,animate_by Which group variables, if any, should be used for faceting and/or
#' animation? (defaults: `NULL`)
#' @param animation_follow Should the animation follow the data (zoom in and out)? (default: `FALSE`)
#' @param category.ids Vector of category IDs to be plotted or leave `NULL` to plot all groups. (default: `NULL`). Only relevant
#' if `data.exposure` is provided.
#' @param category.labels Vector of category labels of same length as `category.ids` or `NULL` to use defaults. (default: `NULL`)
#' Only relevant if `data.exposure` is provided.
#' @param category.colors Vector of colors of same length as category.ids or `NULL` to use defaults. (default: `NULL`)
#' Only relevant if `data.exposure` is provided.
#' @param ... additional arguments handed to geom_polygon.
#'
#' @return ggplot object.
#'
#' @seealso TBD
#' @keywords TBD
#' @rdname plot_expected_categories
#' @export
plot_expected_categories_contour2D.NIW_ideal_adaptor <- function(
    x,
    levels = c(1/2, 2/3, 4/5, 9/10, 19/20),
    data.exposure = NULL,
    data.test = NULL,
    facet_rows_by = NULL, facet_cols_by = NULL, facet_wrap_by = NULL, animate_by = NULL, animation_follow = F,
    category.ids = NULL, category.labels = NULL, category.colors = NULL,
    ...
) {
  facet_rows_by = enquo(facet_rows_by)
  facet_cols_by = enquo(facet_cols_by)
  facet_wrap_by = enquo(facet_wrap_by)
  animate_by = enquo(animate_by)
  check_compatibility_between_NIW_belief_and_data(x, data.exposure, data.test,
                                                  !! facet_rows_by, !! facet_cols_by, !! facet_wrap_by, !! animate_by)
  # Remember groups
  cue.labels = get_cue_labels_from_model(x)
  assert_that(length(cue.labels) == 2, msg = "Expecting exactly two cues for plotting.")

  # Setting aes defaults
  if (is.null(category.ids)) category.ids = levels(x$category)
  if (is.null(category.labels)) category.labels = levels(x$category)
  if (is.null(category.colors)) category.colors = get_default_colors("category", category.ids)

  suppressMessages(
    x %<>%
      mutate(Sigma = get_expected_Sigma_from_S(S, nu)) %>%
      crossing(level = levels) %>%
      mutate(ellipse = pmap(.l = list(Sigma, m, level), ellipse.pmap)) %>%
      # This step is necessary since unnest() can't yet unnest lists of matrices
      # (bug was reported and added as milestone, 11/2019)
      mutate(ellipse = map(ellipse, ~ as_tibble(.x, .name_repair = "unique"))) %>%
      select(-c(kappa, nu, m, S, Sigma, lapse_rate)) %>%
      unnest(ellipse) %>%
      # Get group structure again, as crossing apparently removes it
      group_by(!!! syms(group_vars(x)))
  )

  p <- ggplot(x,
             aes(
               x = .data[[cue.labels[1]]],
               y = .data[[cue.labels[2]]],
               fill = .data$category)) +
    geom_polygon(aes(alpha = 1 - .data$level,
                     group = interaction(
                       .data$category,
                       .data$level,
                       !!! syms(group_vars(x)))),
                 ...) +
    { if (!is.null(data.test))
      add_test_locations_to_2D_plot(data = data.test, cue.labels = cue.labels) } +
    { if (!is.null(data.exposure))
      add_exposure_locations_to_2D_plot(data = data.exposure, cue.labels = cue.labels,
                                        category.ids = category.ids, category.labels = category.labels, category.colors) } +
    scale_x_continuous(cue.labels[1]) +
    scale_y_continuous(cue.labels[2]) +
    scale_fill_manual("Category",
                      breaks = category.ids,
                      labels = category.labels,
                      values = category.colors) +
    scale_alpha_continuous("Cumulative probability",
                           range = c(0.05, .5),
                           breaks = 1 - levels,
                           labels = round(levels, 2))

  p = facet_or_animate(p, !!facet_rows_by, !!facet_cols_by, !! facet_wrap_by, !!animate_by, animation_follow)
  return(p)
}


# this is a temporary hack to deal with the fact that grouped ideal adaptors are just tibbles atm
#' @rdname plot_expected_categories
#' @export
plot_expected_categories_contour2D.tbl_df <- function(...) {
  plot_expected_categories_contour2D.NIW_ideal_adaptor(...)
}


#' @aliases plot_expected_ibbu_stanfit_categories_contour2D
#' @rdname plot_expected_categories
#' @export
plot_expected_categories_contour2D.ideal_adaptor_stanfit <- function(
    model,
    categories = get_category_levels(model),
    groups = get_group_levels(model, include_prior = T),
    cues = get_cue_levels(model),
    plot.test = T,
    plot.exposure = F,
    annotate_inferred_category_means = c("rug", "text"),
    untransform_cues = FALSE,
    levels = plogis(seq(-15, qlogis(.95), length.out = 20)),
    category.colors = get_default_colors("category", categories)
) {
  assert_that(all(annotate_inferred_category_means %in% c("rug", "text")))
  cues <- unique(cues)
  assert_that(length(cues) == 2)

  d.pars <-
    get_expected_category_statistic(
      model,
      categories = categories,
      groups = groups,
      untransform_cues = untransform_cues)

  d.pars %<>%
    rename(x = Sigma.mean, centre = mu.mean) %>%
    crossing(level = levels) %>%
    mutate(ellipse = pmap(.l = list(x, centre, level), ellipse.pmap)) %>%
    unnest(ellipse) %>%
    group_by(across(-ellipse)) %>%
    transmute(cue1 = ellipse[,1], cue2 = ellipse[,2])

  xlim <- range(d.pars$cue1)
  ylim <- range(d.pars$cue2)
  groups_found <- levels(d.pars$group)
  d.pars %>%
    ggplot(
      aes(x = .data$cue1,
          y = .data$cue2,
          fill = .data$category,
          alpha = 1 - .data$level,
          group = paste(.data$category, .data$level))) +
    geom_polygon() +
    { if ("rug" %in% annotate_inferred_category_means)
      geom_rug(
        data = . %>%
          distinct(group, category, centre),
        aes(
          x = map_dbl(centre, ~ .x[1]),
          y = map_dbl(centre, ~ .x[2]),
          color = category),
        inherit.aes = F) } +
    { if ("text" %in% annotate_inferred_category_means)
      list(
        geom_text(
          data = . %>%
            ungroup() %>%
            distinct(group, category, centre),
          aes(
            x = map_dbl(centre, ~ .x[1]),
            label = map(centre, ~ paste(signif(.x[1], 2))),
            color = category),
          y = min(ylim),
          angle = 90,
          hjust = 0,
          inherit.aes = F),
        geom_text(
          data = . %>%
            ungroup() %>%
            distinct(group, category, centre),
          aes(
            y = map_dbl(centre, ~ .x[2]),
            label= map(centre, ~ paste(signif(.x[2], 2))),
            color = category),
          x = min(xlim),
          angle = 0,
          hjust = 0,
          inherit.aes = F)) } +
    scale_fill_manual(
      "Category",
      breaks = categories,
      values = category.colors,
      aesthetics = c("color", "fill")) +
    new_scale_color() +
    # Optionally plot test data
    { if (plot.test)
      add_test_data_to_2D_plot(
        get_test_data(model, groups = setdiff(groups_found, "prior"), .from_staninput = TRUE) %>%
          { if (untransform_cues) get_untransform_function(model)(.) else . } %>%
          ungroup() %>%
          distinct(group, !!! syms(cues)),
        cue.labels = cues) } +
    # Optionally plot exposure data
    { if (plot.exposure)
      add_exposure_summary_to_2D_plot(
        get_exposure_category_statistic.ideal_adaptor_stanfit(
          model,
          categories = levels(d.pars$category),
          groups = setdiff(groups_found, "prior"),
          untransform_cues = untransform_cues)) } +
    scale_x_continuous(cues[1]) +
    scale_y_continuous(cues[2]) +
    scale_color_manual(
      "Category",
      breaks = categories,
      values = lighten(category.colors, amount = .5)) +
    scale_alpha("", range = c(0.1,.9)) +
    facet_wrap(~ group)
}


#' Plot expected univariate (1D) category likelihoods
#'
#' Plot univariate Gaussian categories expected given NIW belief(s). One NIW belief describes the uncertainty about the
#' category statistics of all categories. This includes the m (the mean of category means \eqn{\mu}), S (the scattermatrix),
#' kappa (the strength of the belief in m) and nu (the strength of the belief in S). For the univariate case, m and S are
#' scalars \insertCite{@see @murphy2012 p. 136}{MVBeliefUpdatr}.
#'
#' It is possible to hand more than one NIW belief to this function, and to facet or animate by variables that uniquely
#' identify the different beliefs. For example, one can plot
#' different priors for different talkers (grouping by talker), or different posteriors for different exposure conditions
#' (grouping by exposure condition), the incremental updating of NIW beliefs (grouping by observations), or any combinations
#' of these.
#'
#' @param x An \code{\link{NIW_belief}} or \code{\link{NIW_beliefs}} object.
#' @param levels Levels of the confidence ellipses. (default: .5, .66, .8, .9., and .95)
#' @param data.exposure Optional \code{tibble} or \code{data.frame} that contains exposure data to be plotted. (default: `NULL`)
#' @param data.test Optional \code{tibble} or \code{data.frame} that contains test data to be plotted. (default: `NULL`)
#' @param facet_rows_by,facet_cols_by,facet_wrap_by,animate_by Which group variables, if any, should be used for faceting and/or
#' animation? (defaults: `NULL`)
#' @param animation_follow Should the animation follow the data (zoom in and out)? (default: `FALSE`)
#' @param xlim,ylim Limits for the x- and y-axis.
#' @param category.ids Vector of category IDs to be plotted or leave `NULL` to plot all groups. (default: `NULL`)
#' @param category.labels Vector of group labels of same length as `category.ids` or `NULL` to use defaults. (default: `NULL`)
#' @param category.colors Vector of colors of same length as category.ids or `NULL` to use defaults. (default: `NULL`)
#' @param category.linetypes Vector of linetypes of same length as category.ids or `NULL` to use defaults. (default: `NULL`)
#' Currently being ignored.
#' @param ... additional arguments to geom_line.
#'
#' @return ggplot object.
#'
#' @seealso TBD
#' @keywords TBD
#' @references \insertRef{murphy2012}{MVBeliefUpdatr}
#' @rdname plot_expected_categories
#' @export
plot_expected_categories_density1D.NIW_ideal_adaptor <- function(
    x,
    data.exposure = NULL,
    data.test = NULL,
    facet_rows_by = NULL, facet_cols_by = NULL, facet_wrap_by = NULL, animate_by = NULL, animation_follow = F,
    xlim, ylim = NULL, x.expand = c(0, 0),
    category.ids = NULL, category.labels = NULL, category.colors = NULL, category.linetypes = NULL,
    ...
) {
  facet_rows_by <- enquo(facet_rows_by)
  facet_cols_by <- enquo(facet_cols_by)
  facet_wrap_by <- enquo(facet_wrap_by)
  animate_by <- enquo(animate_by)
  check_compatibility_between_NIW_belief_and_data(x, data.exposure, data.test,
                                                  !! facet_rows_by, !! facet_cols_by, !! facet_wrap_by, !! animate_by)
  # Remember groups
  cue.labels <- get_cue_labels_from_model(x)
  assert_that(length(cue.labels) == 1, msg = "Expecting exactly one cue for plotting.")

  if (is_missing(xlim)) {
    if (!is.null(data.exposure) & !is.null(data.test))
      xlim <- range(range(data.exposure[[cue.labels[1]]]), range(data.test[[cue.labels[1]]])) else
        if (!is.null(data.exposure))
          xlim <- range(data.exposure[[cue.labels[1]]]) else
            if (!is.null(data.test))
              xlim <- range(data.test[[cue.labels[1]]])
  }
  assert_that(!is_missing(xlim), msg = "`xlim` must be specified")

  # Setting aes defaults
  if (is.null(category.ids)) category.ids <- levels(x$category)
  if (is.null(category.labels)) category.labels <- levels(x$category)
  if (is.null(category.colors)) category.colors <- get_default_colors("category", category.ids)
  if (is.null(category.linetypes)) category.linetypes <- rep(1, length(category.ids))

  if (any(!quo_is_null(facet_rows_by),
          !quo_is_null(facet_cols_by),
          !quo_is_null(animate_by))) x %<>% group_by(!! facet_rows_by, !! facet_cols_by, !! facet_wrap_by, !! animate_by,
                                                     .add = TRUE)

  stat_functions <-
    x %>%
    mutate(
      mu = get_expected_mu_from_m(m),
      Sigma = get_expected_Sigma_from_S(S, nu)) %>%
    group_by(category, .add = T) %>%
    group_map(
      .keep = T,
      .f = function(.x)
        stat_function(
          data = .x,
          mapping = aes(color = category),
          fun = dnorm,
          args = list(mean = .x$mu, sd = .x$Sigma^.5),
          ...))

  p <- ggplot(mapping = aes(color = category)) +
    stat_functions +
    { if (!is.null(data.test))
      add_test_locations_to_1D_plot(data = data.test, cue.labels = cue.labels) } +
    { if (!is.null(data.exposure))
      add_exposure_locations_to_1D_plot(data = data.exposure, cue.labels = cue.labels,
                                        category.ids = category.ids, category.labels = category.labels, category.colors) } +
    scale_x_continuous(cue.labels, limits = xlim, expand = x.expand) +
    scale_y_continuous("Density", limits = ylim) +
    scale_color_manual("Category",
                       breaks = category.ids,
                       labels = category.labels,
                       values = category.colors)

  p <- facet_or_animate(p, !!facet_rows_by, !!facet_cols_by, !! facet_wrap_by, !!animate_by, animation_follow)
  return(p)
}

# this is a temporary hack to deal with the fact that grouped ideal adaptors are just tibbles atm
#' @rdname plot_expected_categories
#' @export
plot_expected_categories_density1D.tbl_df <- function(...) {
  plot_expected_categories_density1D.NIW_ideal_adaptor(...)
}


#' @aliases plot_expected_ibbu_stanfit_categories_density1D
#' @rdname plot_expected_categories
#' @export
plot_expected_categories_density1D.ideal_adaptor_stanfit <- function(
    model,
    categories = get_category_levels(model),
    groups = get_group_levels(model, include_prior = T),
    cues = get_cue_levels(model),
    ndraws = NULL,
    plot.test = T,
    plot.exposure = F,
    annotate_inferred_category_means = c("rug", "text"),
    untransform_cues = FALSE,
    category.colors = get_default_colors("category", categories),
    xlim = NULL, resolution = 101
) {
  assert_that(all(annotate_inferred_category_means %in% c("rug", "text")))
  cues <- unique(cues)
  assert_that(length(cues) == 1)

  d.pars <-
    get_draws(
      model,
      categories = categories,
      groups = groups,
      ndraws = ndraws,
      wide = F,
      nest = T,
      untransform_cues = untransform_cues)

  groups_found <- levels(d.pars$group)
  # By default get plot dimensions that are centered around the test data
  xlim <-
    if (is.null(xlim)) {
      get_test_data(model, groups = setdiff(groups_found, "prior"), .from_staninput = TRUE) %>%
        pull(cues[1]) %>%
        { (range(.) - mean(.)) * 1.5 + mean(.) }
    } else xlim

  d.pars %<>%
    crossing(
      cue1 = seq(min(xlim), max(xlim), length.out = resolution)) %>%
    mutate(x = map(cue1, ~ c(.x))) %>%
    mutate(density = pmap_dbl(.l = list(x, m, S, kappa, nu), get_NIW_posterior_predictive.pmap, log = F)) %>%
    # Marginalize over MCMC draws
    group_by(group, category, cue1) %>%
    summarise(density = mean(density))

  d.pars %>%
    ggplot(
      aes(x = .data$cue1,
          color = .data$category,
          fill = .data$category,
          y = .data$density)) +
    geom_line(aes(group = interaction(.data$category, .data$group))) +
    { if ("rug" %in% annotate_inferred_category_means)
      geom_rug(
        data = . %>%
          group_by(group, category) %>%
          summarise(across(c(cue1), mean)),
        aes(x = .data$cue1,
            color = .data$category),
        sides = "b",
        inherit.aes = F) } +
    { if ("text" %in% annotate_inferred_category_means)
      list(
        geom_text(
          data = . %>%
            group_by(group, category) %>%
            summarise(across(c(cue1), mean)),
          aes(
            x = .data$cue1,
            label= signif(.data$cue1, 2),
            color = .data$category),
          y = 0,
          angle = 90,
          hjust = 0,
          inherit.aes = F)) } +
    scale_fill_manual(
      "Category",
      breaks = categories,
      values = category.colors,
      aesthetics = c("color", "fill")) +
    new_scale_color() +
    # Optionally plot test data
    { if (plot.test)
      add_test_data_to_1D_plot(
        get_test_data(model, groups = setdiff(groups_found, "prior"), .from_staninput = TRUE) %>%
          { if (untransform_cues) get_untransform_function(model)(.) else . } %>%
          ungroup() %>%
          distinct(group, !!! syms(cues)),
        cue.labels = cues) } +
    # Optionally plot exposure data
    { if (plot.exposure)
      add_exposure_summary_to_1D_plot(
        get_exposure_category_statistic.ideal_adaptor_stanfit(
          model,
          categories = levels(d.pars$category),
          groups = setdiff(groups_found, "prior"),
          untransform_cues = untransform_cues)) } +
    scale_x_continuous(cues[1]) +
    scale_color_manual(
      "Category",
      breaks = categories,
      values = lighten(category.colors, amount = .5)) +
    coord_fixed(xlim = xlim) +
    facet_wrap(~ .data$group)
}


#' @aliases plot_expected_ibbu_stanfit_categories_density2D
#' @rdname plot_expected_categories
#' @export
plot_expected_categories_density2D.ideal_adaptor_stanfit <- function(
    model,
    categories = get_category_levels(model),
    groups = get_group_levels(model, include_prior = T),
    cues = get_cue_levels(model),
    ndraws = NULL,
    plot.test = T,
    plot.exposure = F,
    annotate_inferred_category_means = c("rug", "text"),
    untransform_cues = FALSE,
    category.colors = get_default_colors("category", categories),
    xlim = NULL, ylim = NULL, resolution = 25
) {
  assert_that(all(annotate_inferred_category_means %in% c("rug", "text")))
  cues <- unique(cues)
  assert_that(length(cues) == 2)

  d.pars <-
    get_draws(
      model,
      categories = categories,
      groups = groups,
      ndraws = ndraws,
      wide = F,
      nest = T,
      untransform_cues = untransform_cues)

  groups_found <- levels(d.pars$group)
  # By default get plot dimensions that are centered around the test data
  xlim <-
    if (is.null(xlim)) {
      get_test_data(model, groups = setdiff(groups_found, "prior"), .from_staninput = TRUE) %>%
        pull(cues[1]) %>%
        { (range(.) - mean(.)) * 1.5 + mean(.) }
    } else xlim
  ylim <-
    if (is.null(ylim)) {
      get_test_data(model, groups = setdiff(groups_found, "prior"), .from_staninput = TRUE) %>%
        pull(cues[2]) %>%
        { (range(.) - mean(.)) * 1.5 + mean(.) }
    } else ylim

  d.pars %<>%
    crossing(
      cue1 = seq(min(xlim), max(xlim), length.out = resolution),
      cue2 = seq(min(ylim), max(ylim), length.out = resolution)) %>%
    mutate(x = map2(cue1, cue2, ~ c(.x, .y))) %>%
    mutate(
      density = pmap_dbl(.l = list(x, m, S, kappa, nu), get_NIW_posterior_predictive.pmap)) %>%
    # Marginalize over MCMC draws
    group_by(group, category, cue1, cue2) %>%
    summarise(density = mean(density))

  d.pars %>%
    ggplot(
      aes(x = .data$cue1,
          y = .data$cue2,
          color = .data$category,
          fill = .data$category,
          z = .data$density)) +
    geom_contour() +
    { if ("rug" %in% annotate_inferred_category_means)
      geom_rug(
        data = . %>%
          group_by(group, category) %>%
          summarise(across(c(cue1, cue2), mean)),
        aes(x = .data$cue1,
            y = .data$cue2,
            color = .data$category),
        inherit.aes = F) } +
    { if ("text" %in% annotate_inferred_category_means)
      list(
        geom_text(
          data = . %>%
            group_by(group, category) %>%
            summarise(across(c(cue1, cue2), mean)),
          aes(
            x = .data$cue1,
            label= signif(.data$cue1, 2),
            color = .data$category),
          y = min(ylim),
          angle = 90,
          hjust = 0,
          inherit.aes = F),
        geom_text(
          data = . %>%
            group_by(group, category) %>%
            summarise(across(c(cue1, cue2), mean)),
          aes(
            y = .data$cue2,
            label= signif(.data$cue2, 2),
            color = .data$category),
          x = min(xlim),
          angle = 0,
          hjust = 0,
          inherit.aes = F)) } +
    scale_fill_manual(
      "Category",
      breaks = categories,
      values = category.colors,
      aesthetics = c("color", "fill")) +
    new_scale_color() +
    # Optionally plot test data
    { if (plot.test)
      add_test_data_to_2D_plot(
        get_test_data(model, groups = setdiff(groups_found, "prior"), .from_staninput = TRUE) %>%
          { if (untransform_cues) get_untransform_function(model)(.) else . } %>%
          ungroup() %>%
          distinct(group, !!! syms(cues)),
        cue.labels = cues) } +
    # Optionally plot exposure data
    { if (plot.exposure)
      add_exposure_summary_to_2D_plot(
        get_exposure_category_statistic.ideal_adaptor_stanfit(
          model,
          categories = levels(d.pars$category),
          groups = setdiff(groups_found, "prior"),
          untransform_cues = untransform_cues)) } +
    scale_x_continuous(cues[1]) +
    scale_y_continuous(cues[2]) +
    scale_color_manual(
      "Category",
      breaks = categories,
      values = lighten(category.colors, amount = .5)) +
    coord_fixed(xlim = xlim, ylim = ylim, ratio = 1) +
    facet_wrap(~ .data$group)
}
