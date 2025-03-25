# Additional plotting functions in plot-expected-categories.R

#' Plot expected categorization function for univariate (1D) categories.
#'
#' Plot categorization function for univariate Gaussian categories expected given NIW parameters.
#'
#' @param target_category The index of the category for which categorization should be shown. (default: `1`)
#' @param xlim,ylim Limits for the x- and y-axis.
#' @param logit Should the categorization function be plotted in logit (`TRUE`) or probabilities (`FALSE`)?
#' (default: `FALSE`)
#' @inheritParams plot_expected_categories_density1D
#'
#' @return ggplot object.
#'
#' @seealso TBD
#' @keywords TBD
#' @rdname plot_expected_categorization_function_1D
#' @export
#'
plot_expected_categorization_function_1D <- function(
  x,
  data.exposure = NULL,
  data.test = NULL,
  target_category = 1,
  logit = F,
  xlim, ylim = NULL, x.expand = c(0, 0),
  facet_rows_by = NULL, facet_cols_by = NULL, facet_wrap_by = NULL, animate_by = NULL, animation_follow = F,
  category.ids = NULL, category.labels = NULL, category.colors = NULL, category.linetypes = NULL,
  ...
) {
  message("TO DO: implement noise and lapse rate handling. (already implemented in underling functions. just not integrated into plotting).")

  facet_rows_by <- enquo(facet_rows_by)
  facet_cols_by <- enquo(facet_cols_by)
  facet_wrap_by <- enquo(facet_wrap_by)
  animate_by <- enquo(animate_by)
  check_compatibility_between_NIW_belief_and_data(x, data.exposure, data.test,
                                                  !! facet_rows_by, !! facet_cols_by, !! facet_wrap_by, !! animate_by)
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
    group_map(
      .keep = T,
      .f = function(.x, .y) {
        cat_function <- get_categorization_function_from_NIW_ideal_adaptor(.x)
        stat_function(
          data = .x,
          fun = cat_function,
          args = list(target_category = target_category, logit = logit), ...) })

  p <-
    ggplot() +
    stat_functions +
    { if (!is.null(data.test))
      add_test_data_to_1D_plot(data = data.test, cue.labels = cue.labels) } +
    { if (!is.null(data.exposure))
      add_exposure_data_to_1D_plot(data = data.exposure, cue.labels = cue.labels,
                                   category.ids = category.ids, category.labels = category.labels, category.colors) } +
    scale_x_continuous(name = cue.labels, limits = xlim, expand = x.expand) +
    scale_y_continuous(name = if (logit)
      paste0("log-odds(resp = ", category.labels[target_category], ")") else
        paste0("p(resp = ", category.labels[target_category], ")")) +
    coord_cartesian(ylim = ylim)

  p <- facet_or_animate(p, !!facet_rows_by, !!facet_cols_by, !! facet_wrap_by, !!animate_by, animation_follow)
  return(p)
}

#' Plot expected categorization function for bivariate (2D) categories.
#'
#' Plot categorization function for bivariate Gaussian categories expected given NIW parameters.
#'
#' @param target_category The index of the category for which categorization should be shown. (default: `1`)
#' @param xlim,ylim Limits for the x- and y-axis.
#' @param resolution How many steps along x and y should be calculated? Note that computational
#' complexity increases quadratically with resolution. (default: 25)
#' @param logit Should the categorization function be plotted in logit (`TRUE`) or probabilities (`FALSE`)?
#' (default: `FALSE`)
#' @inheritParams plot_expected_categories_contour2D
#'
#' @return ggplot object.
#'
#' @seealso TBD
#' @keywords TBD
#' @rdname plot_expected_categorization_function_2D
#' @export
#'
plot_expected_categorization_function_2D <- function(
  x,
  data.exposure = NULL,
  data.test = NULL,
  target_category = 1,
  logit = F,
  xlim, ylim, resolution = 25,
  facet_rows_by = NULL, facet_cols_by = NULL, facet_wrap_by = NULL, animate_by = NULL, animation_follow = F,
  category.ids = NULL, category.labels = NULL, category.colors = NULL,
  ...
) {
  message("TO DO: implement noise and lapse rate handling. (already implemented in underling functions. just not integrated into plotting).")
  facet_rows_by <- enquo(facet_rows_by)
  facet_cols_by <- enquo(facet_cols_by)
  facet_wrap_by <- enquo(facet_wrap_by)
  animate_by <- enquo(animate_by)
  check_compatibility_between_NIW_belief_and_data(x, data.exposure, data.test,
                                                  !! facet_rows_by, !! facet_cols_by, !! facet_wrap_by, !! animate_by)
  cue.labels <- get_cue_labels_from_model(x)
  assert_that(length(cue.labels) == 2, msg = "Expecting exactly two cues for plotting.")
  if (is_missing(xlim)) {
    if (!is.null(data.exposure) & !is.null(data.test))
      xlim <- range(range(data.exposure[[cue.labels[1]]]), range(data.test[[cue.labels[1]]])) else
        if (!is.null(data.exposure))
          xlim <- range(data.exposure[[cue.labels[1]]]) else
            if (!is.null(data.test))
              xlim <- range(data.test[[cue.labels[1]]])
  }
  if (is_missing(ylim)) {
    if (!is.null(data.exposure) & !is.null(data.test))
      ylim <- range(range(data.exposure[[cue.labels[2]]]), range(data.test[[cue.labels[2]]])) else
        if (!is.null(data.exposure))
          ylim <- range(data.exposure[[cue.labels[2]]]) else
            if (!is.null(data.test))
              ylim <- range(data.test[[cue.labels[2]]])
  }
  assert_that(!is_missing(xlim), msg = "`xlim` must be specified")
  assert_that(!is_missing(ylim), msg = "`ylim` must be specified")

  # Setting aes defaults
  if (is.null(category.ids)) category.ids <- levels(x$category)
  if (is.null(category.labels)) category.labels <- levels(x$category)
  if (is.null(category.colors)) category.colors <- get_default_colors("category", category.ids)

  if (any(!quo_is_null(facet_rows_by),
          !quo_is_null(facet_cols_by),
          !quo_is_null(animate_by))) x %<>% group_by(!! facet_rows_by, !! facet_cols_by, !! animate_by,
                                                     .add = TRUE)

  d <- crossing(
    !! sym(cue.labels[1]) := seq(min(xlim), max(xlim), length.out = resolution),
    !! sym(cue.labels[2]) := seq(min(ylim), max(ylim), length.out = resolution))

  x %<>%
    nest() %>%
    mutate(f = map(data, get_categorization_function_from_NIW_ideal_adaptor)) %>%
    # Join in vectored cues
    cross_join(
      d %>%
        transmute(x = pmap(.l = list(!!! syms(cue.labels)), .f = ~ c(...))) %>%
        nest(cues = everything())) %>%
    mutate(
      p_cat = invoke_map(.f = f, .x = cues, target_category = target_category, logit = logit),
      cues = NULL,
      f = NULL) %>%
    # Join separate cues back in
    cross_join(d %>% nest(cues = everything())) %>%
    unnest(c(cues, p_cat))

  p <- ggplot(x,
             mapping = aes(
               x = .data[[cue.labels[1]]],
               y = .data[[cue.labels[2]]])) +
    geom_raster(mapping = aes(fill = if (logit) qlogis(.data$p_cat) else .data$p_cat), ...) +
    # geom_contour(
    #   mapping = aes(z = if (logit) qlogis(.data$p_cat) else .data$p_cat)) +
    { if (!is.null(data.test))
      add_test_locations_to_2D_plot(data = data.test, cue.labels = cue.labels) } +
    { if (!is.null(data.exposure))
      add_exposure_locations_to_2D_plot(data = data.exposure, cue.labels = cue.labels,
                                   category.ids = category.ids, category.labels = category.labels, category.colors) } +
    scale_x_continuous(cue.labels[1]) +
    scale_y_continuous(cue.labels[2]) +
    # For now think about two colors and categories
    scale_fill_gradient2(paste0("p(resp = ", category.labels[target_category], ")"),
                         low = category.colors[1],
                         mid = "white",
                         high = category.colors[2],
                         midpoint = if (logit) 0 else .5) +
    coord_cartesian(xlim = xlim, ylim = ylim)

  p <- facet_or_animate(p, !!facet_rows_by, !!facet_cols_by, !! facet_wrap_by, !!animate_by, animation_follow)
  return(p)
}
