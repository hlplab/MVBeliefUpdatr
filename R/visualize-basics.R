#' @import ggplot2
#' @importFrom ellipse ellipse
#' @importFrom tidybayes mean_hdi
#' @importFrom scales trans_new
NULL


#' Symmetric log transform and function.
#'
#' Makes it possible to use log-stepped scales or coordinate systems even when negative values
#' are included in the data. E.g., in `coord_trans(x = "symlog")`. `symlog`` applies a modified
#' logarithm scale to the specified or current axes that handles negative values while maintaining
#' continuity across zero:
#'
#' y = sign(x) * log10(1 + abs(x) / 10^C )
#'
#' where the scaling constant C determines the resolution of the data around zero. The smallest
#' order of magnitude shown on either side of zero will be 10^ceil(C). If applies as a transform
#' for a ggplot2 coordinate system, C is taken to be 0.
#'
#' Taken from https://www.mathworks.com/matlabcentral/fileexchange/57902-symlog.
#'
#' @references Webber (2012). Measurement Science and Technology .
#' @examples
#' TBD
#' @rdname symlog
#' @export
symlog = function(x, C = 0) sign(x) * log10(1 + abs(x) / 10^C)

#' @rdname symlog
#' @export
inv_symlog = function(x, C = 0) sign(x) * (10^abs(x) * 10^C - 10^C)

#' @rdname symlog
#' @export
symlog_trans = function(){
  scales::trans_new("symlog",
                    transform = function(x) sign(x) * log10(1 + abs(x)),
                    inverse = function(x) sign(x) * (10^abs(x) - 1))
}


#' Get default colors
#'
#' @param var Variable for which default color is requested. Can be either `"category"` or
#' `"group"`.
#'
#' @examples
#' TBD
#' @export
get_default_colors = function(var, n) {
  assert_that(all(var %in% c("category", "group")))
  assert_that(is.count(n))

  if (var == "category")
    c = scales::hue_pal()(n)
  else
    c = scales::brewer_pal(palette = "Set1")(n)

  return(c)
}


#' Get plot limits.
#'
#' Get x and y limits from a ggplot.
#'
#' @param plot A `ggplot` object.
#'
#' @return List with elements x and y, each of which is a vector of two values..
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#'
#' @export
get_plot_limits = function(plot) {
  list(x = ggplot_build(plot)$layout$panel_scales_x[[1]]$range$range,
       y = ggplot_build(obj)$layout$panel_scales_y[[1]]$range$range)
}




#' Get suitable limits for coordinate system based on the MCMC samples of a variable.
#'
#' Useful for, for example, plotting of distribution.
#'
#' @param data A `tibble` or `data.frame` that contains a `measure`.
#' @param measure Name of variable in `data` for which limits are sought.
#' @param hdi.prob Proportion of MCMC samples that are within the limits.
#' @param min,max If min or max are specified, then those limits are returned instead of the HDI-based limits.
#'
#' @return Vector with two values.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#'
#' @export
get_limits = function(data, measure, hdi.prob = .99, min = NULL, max = NULL) {
  data %>%
    mean_hdi((!! rlang::sym(measure)), .width = hdi.prob) %>%
    ungroup() %>%
    summarise(.lower = if (!is.null(min)) min else min(.lower),
              .upper = if (!is.null(max)) max else max(.upper)) %>%
    as.numeric()
}

ellipse.pmap = function(x, centre, level, ...)
  ellipse(x = x, centre = centre, level = level)

add_exposure_data_to_1D_plot = function(
  data,
  cue.labels,
  category.ids,
  category.labels,
  category.colors
) {
  cue.labels[2] = "cue2"
  data %<>% mutate(cue2 = 0)
  add_exposure_data_to_2D_plot(data, cue.labels, category.ids, category.labels, category.colors)
}

add_test_data_to_1D_plot = function(data, cue.labels) {
  cue.labels[2] = "cue2"
  data %<>% mutate(cue2 = 0)
  add_test_data_to_2D_plot(data, cue.labels)
}

add_exposure_data_to_2D_plot = function(
  data,
  cue.labels,
  category.ids,
  category.labels,
  category.colors
) {
  list(
    geom_point(
      data = data,
      mapping = aes(
        x = .data[[cue.labels[1]]],
        y = .data[[cue.labels[2]]],
        shape = .data$category,
        color = .data$category),
      size = 3, alpha = .9),
    scale_shape("Category",
                breaks = category.ids,
                labels = category.labels),
    scale_color_manual("Category",
                       breaks = category.ids,
                       labels = category.labels,
                       values = category.colors))
}

add_test_data_to_2D_plot = function(data, cue.labels) {
  list(
    geom_point(
      data = data,
      mapping = aes(
        x = .data[[cue.labels[1]]],
        y = .data[[cue.labels[2]]]),
      inherit.aes = F,
      color = "black", size = 1))
}

facet_or_animate = function(p, facet_rows_by, facet_cols_by, animate_by, animation_follow) {
  facet_rows_by = enquo(facet_rows_by)
  facet_cols_by = enquo(facet_cols_by)
  animate_by = enquo(animate_by)

  if (!quo_is_null(facet_rows_by) | !quo_is_null(facet_cols_by))
    p = p + facet_grid(
      rows = vars(!! facet_rows_by),
      cols = vars(!! facet_cols_by),
      labeller = label_both)
  if (!quo_is_null(animate_by)) {
    message("Preparing for rendering. This might take a moment.\n")
    p = p +
      labs(title = paste0(as_name(animate_by), ": {closest_state}")) +
      transition_states(!! animate_by,
                        transition_length = 1,
                        state_length = 1) +
      { if (animation_follow) view_follow() } +
      enter_fade() +
      exit_fade()
  }

  return(p)
}

