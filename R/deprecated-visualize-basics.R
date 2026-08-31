#' @include internal-asserts.R
#' @import ggplot2
#' @importFrom ellipse ellipse
#' @importFrom lifecycle deprecate_warn
#' @importFrom purrr map_dbl pmap
#' @importFrom scales trans_new
#' @importFrom tidybayes mean_hdi
NULL

# -----------------------------------------------------------------------------
# deprecated
# -----------------------------------------------------------------------------

#' Deprecated: symlog
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `symlog()` is deprecated. Use [scales::pseudo_log_trans()] instead.
#'
#' @param x Numeric vector to transform.
#' @param C Scaling constant determining resolution around zero. Defaults to 0.
#' @return Transformed numeric vector.
#' @rdname deprecated-functions
#' @export
symlog <- function(x, C = 0) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "symlog()",
    details = "Use scales::pseudo_log_trans() instead."
  )
  sign(x) * log10(1 + abs(x) / 10^C)
}

#' Deprecated: inv_symlog
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `inv_symlog()` is deprecated.
#'
#' @param x Numeric vector to invert.
#' @param C Scaling constant. Defaults to 0.
#' @return Inverse transformed numeric vector.
#' @rdname deprecated-functions
#' @export
inv_symlog <- function(x, C = 0) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "inv_symlog()",
    details = "Use scales::pseudo_log_trans() instead."
  )
  sign(x) * (10^abs(x) * 10^C - 10^C)
}

#' Deprecated: symlog_trans
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `symlog_trans()` is deprecated. Use [scales::pseudo_log_trans()] instead.
#'
#' @return A `scales` transformation object.
#' @rdname deprecated-functions
#' @export
symlog_trans <- function() {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "symlog_trans()",
    with = "scales::pseudo_log_trans()"
  )
  scales::trans_new(
    "symlog",
    transform = function(x) sign(x) * log10(1 + abs(x)),
    inverse = function(x) sign(x) * (10^abs(x) - 1)
  )
}

#' Deprecated: get_default_colors
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_default_colors()` is deprecated.
#'
#' @param var Variable name (`"category"` or `"group"`).
#' @param levels Character vector of factor levels.
#' @return Character vector of color values.
#' @rdname deprecated-functions
#' @export
get_default_colors <- function(var, levels) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "get_default_colors()"
  )
  .assert_that(all(var %in% c("category", "group")))
  .assert_that(is.character(levels))
  n <- length(levels)

  if (var == "category") {
    .assert_that(
      n <= 12,
      msg = "Cannot provide default colors for more than 12 levels."
    )
    palette.colors(n, "Set 3")
  } else {
    if ("prior" %in% levels) {
      .assert_that(
        n - 1 <= 36,
        msg = "Cannot provide default colors for more than 36 levels."
      )
      color <- character(n)
      color[which(levels != "prior")] <- palette.colors(
        n - 1,
        "Polychrome 36"
      )
      color[which(levels == "prior")] <- "darkgray"
      color
    } else {
      .assert_that(
        n <= 36,
        msg = "Cannot provide default colors for more than 36 levels."
      )
      palette.colors(n, "Polychrome 36")
    }
  }
}

#' Deprecated: get_default_shapes
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_default_shapes()` is deprecated.
#'
#' @param var Variable name.
#' @param levels Character vector of factor levels.
#' @return Integer vector of shapes.
#' @rdname deprecated-functions
#' @export
get_default_shapes <- function(var, levels) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "get_default_shapes()"
  )
  .assert_that(all(var %in% c("category", "group")))
  .assert_that(is.character(levels))
  seq_along(levels)
}

#' Deprecated: get_default_linetypes
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_default_linetypes()` is deprecated.
#'
#' @param var Variable name.
#' @param levels Character vector of factor levels.
#' @return Integer vector of linetypes.
#' @rdname deprecated-functions
#' @export
get_default_linetypes <- function(var, levels) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "get_default_linetypes()"
  )
  .assert_that(all(var %in% c("category", "group")))
  .assert_that(is.character(levels))
  seq_along(levels)
}

#' Deprecated: get_plot_limits
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_plot_limits()` is deprecated.
#'
#' @param plot A `ggplot` object.
#' @return List with elements `x` and `y`.
#' @rdname deprecated-functions
#' @export
get_plot_limits <- function(plot) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "get_plot_limits()"
  )
  built <- ggplot2::ggplot_build(plot)
  list(
    x = built$layout$panel_scales_x[[1]]$range$range,
    y = built$layout$panel_scales_y[[1]]$range$range
  )
}

#' Deprecated: get_limits
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_limits()` is deprecated.
#'
#' @param data A `tibble` or `data.frame`.
#' @param measure Variable name.
#' @param by Optional grouping variable.
#' @param hdi.prob HDI probability mass.
#' @param min Optional minimum.
#' @param max Optional maximum.
#' @return Numeric vector of length 2.
#' @rdname deprecated-functions
#' @export
get_limits <- function(
  data,
  measure,
  by = NULL,
  hdi.prob = 0.99,
  min = NULL,
  max = NULL
) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "get_limits()"
  )
  data %>%
    tidybayes::mean_hdi(!!rlang::sym(measure), .width = hdi.prob) %>%
    dplyr::ungroup() %>%
    {
      if (!is.null(by)) dplyr::group_by(., !!rlang::sym(by)) else .
    } %>%
    dplyr::summarise(
      .lower = if (!is.null(min)) min else base::min(.data$.lower),
      .upper = if (!is.null(max)) max else base::max(.data$.upper)
    ) %>%
    as.numeric()
}

#' Deprecated: ellipse.pmap
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `ellipse.pmap()` is deprecated.
#'
#' @param x Covariance matrix.
#' @param centre Center vector.
#' @param level Probability level.
#' @param ... Additional arguments.
#' @return Ellipse coordinates matrix.
#' @rdname deprecated-functions
#' @export
ellipse.pmap <- function(x, centre, level, ...) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "ellipse.pmap()"
  )
  ellipse::ellipse(x = x, centre = centre, level = level, ...)
}

#' Deprecated: add_exposure_data_to_1D_plot
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `add_exposure_data_to_1D_plot()` is deprecated. Use [plot_categories()]
#' instead.
#'
#' @param data Data frame.
#' @param cue.labels Cue column names.
#' @param category.ids Category values.
#' @param category.labels Category labels.
#' @param category.colors Category colors.
#' @return A list of ggplot layers.
#' @rdname deprecated-functions
#' @export
add_exposure_data_to_1D_plot <- function(
  data,
  cue.labels,
  category.ids,
  category.labels,
  category.colors
) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "add_exposure_data_to_1D_plot()",
    with = "plot_categories()"
  )
  cue.labels[2] <- "cue2"
  data$cue2 <- 0
  add_exposure_data_to_2D_plot(
    data,
    cue.labels,
    category.ids,
    category.labels,
    category.colors
  )
}

#' Deprecated: add_test_data_to_1D_plot
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `add_test_data_to_1D_plot()` is deprecated. Use [plot_categories()]
#' instead.
#'
#' @param data Data frame.
#' @param cue.labels Cue column names.
#' @return A list of ggplot layers.
#' @rdname deprecated-functions
#' @export
add_test_data_to_1D_plot <- function(data, cue.labels) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "add_test_data_to_1D_plot()",
    with = "plot_categories()"
  )
  cue.labels[2] <- "cue2"
  data$cue2 <- 0
  add_test_data_to_2D_plot(data, cue.labels)
}

#' Deprecated: add_exposure_data_to_2D_plot
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `add_exposure_data_to_2D_plot()` is deprecated. Use [plot_categories()]
#' instead.
#'
#' @param data Data frame.
#' @param cue.labels Cue column names.
#' @param category.ids Category values.
#' @param category.labels Category labels.
#' @param category.colors Category colors.
#' @return A list of ggplot layers.
#' @rdname deprecated-functions
#' @export
add_exposure_data_to_2D_plot <- function(
  data,
  cue.labels,
  category.ids,
  category.labels,
  category.colors
) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "add_exposure_data_to_2D_plot()",
    with = "plot_categories()"
  )
  list(
    ggplot2::geom_point(
      data = data,
      mapping = ggplot2::aes(
        x = .data[[cue.labels[1]]],
        y = .data[[cue.labels[2]]],
        shape = .data$category,
        color = .data$category
      ),
      size = 3,
      alpha = 0.9
    ),
    ggplot2::scale_shape(
      "Category",
      breaks = category.ids,
      labels = category.labels
    ),
    ggplot2::scale_color_manual(
      "Category",
      breaks = category.ids,
      labels = category.labels,
      values = category.colors
    )
  )
}

#' Deprecated: add_test_data_to_2D_plot
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `add_test_data_to_2D_plot()` is deprecated. Use [plot_categories()]
#' instead.
#'
#' @param data Data frame.
#' @param cue.labels Cue column names.
#' @return A list of ggplot layers.
#' @rdname deprecated-functions
#' @export
add_test_data_to_2D_plot <- function(data, cue.labels) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "add_test_data_to_2D_plot()",
    with = "plot_categories()"
  )
  list(
    ggplot2::geom_point(
      data = data,
      mapping = ggplot2::aes(
        x = .data[[cue.labels[1]]],
        y = .data[[cue.labels[2]]]
      ),
      inherit.aes = FALSE,
      color = "black",
      size = 1,
      alpha = 0.75
    )
  )
}

#' Deprecated: add_exposure_summary_to_1D_plot
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `add_exposure_summary_to_1D_plot()` is deprecated. Use [plot_categories()]
#' instead.
#'
#' @param data Data frame.
#' @return A list of ggplot layers.
#' @rdname deprecated-functions
#' @export
add_exposure_summary_to_1D_plot <- function(data) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "add_exposure_summary_to_1D_plot()",
    with = "plot_categories()"
  )
  data %>%
    dplyr::group_by(.data$category) %>%
    dplyr::summarise(
      mean = list(mean(.data$cue1)),
      sd = list(stats::sd(.data$cue1))
    ) %>%
    dplyr::group_map(
      ~ ggplot2::stat_function(
        fun = function(x) stats::dnorm(x, mean = .x$mean[[1]], sd = .x$sd[[1]]),
        mapping = ggplot2::aes(
          x = .data$cue1,
          color = .data$category
        ),
        linetype = 2,
        inherit.aes = FALSE
      )
    )
}

#' Deprecated: add_exposure_summary_to_2D_plot
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `add_exposure_summary_to_2D_plot()` is deprecated. Use [plot_categories()]
#' instead.
#'
#' @param data Data frame.
#' @param level Probability level.
#' @return A list of ggplot layers.
#' @rdname deprecated-functions
#' @export
add_exposure_summary_to_2D_plot <- function(data, level = 0.95) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "add_exposure_summary_to_2D_plot()",
    with = "plot_categories()"
  )
  list(
    ggplot2::geom_point(
      data = data %>%
        dplyr::mutate(
          cue1 = purrr::map_dbl(.data$mean, ~ .x[1]),
          cue2 = purrr::map_dbl(.data$mean, ~ .x[2])
        ),
      mapping = ggplot2::aes(
        x = .data$cue1,
        y = .data$cue2,
        color = .data$category
      ),
      inherit.aes = FALSE,
      size = 1
    ),
    ggplot2::geom_path(
      data = data %>%
        tidyr::crossing(level = level) %>%
        dplyr::mutate(
          ellipse = purrr::pmap(
            list(.data$cov, .data$mean, .data$level),
            ellipse.pmap
          )
        ) %>%
        tidyr::unnest(.data$ellipse) %>%
        dplyr::group_by(dplyr::across(-c("ellipse"))) %>%
        dplyr::transmute(
          cue1 = .data$ellipse[, 1],
          cue2 = .data$ellipse[, 2]
        ),
      mapping = ggplot2::aes(
        x = .data$cue1,
        y = .data$cue2,
        color = .data$category
      ),
      linetype = 2,
      inherit.aes = FALSE
    )
  )
}

#' Deprecated: facet_or_animate
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `facet_or_animate()` is deprecated.
#'
#' @param p A `ggplot` object.
#' @param facet_rows_by Row facet variable.
#' @param facet_cols_by Col facet variable.
#' @param facet_wrap_by Wrap facet variable.
#' @param animate_by Animation variable.
#' @param animation_follow Logical whether animation follows data.
#' @return A `ggplot` or animated object.
#' @rdname deprecated-functions
#' @export
facet_or_animate <- function(
  p,
  facet_rows_by,
  facet_cols_by,
  facet_wrap_by,
  animate_by,
  animation_follow
) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "facet_or_animate()"
  )
  facet_rows_by <- rlang::enquo(facet_rows_by)
  facet_cols_by <- rlang::enquo(facet_cols_by)
  facet_wrap_by <- rlang::enquo(facet_wrap_by)
  animate_by <- rlang::enquo(animate_by)

  if (!rlang::quo_is_null(facet_rows_by) ||
      !rlang::quo_is_null(facet_cols_by)) {
    p <- p + ggplot2::facet_grid(
      rows = ggplot2::vars(!!facet_rows_by),
      cols = ggplot2::vars(!!facet_cols_by),
      labeller = ggplot2::label_both
    )
  } else if (!rlang::quo_is_null(facet_wrap_by)) {
    p <- p + ggplot2::facet_wrap(
      facets = ggplot2::vars(!!facet_wrap_by),
      labeller = ggplot2::label_both
    )
  }

  if (!rlang::quo_is_null(animate_by)) {
    message("Preparing for rendering. This might take a moment.\n")
    p <- p +
      ggplot2::labs(
        title = paste0(rlang::as_name(animate_by), ": {closest_state}")
      ) +
      gganimate::transition_states(
        !!animate_by,
        transition_length = 1,
        state_length = 1
      ) +
      {
        if (animation_follow) gganimate::view_follow()
      } +
      gganimate::enter_fade() +
      gganimate::exit_fade()
  }

  p
}

#' Deprecated: plot_pairwise_cue_correlation_matrix
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `plot_pairwise_cue_correlation_matrix()` is deprecated. Use
#' `plot_parameters(model, type = "correlation")` instead.
#'
#' @param data Data frame.
#' @param cues Cue columns.
#' @param category Category column.
#' @param category.colors Category colors.
#' @return A `ggplot` object.
#' @rdname deprecated-functions
#' @export
plot_pairwise_cue_correlation_matrix <- function(
  data,
  cues,
  category = category,
  category.colors = seq_along(unique(data$category))
) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "plot_pairwise_cue_correlation_matrix()",
    with = "plot_parameters()"
  )
  .panel_x <- .panel_y <- NULL
  cues <- rlang::enquos(cues)
  category <- rlang::enquo(category)

  data %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = .panel_x,
        y = .panel_y,
        colour = !!category,
        fill = !!category
      )
    ) +
    ggplot2::scale_colour_manual(values = category.colors) +
    ggplot2::scale_fill_manual(values = category.colors) +
    ggplot2::geom_point(alpha = 0.6, shape = 16, size = 1) +
    ggforce::geom_autodensity(alpha = 0.04, position = "identity") +
    ggplot2::stat_ellipse(type = "norm") +
    ggforce::facet_matrix(
      ggplot2::vars(!!!cues),
      layer.diag = 2,
      layer.upper = 3,
      grid.y.diag = FALSE
    )
}
