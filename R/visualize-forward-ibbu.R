#' @import ggplot2 cowplot gganimate transformr av viridis
#' @importFrom ellipse ellipse
#' @importFrom mvtnorm dmvt
#' @importFrom scales trans_new
NULL

#' #' Checking that input to plotting function is compatible
#' check_compatibility_between_NIW_belief_and_data = function(
#'   x,
#'   data.exposure,
#'   data.test,
#'   grouping.var,
#'   panel.group,
#'   animate.group,
#'   return.cues = T
#' ) {
#'   assert_that(!all(panel.group, animate.group))
#'   assert_NIW_belief(x)
#'
#'   cue.labels = get_cue_labels_from_NIW_belief(x)
#'   message(paste("The following variables are assumed to be cues:", paste(cue.labels, collapse = ", ")))
#'
#'   assert_that(!all(!is.null(data.exposure), cue.labels %nin% names(data.exposure)),
#'               msg = "Can't plot exposure data: cue names in exposure data must match those in the NIW belief object.")
#'   assert_that(!all(!is.null(data.exposure), "category" %nin% names(data.exposure)),
#'               msg = "Can't plot exposure data: exposure data does not contain column category.")
#'   assert_that(!all(!is.null(data.exposure), !is.null(grouping.var), grouping.var %nin% names(data.exposure)),
#'               msg = "Can't plot exposure data: if a grouping variable is specified, it must be present in the exposure data.")
#'   assert_that(!all(!is.null(data.test), is.null(data.test)),
#'               msg = "Can't plot test data: No test data provided.")
#'   assert_that(!all(!is.null(data.test), cue.labels %nin% names(data.test)),
#'               msg = "Can't plot test data: cue names in test data must match those in the NIW belief object.")
#'
#'   if (return.cues) return(cue.labels)
#' }
#'

#' Plot expected bivariate (2D) categories.
#'
#' Plot bivariate Gaussian categories expected given NIW belief(s). One NIW belief describes the uncertainty about the
#' category statistics of all categories. This includes the M (the mean of category means \eqn{\mu}), S (the scattermatrix),
#' kappa (the strength of the belief in M) and nu (the strength of the belief in S).
#'
#' It is possible to hand more than one NIW belief to this function, one for each level of a grouping variable. This allows
#' to plot different sets of beliefs, for example, across different panels or as an animation. For example, one can plot
#' different priors for different talkers (grouping by talker), or different posteriors for different exposure conditions
#' (grouping by exposure condition), or the incremental updating of NIW beliefs (grouping by observations).
#'
#' @param x NIW belief object.
#' @param grouping.var Grouping variable that is used to create separate instances of the categories. These instances
#' can be plotted in separate panels (if panel.group is `TRUE`) or used to create animates (if animate.group is `TRUE`).
#' @param panel.group,animate.group Determines whether the grouping variable is for paneling or animation. (both defaults: `FALSE`)
#' @param levels Levels of the confidence ellipses. (default: .5, .66, .8, .9., and .95)
#' @param data.exposure Optional \code{tibble} or \code{data.frame} that contains exposure data to be plotted. (default: `NULL`)
#' @param data.test Optional \code{tibble} or \code{data.frame} that contains test data to be plotted. (default: `NULL`)
#' @param category.ids Vector of category IDs to be plotted or leave `NULL` to plot all groups. (default: `NULL`) It is possible
#' to use \code{\link[tidybayes]{recover_types}} on the stanfit object prior to handing it to this plotting function.
#' @param category.labels Vector of group labels of same length as `category.ids` or `NULL` to use defaults. (default: `NULL`)
#' @param category.colors Vector of colors of same length as category.ids or `NULL` to use defaults. (default: `NULL`)
#' @param category.linetypes Vector of linetypes of same length as category.ids or `NULL` to use defaults. (default: `NULL`)
#' Currently being ignored.
#'
#' @return ggplot object.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @rdname plot_expected_categories_2D
#' @export
#'
plot_expected_categories_contour2D = function(
  x,
  grouping.var = NULL, panel.group = F, animate.group = F,
  levels = c(1/2, 2/3, 4/5, 9/10, 19/20),
  data.exposure = NULL,
  data.test = NULL,
  category.ids = NULL, category.labels = NULL, category.colors = NULL, category.linetypes = NULL
) {
  cue.labels = get_cue_labels_from_NIW_belief(x)
  assert_that(length(cue.labels) == 2, msg = "Expecting exactly two cues for plotting.")

  # Setting aes defaults
  if(is.null(category.ids)) category.ids = levels(x$category)
  if(is.null(category.labels)) category.labels = levels(x$category)
  if(is.null(category.colors)) category.colors = get_default_colors("category", length(category.ids))
  if(is.null(category.linetypes)) category.linetypes = rep(1, length(category.ids))

  x %<>%
    { if(!is.null(grouping.var)) mutate(., !! sym(grouping.var) := factor(!! sym(grouping.var))) else . } %>%
    mutate(Sigma = map2(S, nu, get_Sigma_from_S)) %>%
    crossing(level = levels) %>%
    mutate(ellipse = pmap(.l = list(Sigma, M, level), ellipse.pmap)) %>%
    # This step is necessary since unnest() can't yet unnest lists of matrices
    # (bug was reported and added as milestone, 11/2019)
    mutate(ellipse = map(ellipse, as_tibble)) %>%
    select(-c(kappa, nu, M, S, Sigma, lapse_rate)) %>%
    unnest(ellipse)

  p = ggplot(x,
             aes(
               x = .data[[cue.labels[1]]],
               y = .data[[cue.labels[2]]],
               fill = .data$category)) +
    geom_polygon(aes(alpha = 1 - .data$level,
                     group = paste(.data$category, .data$level))) +
    # Optionally plot test data
    { if (!is.null(data.test))
      geom_point(
        data = data.test,
        mapping = aes(
          x = .data[[cue.labels[1]]],
          y = .data[[cue.labels[2]]]),
        inherit.aes = F,
        color = "black", size = 1) } +
    # Optionally plot exposure data
    { if (!is.null(data.exposure))
      list(
        geom_point(
          data = data.exposure,
          mapping = aes(shape = .data$category, color = .data$category),
          size = 3, alpha = .9),
        scale_shape("Category"),
        scale_color_manual("Category",
                           breaks = category.labels,
                           labels = category.labels,
                           values = category.colors)) } +
    scale_x_continuous(cue.labels[1]) +
    scale_y_continuous(cue.labels[2]) +
    scale_fill_manual("Category",
                      breaks = category.ids,
                      labels = category.labels,
                      values = category.colors) +
    scale_alpha_continuous("Cumulative\nprobability",
                           range = c(0, .3),
                           breaks = round(1 - levels, 2)) +
    theme_bw()

  if (!is.null(grouping.var)) {
    if (panel.group) p = p + facet_wrap(vars(!! sym(grouping.var)),
                                        labeller = label_both)
    else if (animate.group) {
      message("Preparing for rendering. This might take a moment.\n")
      p = p +
        labs(title = paste0(grouping.var, ": {closest_state}")) +
        transition_states(!! sym(grouping.var),
                          transition_length = 1,
                          state_length = 1) +
        enter_fade() +
        exit_fade()
    }
  }

  return(p)
}




#' Checking that input to plotting function is compatible
check_compatibility_between_NIW_belief_and_data = function(
  x,
  data.exposure,
  data.test,
  .group_by,
  facet_rows_by, facet_cols_by, animate_by
) {
  .group_by = enquo(.group_by)
  facet_rows_by = enquo(facet_rows_by)
  facet_cols_by = enquo(facet_cols_by)
  animate_by = enquo(animate_by)

  assert_that(is.NIW_belief(x))

  if (!quo_is_null(.group_by)) {
    assert_that(all(as_name(.group_by) %in% names(x)),
                msg = paste("Grouping variable(s) ", as_name(.group_by), " not found in x."))

    if (!is_empty(groups(x))) message("Overriding grouping structure of x with groups specified in .group_by.")
    x %<>%
      group_by(!! .group_by)
  }

  if (!quo_is_null(facet_rows_by)) {
    assert_that(all(as_name(facet_rows_by) %in% groups(x)),
                msg = paste(as_name(facet_rows_by), " not found in the groups of x."))
    assert_that(!all(!is.null(data.exposure), as_name(facet_rows_by) %nin% names(data.exposure)),
                msg = "Can't plot exposure data: when facet_rows_by is specified, it must be present in the exposure data.")
  }
  if (!quo_is_null(facet_cols_by)) {
    assert_that(all(as_name(facet_cols_by) %in% groups(x)),
                msg = paste(as_name(facet_cols_by), " not found in the groups of x."))
    assert_that(!all(!is.null(data.exposure), as_name(facet_cols_by) %nin% names(data.exposure)),
                msg = "Can't plot exposure data: when facet_cols_by is specified, it must be present in the exposure data.")
  }
  if (!quo_is_null(animate_by)) {
    assert_that(all(as_name(animate_by) %in% groups(x)),
                msg = paste(as_name(animate_by), " not found in the groups of x."))
    assert_that(!all(!is.null(data.exposure), as_name(animate_by) %nin% names(data.exposure)),
                msg = "Can't plot exposure data: when animate_by is specified, it must be present in the exposure data.")
  }

  cue.labels = get_cue_labels_from_NIW_belief(x)
  assert_that(!all(!is.null(data.exposure), cue.labels %nin% names(data.exposure)),
              msg = "Can't plot exposure data: cue names in exposure data must match those in the NIW belief object.")
  assert_that(!all(!is.null(data.exposure), "category" %nin% names(data.exposure)),
              msg = "Can't plot exposure data: exposure data does not contain column category.")
  assert_that(!all(!is.null(data.test), is.null(data.test)),
              msg = "Can't plot test data: No test data provided.")
  assert_that(!all(!is.null(data.test), cue.labels %nin% names(data.test)),
              msg = "Can't plot test data: cue names in test data must match those in the NIW belief object.")

  return(x)
}


#' Plot expected bivariate (2D) categories.
#'
#' Plot bivariate Gaussian categories expected given NIW belief(s). One NIW belief describes the uncertainty about the
#' category statistics of all categories. This includes the M (the mean of category means \eqn{\mu}), S (the scattermatrix),
#' kappa (the strength of the belief in M) and nu (the strength of the belief in S).
#'
#' It is possible to hand more than one NIW belief to this function, as long a the grouping structure is clear. This allows
#' to plot different sets of beliefs, for example, across different panels or as an animation. Grouping structure can be
#' indicated through the argument \code{.group_by} or by providing a grouped NIW belief object. Any of the grouping variables
#' can then be used to facet or animate. For example, one can plot
#' different priors for different talkers (grouping by talker), or different posteriors for different exposure conditions
#' (grouping by exposure condition), the incremental updating of NIW beliefs (grouping by observations), or any combinations
#' of these.
#'
#' @param x NIW belief object.
#' @param .group_by Grouping variable that is used to create separate instances of the categories. These instances
#' can be plotted in separate facets or used to animate the plot.
#' @param facet_rows_by,facet_cols_by,animate_by Which group variables, if any, should be used for faceting and/or
#' animation? (defaults: `NULL`)
#' @param levels Levels of the confidence ellipses. (default: .5, .66, .8, .9., and .95)
#' @param data.exposure Optional \code{tibble} or \code{data.frame} that contains exposure data to be plotted. (default: `NULL`)
#' @param data.test Optional \code{tibble} or \code{data.frame} that contains test data to be plotted. (default: `NULL`)
#' @param category.ids Vector of category IDs to be plotted or leave `NULL` to plot all groups. (default: `NULL`) It is possible
#' to use \code{\link[tidybayes]{recover_types}} on the stanfit object prior to handing it to this plotting function.
#' @param category.labels Vector of group labels of same length as `category.ids` or `NULL` to use defaults. (default: `NULL`)
#' @param category.colors Vector of colors of same length as category.ids or `NULL` to use defaults. (default: `NULL`)
#' @param category.linetypes Vector of linetypes of same length as category.ids or `NULL` to use defaults. (default: `NULL`)
#' Currently being ignored.
#'
#' @return ggplot object.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @rdname plot_expected_categories_2D
#' @export
#'
plot_expected_categories_contour2D_new = function(
  x,
  .group_by = NULL,
  facet_rows_by = NULL, facet_cols_by = NULL, animate_by = NULL,
  levels = c(1/2, 2/3, 4/5, 9/10, 19/20),
  data.exposure = NULL,
  data.test = NULL,
  category.ids = NULL, category.labels = NULL, category.colors = NULL, category.linetypes = NULL
) {
  facet_rows_by = enquo(facet_rows_by)
  facet_cols_by = enquo(facet_cols_by)
  animate_by = enquo(animate_by)
  x = check_compatibility_between_NIW_belief_and_data(x, data.exposure, data.test,
                                                      .group_by,
                                                      !! facet_rows_by, !! facet_cols_by, !! animate_by)
  # Remember groups
  cue.labels = get_cue_labels_from_NIW_belief(x)
  assert_that(length(cue.labels) == 2, msg = "Expecting exactly two cues for plotting.")

  # Setting aes defaults
  if(is.null(category.ids)) category.ids = levels(x$category)
  if(is.null(category.labels)) category.labels = levels(x$category)
  if(is.null(category.colors)) category.colors = get_default_colors("category", length(category.ids))
  if(is.null(category.linetypes)) category.linetypes = rep(1, length(category.ids))

  x %<>%
    mutate(Sigma = map2(S, nu, get_Sigma_from_S)) %>%
    crossing(level = levels) %>%
    mutate(ellipse = pmap(.l = list(Sigma, M, level), ellipse.pmap)) %>%
    # This step is necessary since unnest() can't yet unnest lists of matrices
    # (bug was reported and added as milestone, 11/2019)
    mutate(ellipse = map(ellipse, as_tibble)) %>%
    select(-c(kappa, nu, M, S, Sigma, lapse_rate)) %>%
    unnest(ellipse)

  p = ggplot(x,
             aes(
               x = .data[[cue.labels[1]]],
               y = .data[[cue.labels[2]]],
               fill = .data$category)) +
    geom_polygon(aes(alpha = 1 - .data$level,
                     group = paste(.data$category, .data$level))) +
    # Optionally plot test data
    { if (!is.null(data.test))
      geom_point(
        data = data.test,
        mapping = aes(
          x = .data[[cue.labels[1]]],
          y = .data[[cue.labels[2]]]),
        inherit.aes = F,
        color = "black", size = 1) } +
    # Optionally plot exposure data
    { if (!is.null(data.exposure))
      list(
        geom_point(
          data = data.exposure,
          mapping = aes(shape = .data$category, color = .data$category),
          size = 3, alpha = .9),
        scale_shape("Category"),
        scale_color_manual("Category",
                           breaks = category.labels,
                           labels = category.labels,
                           values = category.colors)) } +
    scale_x_continuous(cue.labels[1]) +
    scale_y_continuous(cue.labels[2]) +
    scale_fill_manual("Category",
                      breaks = category.ids,
                      labels = category.labels,
                      values = category.colors) +
    scale_alpha_continuous("Cumulative\nprobability",
                           range = c(0, .3),
                           breaks = round(1 - levels, 2)) +
    theme_bw()

  if (!quo_is_null(facet_rows_by) | !quo_is_null(facet_cols_by))
    p = p + facet_grid(
      rows = vars(!! facet_rows_by),
      cols = vars(!! facet_cols_by),
      labeller = label_both)
  if (!quo_is_null(animate_by)) {
    message("Preparing for rendering. This might take a moment.\n")
    p = p +
      labs(title = paste0(!! animate_by, ": {closest_state}")) +
      transition_states(!! animate_by,
                        transition_length = 1,
                        state_length = 1) +
      enter_fade() +
      exit_fade()
  }

  return(p)
}


#' Plot expected categorization function for bivariate (2D) categories.
#'
#' Plot categorization function for bivariate Gaussian categories expected given NIW parameters.
#'
#' @param x NIW belief object.
#' @param grouping.var Grouping variable that is used to create separate instances of the categories. These instances
#' can be plotted in separate panels (if panel.group is `TRUE`) or used to create animates (if animate.group is `TRUE`).
#' @param panel.group,animate.group Determines whether the grouping variable is for paneling or animation. (both defaults: `FALSE`)
#' @param data.exposure Optional \code{tibble} or \code{data.frame} that contains exposure data to be plotted. (default: `NULL`)
#' @param data.test Optional \code{tibble} or \code{data.frame} that contains test data to be plotted. (default: `NULL`)
#' @param target_category The index of the category for which categorization should be shown. (default: `1`)
#' @param xlim,ylim Limits for the x- and y-axis.
#' @param resolution How many steps along x and y should be calculated? Note that computational
#' complexity increases quadratically with resolution. (default: 25)
#' @param category.ids Vector of category IDs to be plotted or leave `NULL` to plot all groups. (default: `NULL`) It is possible
#' to use \code{\link[tidybayes]{recover_types}} on the stanfit object prior to handing it to this plotting function.
#' @param category.labels Vector of group labels of same length as `category.ids` or `NULL` to use defaults. (default: `NULL`)
#' @param category.colors Vector of colors of same length as category.ids or `NULL` to use defaults. (default: `NULL`)
#' @param category.linetypes Vector of linetypes of same length as category.ids or `NULL` to use defaults. (default: `NULL`)
#' Currently being ignored.
#'
#' @return ggplot object.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @rdname plot_expected_categorization_function_2D
#' @export
#'
plot_expected_categorization_function_2D = function(
  x,
  grouping.var = NULL, panel.group = F, animate.group = F,
  data.exposure = NULL,
  data.test = NULL,
  target_category = 1,
  xlim, ylim, resolution = 25,
  category.ids = NULL, category.labels = NULL, category.colors = NULL, category.linetypes = NULL
) {
  cue.labels = check_compatibility_between_NIW_belief_and_data(x, data.exposure, data.test,
                                                               grouping.var, panel.group, animate.group,
                                                               return.cues = T)
  assert_that(length(cue.labels) == 2, msg = "Expecting exactly two cues for plotting.")
  assert_that(!missing(xlim), msg = "`xlim` must be specified")
  assert_that(!missing(ylim), msg = "`ylim` must be specified")

  # Setting aes defaults
  if(is.null(category.ids)) category.ids = levels(x$category)
  if(is.null(category.labels)) category.labels = levels(x$category)
  if(is.null(category.colors)) category.colors = get_default_colors("category", length(category.ids))
  if(is.null(category.linetypes)) category.linetypes = rep(1, length(category.ids))

  x %<>%
    { if(!is.null(grouping.var)) mutate(., !! sym(grouping.var) := factor(!! sym(grouping.var))) else . }

  d = crossing(
    !! sym(cue.labels[1]) := seq(min(xlim), max(xlim), length.out = resolution),
    !! sym(cue.labels[2]) := seq(min(ylim), max(ylim), length.out = resolution)
  )

  d %<>%
    cbind(get_posterior_predictives_from_NIW_beliefs(d, x, wide = T, log = T, grouping.var = grouping.var))

  # TO BE DONE: handle case that grouping var might be part OF THE OUTPUT <------------------------- CONTINUE HERE
  log_p = d %>%
    select(starts_with("lpp."))

  d %<>%
    mutate(p_target =
             exp(
               log_p[,target_category] + log(priors[target_category]) -
                 log(rowSums(exp(log_p) * priors))) *
             # Assuming a uniform (unbiased) lapse rate:
             (1 - lapse_rate) + lapse_rate / n.cat)




  p = ggplot(x,
             aes(
               x = .data[[cue.labels[1]]],
               y = .data[[cue.labels[2]]],
               fill = .data$category)) +
    geom_polygon(aes(alpha = 1 - .data$level,
                     group = paste(.data$category, .data$level))) +
    # Optionally plot test data
    { if (!is.null(data.test))
      geom_point(
        data = data.test,
        mapping = aes(
          x = .data[[cue.labels[1]]],
          y = .data[[cue.labels[2]]]),
        inherit.aes = F,
        color = "black", size = 1) } +
    # Optionally plot exposure data
    { if (!is.null(data.exposure))
      list(
        geom_point(
          data = data.exposure,
          mapping = aes(shape = .data$category, color = .data$category),
          size = 3, alpha = .9),
        scale_shape("Category"),
        scale_color_manual("Category",
                           breaks = category.labels,
                           labels = category.labels,
                           values = category.colors)) } +
    scale_x_continuous(cue.labels[1]) +
    scale_y_continuous(cue.labels[2]) +
    scale_fill_manual("Category",
                      breaks = category.ids,
                      labels = category.labels,
                      values = category.colors) +
    scale_alpha_continuous("Cumulative\nprobability",
                           range = c(0, .3),
                           breaks = round(1 - levels, 2)) +
    theme_bw()

  if (!is.null(grouping.var)) {
    if (panel.group) p = p + facet_wrap(vars(!! sym(grouping.var)),
                                        labeller = label_both)
    else if (animate.group) {
      message("Preparing for rendering. This might take a moment.\n")
      p = p +
        labs(title = paste0(grouping.var, ": {closest_state}")) +
        transition_states(!! sym(grouping.var),
                          transition_length = 1,
                          state_length = 1) +
        enter_fade() +
        exit_fade()
    }
  }

  return(p)
}
