#' @import ggplot2 cowplot gganimate transformr av viridis
#' @importFrom ellipse ellipse
#' @importFrom mvtnorm dmvt
#' @importFrom scales trans_new
NULL


#' Checking that input to plotting function is compatible
check_compatibility_between_NIW_belief_and_data = function(
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







#' Plot NIW belief or NIW beliefs object.
#'
#' Plot the parameters of an NIW_belief or NIW_beliefs object.
#'
#' @param x An \code{\link{NIW_belief}} or \code{\link{NIW_beliefs}} object.
#' @param group.colors Vector of fill colors of same length as the number of unique groups in the NIW_belief(s) object, or
#' `NULL` to use defaults. (default: `NULL`)
#' @param facet_rows_by,facet_cols_by,facet_wrap_by,animate_by Which group variables, if any, should be used for faceting and/or
#' animation? (defaults: `NULL`)
#' @param animation_follow Should the animation follow the data (zoom in and out)? (default: `FALSE`)
#'
#' @return ggplot object.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
plot_NIW_belief_parameters = function(
  x,
  group.colors = NULL,
  facet_rows_by = NULL, facet_cols_by = NULL, animate_by = NULL, animation_follow = F
) {
  error("This function is not yet implemented.")

  # Check out this function from another project:
  # plot_VOT_NIW_belief_1D <- function(belief, sigma_max = NULL, prior = NULL) {
  #   mu_sigma <- belief %>%
  #     mutate(
  #       mu = get_expected_mu_from_m(m),
  #       sigma = get_expected_Sigma_from_S(S, nu))
  #
  #   if (is.null(sigma_max)) sigma_max = max(mu_sigma$sigma) * 2
  #   if (!is.null(prior))
  #     prior.mu_sigma <- prior %>%
  #       mutate(
  #         mu = get_expected_mu_from_m(m),
  #         sigma = get_expected_Sigma_from_S(S, nu))
  #
  #   belief %>%
  #     crossing(
  #       mu = seq_range(VOT_range, n = VOT_resolution),
  #       sigma = seq_range(1:sigma_max^.5, n = VOT_resolution)^2) %>%
  #     { if ("Subject" %in% names(.)) group_by(., Subject) else . } %>%
  #     # TO DO: Since mu and sigma can probably be vectors this can be made more efficient by first nesting and then unnesting
  #     # (see what I did for the test_plot function below). The line above this has been in added in anticipation of that
  #     # change (it's currently not required since the density is obtained line by line).
  #     mutate(l = unlist(pmap(
  #       .l = list(mu, m, kappa, sigma, S, nu),
  #       .f = dnorminvwishart))) %>%
  #     ggplot(aes(x = mu, y = sigma, color = category, group = category)) +
  #     { if (is.null(prior))
  #       list(
  #         geom_raster(
  #           data = ~ filter(., category == "/b/"),
  #           mapping = aes(fill = category, alpha = l),
  #           interpolate = T),
  #         geom_raster(
  #           data = ~ filter(., category == "/p/"),
  #           mapping = aes(fill = category, alpha = l),
  #           interpolate = T)) } +
  #     geom_contour(aes(z = l, color = category), breaks = 10^(-10:-3), size = .5) +
  #     { if (!is.null(prior))
  #       geom_contour(
  #         data = prior %>%
  #           crossing(
  #             mu = seq_range(VOT_range, n = VOT_resolution),
  #             sigma = seq_range(1:sigma_max^.5, n = VOT_resolution)^2) %>%
  #           mutate(l = unlist(pmap(
  #             .l = list(mu, m, kappa, sigma, S, nu),
  #             .f = dnorminvwishart))),
  #         aes(z = l, color = category), breaks = 10^(-10:-3), size = .5, alpha = .1) } +
  #     { if (is.null(prior))
  #       geom_point(
  #         data = mu_sigma,
  #         aes(shape = category),
  #         color = "black") } +
  #     { if (!is.null(prior))
  #       list(
  #         geom_point(
  #           data = prior.mu_sigma,
  #           aes(shape = category),
  #           color = "black",
  #           alpha = .5),
  #         geom_segment(
  #           data =
  #             mu_sigma %>%
  #             left_join(
  #               prior.mu_sigma %>%
  #                 rename_at(vars(mu, sigma), ~ paste0("prior_", .x)),
  #               by = "category"),
  #           aes(x = prior_mu, y = prior_sigma, xend = mu, yend = sigma),
  #           arrow = arrow(angle = 15, length = unit(0.1, "inches"), ends = "last", type = "closed"),
  #           color = "black",
  #           size = .5,
  #           alpha = .75)) } +
  #     scale_x_continuous(name = bquote(mu ~ "(msec VOT)")) +
  #     scale_y_sqrt(name = bquote(sigma^2 ~ "(" ~ msec^2 ~ ")"), limits = c(0, sigma_max)) +
  #     scale_color_discrete("Category") +
  #     scale_fill_discrete("Category") +
  #     scale_shape_discrete("Category") +
  #     { if (is.null(prior)) scale_alpha_continuous("density", range = c(0,1), guide = "none") } +
  #     coord_cartesian(expand = 0) +
  #     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  # }
}



#' Plot expected univariate (1D) categories
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
#' @param category.ids Vector of category IDs to be plotted or leave `NULL` to plot all groups. (default: `NULL`) It is possible
#' to use \code{\link[tidybayes]{recover_types}} on the stanfit object prior to handing it to this plotting function.
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
#' @examples
#' TBD
#' @rdname plot_expected_categories
#' @export
plot_expected_categories_density1D = function(
  x,
  data.exposure = NULL,
  data.test = NULL,
  facet_rows_by = NULL, facet_cols_by = NULL, facet_wrap_by = NULL, animate_by = NULL, animation_follow = F,
  xlim, ylim = NULL,
  category.ids = NULL, category.labels = NULL, category.colors = NULL, category.linetypes = NULL,
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
  assert_that(length(cue.labels) == 1, msg = "Expecting exactly one cue for plotting.")

  if (is_missing(xlim)) {
    if (!is.null(data.exposure) & !is.null(data.test))
      xlim = range(range(data.exposure[[cue.labels[1]]]), range(data.test[[cue.labels[1]]])) else
        if (!is.null(data.exposure))
          xlim = range(data.exposure[[cue.labels[1]]]) else
            if (!is.null(data.test))
              xlim = range(data.test[[cue.labels[1]]])
  }
  assert_that(!is_missing(xlim), msg = "`xlim` must be specified")

  # Setting aes defaults
  if(is.null(category.ids)) category.ids = levels(x$category)
  if(is.null(category.labels)) category.labels = levels(x$category)
  if(is.null(category.colors)) category.colors = get_default_colors("category", length(category.ids))
  if(is.null(category.linetypes)) category.linetypes = rep(1, length(category.ids))

  if (any(!quo_is_null(facet_rows_by),
          !quo_is_null(facet_cols_by),
          !quo_is_null(animate_by))) x %<>% group_by(!! facet_rows_by, !! facet_cols_by, !! facet_wrap_by, !! animate_by,
                                                     .add = TRUE)

  list.stat_functions <-
    x %>%
    mutate(
      mu = get_expected_mu_from_m(m),
      Sigma = get_expected_Sigma_from_S(S, nu)) %>%
    group_by(category, .add = T) %>%
    group_map(
      .keep = T,
      .f = function(.x, .y)
        stat_function(
          data = .x,
          mapping = aes(color = category),
          fun = dnorm,
          args = list(mean = .x$mu, sd = .x$Sigma^.5),
          ...))

  p = ggplot(mapping = aes(color = category)) +
    list.stat_functions +
    { if (!is.null(data.test))
      add_test_data_to_1D_plot(data = data.test, cue.labels = cue.labels) } +
    { if (!is.null(data.exposure))
      add_exposure_data_to_1D_plot(data = data.exposure, cue.labels = cue.labels,
                                   category.ids = category.ids, category.labels = category.labels, category.colors) } +
    scale_x_continuous(cue.labels, limits = xlim) +
    scale_y_continuous("Density", limits = ylim) +
    scale_color_manual("Category",
                      breaks = category.ids,
                      labels = category.labels,
                      values = category.colors) +
    theme_bw()

  p = facet_or_animate(p, !!facet_rows_by, !!facet_cols_by, !! facet_wrap_by, !!animate_by, animation_follow)
  return(p)
}

#' Plot expected bivariate (2D) categories
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
#' @param category.ids Vector of category IDs to be plotted or leave `NULL` to plot all groups. (default: `NULL`) It is possible
#' to use \code{\link[tidybayes]{recover_types}} on the stanfit object prior to handing it to this plotting function.
#' @param category.labels Vector of group labels of same length as `category.ids` or `NULL` to use defaults. (default: `NULL`)
#' @param category.colors Vector of colors of same length as category.ids or `NULL` to use defaults. (default: `NULL`)
#' @param category.linetypes Vector of linetypes of same length as category.ids or `NULL` to use defaults. (default: `NULL`)
#' Currently being ignored.
#' @param ... additional arguments handed to geom_polygon.
#'
#' @return ggplot object.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @rdname plot_expected_categories
#' @export
plot_expected_categories_contour2D = function(
  x,
  levels = c(1/2, 2/3, 4/5, 9/10, 19/20),
  data.exposure = NULL,
  data.test = NULL,
  facet_rows_by = NULL, facet_cols_by = NULL, animate_by = NULL, animation_follow = F,
  category.ids = NULL, category.labels = NULL, category.colors = NULL, category.linetypes = NULL,
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
  if(is.null(category.ids)) category.ids = levels(x$category)
  if(is.null(category.labels)) category.labels = levels(x$category)
  if(is.null(category.colors)) category.colors = get_default_colors("category", length(category.ids))
  if(is.null(category.linetypes)) category.linetypes = rep(1, length(category.ids))

  x %<>%
    mutate(Sigma = get_expected_Sigma_from_S(S, nu)) %>%
    crossing(level = levels) %>%
    mutate(ellipse = pmap(.l = list(Sigma, m, level), ellipse.pmap)) %>%
    # This step is necessary since unnest() can't yet unnest lists of matrices
    # (bug was reported and added as milestone, 11/2019)
    mutate(ellipse = map(ellipse, as_tibble)) %>%
    select(-c(kappa, nu, m, S, Sigma, lapse_rate)) %>%
    unnest(ellipse) %>%
    # Get group structure again, as crossing apparently removes it
    group_by(!!! syms(group_vars(x)))

  p = ggplot(x,
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
      add_test_data_to_2D_plot(data = data.test, cue.labels = cue.labels) } +
    { if (!is.null(data.exposure))
      add_exposure_data_to_2D_plot(data = data.exposure, cue.labels = cue.labels,
                                   category.ids = category.ids, category.labels = category.labels, category.colors) } +
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

  p = facet_or_animate(p, !!facet_rows_by, !!facet_cols_by, !! facet_wrap_by, !!animate_by, animation_follow)
  return(p)
}




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
#' @examples
#' TBD
#' @rdname plot_expected_categorization_function_1D
#' @export
#'
plot_expected_categorization_function_1D = function(
  x,
  data.exposure = NULL,
  data.test = NULL,
  target_category = 1,
  logit = F,
  xlim, ylim = NULL,
  facet_rows_by = NULL, facet_cols_by = NULL, facet_wrap_by = NULL, animate_by = NULL, animation_follow = F,
  category.ids = NULL, category.labels = NULL, category.colors = NULL, category.linetypes = NULL,
  ...
) {
  facet_rows_by = enquo(facet_rows_by)
  facet_cols_by = enquo(facet_cols_by)
  facet_wrap_by = enquo(facet_wrap_by)
  animate_by = enquo(animate_by)
  check_compatibility_between_NIW_belief_and_data(x, data.exposure, data.test,
                                                  !! facet_rows_by, !! facet_cols_by, !! facet_wrap_by, !! animate_by)
  cue.labels = get_cue_labels_from_model(x)
  assert_that(length(cue.labels) == 1, msg = "Expecting exactly one cue for plotting.")

  if (is_missing(xlim)) {
    if (!is.null(data.exposure) & !is.null(data.test))
      xlim = range(range(data.exposure[[cue.labels[1]]]), range(data.test[[cue.labels[1]]])) else
        if (!is.null(data.exposure))
          xlim = range(data.exposure[[cue.labels[1]]]) else
            if (!is.null(data.test))
              xlim = range(data.test[[cue.labels[1]]])
  }
  assert_that(!is_missing(xlim), msg = "`xlim` must be specified")

  # Setting aes defaults
  if(is.null(category.ids)) category.ids = levels(x$category)
  if(is.null(category.labels)) category.labels = levels(x$category)
  if(is.null(category.colors)) category.colors = get_default_colors("category", length(category.ids))
  if(is.null(category.linetypes)) category.linetypes = rep(1, length(category.ids))

  if (any(!quo_is_null(facet_rows_by),
          !quo_is_null(facet_cols_by),
          !quo_is_null(animate_by))) x %<>% group_by(!! facet_rows_by, !! facet_cols_by, !! facet_wrap_by, !! animate_by,
                                                     .add = TRUE)

  list.stat_functions <-
    x %>%
      group_map(
        .keep = T,
        .f = function(.x, .y) {
          cat_function <- get_categorization_function_from_NIW_ideal_adaptor(.x, logit = logit)
          stat_function(
            data = .x,
            fun = cat_function, ...) })

  p =
    ggplot() +
    list.stat_functions +
    { if (!is.null(data.test))
      add_test_data_to_1D_plot(data = data.test, cue.labels = cue.labels) } +
    { if (!is.null(data.exposure))
      add_exposure_data_to_1D_plot(data = data.exposure, cue.labels = cue.labels,
                                   category.ids = category.ids, category.labels = category.labels, category.colors) } +
    scale_x_continuous(name = cue.labels, limits = xlim) +
    scale_y_continuous(name = if (logit)
      paste0("log-odds(resp = ", category.labels[target_category], ")") else
        paste0("p(resp = ", category.labels[target_category], ")")) +
    coord_cartesian(ylim = ylim) +
    theme_bw()

  p = facet_or_animate(p, !!facet_rows_by, !!facet_cols_by, !! facet_wrap_by, !!animate_by, animation_follow)
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
#' @examples
#' TBD
#' @rdname plot_expected_categorization_function_2D
#' @export
#'
plot_expected_categorization_function_2D = function(
  x,
  data.exposure = NULL,
  data.test = NULL,
  target_category = 1,
  logit = F,
  xlim, ylim, resolution = 25,
  facet_rows_by = NULL, facet_cols_by = NULL, facet_wrap_by = NULL, animate_by = NULL, animation_follow = F,
  category.ids = NULL, category.labels = NULL, category.colors = NULL, category.linetypes = NULL,
  ...
) {
  facet_rows_by = enquo(facet_rows_by)
  facet_cols_by = enquo(facet_cols_by)
  facet_wrap_by = enquo(facet_wrap_by)
  animate_by = enquo(animate_by)
  check_compatibility_between_NIW_belief_and_data(x, data.exposure, data.test,
                                                  !! facet_rows_by, !! facet_cols_by, !! facet_wrap_by, !! animate_by)
  cue.labels = get_cue_labels_from_model(x)
  assert_that(length(cue.labels) == 2, msg = "Expecting exactly two cues for plotting.")
  if (is_missing(xlim)) {
    if (!is.null(data.exposure) & !is.null(data.test))
      xlim = range(range(data.exposure[[cue.labels[1]]]), range(data.test[[cue.labels[1]]])) else
        if (!is.null(data.exposure))
          xlim = range(data.exposure[[cue.labels[1]]]) else
            if (!is.null(data.test))
              xlim = range(data.test[[cue.labels[1]]])
  }
  if (is_missing(ylim)) {
    if (!is.null(data.exposure) & !is.null(data.test))
      ylim = range(range(data.exposure[[cue.labels[2]]]), range(data.test[[cue.labels[2]]])) else
        if (!is.null(data.exposure))
          ylim = range(data.exposure[[cue.labels[2]]]) else
            if (!is.null(data.test))
              ylim = range(data.test[[cue.labels[2]]])
  }
  assert_that(!is_missing(xlim), msg = "`xlim` must be specified")
  assert_that(!is_missing(ylim), msg = "`ylim` must be specified")

  # Setting aes defaults
  if(is.null(category.ids)) category.ids = levels(x$category)
  if(is.null(category.labels)) category.labels = levels(x$category)
  if(is.null(category.colors)) category.colors = get_default_colors("category", length(category.ids))
  if(is.null(category.linetypes)) category.linetypes = rep(1, length(category.ids))

  if (any(!quo_is_null(facet_rows_by),
          !quo_is_null(facet_cols_by),
          !quo_is_null(animate_by))) x %<>% group_by(!! facet_rows_by, !! facet_cols_by, !! animate_by,
                                                     .add = TRUE)

  d = crossing(
    !! sym(cue.labels[1]) := seq(min(xlim), max(xlim), length.out = resolution),
    !! sym(cue.labels[2]) := seq(min(ylim), max(ylim), length.out = resolution))

  x %<>%
    nest() %>%
    mutate(f = map(data, get_categorization_function_from_NIW_ideal_adaptor, logit = logit)) %>%
    # Join in vectored cues
    left_join(
      d %>%
        transmute(x = pmap(.l = list(!!! syms(cue.labels)), .f = ~ c(...))) %>%
        nest(cues = everything()),
      by = character()) %>%
    mutate(
      p_cat = invoke_map(.f = f, .x = cues, target_category = target_category),
      cues = NULL,
      f = NULL) %>%
    # Join separate cues back in
    left_join(d %>% nest(cues = everything()), by = character()) %>%
    unnest(c(cues, p_cat))

  p = ggplot(x,
             mapping = aes(
               x = .data[[cue.labels[1]]],
               y = .data[[cue.labels[2]]])) +
    geom_raster(mapping = aes(fill = if (logit) qlogis(.data$p_cat) else .data$p_cat), ...) +
    # geom_contour(
    #   mapping = aes(z = if (logit) qlogis(.data$p_cat) else .data$p_cat)) +
    { if (!is.null(data.test))
      add_test_data_to_2D_plot(data = data.test, cue.labels = cue.labels) } +
    { if (!is.null(data.exposure))
      add_exposure_data_to_2D_plot(data = data.exposure, cue.labels = cue.labels,
                                   category.ids = category.ids, category.labels = category.labels, category.colors) } +
    scale_x_continuous(cue.labels[1]) +
    scale_y_continuous(cue.labels[2]) +
    # For now think about two colors and categories
    scale_fill_gradient2(paste0("p(resp = ", category.labels[target_category], ")"),
                         low = category.colors[1],
                         mid = "white",
                         high = category.colors[2],
                         midpoint = if (logit) 0 else .5) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_bw()

  p = facet_or_animate(p, !!facet_rows_by, !!facet_cols_by, !! facet_wrap_by, !!animate_by, animation_follow)
  return(p)
}
