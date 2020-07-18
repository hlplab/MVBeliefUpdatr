#' @import ggplot2 cowplot gganimate transformr av viridis
#' @importFrom ellipse ellipse
#' @importFrom mvtnorm dmvt
#' @importFrom scales trans_new
NULL


#' Plot expected bivariate (2D) categories.
#'
#' Plot bivariate Gaussian categories expected given NIW parameters.
#'
#' @param x NIW belief object.
#' @param grouping.var Grouping variable that is used to create separate instances of the categories. These instances
#' can be plotted in separate panels (if panel.group is `TRUE`) or used to create animates (if animate.group is `TRUE`).
#' @param panel.group,animate.group Determines whether the grouping variable is for paneling or animation. (both defaults: `FALSE`)
#' @param levels Levels of the confidence ellipses. (default: .5, .66, .8, .9., and .95)
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
  plot.exposure = F, plot.test = F,
  category.ids = NULL, category.labels = NULL, category.colors = NULL, category.linetypes = NULL
) {
  assert_that(!all(panel.group, animate.group))
  assert_that(is.NIW_belief(x))
  assert_that(!all(plot.exposure, is.null(data.exposure)),
              msg = "Can't plot exposure data: No exposure data provided.")
  assert_that(!all(plot.exposure, ),
              msg = "Can't plot exposure data: cue names in exposure data must match those in the NIW belief object.")
  assert_that(!all(plot.exposure, !is.null(grouping.var), !(grouping.var %in% names(data.exposure))),
              msg = "Can't plot exposure data: if a grouping variable is specified, it must be present in the exposure data.")
  assert_that(!all(plot.test, is.null(data.test)),
              msg = "Can't plot test data: No test data provided.")

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
    select(-c(kappa, nu, M, S, Sigma, lapse)) %>%
    unnest(ellipse)

  # This is somewhat fragile. E.g., if x has more columns than required. It assumes that the all
  # colums that are not the grouping variable, category, or level MUST be cues.
  cue.names = setdiff(names(x), c(if (!is.null(grouping.var)) grouping.var else NA, "category", "level"))
  message(paste("The following variables are assumed to be cues:", paste(cue.names, collapse = ", ")))
  x %<>%
    rename_at(cue.names,
              function(x) paste0("cue", which(x == cue.names)))

  p = ggplot(x,
             aes(x = cue1, y = cue2,
                 fill = category,
                 alpha = 1-level,
                 group = paste(category, level))) +
    geom_polygon() +
    # Optionally plot test data
    { if (plot.test)
      geom_point(
        data = fit.input$x_test %>%
          rename_at(cue.names,
                    function(x) paste0("cue", which(x == cue.names))),
        mapping = aes(cue1, cue2),
        inherit.aes = F,
        color = "black", size = 1
      )} +
    # Optionally plot exposure data
    { if (plot.exposure)
      geom_point(
        data = data.exposure,
        mapping = aes(cue1, cue2, shape = category, color = category),
        inherit.aes = F, size = 2
      ) +
        geom_path(
          data = crossing(
            group = levels(d$group),
            category = levels(d$category),
            level = .95
          ) %>%
            mutate(
              x = map2(category, group, get_ibbu_exposure_sigma(fit.input, .x, .y)),
              centre = map2(category, group, get_ibbu_exposure_mean(fit.input, .x, .y))
            ) %>%
            mutate(ellipse = pmap(., ellipse.pmap)) %>%
            mutate(ellipse = map(ellipse, as_tibble)) %>%
            unnest(ellipse) %>%
            rename_at(cue.names,
                      function(x) paste0("cue", which(x == cue.names))),
          mapping = aes(cue1, cue2, shape = category, color = category),
          linetype = 2,
          inherit.aes = F)
    } +
    scale_x_continuous(cue.names[1]) +
    scale_y_continuous(cue.names[2]) +
    scale_fill_manual("Category",
                      breaks = category.ids,
                      labels = category.labels,
                      values = category.colors) +
    scale_alpha("", range = c(.1, .5), breaks = levels) +
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



#' Plot expected categorization function for bivariate (2D) categories.
#'
#' Plot categorization function for bivariate Gaussian categories expected given NIW parameters.
#'
#' @param x NIW belief object.
#' @param grouping.var Grouping variable that is used to create separate instances of the categories. These instances
#' can be plotted in separate panels (if panel.group is `TRUE`) or used to create animates (if animate.group is `TRUE`).
#' @param panel.group,animate.group Determines whether the grouping variable is for paneling or animation. (both defaults: `FALSE`)
#' @param levels Levels of the confidence ellipses. (default: 5 steps from close to 0 to 95 percent probability mass)
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
  levels = plogis(seq(-15, qlogis(.95), length.out = 20)),
  # data.exposure = NULL,
  # data.test = NULL,
  # plot.exposure = F, plot.test = F,
  category.ids = NULL, category.labels = NULL, category.colors = NULL, category.linetypes = NULL
) {
  assert_that(!all(panel.group, animate.group))
  assert_that(is.NIW_belief(x))
  # assert_that(!all(is.null(data.exposure), plot.exposure),
  #             msg = "No exposure data provided.")
  # assert_that(!all(is.null(data.test), plot.test),
  #             msg = "No test data provided.")

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
    select(-c(kappa, nu, M, S, Sigma, lapse)) %>%
    unnest(ellipse)

  # This is somewhat fragile. E.g., if x has more columns than required. It assumes that the all
  # colums that are not the grouping variable, category, or level MUST be cues.
  cue.names = setdiff(names(x), c(if (!is.null(grouping.var)) grouping.var else NA, "category", "level"))
  message(paste("The following variables are assumed to be cues:", paste(cue.names, collapse = ", ")))
  x %<>%
    rename_at(cue.names,
              function(x) paste0("cue", which(x == cue.names)))

  p = ggplot(x,
             aes(x = cue1, y = cue2,
                 fill = category,
                 alpha = 1-level,
                 group = paste(category, level))) +
    geom_polygon() +
    scale_x_continuous(cue.names[1]) +
    scale_y_continuous(cue.names[2]) +
    scale_fill_manual("Category",
                      breaks = category.ids,
                      labels = category.labels,
                      values = category.colors) +
    scale_alpha("",
                range = c(0.1,.9)) +
    theme_bw()

  if (!is.null(grouping.var)) {
    if (panel.group) p = p + facet_wrap(vars(!! sym(grouping.var)))
    else if (animate.group) {
      message("Preparing for rendering. This might take a moment.\n")
      p = p +
        labs(title = "Observation: {closest_state}") +
        transition_states(!! sym(grouping.var),
                          transition_length = 1,
                          state_length = 1) +
        enter_fade() +
        exit_fade()
    }
  }

  return(p)
}
