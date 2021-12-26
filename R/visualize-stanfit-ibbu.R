#' @import ggplot2 cowplot gganimate transformr av
#' @importFrom ellipse ellipse
#' @importFrom mvtnorm dmvt
#' @importFrom tidybayes mean_hdi
#' @importFrom ggridges geom_density_ridges
#' @importFrom forcats fct_rev
NULL


#' Plot distribution of IBBU parameters.
#'
#' Plot distribution of post-warmup MCMC samples for all parameters representing the
#' prior and/or posterior beliefs.
#'
#' @param model mv-ibbu-stanfit object.
#' @param which Should parameters for the prior, posterior, or both be plotted? (default: `"both"`)
#' @param ndraws Number of draws to plot (or use to calculate the CIs), or `NULL` if all draws are to be returned. (default: `NULL`)
#' @param untransform_cues Should m_0 and S_0 be transformed back into the original cue space? (default: `TRUE`)
#' @param group.ids Vector of group IDs to be plotted or leave `NULL` to plot all groups. (default: `NULL`) It is possible
#' to use \code{\link[tidybayes]{recover_types}} on the stanfit object prior to handing it to this plotting function.
#' @param group.labels Vector of group labels of same length as group.ids or `NULL` to use defaults. (default: `NULL`)
#' The default labels each categorization function based on whether it is showing prior or posterior categorization,
#' and by its group ID.
#' @param group.colors Vector of fill colors of same length as group.ids or `NULL` to use defaults. (default: `NULL`)
#' @param panel_scaling Should the relative scaling be calculated separately for each panel? If not the scaling is calculated
#' globally. (default: `FALSE`)
#'
#' @return ggplot object.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
#'
plot_ibbu_stanfit_parameters = function(
  model,
  which = "both",
  ndraws = NULL,
  untransform_cues = TRUE,
  group.ids = NULL, group.labels = NULL, group.colors = NULL,
  panel_scaling = F
) {
  d.pars = model %>%
    add_ibbu_stanfit_draws(
      which = which,
      ndraws = ndraws,
      untransform_cues = untransform_cues,
      nest = F)

  if (is.null(group.ids)) group.ids = levels(d.pars$group)
  # Setting aes defaults
  if (which  == "prior") {
    if(is.null(group.labels)) group.labels[1] = "prior"
    if(is.null(group.colors)) group.colors[1] = "darkgray"
  } else if (which == "posterior") {
    if(is.null(group.labels)) group.labels = paste0("posterior (", group.ids, ")")
    if(is.null(group.colors)) group.colors = get_default_colors("group", length(group.ids))
  } else {
    if(is.null(group.labels)) group.labels = ifelse(group.ids == "prior", group.ids, paste0("posterior (", group.ids, ")"))
    if(is.null(group.colors)) group.colors = c("darkgray", get_default_colors("group", length(group.ids) - 1))
  }

  p.m = d.pars %>%
    select(.draw, group, category, cue, m) %>%
    distinct() %T>%
    { get_limits(., "m") ->> x.limits } %>%
    ggplot(aes(
      y = fct_rev(.data$category),
      x = .data$m,
      fill = .data$group)) +
    ggridges::geom_density_ridges(alpha = .5, color = NA,
                                  panel_scaling = panel_scaling, scale = .95,
                                  stat = "density", aes(height = ..density..),
                                  # trim in order to increase resolution and avoid misleading
                                  # overlap with zero for S matrix; if not trimmed, density range
                                  # is estimated for the entire data in each plot
                                  trim = T) +
    geom_vline(xintercept = 0, color = "darkgray") +
    scale_x_continuous("Mean of category means") +
    scale_y_discrete("Category", expand = c(0,0)) +
    scale_fill_manual(
      "Group",
      breaks = group.ids,
      labels = group.labels,
      values = group.colors
    ) +
    coord_cartesian(xlim = x.limits, default = T) +
    facet_grid(~ .data$cue) +
    theme(legend.position = "right")
  legend = cowplot::get_legend(p.m)

  p.m = p.m + theme(legend.position = "none")
  p.S = suppressMessages(
    p.m %+%
      (d.pars %>%
         select(.draw, group, category, cue, cue2, S) %>%
         distinct() %T>%
         { get_limits(., "S") ->> x.limits }) +
      aes(x = .data$S) +
      scale_x_continuous("Scatter matrix",
                         breaks = inv_symlog(
                           seq(
                             ceiling(symlog(min(x.limits))),
                             floor(symlog(max(x.limits))),
                           ))
      ) +
      coord_trans(x = "symlog", xlim = x.limits) +
      facet_grid(.data$cue2 ~ .data$cue))

  p.KN = suppressWarnings(
    suppressMessages(
      p.m %+%
        (d.pars %>%
           select(.draw, group, category, kappa, nu) %>%
           distinct() %>%
           gather(key = "key", value = "value", -c(.draw, group, category)) %T>%
           { get_limits(., "value", min = 1) ->> x.limits } ) +
        aes(x = .data$value) +
        scale_x_continuous("Pseudocounts",
                           breaks = 10^(
                             seq(
                               ceiling(log10(min(x.limits))),
                               floor(log10(max(x.limits)))
                             ))) +
        scale_y_discrete("", expand = c(0,0)) +
        coord_trans(x = "log10", xlim = x.limits) +
        facet_grid(~ .data$key)))

  p.LR =
    d.pars %>%
    select(.draw, lapse_rate) %>%
    distinct() %T>%
    { get_limits(., "lapse_rate") ->> x.limits } %>%
    ggplot(aes(x = .data$lapse_rate)) +
    geom_density(color = NA, fill = "darkgray", alpha = .5,
                 stat = "density") +
    scale_x_continuous("Lapse rate")  +
    scale_y_discrete("") +
    scale_fill_manual(
      "Group",
      breaks = group.ids,
      labels = group.labels,
      values = group.colors
    ) +
    coord_cartesian(xlim = x.limits) +
    theme(legend.position = "none")

  K = length(unique(d.pars$cue))
  p = suppressWarnings(cowplot::plot_grid(
    cowplot::plot_grid(plotlist = list(p.m, p.KN), nrow = 1, rel_widths = c(K,2)),
    cowplot::plot_grid(plotlist = list(
      p.S,
      cowplot::plot_grid(plotlist = list(legend, p.LR),
                         nrow = 2, rel_heights = c(.5, .5))),
      nrow = 1, rel_widths = c(K,1)),
    rel_heights = c(1.5, K), nrow = 2, axis = "lrtb"))
  return(p)
}





#' Plot expected bivariate (2D) categories from MV IBBU stanfit.
#'
#' Plot bivariate Gaussian categories expected given the parameters inferred by incremental Bayesian belief-
#' updating (IBBU). Specifically, the categories are derived by marginalizing over the uncertainty represented
#' by the (post-warmup) MCMC samples. Two methods are available (specified by `type`), which differ in their
#' computational demands and speed.
#'
#' @param model mv-ibbu-stanfit object.
#' @param type Either `"contour"` or `"density"`, specifying the type of plot. Note that the contour plot is *much*
#' faster. It simply gets the expected values of \code{mu} (based on the NIW parameter \code{m}) and \code{Sigma}
#' (based on the NIW parameters \code{S} and \code{nu}) at each MCMC draw, and then averages over
#' all MCMC draws. The plotted categories represent those means of the expected \code{mu} and \code{Sigma}. The
#' density plot instead calculates the posterior predictive for each MCMC draw (i.e, the multivariate Student-T
#' density based on the NIW parameters \code{m, S, kappa, nu}), and then averages those densities. Since this is
#' done for *all* points defined by the data.grid this can be rather computationally expensive and slow.
#' @param plot.test,plot.exposure Should the test and/or exposure stimuli be plotted? (default: `TRUE` for `plot.test`,
#' `FALSE` for `plot.exposure`) The test items are plotted as black points. The exposure mean is plotted as point,
#' and the .95 interval of cue distributions during exposure are plotted as dashed ellipse in the same color as the
#' expected categories.
#' @param summarize Should one expected categories be plotted, marginalizing over MCMC draws (`TRUE`), or should separate
#' expected categories be plotted for each MCMC draw (`FALSE`)? (default: `TRUE`) Currently being ignored.
#' @param ndraws Number of draws to plot (or use to calculate the CIs), or `NULL` if all draws are to be returned.
#' (default: `NULL`) Currently being ignored.
#' @param untransform_cues Should m_0 and S_0 be transformed back into the original cue space? (default: `TRUE`)
#' @param levels Used only if `type` is `"contour"`. levels The cumulative probability levels that should be plotted (using
#' `geom_polygon()`) around the mean. By default the most transparent ellipse still drawn corresponds to .95.
#' @param xlim,ylim For density plots. Limits for the x- and y-axis.
#' @param resolution For density plots. How many steps along x and y should be calculated? Note that computational
#' complexity increases quadratically with resolution. (default: 25)
#' @param category.ids Vector of category IDs to be plotted or leave `NULL` to plot all groups. (default: `NULL`) It is possible
#' to use \code{\link[tidybayes]{recover_types}} on the stanfit object prior to handing it to this plotting function.
#' @param category.labels Vector of group labels of same length as `category.ids` or `NULL` to use defaults. (default: `NULL`)
#' @param category.colors Vector of colors of same length as category.ids or `NULL` to use defaults. (default: `NULL`)
#' @param category.linetypes Vector of linetypes of same length as category.ids or `NULL` to use defaults. (default: `NULL`)
#' Currently being ignored.
#' @param data.grid.xlim,data.grid.ylim,data.grid.resolution Used only if `type` is `"density"`. Limits for x- and y-axis as
#' well as resolution of the data.grid, defining the range over which the posterior predictive (multivariate Student-T density)
#' is calculated. Note that the number of densities to calculate is a *quadratic* function of `data.grid.resolution`. The default
#' for `data.grid.resolution` is 10, corresponding to 100 densities to be calculated for each MCMC draw.
#'
#'
#' @return ggplot object.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @rdname plot_expected_ibbu_stanfit_categories_2D
#' @export
plot_expected_ibbu_stanfit_categories_2D = function(
  model,
  type,
  ...
) {
  assert_that(all(type %in% c("contour", "density"), length(type) == 1))
  if (type == "contour")
    plot_expected_ibbu_stanfit_categories_contour2D(model = model, ...)
  else {
    plot_expected_ibbu_stanfit_categories_density2D(model = model, ...)
  }
}

#' @rdname plot_expected_ibbu_stanfit_categories_2D
#' @export
plot_expected_ibbu_stanfit_categories_contour2D = function(
  model,
  levels = plogis(seq(-15, qlogis(.95), length.out = 20)),
  plot.test = T, plot.exposure = F,
  untransform_cues = TRUE,
  category.ids = NULL, category.labels = NULL, category.colors = NULL, category.linetypes = NULL
) {
  fit.input = get_input_from_stanfit(model)
  assert_that(!all(is.null(fit.input), plot.test))
  d <- get_expected_category_statistic_from_stanfit(model, untransform_cues = untransform_cues)

  # Setting aes defaults
  if(is.null(category.ids)) category.ids = levels(d$category)
  if(is.null(category.labels)) category.labels = levels(d$category)
  if(is.null(category.colors)) category.colors = get_default_colors("category", length(category.ids))
  if(is.null(category.linetypes)) category.linetypes = rep(1, length(category.ids))

  d %<>%
    rename(x = Sigma.mean, centre = mu.mean) %>%
    crossing(level = levels) %>%
    mutate(ellipse = pmap(., ellipse.pmap)) %>%
    # This step is necessary since unnest() can't yet unnest lists of matrices
    # (bug was reported and added as milestone, 11/2019)
    mutate(ellipse = map(ellipse, as_tibble)) %>%
    unnest(ellipse)

  cue.names = setdiff(names(d), c("group", "category", "centre", "x", "level"))
  d %<>%
    rename_at(cue.names,
              function(x) paste0("cue", which(x == cue.names)))
  message("If cues labels can be extracted from stanfit, this code can be improved, and .data[[cue.labels[1]]], etc. can be used.
          See plot_expected_categories for NIW_belief objects.")

  ggplot(d,
         aes(x = .data$cue1,
             y = .data$cue2,
             fill = .data$category,
             alpha = 1 - .data$level,
             group = paste(.data$category, .data$level))) +
    geom_polygon() +
    # Optionally plot test data
    { if (plot.test)
      geom_point(
        data = fit.input$x_test %>%
          rename_at(cue.names,
                    function(x) paste0("cue", which(x == cue.names))),
        mapping = aes(x = .data$cue1, y = .data$cue2),
        inherit.aes = F,
        color = "black", size = 1
    )} +
    # Optionally plot exposure data
    { if (plot.exposure)
      geom_point(
        data =
          get_ibbu_stanfit_exposure_mean(
            fit.input,
            category = levels(d$category),
            group = levels(d$group)) %>%
          rename_at(cue.names,
                    function(x) paste0("cue", which(x == cue.names))),
        mapping = aes(
          x = .data$cue1,
          y = .data$cue2,
          shape = .data$category,
          color = .data$category),
        inherit.aes = F, size = 2
      ) +
        geom_path(
          data = crossing(
            group = levels(d$group),
            category = levels(d$category),
            level = .95
          ) %>%
            mutate(
              x = map2(category, group, get_ibbu_stanfit_exposure_sigma(fit.input, .x, .y)),
              centre = map2(category, group, get_ibbu_stanfit_exposure_mean(fit.input, .x, .y))
            ) %>%
            mutate(ellipse = pmap(., ellipse.pmap)) %>%
            mutate(ellipse = map(ellipse, as_tibble)) %>%
            unnest(ellipse) %>%
            rename_at(cue.names,
                      function(x) paste0("cue", which(x == cue.names))),
          mapping = aes(
            x = .data$cue1,
            y = .data$cue2,
            shape = .data$category,
            color = .data$category),
          linetype = 2,
          inherit.aes = F)
    } +
    scale_x_continuous(cue.names[1]) +
    scale_y_continuous(cue.names[2]) +
    scale_fill_manual("Category",
                      breaks = category.ids,
                      labels = category.labels,
                      values = category.colors) +
    scale_alpha("",
                range = c(0.1,.9)) +
    facet_wrap(~ group)
}


#' @rdname plot_expected_ibbu_stanfit_categories_2D
#' @export
plot_expected_ibbu_stanfit_categories_density2D = function(
  model,
  plot.test = T, plot.exposure = F,
  untransform_cues = TRUE,
  category.ids = NULL, category.labels = NULL, category.colors = NULL, category.linetypes = NULL,
  xlim, ylim, resolution = 25
) {
  fit.input = get_input_from_stanfit(model)
  assert_that(is.NIW_ideal_adaptor_stanfit(model) | is.NIW_ideal_adaptor_MCMC(model))
  assert_that(!all(is.null(fit.input), plot.test))

  if (is.NIW_ideal_adaptor_stanfit(model))
    d = add_ibbu_stanfit_draws(model, which = which, wide = F, nest = T, untransform_cues = untransform_cues)
  else
    d = model

  # Setting aes defaults
  if(is.null(category.ids)) category.ids = levels(d$category)
  if(is.null(category.labels)) category.labels = levels(d$category)
  if(is.null(category.colors)) category.colors = get_default_colors("category", length(category.ids))
  if(is.null(category.linetypes)) category.linetypes = rep(1, length(category.ids))

  cue.names = row.names(d$m[[1]])
  d %<>%
    crossing(
      cue1 = seq(min(xlim), max(xlim), length.out = resolution),
      cue2 = seq(min(ylim), max(ylim), length.out = resolution)) %>%
    mutate(x = map2(cue1, cue2, ~ c(.x, .y))) %>%
    mutate(
      density = pmap(., get_NIW_posterior_predictive.pmap),
      density = unlist(density)
    ) %>%
    # Marginalize over MCMC draws
    group_by(group, category, cue1, cue2) %>%
    summarise(density = mean(density))

  ggplot(d,
         aes(x = .data$cue1,
             y = .data$cue2,
             color = .data$category,
             fill = .data$category,
             z = .data$density)) +
    geom_contour() +
    # Optionally plot test data
    { if (plot.test)
      geom_point(
        data = fit.input$x_test %>%
          rename_at(cue.names,
                    function(x) paste0("cue", which(x == cue.names))),
        mapping = aes(
          x = .data$cue1,
          y = .data$cue2),
        inherit.aes = F,
        color = "black", size = 1
      )} +
    # Optionally plot exposure data
    { if (plot.exposure)
      geom_point(
        data = get_ibbu_stanfit_exposure_mean(fit.input,
                                              category = levels(d$category),
                                              group = levels(d$group)) %>%
          rename_at(cue.names,
                    function(x) paste0("cue", which(x == cue.names))),
        mapping = aes(
          x = .data$cue1,
          y = .data$cue2,
          shape = .data$category,
          color = .data$category),
        inherit.aes = F, size = 2) +
        geom_path(
          data = crossing(
            group = levels(d$group),
            category = levels(d$category),
            level = .95) %>%
            mutate(
              x = map2(category, group, get_ibbu_stanfit_exposure_sigma(fit.input, .x, .y)),
              centre = map2(category, group, get_ibbu_stanfit_exposure_mean(fit.input, .x, .y))) %>%
            mutate(ellipse = pmap(., ellipse.pmap)) %>%
            mutate(ellipse = map(ellipse, as_tibble)) %>%
            unnest(ellipse) %>%
            rename_at(cue.names,
                      function(x) paste0("cue", which(x == cue.names))),
          mapping = aes(
            x = .data$cue1,
            y = .data$cue2,
            shape = .data$category,
            color = .data$category),
          linetype = 2,
          inherit.aes = F)
    } +
    scale_x_continuous(cue.names[1]) +
    scale_y_continuous(cue.names[2]) +
    scale_color_manual("Category",
                       breaks = category.ids,
                       labels = category.labels,
                       values = category.colors) +
    coord_fixed(xlim = xlim, ylim = ylim, ratio = 1) +
    facet_wrap(~ .data$group)
}





#' Plot prior and posterior categorization of test tokens.
#'
#' Categorize test tokens by prior and/or posterior beliefs, and plot the resulting categorization function along
#' a one-dimensional continuum (regardless of the dimensionality of the cue space in which categorization takes
#' place). This provides the type of categorization plot typical for, for example, perceptual recalibration or
#' phonetic tuning studies.
#'
#' Tokens are sorted based on the increasing probability of a \code{target_category} response for the condition
#' (\code{group}, e.g., prior or a specific exposure group) specified in \code{sort.by}. By default both the mean
#' categorization and confidence intervals are plotted.
#' If `summarize=TRUE`, the function marginalizes over all posterior samples. The number of samples
#' is determined by ndraws. If ndraws is NULL, all samples are used. Otherwise ndraws random
#' samples will be used. If `summarize=FALSE`, separate categorization plots for all ndraws
#' individual samples will be plotted in separate panels.
#'
#' @param model \code{\link{NIW_ideal_adaptor_stanfit}} object.
#' @param data.test Optionally, a \code{tibble} or \code{data.frame} with test data.
#' If `NULL` the input will be extracted from fit. (default: `NULL`).
#' @param which Should categorization for the prior, posterior, or both be plotted? (default: `"both"`)
#' @param summarize Should one categorization function (optionally with CIs) be plotted (`TRUE`) or should separate
#' unique categorization function be plotted for each MCMC draw (`FALSE`)? (default: `TRUE`)
#' @param ndraws Number of draws to plot (or use to calculate the CIs), or `NULL` if all draws are to be returned.
#' (default: `NULL`)
#' @param confidence.intervals The two confidence intervals that should be plotted (using `geom_ribbon`) around the mean.
#' (default: `c(.66, .95)`)
#' @param target_category The index of the category for which categorization should be shown. (default: `1`)
#' @param panel.group Should the groups be plotted in separate panels? (default: `FALSE`)
#' @param group.ids Vector of group IDs to be plotted or leave `NULL` to plot all groups. (default: `NULL`)
#' @param group.labels Vector of group labels of same length as `group.ids` or `NULL` to use defaults. (default: `NULL`)
#' The default labels each categorization function based on whether it is showing prior or posterior categorization,
#' and by its group ID.
#' @param group.colors Vector of colors of same length as group.ids or `NULL` to use defaults. (default: `NULL`)
#' @param group.linetypes Vector of linetypes of same length as group.ids or `NULL` to use defaults. (default: `NULL`)
#' @param sort_by Which group, if any, should the x-axis be sorted by (in increasing order of posterior probability
#' from left to right). Set to 0 for sorting by prior (default). Set to `NULL` if no sorting is desired.
#'
#' @return ggplot object.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
plot_ibbu_stanfit_test_categorization = function(
  model,
  data.test = NULL,
  which = "both",
  summarize = T,
  ndraws = NULL,
  untransform_cues = TRUE,
  confidence.intervals = c(.66, .95),
  target_category = 1,
  panel.group = FALSE,
  group.ids = NULL, group.labels = NULL, group.colors = NULL, group.linetypes = NULL,
  sort.by = "prior"
) {
  assert_NIW_ideal_adaptor_stanfit(model)
  if (is.null(data.test)) {
    data.test = get_test_data_from_stanfit(model)
  }
  assert_that(is.flag(summarize))
  assert_that(is.null(confidence.intervals) |
                all(is.numeric(confidence.intervals),
                    length(confidence.intervals) == 2,
                    all(between(confidence.intervals, 0, 1))),
              msg = "Confidence intervals must be NULL (if not CIs are desired) or a vector of two probabilities.")
  assert_that(is.null(group.labels) | is.character(group.labels))
  assert_that(is.null(group.linetypes) | is.numeric(group.linetypes))
  assert_that(is.null(sort.by) | length(sort.by) == 1)


  # Set confidence intervals
  if (!is.null(confidence.intervals)) {
    ci.offset = (1 - confidence.intervals) / 2
    confidence.intervals = c(ci.offset, 1-ci.offset)
  }
  confidence.intervals = sort(confidence.intervals)

  # Get prior and posterior parameters
  d.pars =
    add_ibbu_stanfit_draws(
      model,
      which = which,
      summarize = F,
      wide = F,
      ndraws = ndraws,
      untransform_cues = untransform_cues)

  # Now set ndraws to the number of MCMC samples
  ndraws = if (is.null(ndraws)) get_number_of_draws(model) else ndraws
  if (ndraws > 500)
    message(paste("Marginalizing over", ndraws, "MCMC samples. This might take some time.\n"))

  # If group.ids are NULL set them to the levels of groups found in the extraction
  # of posteriors from model
  if (is.null(group.ids))  group.ids = levels(d.pars$group)
  assert_that(all(group.ids %in% unique(d.pars$group)),
              msg = "Some group.ids were not found in the stanfit object.")
  assert_that(sort.by %in% group.ids,
              msg = "Sort.by must be NULL or one of the group IDs (group.ids).")

  # Setting aes defaults
  if(is.null(group.labels)) group.labels = paste0("posterior (", group.ids[-1], ")")
  if(is.null(group.colors)) group.colors = get_default_colors("group", length(group.ids) - 1)
  if(is.null(group.linetypes)) group.linetypes = rep(1, length(group.ids) - 1)
  # If no specific color for prior was specified
  if(length(group.labels) < length(group.ids)) group.labels = c("prior", group.labels)
  if(length(group.colors) < length(group.ids)) group.colors = c("darkgray", group.colors)
  if(length(group.linetypes) < length(group.ids))  group.linetypes = c(3, group.linetypes)

  # Exclude all groups that are not in group.ids
  d.pars %<>%
    filter(group %in% group.ids)

  # Prepare test_data
  message("Using IBBU stanfit input to extract information about test data.")
  cue.labels = get_cue_labels_from_model(d.pars)
  test_data = data.test %>%
    distinct() %>%
    { if (untransform_cues) get_untransform_function_from_stanfit(model)(.) else . } %>%
    # CHECK: Could be replaced by make_vector_column
    transmute(x = pmap(.l = list(!!! syms(cue.labels)), .f = ~ c(...))) %>%
    nest(cues = x) %>%
    crossing(group = levels(d.pars$group))

  d.pars %<>%
    # Write a categorization function for each draw
    group_by(group, .draw) %>%
    do(f = get_categorization_function_from_grouped_ibbu_stanfit_draws(., logit = T)) %>%
    left_join(test_data) %>%
    group_by(group, .draw) %>%
    mutate(
      p_cat = invoke_map(.f = f, .x = cues, target_category = target_category),
      f = NULL) %>%
    unnest(c(cues, p_cat))

  if (summarize) {
    d.pars %<>%
      # For each unique group and test token obtain the CIs and the mean.
      group_by(group, x) %>%
      summarise_at(
        "p_cat",
        .funs = list(
          # na.rm = T excludes cases that might result from estimated probabilities of 0 and 1 (infinities in log-odds)
          y.outer.min = function(x) plogis(quantile(x, confidence.intervals[1], na.rm = T)),
          y.outer.max = function(x) plogis(quantile(x, confidence.intervals[4], na.rm = T)),
          y.inner.min = function(x) plogis(quantile(x, confidence.intervals[2], na.rm = T)),
          y.inner.max = function(x) plogis(quantile(x, confidence.intervals[3], na.rm = T)),
          p_cat = function(x) plogis(mean(x, na.rm = T)))
      )
  } else {
    d.pars %<>%
      mutate(p_cat = plogis(p_cat))
  }

  d.pars %<>%
    # Get cues as character strings (just in case)
    mutate(
      token.cues = map(x, ~paste(.x, collapse = ",\n")),
      x = NULL) %>%
    ungroup()

  # If sort.by is specified, sort levels of x-axis by that group.
  if (!is.null(sort.by)) {
    sort.levels = d.pars %>%
      filter(group == sort.by) %>%
      group_by(token.cues) %>%
      summarise(p_cat = mean(p_cat)) %>%
      arrange(p_cat) %>% pull(token.cues)

    d.pars %<>%
      mutate(
        token.cues = factor(token.cues,
                            levels = sort.levels),
        token = factor(as.numeric(token.cues), levels = 1:length(levels(token.cues))))
  }

  if (is.null(get_category_levels(model)))
    category1 = "category 1" else category1 = get_category_levels(model, 1)

  p = d.pars %>%
    ungroup() %>%
    mutate(group = factor(group, levels = group.ids)) %>%
    ggplot(aes(
      x = .data$token,
      y = .data$p_cat,
      color = .data$group,
      linetype = .data$group)) +
    scale_x_discrete("Test token",
                     breaks = levels(.data$token),
                     labels = paste0(levels(.data$token), "\n",
                                     levels(.data$token.cues))) +
    scale_y_continuous(
      paste0("Predicted proportion of ", category1, "-responses"),
      limits = c(0,1)
    ) +
    scale_color_manual(
      "Group",
      breaks = group.ids,
      labels = group.labels,
      values = group.colors
    ) +
    scale_linetype_manual(
      "Group",
      breaks = group.ids,
      labels = group.labels,
      values = group.linetypes
    )

  if (summarize & !is.null(confidence.intervals)) {
    p = p +
      geom_ribbon(
        aes(
          x = as.numeric(.data$token),
          ymin = .data$y.outer.min,
          ymax = .data$y.outer.max,
          shape = .data$group,
          fill = .data$group),
        color = NA, alpha = .1
      ) +
      geom_ribbon(
        aes(
          x = as.numeric(.data$token),
          ymin = .data$y.inner.min,
          ymax = .data$y.inner.max,
          fill = .data$group),
        color = NA, alpha = .3
      ) +
      scale_fill_manual(
        "Group",
        breaks = group.ids,
        labels = group.labels,
        values = group.colors
      ) +
      scale_shape_discrete(
        "Group",
        breaks = group.ids,
        labels = group.labels
      )

    # Place information about confidence intervals on plot.
    p = p +
      ggtitle(paste0((confidence.intervals[4]-confidence.intervals[1]) * 100,
                     "% and ",
                     (confidence.intervals[3]-confidence.intervals[2]) * 100,
                     "% CIs\nbased on ",
                     ndraws,
                     " posterior samples.")
      )
  }

  p = p +
    geom_point(alpha = .9) +
    geom_line(size = 1, alpha = .9, aes(x = as.numeric(.data$token)))

  if (!summarize & panel.group) p = p + facet_grid(group ~ .draw) else
    if (panel.group) p = p + facet_wrap(~ group) else
      if (!summarize) p = p + facet_wrap(~ .draw)

  return(p)
}
