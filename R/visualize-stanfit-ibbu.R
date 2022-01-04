#' @import ggplot2 cowplot gganimate transformr av
#' @importFrom ellipse ellipse
#' @importFrom mvtnorm dmvt
#' @importFrom tidybayes mean_hdi
#' @importFrom ggridges geom_density_ridges
#' @importFrom ggforce facet_matrix geom_autodensity
#' @importFrom forcats fct_rev
NULL


#' Plot distribution of IBBU parameters.
#'
#' Plot distribution of post-warmup MCMC samples for all parameters representing the
#' prior and/or posterior beliefs.
#'
#' @param model mv-ibbu-stanfit object.
#' @param categories,groups,cues Character vector of categories/groups/cues to be plotted. Typically, the levels of these factors
#' are automatically added to the fit during the creation of the fit. If necessary, however, it is possible to use
#' \code{\link[tidybayes]{recover_types}} on the stanfit object to add or change these levels later.
#' (default: all categories/groups/cues will be plotted)
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
#'
#' @rdname plot_ibbu_stanfit_parameters
#' @export
plot_ibbu_stanfit_parameters = function(
  model,
  categories = get_category_levels_from_stanfit(model),
  groups = get_group_levels_from_stanfit(model, include_prior = T),
  cues = get_cue_levels_from_stanfit(model),
  which = "both",
  ndraws = NULL,
  untransform_cues = TRUE,
  group.ids = NULL, group.labels = NULL, group.colors = NULL,
  panel_scaling = F
) {
  d.pars <-
    model %>%
    add_ibbu_stanfit_draws(
      which = which,
      ndraws = ndraws,
      untransform_cues = untransform_cues,
      nest = F)

  if (is.null(group.ids)) group.ids = levels(d.pars$group)
  # Setting aes defaults
  if (which == "prior") {
    if(is.null(group.labels)) group.labels[1] = "prior"
    if(is.null(group.colors)) group.colors[1] = "darkgray"
  } else if (which == "posterior") {
    if (is.null(group.labels)) group.labels = paste0("posterior (", group.ids, ")")
    if (is.null(group.colors)) group.colors = get_default_colors("group", group.ids)
  } else {
    if (is.null(group.labels)) group.labels = ifelse(group.ids == "prior", group.ids, paste0("posterior (", group.ids, ")"))
    if (is.null(group.colors)) group.colors = get_default_colors("group", group.ids)
  }

  p.m <- d.pars %>%
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
    # # Add empty data points to control scale: https://stackoverflow.com/questions/51735481/ggplot2-change-axis-limits-for-each-individual-facet-panel
    # geom_blank() +
    scale_x_continuous("Mean of category means") +
    scale_y_discrete("Category", expand = expansion(mult = c(0 , 0.1))) +
    scale_fill_manual(
      "Group",
      breaks = groups,
      values = group.colors) +
    coord_cartesian(default = T) +
    facet_grid(~ .data$cue, scales = "free_x") +
    theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))
  legend = cowplot::get_legend(p.m)

  p.m <- p.m + theme(legend.position = "none")
  p.S <- suppressMessages(
    p.m %+%
      (d.pars %>%
         select(.draw, group, category, cue, cue2, S) %>%
         distinct() %T>%
         { get_limits(., "S") ->> x.limits }) +
      aes(x = .data$S) +
      scale_x_continuous(
        "Scatter matrix",
        breaks = inv_symlog(
          seq(
            ceiling(symlog(min(x.limits))),
            floor(symlog(max(x.limits)))))) +
      coord_trans(x = "symlog") +
      facet_grid(.data$cue2 ~ .data$cue, scales = "free_x")) +
    theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

  p.KN <- suppressWarnings(
    suppressMessages(
      p.m %+%
        (d.pars %>%
           select(.draw, group, category, kappa, nu) %>%
           distinct() %>%
           gather(key = "key", value = "value", -c(.draw, group, category)) %T>%
           { get_limits(., "value", min = 1) ->> x.limits } ) +
        aes(x = .data$value) +
        scale_x_continuous(
          "Pseudocounts",
          breaks = 10^(
            seq(
              ceiling(log10(min(x.limits))),
              floor(log10(max(x.limits)))))) +
        scale_y_discrete("", expand = expansion(mult = c(0 , 0.1))) +
        coord_trans(x = "log10", xlim = x.limits) +  # <------------------------ xlim = x.limits for now put back in since plots was empty without it.
        facet_grid(~ .data$key, scales = "free_x"))) +
    theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

  p.LR <-
    d.pars %>%
    select(.draw, lapse_rate) %>%
    distinct() %T>%
    { get_limits(., "lapse_rate") ->> x.limits } %>%
    ggplot(aes(x = .data$lapse_rate)) +
    geom_density(color = NA, fill = "darkgray", alpha = .5,
                 stat = "density") +
    scale_x_continuous("Lapse rate")  +
    scale_y_discrete("", expand = expansion(mult = c(0 , 0.1))) +
    scale_fill_manual(
      "Group",
      breaks = groups,
      values = group.colors
    ) +
    coord_cartesian(xlim = x.limits) +
    theme(legend.position = "none")

  K <- length(unique(d.pars$cue))
  p <- suppressWarnings(cowplot::plot_grid(
    cowplot::plot_grid(plotlist = list(p.m, p.KN), nrow = 1, rel_widths = c(K, 2), axis = "btlr", align = "hv"),
    cowplot::plot_grid(plotlist = list(
      p.S,
      cowplot::plot_grid(plotlist = list(legend, p.LR),
                         nrow = 2, rel_heights = c(.5, .5))),
      nrow = 1, rel_widths = c(K,1)),
    rel_heights = c(1.5, K), nrow = 2, axis = "btlr", align = "hv"))

  return(p)
}




#' Plot correlations between IBBU parameters.
#'
#' Plot correlations between post-warmup MCMC samples for all parameters representing the prior and/or posterior beliefs.
#'
#' @param model mv-ibbu-stanfit object.
#' @param categories,groups,cues Character vector of categories/groups/cues to be plotted. Typically, the levels of these factors
#' are automatically added to the fit during the creation of the fit. If necessary, however, it is possible to use
#' \code{\link[tidybayes]{recover_types}} on the stanfit object to add or change these levels later.
#' (default: all categories and cues; `"prior"` for group since those are the only free parameters)
#' @param ndraws Number of draws to plot (or use to calculate the CIs), or `NULL` if all draws are to be returned. (default: `NULL`)
#' @param untransform_cues Should m_0 and S_0 be transformed back into the original cue space? (default: `TRUE`)
#' @param category.colors Vector of fill colors of same length as category or `NULL` to use defaults. (default: `NULL`)
#'
#' @return ggplot object.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#'
#' @export
plot_ibbu_stanfit_parameter_correlations = function(
  model,
  categories = get_category_levels_from_stanfit(model),
  groups = "prior",
  cues = get_cue_levels_from_stanfit(model),
  ndraws = NULL,
  untransform_cues = TRUE,
  category.colors = get_default_colors("category", categories)
) {
  d.pars <-
    model %>%
    add_ibbu_stanfit_draws(
      categories = categories,
      groups = groups,
      ndraws = ndraws,
      untransform_cues = untransform_cues,
      nest = T)

  d.pars %<>%
    mutate(
      S_tau = map(S, cov2tau),
      S_rho = map(S, cov2cor)) %>%
    select(-S) %>%
    unnest(c(m, S_tau, S_rho)) %>%
    group_by(across(-c(m, S_tau, S_rho))) %>%
    mutate(cue1 = .env$cues) %>%
    group_by(across(-c(S_rho))) %>%
    transmute(!! sym(cues[1]) := S_rho[,1], !! sym(cues[2]) := S_rho[,2]) %>%
    pivot_longer(cols = .env$cues, values_to = "S_rho", names_to = "cue2") %>%
    ungroup() %>%
    select(cue1, cue2, everything()) %>%
    filter(cue1 != cue2)

  # Removing redundant (duplicate) correlation information
  d.pars %<>%
    select(-c(cue2, S_rho)) %>%
    pivot_wider(names_from = "cue1", values_from = c("m", "S_tau")) %>%
    left_join(
      d.pars %>%
        group_by(.chain, .iteration, .draw, group, category) %>%
        mutate(combination = map2(.data$cue1, .data$cue2, ~paste(sort(c(.x, .y)), collapse = "_")) %>% unlist()) %>%
        distinct(combination, .keep_all = T) %>%
        select(-c(cue1, cue2, m, S_tau)) %>%
        pivot_wider(names_from = "combination", values_from = "S_rho", names_prefix = "S_rho_"),
      by = c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "lapse_rate"))

  p <- d.pars %>%
    ggplot(aes(x = .panel_x, y = .panel_y)) +
    geom_point(alpha = 0.2, shape = 16, size = 0.8, aes(color = category)) +
    geom_smooth(alpha = 0.2, aes(color = category)) +
    geom_autodensity(aes(fill = category), alpha = .5, position = position_identity()) +
    geom_density2d(aes(color = category), contour_var = "ndensity") +
    scale_color_manual("Category",
                      breaks = categories,
                      values = category.colors,
                      aesthetics = c("color", "fill")) +
    facet_matrix(
      vars(starts_with("kappa"), starts_with("nu"), starts_with("m_"), starts_with("S_")),
      layer.lower = c(1,2), layer.diag = 3, layer.upper = 4) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}






#' Plot expected bivariate (2D) categories from MV IBBU stanfit.
#'
#' Plot bivariate Gaussian categories expected given the parameters inferred by incremental Bayesian belief-
#' updating (IBBU). Specifically, the categories are derived by marginalizing over the uncertainty represented
#' by the (post-warmup) MCMC samples. Two methods are available (specified by `type`), which differ in their
#' computational demands and speed.
#'
#' @param model An mv-ibbu-stanfit object.
#' @param type Either `"contour"` or `"density"`, specifying the type of plot. Note that the contour plot is *much*
#' faster. It simply gets the expected values of \code{mu} (based on the NIW parameter \code{m}) and \code{Sigma}
#' (based on the NIW parameters \code{S} and \code{nu}) at each MCMC draw, and then averages over
#' all MCMC draws. The plotted categories represent those means of the expected \code{mu} and \code{Sigma}. The
#' density plot instead calculates the posterior predictive for each MCMC draw (i.e, the multivariate Student-T
#' density based on the NIW parameters \code{m, S, kappa, nu}), and then averages those densities. Since this is
#' done for *all* points defined by the data.grid this can be rather computationally expensive and slow.
#' @param categories,groups Character vector of categories/groups to be plotted. Typically, the levels of these factors
#' are automatically added to the fit during the creation of the fit. If necessary, however, it is possible to use
#' \code{\link[tidybayes]{recover_types}} on the stanfit object to add or change these levels later.
#' (default: all categories/groups will be plotted)
#' @param plot.test,plot.exposure Should the test and/or exposure stimuli be plotted? (default: `TRUE` for `plot.test`,
#' `FALSE` for `plot.exposure`) The test items are plotted as black points. The exposure mean is plotted as point,
#' and the .95 interval of cue distributions during exposure are plotted as dashed ellipse in the same color as the
#' expected categories.
#' @param annotate_inferred_category_means Character vector indicating whether the location and value of the mean be
#' indicated through data rugs (`"rug"`) and/or text labels (`"text"`)? Set to NULL to ignore. (default: `c("rug", "text")`)
#' @param untransform_cues Should m_0 and S_0 be transformed back into the original cue space? (default: `TRUE`)
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
  categories = get_category_levels_from_stanfit(model),
  groups = get_group_levels_from_stanfit(model, include_prior = T),
  plot.test = T,
  plot.exposure = F,
  annotate_inferred_category_means = c("rug", "text"),
  untransform_cues = TRUE,
  levels = plogis(seq(-15, qlogis(.95), length.out = 20)),
  category.colors = get_default_colors("category", categories)
) {
  assert_that(all(annotate_inferred_category_means %in% c("rug", "text")))

  d <-
    get_expected_category_statistic_from_stanfit(
      model,
      categories = categories,
      groups = groups,
      untransform_cues = untransform_cues)

  cue.names <- get_cue_levels_from_stanfit(model)
  d %<>%
    rename(x = Sigma.mean, centre = mu.mean) %>%
    crossing(level = levels) %>%
    mutate(ellipse = pmap(.l = list(x, centre, level), ellipse.pmap)) %>%
    unnest(ellipse) %>%
    group_by(across(-ellipse)) %>%
    transmute(cue1 = ellipse[,1], cue2 = ellipse[,2])

  min.cue1 <- min(d$cue1)
  min.cue2 <- min(d$cue2)
  d %>%
    ggplot(
    aes(x = .data$cue1,
        y = .data$cue2,
        fill = .data$category,
        alpha = 1 - .data$level,
        group = paste(.data$category, .data$level))) +
    geom_polygon() +
    # Optionally plot test data
    { if (plot.test)
      geom_point(
        data =
          get_test_data_from_stanfit(model) %>%
          { if (untransform_cues) get_untransform_function_from_stanfit(model)(.) else . } %>%
          rename_at(cue.names,
                    function(x) paste0("cue", which(x == cue.names))) %>%
          ungroup() %>%
          distinct(group, cue1, cue2),
        mapping = aes(x = .data$cue1, y = .data$cue2),
        inherit.aes = F,
        color = "black", size = 1, alpha = .75
    )} +
    # Optionally plot exposure data
    { if (plot.exposure)
      list(
        geom_point(
          data =
            get_exposure_mean_from_stanfit(
              model,
              categories = levels(d$category),
              groups = setdiff(levels(d$group), "prior"),
              untransform_cues = untransform_cues) %>%
            mutate(cue1 = unlist(map(.data$mean, ~ .x[1])), cue2 = unlist(map(.data$mean, ~ .x[2]))),
          mapping = aes(
            x = .data$cue1,
            y = .data$cue2,
            color = .data$category),
          inherit.aes = F, size = 2),
        geom_path(
          data =
            get_exposure_statistic_from_stanfit(
              model,
              categories = levels(d$category),
              groups = setdiff(levels(d$group), "prior"),
              untransform_cues = untransform_cues) %>%
            rename(x = cov, centre = mean) %>%
            crossing(level = .95) %>%
            mutate(ellipse = pmap(.l = list(x, centre, level), ellipse.pmap)) %>%
            unnest(ellipse) %>%
            group_by(across(-ellipse)) %>%
            transmute(cue1 = ellipse[,1], cue2 = ellipse[,2]),
          mapping = aes(
            x = .data$cue1,
            y = .data$cue2,
            color = .data$category),
          linetype = 2,
          inherit.aes = F)) } +
    { if ("rug" %in% annotate_inferred_category_means)
      geom_rug(
        data = . %>%
          distinct(group, category, centre),
        aes(
          x = map(centre, ~ .x[1]) %>% unlist(),
          y = map(centre, ~ .x[2]) %>% unlist(),
          color = category),
        inherit.aes = F) } +
    { if ("text" %in% annotate_inferred_category_means)
      list(
        geom_text(
          data = . %>%
            ungroup() %>%
            distinct(group, category, centre),
          aes(
            x = map(centre, ~ .x[1]) %>% unlist(),
            label= map(centre, ~ paste(signif(.x[1], 2))),
            color = category),
          y = min.cue2,
          angle = 90,
          hjust = 0,
          inherit.aes = F),
        geom_text(
          data = . %>%
            ungroup() %>%
            distinct(group, category, centre),
          aes(
            y = map(centre, ~ .x[2]) %>% unlist(),
            label= map(centre, ~ paste(signif(.x[2], 2))),
            color = category),
          x = min.cue1,
          angle = 0,
          hjust = 0,
          inherit.aes = F)) } +
    scale_x_continuous(cue.names[1]) +
    scale_y_continuous(cue.names[2]) +
    scale_fill_manual(
      "Category",
      breaks = categories,
      values = category.colors) +
    scale_color_manual(
      "Category",
      breaks = categories,
      values = lighten(category.colors, amount = .1)) +
    scale_alpha("", range = c(0.1,.9)) +
    facet_wrap(~ group)
}


#' @rdname plot_expected_ibbu_stanfit_categories_2D
#' @export
plot_expected_ibbu_stanfit_categories_density2D = function(
  model,
  categories = get_category_levels_from_stanfit(model),
  groups = get_group_levels_from_stanfit(model, include_prior = T),
  plot.test = T,
  plot.exposure = F,
  annotate_inferred_category_means = c("rug", "text"),
  untransform_cues = TRUE,
  category.colors = get_default_colors("category", categories),
  xlim, ylim, resolution = 25
) {
  assert_that(all(annotate_inferred_category_means %in% c("rug", "text")))

  d <-
    add_ibbu_stanfit_draws(
      model,
      categories = categories,
      groups = groups,
      wide = F,
      nest = T,
      untransform_cues = untransform_cues)

  cue.names = row.names(d$m[[1]])
  d %<>%
    crossing(
      cue1 = seq(min(xlim), max(xlim), length.out = resolution),
      cue2 = seq(min(ylim), max(ylim), length.out = resolution)) %>%
    mutate(x = map2(cue1, cue2, ~ c(.x, .y))) %>%
    mutate(
      density = pmap(., get_NIW_posterior_predictive.pmap),
      density = unlist(density)) %>%
    # Marginalize over MCMC draws
    group_by(group, category, cue1, cue2) %>%
    summarise(density = mean(density))

  min.cue1 <- min(d$cue1)
  min.cue2 <- min(d$cue2)
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
        data =
          get_test_data_from_stanfit(model) %>%
          { if (untransform_cues) get_untransform_function_from_stanfit(model)(.) else . } %>%
          rename_at(cue.names,
                    function(x) paste0("cue", which(x == cue.names))) %>%
          ungroup() %>%
          distinct(group, cue1, cue2),
        mapping = aes(
          x = .data$cue1,
          y = .data$cue2),
        inherit.aes = F,
        color = "black", size = 1
      )} +
    # Optionally plot exposure data
    { if (plot.exposure)
      list(
        geom_point(
        data =
          get_exposure_mean_from_stanfit(
            model,
            categories = levels(d$category),
            groups = setdiff(levels(d$group), "prior"),
            untransform_cues = untransform_cues) %>%
          mutate(cue1 = unlist(map(mean, ~.x[1])), cue2 = unlist(map(mean, ~.x[2]))),
        mapping = aes(
          x = .data$cue1,
          y = .data$cue2,
          color = .data$category),
        inherit.aes = F, size = 2),
      geom_path(
        data =
          crossing(
           category = levels(d$category),
           group = setdiff(levels(d$group), "prior"),
           level = .95) %>%
          mutate(
            x = map2(category, group, ~ get_exposure_ss_from_stanfit(model, categories = .x, groups = .y, untransform_cues = untransform_cues)),
            centre = map2(category, group, ~ get_exposure_mean_from_stanfit(model, categories = .x, groups = .y, untransform_cues = untransform_cues))) %>%
          mutate(ellipse = pmap(.l = list(x, centre, level), ellipse.pmap)) %>%
          unnest(ellipse) %>%
          group_by(across(-ellipse)) %>%
          transmute(cue1 = ellipse[,1], cue2 = ellipse[,2]),
        mapping = aes(
          x = .data$cue1,
          y = .data$cue2,
          color = .data$category),
        linetype = 2,
        inherit.aes = F)) } +
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
          y = min.cue2,
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
          x = min.cue1,
          angle = 0,
          hjust = 0,
          inherit.aes = F)) } +
    scale_x_continuous(cue.names[1]) +
    scale_y_continuous(cue.names[2]) +
    scale_fill_manual(
      "Category",
      breaks = categories,
      values = category.colors) +
    scale_color_manual(
      "Category",
      breaks = categories,
      values = lighten(category.colors, amount = .1)) +
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
#' (\code{group}, e.g., prior or a specific exposure group) specified in \code{sort_by}. By default both the mean
#' categorization and confidence intervals are plotted.
#' If `summarize=TRUE`, the function marginalizes over all posterior samples. The number of samples
#' is determined by ndraws. If ndraws is NULL, all samples are used. Otherwise ndraws random
#' samples will be used. If `summarize=FALSE`, separate categorization plots for all ndraws
#' individual samples will be plotted in separate panels.
#'
#' @param model \code{\link{NIW_ideal_adaptor_stanfit}} object.
#' @param data.test Optionally, a \code{tibble} or \code{data.frame} with test data.
#' If `NULL` the input will be extracted from fit. (default: `NULL`).
#' @param groups Character vector of groups to be plotted. Typically, the levels of these factors
#' are automatically added to the fit during the creation of the fit. If necessary, however, it is possible to use
#' \code{\link[tidybayes]{recover_types}} on the stanfit object to add or change these levels later.
#' (default: all categories/groups will be plotted)
#' @param summarize Should one categorization function (optionally with CIs) be plotted (`TRUE`) or should separate
#' unique categorization function be plotted for each MCMC draw (`FALSE`)? (default: `TRUE`)
#' @param ndraws Number of draws to plot (or use to calculate the CIs), or `NULL` if all draws are to be returned.
#' (default: `NULL`)
#' @param confidence.intervals The two confidence intervals that should be plotted (using `geom_ribbon`) around the mean.
#' (default: `c(.66, .95)`)
#' @param target_category The index of the category for which categorization should be shown. (default: `1`)
#' @param panel.group Should the groups be plotted in separate panels? (default: `FALSE`)
#' @param group.colors,group.linetypes Vector of colors and linetypes of same length as `groups` or `NULL` to use defaults.
#' @param category.colors Vector of colors and linetypes of same length as `categories` or `NULL` to use defaults. Only
#' relevant when `plot_in_cue_space = TRUE`.
#' @param all_test_locations Should predictions be shown for all combinations of test locations and group, or should only
#' combinations be shown that actually occurred in the data? (default: `FALSE`)
#' @param plot_in_cue_space Should predictions be plotted in the cue space? If not, test tokens are essentially treated
#' as factors and sorted along the x-axis based on `sort_by`. (default: `TRUE`)
#' @param sort_by Which group, if any, should the x-axis be sorted by (in increasing order of posterior probability
#' from left to right). Set to 0 for sorting by prior (default). Set to `NULL` if no sorting is desired. (default: `"prior"`)
#'
#' @return ggplot object.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @rdname plot_ibbu_stanfit_test_categorization
#' @export
plot_ibbu_stanfit_test_categorization = function(
  model,
  data.test = NULL,
  groups = get_group_levels_from_stanfit(model, include_prior = TRUE),
  summarize = T,
  ndraws = NULL,
  untransform_cues = TRUE,
  confidence.intervals = c(.66, .95),
  target_category = 1,
  panel.group = if (plot_in_cue_space) TRUE else FALSE,
  group.colors = get_default_colors("group", groups),
  group.linetypes = get_default_linetypes("group", groups),
  category.colors = get_default_colors("category", categories),
  all_test_locations = TRUE,
  plot_in_cue_space = FALSE,
  sort_by = if (plot_in_cue_space) NULL else "prior"
) {
  if (is.null(data.test)) data.test <- get_test_data_from_stanfit(model)
  assert_that(is.flag(summarize))
  assert_that(is.null(confidence.intervals) |
                all(is.numeric(confidence.intervals),
                    length(confidence.intervals) == 2,
                    all(between(confidence.intervals, 0, 1))),
              msg = "Confidence intervals must be NULL (if not CIs are desired) or a vector of two probabilities.")
  assert_that(is.null(sort_by) | length(sort_by) == 1)


  # Set confidence intervals
  if (!is.null(confidence.intervals)) {
    ci.offset <- (1 - confidence.intervals) / 2
    confidence.intervals <- c(ci.offset, 1-ci.offset)
  }
  confidence.intervals = sort(confidence.intervals)

  # Get prior and posterior parameters
  ndraws <- if (is.null(ndraws)) get_number_of_draws(model) else ndraws
  d.pars <-
    add_ibbu_stanfit_draws(
      model,
      groups = groups,
      summarize = F,
      wide = F,
      ndraws = ndraws,
      untransform_cues = untransform_cues)

  if (ndraws > 500)
    message(paste("Marginalizing over", ndraws, "MCMC samples. This might take some time.\n"))

  if (!is.null(sort_by))
    assert_that(sort_by %in% groups,
                msg = paste("sort_by must be NULL or one of the groups:",
                      paste(groups, collapse = ", ")))

  # Exclude all groups that are not in group.id
  d.pars %<>%
    filter(group %in% .env$groups)

  # Prepare test_data
  cue.labels <- get_cue_levels_from_stanfit(model)
  if (all_test_locations) {
    test_data <-
      data.test %>%
      distinct(!!! syms(cue.labels)) %>%
      { if (untransform_cues) get_untransform_function_from_stanfit(model)(.) else . } %>%
      make_vector_column(cols = cue.labels, vector_col = "x", .keep = "all") %>%
      nest(cues_joint = x, cues_separate = cue.labels) %>%
      crossing(group = levels(d.pars$group))
  } else {
    test_data <-
      data.test %>%
      distinct(!!! syms(cue.labels), group.id, group) %>%
      { if (untransform_cues) get_untransform_function_from_stanfit(model)(.) else . } %>%
      make_vector_column(cols = cue.labels, vector_col = "x", .keep = "all") %>%
      group_by(group.id, group) %>%
      nest(cues_joint = x, cues_separate = cue.labels)
  }

  d.pars %<>%
    group_by(group, .draw) %>%
    do(f = get_categorization_function_from_grouped_ibbu_stanfit_draws(., logit = T)) %>%
    right_join(test_data) %>%
    group_by(group, .draw) %>%
    mutate(p_cat = invoke_map(.f = f, .x = cues_joint, target_category = target_category)) %>%
    select(-f) %>%
    unnest(c(cues_joint, cues_separate, p_cat))

  if (summarize) {
    d.pars %<>%
      # For each unique group and test token obtain the CIs and the mean.
      group_by(group, x, !!! syms(cue.labels)) %>%
      summarise_at(
        "p_cat",
        .funs = list(
          # na.rm = T excludes cases that might result from estimated probabilities of 0 and 1 (infinities in log-odds)
          y.outer.min = function(x) plogis(quantile(x, confidence.intervals[1], na.rm = T)),
          y.outer.max = function(x) plogis(quantile(x, confidence.intervals[4], na.rm = T)),
          y.inner.min = function(x) plogis(quantile(x, confidence.intervals[2], na.rm = T)),
          y.inner.max = function(x) plogis(quantile(x, confidence.intervals[3], na.rm = T)),
          p_cat = function(x) plogis(mean(x, na.rm = T))))
  } else {
    d.pars %<>%
      mutate(p_cat = plogis(p_cat))
  }

  d.pars %<>%
    # Get cues as character strings (just in case)
    mutate(token.cues = map(x, ~paste(.x, collapse = ",\n"))) %>%
    ungroup() %>%
    select(-c(x))

  target_category_label <- if (is.null(get_category_levels_from_stanfit(model))) paste("category", target_category) else get_category_levels_from_stanfit(model, target_category)

  if (plot_in_cue_space) {
    if (length(get_cue_levels_from_stanfit(model)) > 2) stop("plot_in_cue_space = T not yet implemented for more than two cues.")

    p <-
      d.pars %>%
      mutate(group = factor(group, levels = .env$groups)) %>%
      ggplot(
        aes(
        x = !! sym(cue.labels[1]),
        y = !! sym(cue.labels[2]),
        fill = .data$p_cat)) +
      geom_raster() +
      scale_x_continuous(cue.labels[1]) +
      scale_y_continuous(cue.labels[2]) +
      scale_fill_gradient2(
        paste0("Predicted proportion\n of", target_category_label, "-responses"),
        midpoint = .5,
        high = category.colors[target_category],
        mid = "white",
        low = if (length(category.colors[which(get_category_levels_from_stanfit(model) != target_category_label)]) > 1) "gray" else category.colors[which(get_category_levels_from_stanfit(model) != target_category_label)],
        limits = c(0,1)) +
      coord_cartesian(expand = F)
  } else {
    # If sort_by is specified, sort levels of x-axis by that group.
    if (!is.null(sort_by)) {
      sort.levels <- d.pars %>%
        filter(group == sort_by) %>%
        group_by(token.cues) %>%
        summarise(p_cat = mean(p_cat)) %>%
        arrange(p_cat) %>% pull(token.cues)

      d.pars %<>%
        mutate(
          token.cues = factor(token.cues,
                              levels = sort.levels),
          token = factor(as.numeric(token.cues), levels = 1:length(levels(token.cues))))
    }

    p <-
      d.pars %>%
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
        paste0("Predicted proportion of ", target_category_label, "-responses"),
        limits = c(0,1)) +
      scale_color_manual(
        "Group",
        breaks = .env$groups,
        values = group.colors) +
      scale_linetype_manual(
        "Group",
        breaks = .env$groups,
        values = group.linetypes)

    if (summarize & !is.null(confidence.intervals)) {
      p <- p +
        geom_ribbon(
          aes(
            x = as.numeric(.data$token),
            ymin = .data$y.outer.min,
            ymax = .data$y.outer.max,
            fill = .data$group),
          color = NA, alpha = .1) +
        geom_ribbon(
          aes(
            x = as.numeric(.data$token),
            ymin = .data$y.inner.min,
            ymax = .data$y.inner.max,
            fill = .data$group),
          color = NA, alpha = .3) +
        scale_fill_manual(
          "Group",
          breaks = groups,
          values = group.colors)

      # Place information about confidence intervals on plot.
      p <- p +
        ggtitle(paste0((confidence.intervals[4]-confidence.intervals[1]) * 100,
                       "% and ",
                       (confidence.intervals[3]-confidence.intervals[2]) * 100,
                       "% CIs\nbased on ", ndraws, " posterior samples."))
    }

    p <- p +
      geom_point(alpha = .9) +
      geom_line(size = 1, alpha = .9, aes(x = as.numeric(.data$token)))
  }
  if (!summarize & panel.group) p <- p + facet_grid(group ~ .draw) else
    if (panel.group) p <- p + facet_wrap(~ group) else
      if (!summarize) p <- p + facet_wrap(~ .draw)

  return(p)
}
