#' @import ggplot2 cowplot gganimate ggnewscale transformr av
#' @importFrom colorspace lighten darken
#' @importFrom ellipse ellipse
#' @importFrom mvtnorm dmvt
#' @importFrom tidybayes mean_hdi
#' @importFrom ggridges geom_density_ridges
#' @importFrom ggforce facet_matrix geom_autodensity
#' @importFrom forcats fct_rev
NULL


#' Plot distribution of parameters of ideal adaptor stanfit
#'
#' Plot distribution of posterior MCMC samples for all parameters representing the
#' prior and/or posterior beliefs.
#'
#' @aliases plot_ibbu_stanfit_parameters
#'
#' @param model An ideal adaptor stanfit object.
#' @param categories,groups,cues Character vector of categories/groups/cues to be plotted. Typically, the levels of these factors
#' are automatically added to the fit during the creation of the fit. If necessary, however, it is possible to use
#' \code{\link[tidybayes]{recover_types}} on the stanfit object to add or change these levels later.
#' (default: all categories/groups/cues will be plotted)
#' @param which Should parameters for the prior, posterior, or both be plotted? (default: `"both"`)
#' @param ndraws Number of draws to plot (or use to calculate the CIs), or `NULL` if all draws are to be returned. (default: `NULL`)
#' @param untransform_cues Should m_0 and S_0 be transformed back into the original cue space? (default: `TRUE`)
#' @param group.colors Vector of fill colors of same length as group.ids or `NULL` to use defaults. (default: `NULL`)
#' @param panel_scaling Should the relative scaling be calculated separately for each panel? If not the scaling is calculated
#' globally. (default: `FALSE`)
#'
#' @return ggplot object.
#'
#' @seealso TBD
#' @keywords TBD
#'
#' @importFrom ggridges geom_density_ridges
plot_parameters <- function(fit, ...) {
  UseMethod("plot_parameters")
}

#' @rdname plot_parameters
#' @export
plot_parameters.NIW_ideal_adaptor_stanfit <- function(
  model,
  categories = get_category_levels(model),
  groups = get_group_levels(model, include_prior = T),
  cues = get_cue_levels(model),
  ndraws = NULL,
  untransform_cues = TRUE,
  panel_scaling = F,
  group.colors = get_default_colors("group", groups)
) {
  # Binding variables that RMD Check gets confused about otherwise
  # (since they are in non-standard evaluations)
  .draw <- group <- category <- cue <- cue2 <- kappa <- nu <- m <- S <- lapse_rate <- x.limits <- NULL

  d.pars <-
    model %>%
    get_draws(
      groups = groups,
      ndraws = ndraws,
      untransform_cues = untransform_cues,
      nest = F)

  p.m <-
    d.pars %>%
    select(.draw, group, category, cue, m) %>%
    distinct() %T>%
    { get_limits(., "m") ->> x.limits } %>%
    ggplot(aes(
      y = fct_rev(.data$category),
      x = .data$m,
      fill = .data$group)) +
    geom_density_ridges(alpha = .5, color = NA,
                                  panel_scaling = panel_scaling, scale = .95,
                                  stat = "density", aes(height = after_stat(density)),
                                  # trim in order to increase resolution and avoid misleading
                                  # overlap with zero for S matrix; if not trimmed, density range
                                  # is estimated for the entire data in each plot
                                  trim = T) +
    geom_vline(xintercept = 0, color = "darkgray") +
    # # Add empty data points to control scale: https://stackoverflow.com/questions/51735481/ggplot2-change-axis-limits-for-each-individual-facet-panel
    # geom_blank() +
    scale_x_continuous("Mean m of category means") +
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
        "Scatter matrix S",
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
        coord_trans(x = "log10", xlim = x.limits) +  # lim = x.limits for now put back in since plots was empty without it.
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
  p <- suppressWarnings(
    cowplot::plot_grid(
      cowplot::plot_grid(plotlist = list(p.m, p.KN), nrow = 1, rel_widths = c(K, 2), axis = "btlr", align = "hv"),
      cowplot::plot_grid(plotlist = list(
        p.S,
        cowplot::plot_grid(plotlist = list(legend, p.LR),
                           nrow = 2, rel_heights = c(.5, .5))),
        nrow = 1, rel_widths = c(K,1)),
      rel_heights = c(1.5, K), nrow = 2, axis = "btlr", align = "hv"))

  return(p)
}




#' Plot correlations between parameters of ideal adaptor stanfit
#'
#' Plot correlations between posterior MCMC samples for all parameters representing the prior and/or posterior beliefs.
#'
#' @aliases plot_ibbu_stanfit_parameter_correlations
#'
#' @param model An ideal adaptor stanfit object.
#' @param categories,groups,cues,pars Character vector of categories, groups, cues, and/or parameters to be plotted. To
#' subset to specific parameters, note the following naming conventions:
#'
#' (1) all `m` parameters are named `m_{cue_name}`
#' (2) the `S` parameter is decomposed into a vector of standard deviations `tau` and a correlation matrix `rho`.
#' (3) all `tau` parameters are named `tau_{cue_name}`
#' (4) all `rho` parameters are named `rho_{cue_name1}__x__{cue_name2}`
#'
#' (default: all categories, cues, and parameters in the model; `"prior"` for group since those are the only free parameters).
#' @param ndraws Number of draws to plot (or use to calculate the CIs), or `NULL` if all draws are to be returned. (default: `NULL`)
#' @param untransform_cues Should m_0 and S_0 be transformed back into the original cue space? (default: `TRUE`)
#' @param category.colors Vector of fill colors of same length as category or `NULL` to use defaults. (default: `NULL`)
#'
#' @return ggplot object.
#'
#' @details
#' Typically, the categories, groups, and cues are automatically added to the fit during the creation of the fit. If necessary,
#' however, it is possible to use \code{\link[tidybayes]{recover_types}} on the stanfit object to add or change these levels later.
#'
#'
#' @seealso TBD
#' @keywords TBD
#'
#' @importFrom ggforce facet_matrix geom_autodensity
#' @importFrom ggnewscale new_scale
#' @importFrom colorspace lighten
plot_parameter_correlations <- function(fit, ...) {
  UseMethod("plot_parameter_correlations")
}

#' @export plot_parameter_correlations
plot_parameter_correlations.NIW_ideal_adaptor_stanfit <- function(
  model,
  categories = get_category_levels(model),
  groups = "prior",
  cues = get_cue_levels(model),
  pars = NULL,
  ndraws = NULL,
  untransform_cues = TRUE,
  category.colors = get_default_colors("category", categories)
) {
  assert_that(is.null(pars) || is.character(pars))

  d.pars <-
    model %>%
    get_draws(
      categories = categories,
      groups = groups,
      ndraws = ndraws,
      untransform_cues = untransform_cues,
      nest = T) %>%
    group_by(across(-S)) %>%
    transmute(
      tau = map(S, cov2tau),
      Rho = map(S, cov2cor))

  if (length(cues) > 1) {
    d.Rho <- d.pars %>% ungroup %>% select(-c(m, tau))
    for (c in 2:length(cues)) {
      for (c2 in 1:(c - 1)) {
        d.Rho %<>%
          mutate(!!sym(paste0("rho_", cues[c], "__x__", cues[c2])) := map_dbl(Rho, ~ .x[c, c2]))
      }
    }
  }

  d.pars %<>%
    select(-Rho) %>%
    unnest(c(m, tau)) %>%
    group_by(across(-c(m, tau))) %>%
    mutate(cue1 = .env$cues) %>%
    pivot_wider(names_from = "cue1", values_from = c("m", "tau")) %>%
    # join in info about Rho if there is any
    { if (length(cues) > 1) {
      left_join(
        .,
        d.Rho %>% select(-Rho),
        by = join_by(.chain, .iteration, .draw, group, category, kappa, nu, lapse_rate))
    } else . } %>%
    ungroup()

  if (!is.null(pars)) {
    assert_that(
      all(pars %in% colnames(d.pars)),
      msg = paste0("The following parameter(s) could not be found in the stanfit object: ", paste(setdiff(pars, colnames(d.pars)), collapse = ", "), "."))

    d.pars %<>% select(.chain, .iteration, .draw, group, category, !!! syms(pars))
  }

  correlation_plot <-
    ggplot(data = NULL, aes(x = .panel_x, y = .panel_y)) +
    geom_density2d(aes(color = category), contour_var = "ndensity", alpha = .5) +
    geom_autodensity(aes(fill = category), alpha = .5, position = position_identity()) +
    geom_point(alpha = 0.1, shape = 16, size = 0.8, aes(color = category)) +
    scale_color_manual("Category",
                       breaks = categories,
                       values = category.colors,
                       aesthetics = c("color", "fill")) +
    new_scale_color() +
    geom_smooth(alpha = 0.5, aes(color = category)) +
    scale_color_manual("Category",
                       breaks = categories,
                       values = lighten(category.colors, amount = .5)) +
    facet_matrix(
      vars(starts_with("kappa"), starts_with("nu"), starts_with("m_"), starts_with("tau_"), starts_with("rho_")),
      layer.lower = c(3,4), layer.diag = 2, layer.upper = 1) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

  if (length(groups) > 1) {
    p <- list()
    for (g in groups) {
      p <- append(p, correlation_plot %+% (d.pars %>% filter(group == g)))
    }
    p <- plot_grid(plotlist = p)
  } else p <- correlation_plot %+% d.pars

  return(p)
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
#' @details
#' Typically, the categories, groups, and cues
#' are automatically added to the fit during the creation of the fit. If necessary, however, it is possible to use
#' \code{\link[tidybayes]{recover_types}} on the stanfit object to add or change these levels later.
#'
#'
#' @seealso TBD
#' @keywords TBD
#'
#' @importFrom purrr map_dbl
plot_expected_categories <- function(fit, ...) {
  UseMethod("plot_expected_categories")
}

plot_expected_categories_contour <- function(fit, ...) {
  UseMethod("plot_expected_categories_contour")
}

plot_expected_categories_density <- function(fit, ...) {
  UseMethod("plot_expected_categories_density")
}

plot_expected_categories_contour2D <- function(fit, ...) {
  UseMethod("plot_expected_categories_contour2D")
}

plot_expected_categories_density1D <- function(fit, ...) {
  UseMethod("plot_expected_categories_density1D")
}

plot_expected_categories_density2D <- function(fit, ...) {
  UseMethod("plot_expected_categories_density2D")
}

#' @rdname plot_expected_categories
#' @export
plot_expected_categories.NIW_ideal_adaptor_stanfit <- function(
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
plot_expected_categories_contour.NIW_ideal_adaptor_stanfit <- function(
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
plot_expected_categories_density.NIW_ideal_adaptor_stanfit <- function(
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

#' @aliases plot_expected_ibbu_stanfit_categories_contour2D
#' @rdname plot_expected_categories
#' @export
plot_expected_categories_contour2D.NIW_ideal_adaptor_stanfit <- function(
  model,
  categories = get_category_levels(model),
  groups = get_group_levels(model, include_prior = T),
  cues = get_cue_levels(model),
  plot.test = T,
  plot.exposure = F,
  annotate_inferred_category_means = c("rug", "text"),
  untransform_cues = TRUE,
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
        get_test_data(model, groups = setdiff(groups_found, "prior")) %>%
          { if (untransform_cues) get_untransform_function_from_stanfit(model)(.) else . } %>%
          ungroup() %>%
          distinct(group, !!! syms(cues)),
        cue.labels = cues) } +
    # Optionally plot exposure data
    { if (plot.exposure)
      add_exposure_summary_to_2D_plot(
        get_exposure_category_statistic.NIW_ideal_adaptor_stanfit(
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


#' @aliases plot_expected_ibbu_stanfit_categories_density1D
#' @rdname plot_expected_categories
#' @export
plot_expected_categories_density1D.NIW_ideal_adaptor_stanfit <- function(
    model,
    categories = get_category_levels(model),
    groups = get_group_levels(model, include_prior = T),
    cues = get_cue_levels(model),
    ndraws = NULL,
    plot.test = T,
    plot.exposure = F,
    annotate_inferred_category_means = c("rug", "text"),
    untransform_cues = TRUE,
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
      get_test_data(model, groups = setdiff(groups_found, "prior")) %>%
        pull(cues[1]) %>%
        { (range(.) - mean(.) * 1.5) + mean(.) }
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
        get_test_data(model, groups = setdiff(groups_found, "prior")) %>%
          { if (untransform_cues) get_untransform_function_from_stanfit(model)(.) else . } %>%
          ungroup() %>%
          distinct(group, !!! syms(cues)),
        cue.labels = cues) } +
    # Optionally plot exposure data
    { if (plot.exposure)
      add_exposure_summary_to_1D_plot(
        get_exposure_category_statistic.NIW_ideal_adaptor_stanfit(
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
plot_expected_categories_density2D.NIW_ideal_adaptor_stanfit <- function(
  model,
  categories = get_category_levels(model),
  groups = get_group_levels(model, include_prior = T),
  cues = get_cue_levels(model),
  ndraws = NULL,
  plot.test = T,
  plot.exposure = F,
  annotate_inferred_category_means = c("rug", "text"),
  untransform_cues = TRUE,
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
      get_test_data(model, groups = setdiff(groups_found, "prior")) %>%
        pull(cues[1]) %>%
        { (range(.) - mean(.) * 1.5) + mean(.) }
    } else xlim
  ylim <-
    if (is.null(ylim)) {
    get_test_data(model, groups = setdiff(groups_found, "prior")) %>%
      pull(cues[2]) %>%
      { (range(.) - mean(.) * 1.5) + mean(.) }
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
        get_test_data(model, groups = setdiff(groups_found, "prior")) %>%
          { if (untransform_cues) get_untransform_function_from_stanfit(model)(.) else . } %>%
          ungroup() %>%
          distinct(group, !!! syms(cues)),
        cue.labels = cues) } +
    # Optionally plot exposure data
    { if (plot.exposure)
      add_exposure_summary_to_2D_plot(
        get_exposure_category_statistic.NIW_ideal_adaptor_stanfit(
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

#' Plot prior and posterior categorization of test tokens for an NIW ideal adaptor stanfit
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
#' @aliases plot_ibbu_stanfit_test_categorization
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
#' @param group.colors,group.shapes,group.linetypes Vector of colors, shapes, and linetypes of same length as `groups` or `NULL` to use defaults.
#' @param category.colors Vector of colors and linetypes of same length as `categories` or `NULL` to use defaults. Only
#' relevant when `plot_in_cue_space = TRUE`.
#' @param all_test_locations Should predictions be shown for all combinations of test locations and group, or should only
#' combinations be shown that actually occurred in the data? (default: `FALSE`)
#' @param plot_in_cue_space Currently only available if the model has one or two cues. Should predictions be plotted in the cue space?
#' If not, test tokens are treated as factors and sorted along the x-axis based on `sort_by`. (default: `TRUE`)
#' @param untransform_cues Should the cues be untransformed before plotting? This should only have visual consequences
#' if `plot_in_cue_space = T`. (default: `TRUE`)
#' @param sort_by Which group, if any, should the x-axis be sorted by (in increasing order of posterior probability
#' from left to right). Set to 0 for sorting by prior (default). Set to `NULL` if no sorting is desired. (default: `"prior"`)
#'
#' @return ggplot object.
#'
#' @seealso TBD
#' @keywords TBD
#'
#' @importFrom dplyr do right_join
#' @importFrom purrr invoke_map
plot_expected_categorization <- function(fit, ...) {
  UseMethod("plot_expected_categorization")
}

#' @rdname plot_expected_categorization
#' @export
plot_expected_categorization.NIW_ideal_adaptor_stanfit <- function(
  model,
  data.test = NULL,
  groups = get_group_levels(model, include_prior = TRUE),
  summarize = T,
  ndraws = NULL,
  confidence.intervals = c(.66, .95),
  target_category = 1,
  panel.group = if (plot_in_cue_space) TRUE else FALSE,
  group.colors = get_default_colors("group", groups),
  group.shapes = get_default_shapes("group", groups),
  group.linetypes = get_default_linetypes("group", groups),
  category.colors = get_default_colors("category", get_category_levels(model)),
  all_test_locations = TRUE,
  plot_in_cue_space = FALSE,
  untransform_cues = TRUE,
  sort_by = if (plot_in_cue_space) NULL else "prior"
) {
  if (is.null(data.test)) data.test <- get_test_data(model)
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
  d.pars <-
    get_draws(
      model,
      groups = groups,
      summarize = F,
      wide = F,
      ndraws = ndraws,
      untransform_cues = untransform_cues)

  # ndraws <- if (is.null(ndraws)) get_number_of_draws(model) else ndraws
  # if (ndraws > 500)
  #   message(paste("Marginalizing over", ndraws, "MCMC samples. This might take some time.\n"))

  if (!is.null(sort_by))
    assert_that(sort_by %in% groups,
                msg = paste("sort_by must be NULL or one of the groups:",
                      paste(groups, collapse = ", ")))

  d.pars %<>%
    filter(group %in% .env$groups)

  # Prepare test_data
  cue.labels <- get_cue_levels(model)
  if (all_test_locations) {
    test_data <-
      data.test %>%
      distinct(!!! syms(cue.labels)) %>%
      { if (untransform_cues) get_untransform_function_from_stanfit(model)(.) else . } %>%
      make_vector_column(cols = cue.labels, vector_col = "x", .keep = "all") %>%
      nest(cues_joint = x, cues_separate = .env$cue.labels) %>%
      crossing(group = levels(d.pars$group))
  } else {
    test_data <-
      data.test %>%
      distinct(!!! syms(cue.labels), group) %>%
      { if (untransform_cues) get_untransform_function_from_stanfit(model)(.) else . } %>%
      make_vector_column(cols = cue.labels, vector_col = "x", .keep = "all") %>%
      group_by(group) %>%
      nest(cues_joint = x, cues_separate = .env$cue.labels)
  }

  # There are some issues with using exec instead of invoke_map in a mutate context since
  # !!! takes precedence over the definition of the function in the map statement
  # (see https://stackoverflow.com/questions/64356232/using-mutate-with-map2-and-exec-instead-of-invoke-map)
  # This works around that issue
  .fun <- function(fn, args) exec(fn, !!!args, target_category = target_category, logit = T)
  d.pars %<>%
    group_by(group, .draw) %>%
    do(f = get_categorization_function_from_grouped_ibbu_stanfit_draws(.)) %>%
    right_join(test_data, by = "group") %>%
    group_by(group, .draw) %>%
    mutate(p_cat = map2(f, cues_joint, .fun)) %>%
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

  target_category_label <- if (is.null(get_category_levels(model))) paste("category", target_category) else get_category_levels(model, target_category)

  if (plot_in_cue_space) {
    if (length(get_cue_levels(model)) == 2) stop("plot_in_cue_space = T not yet implemented for more than two cues (and makes no sense for a single cue).")

    p <-
      d.pars %>%
      mutate(group = factor(group, levels = .env$groups)) %>%
      ggplot(
        aes(
        x = !! sym(cue.labels[1]),
        y = !! sym(cue.labels[2]),
        fill = .data$p_cat)) +
      geom_raster(interpolate = FALSE) +
      scale_x_continuous(cue.labels[1]) +
      scale_y_continuous(cue.labels[2]) +
      scale_fill_gradient2(
        paste0("Predicted proportion\n of ", target_category_label, "-responses"),
        midpoint = .5,
        high = category.colors[target_category],
        mid = "white",
        low = if (length(category.colors[which(get_category_levels(model) != target_category_label)]) > 1) "gray" else category.colors[which(get_category_levels(model) != target_category_label)],
        limits = c(0,1)) +
      coord_cartesian(expand = F)
  } else {
    # If sort_by is specified, sort levels of x-axis by that group.
    if (!is.null(sort_by)) {
      sort.levels <-
        d.pars %>%
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
      ggplot(
        aes(
          x = .data$token,
          y = .data$p_cat,
          color = .data$group,
          linetype = .data$group)) +
      scale_x_discrete(
        "Test token",
        breaks = levels(d.pars$token),
        labels = paste0(levels(d.pars$token), "\n",
                        levels(d.pars$token.cues))) +
      scale_y_continuous(
        paste0("Predicted proportion of ", target_category_label, "-responses"),
        limits = c(0,1)) +
      scale_color_manual(
        "Group",
        breaks = groups,
        values = group.colors) +
      scale_linetype_manual(
        "Group",
        breaks = groups,
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
