# Additional plotting functions in plot-expected-categories.R

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
#' [tidybayes::recover_types] on the stanfit object to add or change these levels later.
#' (default: all categories/groups/cues will be plotted)
#' @param which Should parameters for the prior, posterior, or both be plotted? (default: `"both"`)
#' @param ndraws Number of draws to plot (or use to calculate the CIs), or `NULL` if all draws are to be returned. (default: `NULL`)
#' @param untransform_cues DEPRECATED. Should m_0 and S_0 be transformed back into the original cue space? (default: `FALSE`)
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
#' @export
plot_parameters <- function(fit, ...) {
  UseMethod("plot_parameters")
}

#' @rdname plot_parameters
#' @export
plot_parameters.ideal_adaptor_stanfit <- function(
  model,
  categories = get_category_levels(model),
  groups = get_group_levels(model, include_prior = T),
  cues = get_cue_levels(model),
  ndraws = NULL,
  untransform_cues = FALSE,
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
  legend <- suppressWarnings(cowplot::get_legend(p.m))

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
#' @param untransform_cues DEPRECATED. Should m_0 and S_0 be transformed back into the original cue space? (default: `FALSE`)
#' @param category.colors Vector of fill colors of same length as category or `NULL` to use defaults. (default: `NULL`)
#'
#' @return ggplot object.
#'
#' @details
#' Typically, the categories, groups, and cues are automatically added to the fit during the creation of the fit. If necessary,
#' however, it is possible to use [tidybayes::recover_types] on the stanfit object to add or change these levels later.
#'
#'
#' @seealso TBD
#' @keywords TBD
#'
#' @importFrom ggforce facet_matrix geom_autodensity
#' @importFrom ggnewscale new_scale
#' @importFrom colorspace lighten
#' @export
plot_parameter_correlations <- function(model, ...) {
  UseMethod("plot_parameter_correlations")
}

#' @export
#' @rdname plot_parameter_correlations
plot_parameter_correlations.ideal_adaptor_stanfit <- function(
  model,
  categories = get_category_levels(model),
  groups = "prior",
  cues = get_cue_levels(model),
  pars = NULL,
  ndraws = NULL,
  untransform_cues = FALSE,
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
#' @param model \code{\link{ideal_adaptor_stanfit}} object.
#' @param data.test Optionally, a \code{tibble} or \code{data.frame} with test data.
#'   If `NULL` the input will be extracted from fit. (default: `NULL`).
#' @param groups Character vector of groups to be plotted. Typically, the levels of these factors
#'   are automatically added to the fit during the creation of the fit. If necessary, however, it is possible to use
#'   [tidybayes::recover_types] on the stanfit object to add or change these levels later.
#'   (default: all categories/groups will be plotted)
#' @param summarize Should one categorization function (optionally with CIs) be plotted (`TRUE`) or should separate
#'   unique categorization function be plotted for each MCMC draw (`FALSE`)? (default: `TRUE`)
#' @param ndraws Number of draws to plot (or use to calculate the CIs), or `NULL` if all draws are to be returned.
#'   (default: `NULL`)
#' @param confidence.intervals The two confidence intervals that should be plotted (using `geom_ribbon`) around the mean.
#'   (default: `c(.66, .95)`)
#' @param target_category The index of the category for which categorization should be shown. (default: `1`)
#' @param panel.group Should the groups be plotted in separate panels? (default: `FALSE`)
#' @param group.colors,group.shapes,group.linetypes Vector of colors, shapes, and linetypes of same length as `groups` or `NULL` to use defaults.
#' @param category.colors Vector of colors and linetypes of same length as `categories` or `NULL` to use defaults. Only
#'   relevant when `plot_in_cue_space = TRUE`.
#' @param all_test_locations Should predictions be shown for all combinations of test locations and group, or should only
#'   combinations be shown that actually occurred in the data? (default: `FALSE`)
#' @param plot_in_cue_space Currently only available if the model has one or two cues. Should predictions be plotted in the cue space?
#'   If not, test tokens are treated as factors and sorted along the x-axis based on `sort_by`. (default: `TRUE`)
#' @param plot_test_data Should the test data be plotted? If `plot_in_cue_space = TRUE`, then test data will be shown as points
#'   on top of the raster. If not, then pointranges will be shown. (default: `TRUE`)
#' @param sort_by Which group, if any, should the x-axis be sorted by (in increasing order of posterior probability
#'   from left to right). Set to 0 for sorting by prior (default). Set to `NULL` if no sorting is desired. (default: `"prior"`)
#' @param untransform_cues DEPRECATED. Should the cues be untransformed before plotting? This should only have visual consequences
#'   if `plot_in_cue_space = T`. (default: `FALSE`)
#'
#' @return ggplot object.
#'
#' @seealso TBD
#' @keywords TBD
#'
#' @importFrom dplyr do right_join
#' @importFrom purrr invoke_map
#' @export
plot_expected_categorization <- function(model, ...) {
  UseMethod("plot_expected_categorization")
}

#' @rdname plot_expected_categorization
#' @export
plot_expected_categorization.ideal_adaptor_stanfit <- function(
  model,
  data.test = NULL,
  # handed to get_categorization_function:
  groups = get_group_levels(model, include_prior = T),
  lapse_treatment = c("no_lapses", "sample", "marginalize")[3],
  ndraws = NULL,
  # used by function returned by get_categorization_function:
  target_category = 1,
  logit = F,
  confidence.intervals = c(.66, .95),
  summarize = T,
  panel.group = if (plot_in_cue_space) TRUE else FALSE,
  group.colors = get_default_colors("group", groups),
  group.shapes = get_default_shapes("group", groups),
  group.linetypes = get_default_linetypes("group", groups),
  category.colors = get_default_colors("category", get_category_levels(model)),
  all_test_locations = TRUE,
  plot_in_cue_space = FALSE,
  plot_test_data = TRUE,
  sort_by = if ("prior" %in% groups) "prior" else NULL,
  # deprecated
  untransform_cues = FALSE
) {
  if (is.null(data.test)) data.test <- get_test_data(model, .from_staninput = T)
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
  d.pars <- get_categorization_function(model, groups = groups, lapse_treatment = lapse_treatment, ndraws = ndraws)

  if (!is.null(sort_by))
    assert_that(sort_by %in% groups,
                msg = paste("sort_by must be NULL or one of the groups:",
                      paste(groups, collapse = ", ")))

  # Prepare test_data
  cue.labels <- get_cue_levels(model)
  if (all_test_locations) {
    test_data <-
      data.test %>%
      distinct(!!! syms(cue.labels)) %>%
      { if (untransform_cues) get_untransform_function(model)(.) else . } %>%
      make_vector_column(cols = cue.labels, vector_col = "x", .keep = "all") %>%
      nest(cues_joint = x, cues_separate = .env$cue.labels) %>%
      crossing(group = levels(d.pars$group))
  } else {
    test_data <-
      data.test %>%
      distinct(!!! syms(cue.labels), group) %>%
      { if (untransform_cues) get_untransform_function(model)(.) else . } %>%
      make_vector_column(cols = cue.labels, vector_col = "x", .keep = "all") %>%
      group_by(group) %>%
      nest(cues_joint = x, cues_separate = .env$cue.labels)
  }

  # There are some issues with using exec instead of invoke_map in a mutate context since
  # !!! takes precedence over the definition of the function in the map statement
  # (see https://stackoverflow.com/questions/64356232/using-mutate-with-map2-and-exec-instead-of-invoke-map)
  # This works around that issue.
  #
  # Independently, note that we're first obtaining estimates in logits so that summarization
  # (if requested) is done over logits before we transform into probabilities.
  .fun <- function(fn, args) exec(fn, !!!args, target_category = target_category, logit = T)
  d.pars %<>%
    right_join(test_data, by = "group") %>%
    group_by(group, .draw) %>%
    mutate(p_cat = map2(f, cues_joint, .fun)) %>%
    select(-f) %>%
    unnest(c(cues_joint, cues_separate, p_cat))

  if (summarize) {
    f <- if (logit) I else plogis
    d.pars %<>%
      # For each unique group and test token obtain the CIs and the mean.
      group_by(group, x, !!! syms(cue.labels)) %>%
      summarise(
        across(
          p_cat,
          list(
            # na.rm = T excludes cases that might result from estimated probabilities of 0 and 1 (infinities in log-odds)
            y.outer.min = function(x) f(quantile(x, confidence.intervals[1], na.rm = T)),
            y.outer.max = function(x) f(quantile(x, confidence.intervals[4], na.rm = T)),
            y.inner.min = function(x) f(quantile(x, confidence.intervals[2], na.rm = T)),
            y.inner.max = function(x) f(quantile(x, confidence.intervals[3], na.rm = T)),
            p_cat = function(x) f(mean(x, na.rm = T))),
          .names = "{.fn}"))
  } else if (!logit) {
    d.pars %<>%
      mutate(p_cat = plogis(p_cat))
  }

  d.pars %<>% ungroup()
  target_category_label <- if (is.null(get_category_levels(model))) paste("category", target_category) else get_category_levels(model, target_category)

  if (plot_in_cue_space) {
    if (length(get_cue_levels(model)) != 2) stop("plot_in_cue_space = T not yet implemented for more than two cues (and makes no sense for a single cue).")

    p <-
      d.pars %>%
      mutate(group = factor(group, levels = .env$groups)) %>%
      ggplot(
        aes(
          x = map_dbl(.data$x, ~ .x[1]),
          y = map_dbl(.data$x, ~ .x[2]))) +
      geom_raster(aes(fill = .data$p_cat), interpolate = FALSE) +
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

    if (plot_test_data) p <- p + geom_point(data = test_data %>% unnest(cues_joint), color = "black")
  } else {
    # Get cues as rounded character strings (for the x-axis)
    d.pars %<>% mutate(token.cues = map(x, ~ paste(signif(.x), collapse = ",\n"))) %>% select(-c(x))

    # If sort_by is specified, sort levels of x-axis by that group.
    if (!is.null(sort_by)) {
      sort.levels <-
        d.pars %>%
        filter(group == sort_by) %>%
        group_by(token.cues) %>%
        summarise(p_cat = mean(p_cat)) %>%
        arrange(p_cat) %>%
        pull(token.cues)
    } else {
      sort.levels <- unique(d.pars$token.cues)
    }

    # Get cues as rounded character strings (for the x-axis)
    d.pars %<>%
      mutate(
        token.cues = factor(token.cues, levels = sort.levels),
        token = factor(as.numeric(token.cues)))

    p <-
      d.pars %>%
      ggplot(
        aes(
          x = .data$token.cues,
          y = .data$p_cat,
          color = .data$group,
          linetype = .data$group)) +
      scale_x_discrete("Test token") +
      # scale_x_discrete(
      #   "Test token",
      #   breaks = levels(d.pars$token),
      #   labels = paste0(levels(d.pars$token), "\n",
      #                   levels(d.pars$token.cues))) +
      scale_y_continuous(
        paste0("Predicted ", if (logit) "log-odds" else "proportion", " of ", target_category_label, "-responses"),
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

      if (plot_test_data) {
        message("plot_test_data not yet implemented when plot_in_cue_space = FALSE. Not plotting test data.")

        # Need access to participant information for adequate CIs.
        # p <-
        #   p +
        #   stat_summary(
        #     data =
        #       test_data %>%
        #       group_by(token.cues, group, ParticipantID) %>%
        #       summarise(p_cat = mean(Response == "SH")),
        #     fun.data = mean_cl_boot,
        #     geom = "pointrange", position = position_dodge(.25))
      }

      # Place information about confidence intervals on plot.
      p <- p +
        ggtitle(paste0((confidence.intervals[4]-confidence.intervals[1]) * 100,
                       "% and ",
                       (confidence.intervals[3]-confidence.intervals[2]) * 100,
                       "% CIs\nbased on ", ndraws, " posterior samples."))
    }

    p <-
      p +
      geom_point(alpha = .9) +
      geom_line(size = 1, alpha = .9, aes(x = as.numeric(.data$token)))
  }
  if (!summarize & panel.group) p <- p + facet_grid(group ~ .draw) else
    if (panel.group) p <- p + facet_wrap(~ group) else
      if (!summarize) p <- p + facet_wrap(~ .draw)

  return(p)
}
