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
#' @param fit mv-ibbu-stanfit object.
#' @param which Should parameters for the prior, posterior, or both be plotted? (default: `"both"`)
#' @param n.draws Number of draws to plot (or use to calculate the CIs), or `NULL` if all draws are to be returned. (default: `NULL`)
#' @param group.ids Vector of group IDs to be plotted or leave `NULL` to plot all groups. (default: `NULL`) It is possible
#' to use \code{\link[tidybayes]{recover_types}} on the stanfit object prior to handing it to this plotting function.
#' @param group.labels Vector of group labels of same length as group.ids or `NULL` to use defaults. (default: `NULL`)
#' The default labels each categorization function based on whether it is showing prior or posterior categorization,
#' and by its group ID.
#' @param group.colors Vector of fill colors of same length as group.ids or `NULL` to use defaults. (default: `NULL`)
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
  fit,
  which = "both",
  n.draws = NULL,
  group.ids = NULL, group.labels = NULL, group.colors = NULL,
  panel_scaling = F
) {
  # If n.draws is specified, get the IDs of the specific (randomly drawn n.draws) samples
  if (!is.null(n.draws)) draws = get_random_draw_indices(fit, n.draws)

  d.pars = fit %>%
    add_ibbu_stanfit_draws(
      which = which,
      draws = if (!is.null(n.draws)) draws else NULL,
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

  p.M = d.pars %>%
    select(.draw, group, category, cue, M) %>%
    distinct() %T>%
    { get_limits(., "M") ->> x.limits } %>%
    ggplot(aes(
      y = fct_rev(.data$category),
      x = .data$M,
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
    theme_bw() + theme(legend.position = "right")
  legend = cowplot::get_legend(p.M)

  p.M = p.M + theme(legend.position = "none")
  p.S = suppressMessages(
    p.M %+%
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
      p.M %+%
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
    theme_bw() + theme(legend.position = "none")

  K = length(unique(d.pars$cue))
  p = suppressWarnings(cowplot::plot_grid(
    cowplot::plot_grid(plotlist = list(p.M, p.KN), nrow = 1, rel_widths = c(K,2)),
    cowplot::plot_grid(plotlist = list(
      p.S,
      cowplot::plot_grid(plotlist = list(legend, p.LR),
                         nrow = 2, rel_heights = c(.5, .5))),
      nrow = 1, rel_widths = c(K,1)),
    rel_heights = c(1.5, K), nrow = 2, axis = "lrtb"))
  return(p)
}

#' Get categorization function
#'
#' Returns a categorization function for the first category, based on a set of parameters for the Normal-inverse-wishart (NIW)
#' distribtuion. Ms, Ss, kappas, nus, and priors are assumed to be of the same length and sorted the same way, so that the first
#' element of Ms is corresponding to the same category as the first element of Ss, kappas, nus, and priors, etc.
#'
#' @param Ms List of IBBU-inferred means describing the multivariate normal distribution over category means.
#' @param Ss List of IBBU-inferred scatter matrices describing the inverse Wishart distribution over category
#' covariance matrices.
#' @param kappas List of IBBU-inferred kappas describing the strength of the beliefs into the distribution over catgory means.
#' @param nus List of IBBU-inferred nus describing the strength of the beliefs into the distribution over catgory covariance matrices.
#' @param lapse_rate An IBBU-inferred lapse rate for the categorization responses.
#' @param priors Vector of categories' prior probabilities. (default: uniform prior over categories)
#' @param n.cat Number of categories. Is inferred from the input, but can be set manually.
#' @param logit Should the function that is returned return log-odds (TRUE) or probabilities (FALSE)? (default: TRUE)
#'
#' @return A function that takes as input cue values and returns posterior probabilities of the first category,
#' based on the posterior predictive of the cues given the (IBBU-derived parameters for the) categories' M, S,
#' kappa, nu, and prior, as well as the lapse rate.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
get_categorization_function = function(
  Ms, Ss, kappas, nus, lapse_rate,
  priors = rep(1 / n.cat, n.cat),
  n.cat = length(Ms),
  logit = FALSE
) {
  assert_that(are_equal(length(Ms), length(Ss)),
              are_equal(length(Ms), length(priors)),
              are_equal(length(Ms), length(kappas)),
              are_equal(length(Ms), length(nus)),
              msg = "The number of Ms, Ss, kappas, nus, and priors must be identical.")
  assert_that(between(lapse_rate, 0, 1))

  # Get dimensions of multivariate category
  D = length(Ms[[1]])
  assert_that(nus[[1]] >= D,
    msg = "Nu must be at least K (number of dimensions of the multivariate Gaussian category).")

  f <- function(x) {
    log_p = array()
    for (cat in 1:n.cat) {
      log_p[cat] = get_posterior_predictive(x, Ms[[cat]], Ss[[cat]], kappas[[cat]], nus[[cat]], log = T)
    }

    log_p1 = exp(
      log_p[1] + log(priors[1]) -
        log(sum(exp(log_p) * priors
        ))) * (1 - lapse_rate) + lapse_rate / n.cat

    if (logit)
      return(qlogis(log_p1))
    else
      return(log_p1)
  }

  return(f)
}



#' Get categorization function from grouped IBBU draws
#'
#' Convenience function intended for internal use.
#' @noRd
get_categorization_function_from_grouped_ibbu_stanfit_draws = function(fit, ...) {
  get_categorization_function(
    Ms = fit$M,
    Ss = fit$S,
    kappas = fit$kappa,
    nus = fit$nu,
    lapse_rate = unique(unlist(fit$lapse_rate)),
    ...
  )
}



#' Plot prior and posterior categorization of test tokens.
#'
#' Plot both prior and posterior categorization functions, as well as their confidence intervals.
#' If `summarize=TRUE`, the function marginalizes over all posterior samples. The number of samples
#' is determined by n.draws. If n.draws is NULL, all samples are used. Otherwise n.draws random
#' samples will be used. If `summarize=FALSE`, separate categorization plots for all n.draws
#' individual samples will be plotted in separate panels.
#'
#' @param fit mv-ibbu-stanfit object.
#' @param fit.input Input to the mv-ibbu-stanfit object.
#' @param which Should categorization for the prior, posterior, or both be plotted? (default: `"both"`)
#' @param summarize Should one categorization function (optionally with CIs) be plotted (`TRUE`) or should separate
#' unique categorization function be plotted for each MCMC draw (`FALSE`)? (default: `TRUE`)
#' @param n.draws Number of draws to plot (or use to calculate the CIs), or `NULL` if all draws are to be returned.
#' (default: `NULL`)
#' @param confidence.intervals The two confidence intervals that should be plotted (using `geom_ribbon`) around the mean.
#' (default: `c(.66, .95)`)
#' @param group.ids Vector of group IDs to be plotted or leave `NULL` to plot all groups. (default: `NULL`) It is possible
#' to use \code{\link[tidybayes]{recover_types}} on the stanfit object prior to handing it to this plotting function.
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
#'
plot_stanfit_test_categorization = function(
  fit,
  fit.input,
  which = "both",
  summarize = T,
  n.draws = NULL,
  confidence.intervals = c(.66, .95),
  group.ids = NULL, group.labels = NULL, group.colors = NULL, group.linetypes = NULL,
  sort.by = "prior"
) {
  warning("This function is not yet implemented!")
  assert_that(is.mv_ibbu_stanfit(fit))
  assert_that(!is.null(fit.input))
  assert_that(is.flag(summarize))
  assert_that(is.null(n.draws) | is.count(n.draws))
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

  # If n.draws is specified, get the IDs of the specific (randomly drawn n.draws) samples
  if (!is.null(n.draws)) draws = get_random_draw_indices(fit, n.draws)

  # Get prior and posterior parameters
  d.pars =
    add_ibbu_stanfit_draws(fit,
                   which = which,
                   summarize = F,
                   wide = F,
                   draws = if (!is.null(n.draws)) draws else NULL)

  # Now set n.draws to the number of MCMC samples
  n.draws = if (is.null(n.draws)) get_number_of_draws(fit) else n.draws
  if (n.draws > 500)
    message(paste("Marginalizing over", n.draws, "MCMC samples. This might take some time.\n"))

  # If group.ids are NULL set them to the levels of groups found in the extraction
  # of posteriors from fit
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
  test_data = fit.input$x_test %>%
    distinct() %>%
    transmute(x = pmap(.l = list(!!! syms(cues)), .f = ~ c(...))) %>%
    mutate(token = 1:length(x))

  # Store the number of test tokens since it's reused a number of times
  # THOUGH THERE MIGHT BE MORE ELEGANT SOLUTIONS TO THOSE LINES (SEE BELOW).
  # If you remove this line, make sure all dependencies are dealt with.
  n.tokens = nrow(test_data)

  # THIS PART (RATHER THAN THE SUMMARY BELOW) SEEMS TO BE THE SLOW PART.
  d.pars %<>%
    # Write a categorization function for each draw
    group_by(group, .draw) %>%
    do(f = get_categorization_function_from_grouped_ibbu_stanfit_draws(., logit = T)) %>%
    ungroup() %>%
    # Make as many copies of the data as there are test token (types) and label
    # each row for the test token (so that each data point can be categorized).
    # THERE MUST BE A MORE ELEGANT SOLUTION TO THIS.
    slice(rep(1:n(), each = n.tokens)) %>%
    mutate(token = rep(1:n.tokens, nrow(.) / n.tokens)) %>%
    # Join the test data with the cue information and then apply the categorization
    # function from each row (derived from the prior parameters of that draw) to the
    # token in that row.
    left_join(test_data) %>%
    mutate(
      token.cues = map(x, ~paste(.x, collapse = ",")),
      probability_cat1 = map(x, f),
      x = NULL
    ) %>%
    mutate_all(.funs = unlist) %>%
    ungroup()

  if (summarize) {
    d.pars %<>%
      select(-.draw) %>%
      # For each unique group and test token obtain the CIs and the mean.
      group_by(group, token, token.cues) %>%
      summarise_all(.funs = list(
        # na.rm = T excludes cases that might result from estimated probabilities of 0 and 1 (infinities in log-odds)
        y.outer.min = function(x) plogis(quantile(x, confidence.intervals[1], na.rm = T)),
        y.outer.max = function(x) plogis(quantile(x, confidence.intervals[4], na.rm = T)),
        y.inner.min = function(x) plogis(quantile(x, confidence.intervals[2], na.rm = T)),
        y.inner.max = function(x) plogis(quantile(x, confidence.intervals[3], na.rm = T)),
        probability_cat1 = function(x) plogis(mean(x, na.rm = T)))
      )
  } else {
    d.pars %<>%
      mutate(probability_cat1 = plogis(probability_cat1))
  }

  # If sort.by is specified, sort levels of x-axis by that group.
  if (!is.null(sort.by)) {
    token.levels = unique(test_data$token)[order((d.pars %>% filter(group == sort.by) %>%
                                                    select(group, token, probability_cat1) %>%
                                                    distinct())[["probability_cat1"]])]
    d.pars %<>%
      ungroup() %>%
      mutate(token = factor(token, levels = token.levels))
  }

  if (is.null(get_category_levels(fit)))
    category1 = "category 1" else category1 = get_category_levels(fit, 1)

  p = d.pars %>%
    ungroup() %>%
    mutate(group = factor(group, levels = group.ids)) %>%
    ggplot(aes(
      x = .data$token,
      y = .data$probability_cat1,
      color = .data$group,
      linetype = .data$group)) +
    scale_x_discrete("Test token") +
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
      )

    # Place information about confidence intervals on plot.
    p = p +
      annotate(geom = "text",
               x = mean(as.numeric(as.character(d.pars$token)), na.rm = T),
               y = 1,
               label = paste0((confidence.intervals[4]-confidence.intervals[1]) * 100,
                              "% and ",
                              (confidence.intervals[3]-confidence.intervals[2]) * 100,
                              "% CIs based on ",
                              n.draws,
                              " posterior samples.")
      )
  }

  p = p +
    geom_point(alpha = .9) +
    geom_line(size = 1, alpha = .9, aes(x = as.numeric(.data$token))) +
    theme_bw()

  if (!summarize) p = p + facet_wrap(~ .draw)

  return(p)
}



#' Plot expected bivariate (2D) categories from MV IBBU stanfit.
#'
#' Plot bivariate Gaussian categories expected given the parameters inferred by incremental Bayesian belief-
#' updating (IBBU). Specifically, the categories are derived by marginalizing over the uncertainty represented
#' by the (post-warmup) MCMC samples. Two methods are available (specified by `type`), which differ in their
#' computational demands and speed.
#'
#' @param fit mv-ibbu-stanfit object.
#' @param fit.input Optionally, the input to the mv-ibbu-stanfit object, in which case the test tokens will also be plotted,
#' using `geom_point()`.
#' @param type Either `"contour"` or `"density"`, specifying the type of plot. Note that the contour plot is *much*
#' faster. It simply gets the expected values of \code{mu} (based on the NIW parameter \code{M}) and \code{Sigma}
#' (based on the NIW parameters \code{S} and \code{nu}) at each MCMC draw, and then averages over
#' all MCMC draws. The plotted categories represent those means of the expected \code{mu} and \code{Sigma}. The
#' density plot instead calculates the posterior predictive for each MCMC draw (i.e, the multivariate Student-T
#' density based on the NIW parameters \code{M, S, kappa, nu}), and then averages those densities. Since this is
#' done for *all* points defined by the data.grid this can be rather computationally expensive and slow.
#' @param plot.test,plot.exposure Should the test and/or exposure stimuli be plotted? (default: `TRUE` for `plot.test`,
#' `FALSE` for `plot.exposure`) The test items are plotted as black points. The exposure mean is plotted as point,
#' and the .95 interval of cue distributions during exposure are plotted as dashed ellipse in the same color as the
#' expected categories.
#' @param summarize Should one expected categories be plotted, marginalizing over MCMC draws (`TRUE`), or should separate
#' expected categories be plotted for each MCMC draw (`FALSE`)? (default: `TRUE`) Currently being ignored.
#' @param n.draws Number of draws to plot (or use to calculate the CIs), or `NULL` if all draws are to be returned.
#' (default: `NULL`) Currently being ignored.
#' @param levels Used only if `type` is `"contour"`. levels The cumulative probability levels that should be plotted (using
#' `geom_polygon()`) around the mean. By default the most transparent ellipse still drawn corresponds to .95.
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
  x,
  fit.input = NULL,
  type,
  ...
) {
  assert_that(all(type %in% c("contour", "density"), length(type) == 1))
  if (type == "contour")
    plot_expected_ibbu_stanfit_categories_contour2D(x = x, fit.input = fit.input, ...)
  else {
    plot_expected_ibbu_stanfit_categories_density2D(x = x, fit.input = fit.input, ...)
  }
}

#' @rdname plot_expected_ibbu_stanfit_categories_2D
#' @export
plot_expected_ibbu_stanfit_categories_contour2D = function(
  x,
  fit.input = NULL, # should change in the future
  levels = plogis(seq(-15, qlogis(.95), length.out = 20)),
  plot.test = T, plot.exposure = F,
  category.ids = NULL, category.labels = NULL, category.colors = NULL, category.linetypes = NULL
) {
  assert_that(is.mv_ibbu_stanfit(x) | is.mv_ibbu_MCMC(x))
  assert_that(!all(is.null(fit.input), plot.test))

  d = get_expected_category_statistic(x)

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
    facet_wrap(~ group) +
    theme_bw()
}


#' @rdname plot_expected_ibbu_stanfit_categories_2D
#' @export
plot_expected_ibbu_stanfit_categories_density2D = function(
  x,
  fit.input = NULL, # should change in the future
  plot.test = T, plot.exposure = F,
  category.ids = NULL, category.labels = NULL, category.colors = NULL, category.linetypes = NULL,
  xlim = c(-10, 10), ylim = c(-10, 10), resolution = 10
) {
  assert_that(is.mv_ibbu_stanfit(x) | is.mv_ibbu_MCMC(x))
  assert_that(!all(is.null(fit.input), plot.test))

  if (is.mv_ibbu_stanfit(x))
    d = add_ibbu_stanfit_draws(x, which = which, wide = F, nest = T)
  else
    d = x

  # Setting aes defaults
  if(is.null(category.ids)) category.ids = levels(d$category)
  if(is.null(category.labels)) category.labels = levels(d$category)
  if(is.null(category.colors)) category.colors = get_default_colors("category", length(category.ids))
  if(is.null(category.linetypes)) category.linetypes = rep(1, length(category.ids))

  cue.names = row.names(d$M[[1]])
  d %<>%
    crossing(
      cue1 = seq(min(xlim), max(xlim), length.out = resolution),
      cue2 = seq(min(ylim), max(ylim), length.out = resolution)) %>%
    mutate(x = map2(cue1, cue2, ~ c(.x, .y))) %>%
    mutate(
      density = pmap(., get_posterior_predictive.pmap),
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
    scale_color_manual("Category",
                       breaks = category.ids,
                       labels = category.labels,
                       values = category.colors) +
    coord_fixed(xlim = xlim, ylim = ylim, ratio = 1) +
    facet_wrap(~ .data$group) +
    theme_bw()
}

