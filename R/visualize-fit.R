#' @import assertthat purrr ggplot2
#' @importFrom mvtnorm dmvt
#' @importFrom magrittr %<>%
#' @importFrom dplyr %>% mutate summarise
#' @importFrom rlang !! !!! sym syms expr
NULL

#' Plot prior and posterior categorization of test tokens.
#'
#' Plot both prior and posterior categorization functions, as well as their confidence intervals.
#' If summarize=TRUE, the function marginalizes over all posterior samples. The number of samples
#' is determined by n.draws. If n.draws is NULL, all samples are used. Otherwise n.draws random
#' samples will be used. If summarize=FALSE, separate categorization plots for all n.draws
#' individual samples will be plotted in separate panels.
#'
#' @param fit mv-ibbu-stanfit object.
#' @param fit.input Input to the mv-ibbu-stanfit object.
#' @param n.draws Number of draws to plot (or use to calculate the CIs), or NULL if all draws are to be returned. (default: NULL)
#' @param summarize Should one categorization function (optionally with CIs) be plotted (TRUE) or should separate
#' unique categorization function be plotted for each MCMC draw (FALSE)? (default: FALSE)
#' @param group.ids Vector of group IDs to be plotted or leave NULL to plot all groups. (default: NULL) It is possible
#' to use \code{\link[tidybayes]{recover_types}} on the stanfit object prior to handing it to this plotting function.
#' @param group.labels Vector of group labels of same length as group.ids or NULL to use defaults. (default: NULL)
#' The defaultlabels each categorization function based on whether it is showing prior or posterior categorization,
#' and by its group ID.
#' @param group.colors Vector of colors of same length as group.ids or NULL to use defaults. (default: NULL)
#' @param group.linetypes Vector of linetypes of same length as group.ids or NULL to use defaults. (default: NULL)
#' @param sort_by Which group, if any, should the x-axis be sorted by (in increasing order of posterior probability
#' from left to right). Set to 0 for sorting by prior (default). Set to NULL if no sorting is desired.
#'
#' @return ggplot.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
#'
plot_ibbu_test_categorization = function(
  fit,
  fit.input,
  summarize = T,
  n.draws = NULL,
  confidence.interval = c(.025, .25, .75, .975),
  group.ids = NULL, group.labels = NULL, group.colors = NULL, group.linetypes = NULL,
  sort.by = 0
) {
  assert_that(is.mv_ibbu_stanfit(fit))
  assert_that(is.flag(summarize))
  assert_that(is.null(n.draws) | is.count(n.draws))
  assert_that(is.null(confidence.intervals) |
                all(is.numeric(confidence.intervals),
                    length(confidence.interval) == 4,
                    all(between(confidence.interval, 0, 1))),
              msg = "Confidence intervals must be NULL (if not CIs are desired) or a vector of four probabilities.")
  assert_that(is.null(group.ids) | is.numeric(group.ids))
  assert_that(is.null(group.labels) | is.character(group.labels))
  assert_that(is.null(group.linetypes) | is.numeric(group.linetypes))
  assert_that(is.null(sort.by) | length(sort.by) == 1)

  # If n.draws is specified, get the IDs of the specific (randomly drawn n.draws) samples
  if (!is.null(n.draws))
    draws = sample(1:nrow(as.data.frame(fit, pars = "lapse_rate")), size = n.draws)

  # Get prior and posterior parameters
  d.pars =
    add_ibbu_draws(fit,
                   which = c("both"),
                   summarize = F,
                   wide = F,
                   draws = if (!is.null(n.draws)) draws else NULL)

  # If group.ids are NULL set them to the levels of groups found in the extraction
  # of posteriors from fit
  if (is.null(group.ids))  group.ids = unique(d.pars$group)
  assert_that(all(group.ids %in% unique(d.pars$group)),
              msg = "Some group.ids were not found in the stanfit object.")
  assert_that(sort.by %in% group.ids,
              msg = "Sort.by must be NULL or one of the group IDs (group.ids).")

  # Setting aes defaults
  if(is.null(group.labels)) group.labels = paste0("posterior (", group.ids[-1], ")")
  if(is.null(group.colors)) group.colors = rep("black", length(group.ids) - 1)
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
    rowwise() %>%
    transmute(x = list(c(!!! rlang::syms(cues)))) %>%
    ungroup() %>%
    mutate(token = 1:length(x))

  # Store the number of test tokens since it's reused a number of times
  # THOUGH THERE MIGHT BE MORE ELEGANT SOLUTIONS TO THOSE LINES (SEE BELOW).
  # If you remove this line, make sure all dependencies are dealt with.
  n.tokens = nrow(test_data)
  n.samples = if (is.null(n.draws)) get_number_of_draws(fit) else n.draws
  # intentionally NOT named n.draws, as this is meant to also capture the case when all posterior samples are used.

  if (n.samples > 500)
    cat(paste("You are marginalizing over", n.samples, "samples. This might take some time.\n"))

  # THIS PART (RATHER THAN THE SUMMARY BELOW) SEEMS TO BE THE SLOW PART.
  d.pars %<>%
    # Write a categorization function for each draw
    group_by(group, .draw) %>%
    do(f = get_categorization_function( #  get_categorization_function_from_ibbu_draws(., logit = T))
      Ms = .$M,
      Ss = .$S,
      kappas = .$kappa,
      nus = .$nu,
      lapse_rate = unique(unlist(.$lapse_rate)),
      logit = T
    )) %>%
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
    rowwise() %>%
    do(
      .draw = .$.draw,
      group = .$group,
      token = .$token,
      token.cues = paste(.$x, collapse = ","),
      probability_cat1 = .$f(.$x)
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
        y.outer.min = function(x) plogis(quantile(x, confidence.interval[1], na.rm = T)),
        y.outer.max = function(x) plogis(quantile(x, confidence.interval[4], na.rm = T)),
        y.inner.min = function(x) plogis(quantile(x, confidence.interval[2], na.rm = T)),
        y.inner.max = function(x) plogis(quantile(x, confidence.interval[3], na.rm = T)),
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

  if (is.null(attr(fit, "constructors")$category(1))) {
    category1 = "category 1"
    warning(paste0(class_name, " object does not contain type information about the categories.
                   Consider applying recover_types() from the tidybayes package first."))
  } else category1 = attr(fit, "constructors")$category(1)

  p = d.pars %>%
    ungroup() %>%
    mutate(group = factor(group)) %>%
    ggplot(aes(x = token, y = probability_cat1, color = group, linetype = group)) +
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

  if (summarize & !is.null(confidence.interval)) {
    p = p +
      geom_ribbon(
        aes(x = as.numeric(token), ymin = y.outer.min, ymax = y.outer.max, fill = group),
        color = NA, alpha = .1
      ) +
      geom_ribbon(
        aes(x = as.numeric(token), ymin = y.inner.min, ymax = y.inner.max, fill = group),
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
               label = paste0((confidence.interval[4]-confidence.interval[1]) * 100,
                              "% and ",
                              (confidence.interval[3]-confidence.interval[2]) * 100,
                              "% CIs based on ",
                              n.samples,
                              " posterior samples.")
      )
  }

  p = p +
    geom_point(alpha = .9) +
    geom_line(size = 1, alpha = .9, aes(x = as.numeric(token)))

  if (!summarize) p = p + facet_wrap(~ .draw)

  return(p)
}



get_categorization_function_from_ibbu_draws = function(fit, ...) {
  get_categorization_function(
    Ms = fit$M,
    Ss = fit$S,
    kappas = fit$kappa,
    nus = fit$nu,
    lapse_rate = unique(unlist(fit$lapse_rate)),
    ...
  )
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
#'
get_categorization_function = function(
  Ms, Ss, kappas, nus, lapse_rate,
  priors = 1 / n.cat,
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
  K = length(Ms[[1]])
  assert_that(nus[[1]] < K,
    msg = "Nu must be at least K (number of dimensions of the multivariate Gaussian category).")

  f <- function(x) {
    log_p = array()
    for (cat in 1:n.cat) {
      log_p[cat] = mvtnorm::dmvt(x,
                                 delta = Ms[[cat]],
                                 sigma = Ss[[cat]] * (kappas[[cat]] + 1) / (kappas[[cat]] * (nus[[cat]] - K + 1)),
                                 df = nus[[cat]] - K + 1,
                                 log = TRUE)
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

