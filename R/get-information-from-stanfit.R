#' @importFrom tidybayes spread_draws recover_types
#' @importFrom tidyr spread gather unite
NULL


#' Get parameter names
#'
#' Get the names for all parameters in `fit`.
#' @export
get_params = function(fit) {
  return(fit@model_pars)
}

#' Get number of post-warmup MCMC samples from stanfit
#'
#' Get the total number of post-warmup MCMC samples in `fit`.
#' @export
get_number_of_draws = function(fit) {
  return(length(fit@sim$samples[[1]][[1]]))
}

#' Get indices for random MCMC draws from stanfit
#'
#' Returns `n.draws` indices for random post-warmup MCMC draws (without replacement) from
#' `fit`.
#'
#' @param fit mv_ibbu_stanfit object.
#' @param n.draws Number of indices to be returned. Can't be larger than total number of
#' post-warmup samples across all MCMC chains in `fit`.
#'
#' @return A numeric vector.
get_random_draw_indices = function(fit, n.draws)
{
  n.all.draws = get_number_of_draws(fit)
  assert_that(n.draws <= n.all.draws)

  draws = sample(1:n.all.draws, size = n.draws)
  return(draws)
}


#' Get or restore the original group or category levels.
#'
#' Checks if information is available about the original values and order of the factor levels
#' for the category variable (for which beliefs about means and covariances are inferred) or
#' group variable (e.g., subject or exposure group), respectively. If available,
#' that information is returned. `get_category_levels()` and `get_group_levels()` are
#' convenience functions, calling `get_original_levels()`.
#'
#' @param fit mv_ibbu_stanfit object.
#' @param variable Either "category" or "group".
#' @param indeces A vector of category or group indices that should be turned into the original
#' category levels, or `NULL` if only the unique levels in their original order (as vector of characters)
#' should be returned. (default: `NULL`)
#'
#' @return If no category or group indices are provided, the levels of the category/group are returned (in the
#' original order). Otherwise a vector of the same length as \code{indices}
#'
#' @seealso \code{\link[tidybayes]{recover_types}} from tidybayes, \code{\link{get_constructor}}
#' @keywords TBD
#' @examples
#' TBD
#' @rdname get_original_levels
#' @export
get_original_levels = function(fit, variable = c("category", "group"), indices = NULL) {
  assert_that(is.null(indices) | all(indices > 0))
  f = get_constructor(fit, variable)

  if (is.null(indices)) return(levels(f(c()))) else return(f(indices))
}


#' @rdname get_original_levels
#' @export
get_category_levels = function(fit, indices = NULL) {
  return(get_original_levels(fit, "category", indices))
}

#' @rdname get_original_levels
#' @export
get_group_levels = function(fit, indices = NULL) {
  return(get_original_levels(fit, "group", indices))
}


#' Get tidybayes constructor from IBBU stanfit
#'
#' Gets the tidybayes constructor function from the stanfit object. `get_category_constructor()` and
#' `get_group_constructor()` are convenience functions, calling `get_constructor()`. See \code{
#' \link[tidybayes]{recover_types}}. If variable is
#'
#' @param fit mv_ibbu_stanfit object.
#' @param variable Either "category" or "group". If set to `NULL` then a list of all constructors is
#' returned. That list is `NULL` if not tidybayes constructors are found in fit. (default: c("category", "group"))
#'
#' @return A constructor function, a list of constructor functions, or `NULL`. If a specific constructor
#' function is requested but not found, a warning is shown.
#'
#' @seealso \code{\link[tidybayes]{recover_types}} from tidybayes, \code{\link{get_original_levels}}
#' @keywords TBD
#' @examples
#' TBD
#' @rdname get_constructor
#' @export
get_constructor = function(fit, variable = NULL) {
  assert_NIW_ibbu_stanfit(fit)
  if (is.null(variable)) return(attr(fit, "tidybayes_constructors"))

  assert_that(variable %in% c("category", "group"), msg = "Variable name must be one of category or group.")

  if (is.null(attr(fit, "tidybayes_constructors")[[rlang::sym(variable)]])) {
    warning(paste0(class_name, " object does not contain type information about the ", variable,
                   " variable. Consider applying tidybayes::recover_types() first."))
    return(NULL)
  }

  f = attr(fit, "tidybayes_constructors")[[rlang::sym(variable)]]

  return(f)
}


#' @rdname get_constructor
#' @export
get_category_constructor = function(fit) {
  return(get_constructor(fit, "category"))
}

#' @rdname get_constructor
#' @export
get_group_constructor = function(fit) {
  return(get_constructor(fit, "group"))
}




#' Get category mean mu or covariance matrix sigma from an NIW belief MCMC object.
#'
#' Returns the category means mu and/or category covariance matrix Sigma for the exposure data from an IBBU
#' stanfit or NIW belief MCMC object.
#'
#' @param x An mv_ibbu_stanfit or NIW belief MCMC object.
#' @param category Character vector with categories (or category) for which category statistics are to be
#' returned.  If `NULL` then all categories are included. (default: `NULL`)
#' @param group Character vector with groups (or group) for which category statistics are to be
#' returned. If `NULL` then all groups are included. (default: `NULL`)
#' @param statistic Which category statistic should be returned? `mu` for category mean or `Sigma` for category
#' covariance matrix, or `c("mu", "Sigma")` for both. (default: both)
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @rdname get_category_statistic
#' @export
get_category_statistic = function(x, grouping.vars = NULL,
                                  statistic = c("mu", "Sigma")) {
  assert_that(all(statistic %in% c("mu", "Sigma")))
  stop("get_category_statistics not yet implemented!")

  # More here ######################################
  # Make general so as to extract mu and sigma for any combination of grouping variables
  # Catch case when grouping variables are not specified (NULL)

  # FIX FIX FIX FIX If just one category and group was requested, just return that object, rather
  # than the tibble
  if (nrow(x) == 1) x = x[,paste0(statistic, ".mean")][[1]][[1]]
  return(x)
}


#' Get category mean mu or covariance matrix sigma of exposure data for IBBU
#'
#' Returns the category means mu and/or category covariance matrix Sigma for the exposure data for an incremental
#' Bayesian belief-updating (IBBU) model from an IBBU stanfit or NIW belief MCMC object.
#'
#' @param x An mv_ibbu_stanfit or NIW belief MCMC object.
#' @param category Character vector with categories (or category) for which category statistics are to be
#' returned.  If `NULL` then all categories are included. (default: `NULL`)
#' @param group Character vector with groups (or group) for which category statistics are to be
#' returned. If `NULL` then all groups are included. (default: `NULL`)
#' @param statistic Which category statistic should be returned? `mu` for category mean or `Sigma` for category
#' covariance matrix, or `c("mu", "Sigma")` for both. (default: both)
#'
#' @return If just one group and category was requested, a vector (for the mean) or matrix (for the covariance
#' matrix). If more than one group or category was requested, a tibble with one row for each unique combination
#' of group and category.

#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @rdname get_ibbu_stanfit_exposure_category_statistic
#' @export
get_ibbu_stanfit_exposure_category_statistic = function(x, category = NULL, group = NULL,
                                  statistic = c("mu", "Sigma")) {
  assert_that(is.NIW_ibbu_input(x) | is.NIW_ibbu_stanfit(x))
  stop("get_ibbu_stanfit_exposure_statistics not yet implemented!")

  x = get_ibbu_stanfit_input(x)

  # More here. ######################################
  # filter out group "prior"
  # deal with cases for which there is no exposure data
  # Assume that all cues are used

  return(get_ibbu_stanfit_category_statistic(x, grouping.vars = c("category", "group"), statistic))
}

#' @rdname get_ibbu_stanfit_exposure_category_statistic
#' @export
get_ibbu_stanfit_exposure_mean = function(x, category, group) {
  return(get_ibbu_stanfit_exposure_category_statistic(x, category, group, statistic = "mu"))
}

#' @rdname get_ibbu_stanfit_exposure_category_statistic
#' @export
get_ibbu_stanfit_exposure_Sigma = function(x, category, group) {
  return(get_ibbu_stanfit_exposure_category_statistic(x, category, group, statistic = "Sigma"))
}



#' @rdname get_categorization_function
#' @export
get_categorization_function_from_grouped_ibbu_stanfit_draws = function(fit, ...) {
  get_categorization_function(
    ms = fit$m,
    Ss = fit$S,
    kappas = fit$kappa,
    nus = fit$nu,
    lapse_rate = unique(unlist(fit$lapse_rate)),
    ...
  )
}

#' Get the input data from an MV IBBU stanfit object.
#'
#' Returns the inputs handed to \code{stan} or \code{sampling} during the creation of the \code{stanfit}
#' object.
#'
#' @param x An mv_ibbu_stanfit object.
#'
#' @return A list with element names and structure determined by the type of MV IBBU model.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @rdname get_ibbu_input
#' @export
get_ibbu_stanfit_input = function(x) {
  assert_that(is.NIW_ibbu_input(x) | is.NIW_ibbu_stanfit(x))

  if (is.NIW_ibbu_input(x)) return(x) else {
    stop("Extraction of input data from MV IBBU stanfit not yet implemented!")
  }

  return(x)
}



#' Add MCMC draws of IBBU parameters to a tibble.
#'
#' Add MCMC draws of all parameters from incremental Bayesian belief-updating (IBBU) to a tibble. Both wide
#' (`wide=TRUE`) or long format (`wide=FALSE`) can be chosen as output. By default all post-warmup draws are
#' returned, but if `summarize=TRUE` then just the mean of each parameter is returned instead.
#'
#' By default, the category means and scatter matrices are nested, rather than each of their elements being
#' stored separately (`nest=TRUE`).
#'
#' @param fit mv-ibbu-stanfit object.
#' @param which Should parameters for the prior, posterior, or both be added? (default: `"posterior"`)
#' @param draws Vector with specific draw(s) to be returned, or `NULL` if all draws are to be returned. (default: `NULL`)
#' @param summarize Should the mean of the draws be returned instead of all of the draws? (default: `FALSE`)
#' @param wide Should all parameters be returned in one row? (default: `FALSE`)
#' @param nest Should the category mean vectors and scatter matrices be nested into one cell each, or should each element
#' be stored in a separate cell? (default: `TRUE`)
#' @param category Name of the category variable. (default: "category")
#' @param group Name of the group variable. (default: "group")
#'
#' @return tibble with post-warmup (posterior) MCMC draws of the prior/posterior parameters of the IBBU model
#' (\code{kappa, nu, m, S, lapse_rate}). \code{kappa} and \code{nu} are the pseudocounts that determine the strength of the beliefs
#' into the mean and covariance matrix, respectively. \code{m} is the mean of the multivariate normal distribution over category
#' means mu. \code{S} is the scatter matrix that determines both the covariance of the category means mu, and the
#' Inverse Wishart distribution over category covariance matrices Sigma.
#'
#' The expected value of the category mean mu is \code{m}. The expected value of the category covariance matrix Sigma
#' is \code{S / (nu - D - 1)}, where \code{D} is the dimension of the multivariate Gaussian category. For details,
#' \insertCite{@see @murphy2012 p. 134;textual}{MVBeliefUpdatr}.
#'
#' @seealso TBD
#' @keywords TBD
#' @references \insertRef{murphy2012}{MVBeliefUpdatr}
#' @examples
#' TBD
#' @export
add_ibbu_stanfit_draws = function(
  fit,
  which = "posterior",
  draws = NULL,
  summarize = FALSE,
  wide = FALSE,
  nest = TRUE,
  category = "category",
  group = "group"
) {
  assert_that(which %in% c("prior", "posterior", "both"))
  assert_that(any(is.null(draws), all(draws > 0)))
  assert_that(is.flag(summarize))
  assert_that(is.flag(wide))
  assert_that(!all(wide, !nest),
              msg = "Wide format is currently not implemented without nesting.")

  if (which == "both") {
    d.prior = add_ibbu_stanfit_draws(fit = fit, which = "prior",
                             draws = draws,
                             summarize = summarize, wide = wide, nest = nest)
    d.posterior = add_ibbu_stanfit_draws(fit = fit, which = "posterior",
                   draws = draws,
                   summarize = summarize, wide = wide, nest = nest)
    d.pars = rbind(d.prior, d.posterior) %>%
      mutate(!! rlang::sym(group) :=
               factor(!! rlang::sym(group),
                      levels = c(with(d.prior, levels(!! rlang::sym(group))),
                                 with(d.posterior, levels(!! rlang::sym(group))))))

    return(d.pars)
  } else {
    assert_NIW_ibbu_stanfit(fit)

    # Parameters' names depend on whether prior or posterior is to be extracted.
    postfix = if (which == "prior") "_0" else "_n"
    kappa = paste0("kappa", postfix)
    nu = paste0("nu", postfix)
    m = paste0("m", postfix)
    S = paste0("S", postfix)

    # Variables by which parameters are indexed
    pars.index = if (which == "prior") category else c(category, group)

    # Should m and S be nested?
    if (!nest) {
      if (which == "prior")
        d.pars = fit %>%
          spread_draws(
            !! rlang::sym(kappa),
            !! rlang::sym(nu),
            (!! rlang::sym(m))[!!! rlang::syms(pars.index), cue],
            (!! rlang::sym(S))[!!! rlang::syms(pars.index), cue, cue2],
            lapse_rate,
            n = draws
          )
      else
        d.pars = fit %>%
          spread_draws(
            (!! rlang::sym(kappa))[!!! rlang::syms(pars.index)],
            (!! rlang::sym(nu))[!!! rlang::syms(pars.index)],
            (!! rlang::sym(m))[!!! rlang::syms(pars.index), cue],
            (!! rlang::sym(S))[!!! rlang::syms(pars.index), cue, cue2],
            lapse_rate,
            n = draws
          )

      if (any(c("mu_0", "mu_n", "sigma_0", "sigma_n") %in% names(d.pars))) {
        message("This seems to be an old model. Some of the NIX/NIW parameters are called mu_* or sigma_*. Renaming
                them to m and S.")
        d.pars %<>%
          rename(m = !! rlang::sym(m), S = !! rlang::sym(S)) %>%
          rename_at(vars(ends_with(postfix)), ~ sub(postfix, "", .))
      }

    # If nesting is the goal:
    } else {
      # Get kappa and nu
      if (which == "prior") {
        kappa_nu =
          fit %>%
          tidybayes::spread_draws(
            !! rlang::sym(kappa),
            !! rlang::sym(nu)
          )
      } else {
        kappa_nu =
          fit %>%
          tidybayes::spread_draws(
            (!! rlang::sym(kappa))[!!! rlang::syms(pars.index)],
            (!! rlang::sym(nu))[!!! rlang::syms(pars.index)]
          ) %>%
          group_by(!!! rlang::syms(pars.index))
      }

      # Get lapse rate and join it with kappa and nu.
      d.pars = fit %>%
        tidybayes::spread_draws(lapse_rate) %>%
        { if (!is.null(draws)) filter(., .draw %in% draws) else . } %>%
        { if (summarize)
          dplyr::summarise(.,
                           .chain = "all", .iteration = "all", .draw = "all",
                           lapse_rate = mean(lapse_rate)
          ) else . } %>%
        left_join(kappa_nu %>%
        { if (!is.null(draws)) filter(., .draw %in% draws) else . } %>%
          rename(
            kappa = !! rlang::sym(kappa),
            nu = !! rlang::sym(nu)
          ) %>%
          { if (summarize)
            dplyr::summarise(.,
                             .chain = "all", .iteration = "all", .draw = "all",
                             kappa = mean(kappa),
                             nu = mean(nu)
            ) else . } %>%
          ungroup()
        ) %>%
        # Join in m
        left_join(
          fit %>%
            tidybayes::spread_draws((!! rlang::sym(m))[!!! rlang::syms(pars.index), cue]) %>%
            { if (!is.null(draws)) filter(., .draw %in% draws) else . } %>%
            { if (summarize)
              group_by(., !!! rlang::syms(pars.index), cue) %>%
                # Obtain expected mean of category mean, m_0 or m_N
                dplyr::summarise(., !! rlang::sym(m) := mean(!! rlang::sym(m))
                ) %>%
                mutate(.,
                       .chain = "all", .iteration = "all", .draw = "all"
                ) else . } %>%
            group_by(.chain, .iteration, .draw, !!! rlang::syms(pars.index)) %>%
            summarise(m = list(
              matrix((!! rlang::sym(m)),
                     dimnames = list(unique(cue), NULL),
                     byrow = T,
                     nrow = length((!! rlang::sym(m))))))
        ) %>%
        # Join in S
        left_join(
          fit %>%
            tidybayes::spread_draws((!! rlang::sym(S))[!!! rlang::syms(pars.index), cue, cue2]) %>%
            { if (!is.null(draws)) filter(., .draw %in% draws) else . } %>%
            { if (summarize)
              group_by(., !!! rlang::syms(pars.index), cue, cue2) %>%
                # Obtain expected scatter matrix S_0, or S_N
                dplyr::summarise(., !! rlang::sym(S) := mean(!! rlang::sym(S))
                ) %>%
                mutate(.,
                       .chain = "all", .iteration = "all", .draw = "all"
                ) else . } %>%
            group_by(.chain, .iteration, .draw, !!! rlang::syms(pars.index)) %>%
            summarise(S =
                        list(
                          matrix((!! rlang::sym(S)),
                                 dimnames = list(unique(cue), unique(cue2)),
                                 byrow = T,
                                 nrow = sqrt(length((!! rlang::sym(S)))))))
        )
    }

    # Make sure order of variables is identical for prior or posterior (facilitates processing of the
    # output of this function). If group is prior, then add that variable to d.pars first.
    if (which == "prior") d.pars %<>% mutate((!! rlang::sym(group)) := "prior")
    d.pars %<>% select(.chain, .iteration, .draw,
                       !! rlang::sym(group), !! rlang::sym(category),
                       # Using starts_with in order to capture case in which variables are *not* nested
                       starts_with("kappa"), starts_with("nu"), starts_with("m"), starts_with("S"),
                       lapse_rate)

    # Clean-up
    d.pars %<>%
      ungroup() %>%
      # Make sure that group and category are factors (even if group or category ids are just numbers)
      mutate_at(vars(!! rlang::sym(group), !! rlang::sym(category)), factor)

    if (wide) {
      if (which == "prior")
        d.pars %<>% gather(variable, value, c(m, S))
      else
        d.pars %<>% gather(variable, value, c(kappa, nu, m, S))

      d.pars %<>%
        unite(temp, !!! rlang::syms(pars.index), variable) %>%
        spread(temp, value)
    }

    return(d.pars)
  }
}





#' Get expected category mean mu or covariance matrix sigma
#'
#' Returns the expected value of posterior marginal distribution over category means mu and/or
#' category covariance matrix Sigma, marginalized over all MCMC samples.
#'
#' Each MCMC samples' expected value for the category mean \code{E[mu] = m_n}
#' (i.e, the posterior/updated mean of the multivariate Normal over category means \code{mu}).
#' Marginalizing across all MCMC samples (representing uncertainty in the true value of
#' \code{m_n}), we get \code{E[E[mu]] = mean(m_n)}.
#'
#' Each MCMC samples' expected value for the category covariance matrix
#' \code{E[Sigma] = S_n / (nu_n - D - 1)}, where \code{S_n} is the posterior/updated scatter matrix,
#' \code{nu_n} is the posterior/updated pseudocount representing the strength of the posterior/updated
#' beliefs over category covariance matrices sigma (i.e., the inverse-Wishart), and \code{D} is
#' the dimension of the multivariate Normal. Marginalizing across all MCMC samples
#' (representing uncertainty in the true value of \code{S_n}), we get
#' \code{E[E[Sigma]] = mean(S_n / (nu_n - D - 1))}.
#'
#' @param x An mv_ibbu_stanfit or NIW_belief_MCMC object.
#' @param category Character vector with categories (or category) for which category statistics are to be
#' returned.  If `NULL` then all categories are included. (default: `NULL`)
#' @param group Character vector with groups (or group) for which category statistics are to be
#' returned. If `NULL` then all groups are included. (default: `NULL`)
#' @param statistic Which category statistic should be returned? `mu` for category mean or `Sigma` for category
#' covariance matrix, or `c("mu", "Sigma")` for both. (default: both)
#'
#' @return If just one group and category was requested, a vector (for the mean) or matrix (for the covariance
#' matrix). If more than one group or category was requested, a tibble with one row for each unique combination
#' of group and category.
#'
#' @seealso TBD
#' @keywords TBD
#' @references \insertRef{murphy2012}{MVBeliefUpdatr}
#' @examples
#' TBD
#' @rdname get_expected_category_statistic
#' @export
get_expected_category_statistic = function(x, category = NULL, group = NULL,
                                           statistic = c("mu", "Sigma")) {
  assert_that(all(statistic %in% c("mu", "Sigma")))
  assert_that(is.NIW_ibbu_stanfit(x) | is.NIW_belief_MCMC(x, is.nested = T, is.long = T))
  if (is.NIW_ibbu_stanfit(x))
    x = add_ibbu_stanfit_draws(x, which = "both", wide = F, nest = T)

  assert_that(any(is.null(category), is.character(category), is.numeric(category)))
  assert_that(any(is.null(group), is.character(group), is.numeric(group)))
  if (is.null(category)) category = unique(x$category)
  if (is.null(group)) group = unique(x$group)

  x %<>%
    filter(group %in% !! group, category %in% !! category) %>%
    mutate(Sigma = map2(S, nu, get_Sigma_from_S)) %>%
    group_by(group, category) %>%
    summarise(
      mu.mean = list(m %>% reduce(`+`) / length(m)),
      Sigma.mean = list(Sigma %>% reduce(`+`) / length(Sigma))) %>%
    select(group, category, !!! rlang::syms(paste0(statistic, ".mean")))

  if (!all(sort(unique(as.character(x$group))) == sort(as.character(group))))
    warning("Not all groups were found in x.")

  # If just one category and group was requested, just return that object, rather
  # than the tibble
  if (nrow(x) == 1) x = x[,paste0(statistic, ".mean")][[1]][[1]]
  return(x)
}

#' @rdname str(iistic
#' @export
get_expected_mu = function(x, category, group) {
  return(get_expected_category_statistic(x, category, group, statistic = "mu"))
}

#' @rdname get_expected_category_statistic
#' @export
get_expected_sigma = function(x, category, group) {
  return(get_expected_category_statistic(x, category, group, statistic = "Sigma"))
}
