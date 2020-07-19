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

#' Get number of post-warmup MCMC samples
#'
#' Get the total number of post-warmup MCMC samples in `fit`.
#' @export
get_number_of_draws = function(fit) {
  return(length(fit@sim$samples[[1]][[1]]))
}

#' Get indices for random MCMC draws
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
  f = get_constructor(variable)

  if (is.null(indices)) return(levels(f(c()))) else return(f(indices))
}


#' @rdname get_original_levels
#' @export
get_category_levels = function(fit, indices = NULL) {
  assert_that(is.mv_ibbu_stanfit(fit))
  assert_that(is.null(indices) | all(indices > 0))

  return(get_original_levels(fit, "category", indices))
}

#' @rdname get_original_levels
#' @export
get_group_levels = function(fit, indices = NULL) {
  assert_that(is.mv_ibbu_stanfit(fit))
  assert_that(is.null(indices) | all(indices > 0))

  return(get_original_levels(fit, "group", indices))
}



#' Get tidybayes constructor from IBBU stanfit
#'
#' Gets the tidybayes constructor function from the stanfit object. `get_category_constructor()` and
#' `get_group_constructor()` are convenience functions, calling `get_constructor()`. See \code{
#' \link[tidybayes]{recover_types}}.
#'
#' @param fit mv_ibbu_stanfit object.
#' @param variable Either "category" or "group".
#'
#' @return A constructor function.
#'
#' @seealso \code{\link[tidybayes]{recover_types}}, \code{\link{get_original_levels}}
#' @keywords TBD
#' @examples
#' TBD
#' @rdname get_constructor
#' @export
get_constructor = function(fit, variable = c("category", "group")) {
  assert_that(variable %in% c("category", "group"))

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
  assert_that(is.mv_ibbu_stanfit(fit))

  return(get_constructor(fit, "category"))
}

#' @rdname get_constructor
#' @export
get_group_constructor = function(fit) {
  assert_that(is.mv_ibbu_stanfit(fit))

  return(get_constructor(fit, "group"))
}


#' Get category mean mu or covariance matrix sigma from a tibble of cues.
#'
#' Returns the category means mu and/or category covariance matrix Sigma for the exposure data.
#'
#' @param x An mv_ibbu_stanfit or mv_ibbu_MCMC object.
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
#' The category means mu and/or category covariance matrix Sigma for the exposure data for an incremental
#' Bayesian belief-updating (IBBU) model.
#'
#' @param x An mv_ibbu_stanfit or mv_ibbu_MCMC object.
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
  assert_that(is.mv_ibbu_input(x) | is.mv_ibbu_stanfit(x))
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
  assert_that(is.mv_ibbu_input(x) | is.mv_ibbu_stanfit(x))

  if (is.mv_ibbu_input(x)) return(x) else {
    stop("Extraction of input data from MV IBBU stanfit not yet implemented!")
  }

  return(x)
}



#' Add MCMC draws of IBBU parameters to a tibble.
#'
#' Add MCMC draws of all parameters from incremental Bayesian belief-updating (IBBU) to a tibble. Both wide
#' (`wide=TRUE`) or long format (`wide=FALSE`) can be chosen as output. By default all post-warmup draws are
#' returned, but if `summarize=TRUE` then just the mean of each parameter is returned instead. Users can
#' optionally provide a data frame with the unique values of the `group`, `category`, and cue variables (`cue`, `cue2`;
#' for example through \code{\link[tidyr]{crossing}}). If provided, this information will be used to recover
#' the types in the stanfit object, adding it to the resulting tibble with MCMC draws.
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
#' (\code{kappa, nu, M, S, lapse_rate}). \code{kappa} and \code{nu} are the pseudocounts that determine the strength of the beliefs
#' into the mean and covariance matrix, respectively. \code{M} is the mean of the multivariate normal distribution over category
#' means mu. \code{S} is the scatter matrix that determines both the covariance of the category means mu, and the
#' Inverse Wishart distribution over category covariance matrices Sigma.
#'
#' The expected value of the category mean mu is \code{M}. The expected value of the category covariance matrix Sigma
#' is \code{S / (nu - D - 1)}, where \code{D} is the dimension of the multivariate Gaussian category. For details, see
#' Murphy (2012, p. 134).
#'
#' @seealso TBD
#' @keywords TBD
#' @references Murphy, K. P. (2012). Machine learning: a probabilistic perspective. MIT press.
#' @examples
#' TBD
#' @export
#'
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
    assert_that(is.mv_ibbu_stanfit(fit))

    # Parameters' names depend on whether prior or posterior is to be extracted.
    postfix = if (which == "prior") "_0" else "_n"
    kappa = paste0("kappa", postfix)
    nu = paste0("nu", postfix)
    M = paste0("mu", postfix)
    S = paste0("sigma", postfix)

    # Variables by which parameters are indexed
    pars.index = if (which == "prior") category else c(category, group)

    # Should M and S be nested?
    if (!nest) {
      if (which == "prior")
        d.pars = fit %>%
          spread_draws(
            !! rlang::sym(kappa),
            !! rlang::sym(nu),
            (!! rlang::sym(M))[!!! rlang::syms(pars.index), cue],
            (!! rlang::sym(S))[!!! rlang::syms(pars.index), cue, cue2],
            lapse_rate,
            n = draws
          )
      else
        d.pars = fit %>%
          spread_draws(
            (!! rlang::sym(kappa))[!!! rlang::syms(pars.index)],
            (!! rlang::sym(nu))[!!! rlang::syms(pars.index)],
            (!! rlang::sym(M))[!!! rlang::syms(pars.index), cue],
            (!! rlang::sym(S))[!!! rlang::syms(pars.index), cue, cue2],
            lapse_rate,
            n = draws
          )

      message("Currently mv_ibbu_stanfits have mu_{0,n} and sigma_{0,n} as parameter names. These are actually M_{0,n} and S_{0,n}.\n
              add_ibbu_stanfit_draws() renames these parameters to M and S.")
      # See also naming of parameters at beginning of this else block (no other part of the code needs to change.)
      d.pars %<>%
        rename(M = !! rlang::sym(M), S = !! rlang::sym(S)) %>%
        rename_at(vars(ends_with(postfix)), ~ sub(postfix, "", .))

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
        # Join in M
        left_join(
          fit %>%
            tidybayes::spread_draws((!! rlang::sym(M))[!!! rlang::syms(pars.index), cue]) %>%
            { if (!is.null(draws)) filter(., .draw %in% draws) else . } %>%
            { if (summarize)
              group_by(., !!! rlang::syms(pars.index), cue) %>%
                # Obtain expected mean of category mean, M_0 or M_N
                dplyr::summarise(., !! rlang::sym(M) := mean(!! rlang::sym(M))
                ) %>%
                mutate(.,
                       .chain = "all", .iteration = "all", .draw = "all"
                ) else . } %>%
            group_by(.chain, .iteration, .draw, !!! rlang::syms(pars.index)) %>%
            summarise(M = list(
              matrix((!! rlang::sym(M)),
                     dimnames = list(unique(cue), NULL),
                     byrow = T,
                     nrow = length((!! rlang::sym(M))))))
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
                       starts_with("kappa"), starts_with("nu"), starts_with("M"), starts_with("S"),
                       lapse_rate)

    if (wide) {
      if (which == "prior")
        d.pars %<>% gather(variable, value, c(M, S))
      else
        d.pars %<>% gather(variable, value, c(kappa, nu, M, S))

      d.pars %<>%
        unite(temp, !!! rlang::syms(pars.index), variable) %>%
        spread(temp, value)
    }

    # Final clean-up
    d.pars %<>%
      ungroup() %>%
      # Make sure that group and category are factors (even if group or category ids are just numbers)
      mutate_at(vars(!! rlang::sym(group), !! rlang::sym(category)), factor)

    return(d.pars)
  }
}






