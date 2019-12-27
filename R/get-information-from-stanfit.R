#' @import assertthat purrr
#' @importFrom magrittr %<>%
#' @importFrom tidybayes spread_draws recover_types
#' @importFrom dplyr %>% select filter mutate summarise left_join rename group_by ungroup
#' @importFrom tidyr spread gather unite
#' @importFrom tibble tibble is_tibble
#' @importFrom rlang !! !!! sym syms expr
NULL

get_params = function(fit) {
  return(fit@model_pars)
}

get_number_of_draws = function(fit) {
  return(length(fit@sim$samples[[1]][[1]]))
}

#' Get or restore the original group or category levels.
#'
#' Checks if information is available about the original values and order of the factor levels
#' for the category variable (for which beliefs about means and covariances are inferred) or
#' group variable (e.g., subject or exposure group), respectively. If available,
#' that information is returned. \code{get_category_levels()} and \code{get_group_levels()} are
#' convenience functions, calling \code{get_original_levels()}.
#'
#' @param fit mv-ibbu-stanfit object.
#' @param variable Either \code{"category"} or \code{"group"}.
#' @param indeces A vector of category or group indices that should be turned into the original
#' category levels, or NULL if only the unique levels in their original order (as vector of characters)
#' should be returned. (default: NULL)
#'
#' @return If no category or group indices are provided, the levels of the category/group are returned (in the
#' original order). Otherwise a vector of the same length as \code{indices}
#'
#' @seealso
#' @keywords TBD
#' @examples
#' TBD
#' @rdname get_original_levels
#' @export
get_original_levels = function(fit, variable = c("category", "group"), indices = NULL) {
  assert_that(variable %in% c("category", "group"))

  if (is.null(attr(fit, "constructors")[[rlang::sym(variable)]])) {
    warning(paste0(class_name, " object does not contain type information about the ", variable,
                   " variable. Consider applying recover_types() from the tidybayes package first."))
    return(NULL)
  }

  f = attr(fit, "constructors")[[rlang::sym(variable)]]

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


#' Add MCMC draws of IBBU parameters to a tibble.
#'
#' Add MCMC draws of all parameters from incremental Bayesian belief-updating (IBBU) to a tibble. Both wide
#' (\code{wide=TRUE}) or long format (\code{wide=FALSE}) can be chosen as output. By default all post-warmup draws are
#' returned, but if \code{summarize=TRUE} then just the mean of each parameter is returned instead. Users can
#' optionally provide a data frame with the unique values of the group, category, and cue variables (cue, cue2;
#' for example through \code{\link[tidyr]{crossing}}). If provided, this information will be used to recover
#' the types in the stanfit object, adding it to the resulting tibble with MCMC draws.
#'
#' By default, the category means and scatter matrices are nested, rather than each of their elements being
#' stored separately (\code{nest=TRUE}).
#'
#' @param fit mv-ibbu-stanfit object.
#' @param which Should parameters for the prior, posterior, or both be added? (default: posterior)
#' @param draws Vector with specific draw(s) to be returned, or NULL if all draws are to be returned. (default: NULL)
#' @param summarize Should the mean of the draws be returned instead of all of the draws? (default: FALSE)
#' @param wide Should all parameters be returned in one row? (default: FALSE)
#' @param nest Should the category mean vectors and scatter matrices be nested into one cell each, or should each element
#' be stored in a separate cell? (defaul: TRUE)
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
add_ibbu_draws = function(
  fit,
  which = c("prior", "posterior", "both")[2],
  draws = NULL,
  summarize = FALSE,
  wide = FALSE,
  nest = TRUE
) {
  assert_that(is.mv_ibbu_stanfit(fit))
  assert_that(which %in% c("prior", "posterior", "both"))
  assert_that(any(is.null(draws), all(draws > 0)))
  assert_that(is.flag(summarize))
  assert_that(is.flag(wide))
  assert_that(!all(wide, !nest),
              msg = "Wide format is currently not implemented without nesting.")

  category = "category"
  group = "group"

  if (which == "both") {
    d.prior = add_ibbu_draws(fit = fit, which = "prior",
                             draws = if(!is.null(draws)) draws else NULL,
                             summarize = summarize, wide = wide)
    d.posterior = add_ibbu_draws(fit = fit, which = "posterior",
                   draws = if(!is.null(draws)) draws else NULL,
                   summarize = summarize, wide = wide)
    d.pars = rbind(d.prior, d.posterior) %>%
      mutate(!! rlang::sym(group) :=
               factor(!! rlang::sym(group),
                      levels = c(with(d.prior, levels(!! rlang::sym(group))),
                                 with(d.posterior, levels(!! rlang::sym(group))))))

    return(d.pars)
  } else {
    # Parameters' names depend on whether prior or posterior is to be extracted.
    postfix = if (which == "prior") "_0" else "_n"
    kappa = paste0("kappa", postfix)
    nu = paste0("nu", postfix)
    mu = paste0("mu", postfix)
    sigma = paste0("sigma", postfix)

    # Variables by which parameters are indexed
    pars.index = if (which == "prior") category else c(category, group)

    # Should M and S be nested?
    if (!nest) {
      d.pars = fit %>%
        spread_draws(
          { if (which == "prior") !! rlang::sym(kappa) else !! rlang::sym(kappa)[!!! rlang::syms(pars.index)]},
          { if (which == "prior") !! rlang::sym(nu) else !! rlang::sym(nu)[!!! rlang::syms(pars.index)]},
          !! rlang::sym(mu)[(!!! rlang::syms(pars.index)), cue],
          !! rlang::sym(sigma)[(!!! rlang::syms(pars.index)), cue, cue2],
          lapse_rate
        ) %>%
        rename_at(vars(starts_with("mu")), funs(gsub("^mu_(0|n)", "M", ., perl = T))) %>%
        rename_at(vars(starts_with("sigma")), funs(gsub("^mu_(0|n)", "S", ., perl = T)))
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
        # Join in mu
        left_join(
          fit %>%
            tidybayes::spread_draws((!! rlang::sym(mu))[!!! rlang::syms(pars.index), cue]) %>%
            { if (!is.null(draws)) filter(., .draw %in% draws) else . } %>%
            { if (summarize)
              group_by(., !!! rlang::syms(pars.index), cue) %>%
                # Obtain expected mean M_0 or M_N (this is not mu, although we use that name here)
                dplyr::summarise(., !! rlang::sym(mu) := mean(!! rlang::sym(mu))
                ) %>%
                mutate(.,
                       .chain = "all", .iteration = "all", .draw = "all"
                ) else . } %>%
            group_by(.chain, .iteration, .draw, !!! rlang::syms(pars.index)) %>%
            summarise(M = list(
              matrix((!! rlang::sym(mu)),
                     dimnames = list(unique(cue), NULL),
                     byrow = T,
                     nrow = length((!! rlang::sym(mu))))))
        ) %>%
        # Join in sigma
        left_join(
          fit %>%
            tidybayes::spread_draws((!! rlang::sym(sigma))[!!! rlang::syms(pars.index), cue, cue2]) %>%
            { if (!is.null(draws)) filter(., .draw %in% draws) else . } %>%
            { if (summarize)
              group_by(., !!! rlang::syms(pars.index), cue, cue2) %>%
                # Obtain expected co-variance matrix S_0 or S_N (this is not sigma, although we use that name here)
                dplyr::summarise(., !! rlang::sym(sigma) := mean(!! rlang::sym(sigma))
                ) %>%
                mutate(.,
                       .chain = "all", .iteration = "all", .draw = "all"
                ) else . } %>%
            group_by(.chain, .iteration, .draw, !!! rlang::syms(pars.index)) %>%
            summarise(S =
                        list(
                          matrix((!! rlang::sym(sigma)),
                                 dimnames = list(unique(cue), unique(cue2)),
                                 byrow = T,
                                 nrow = sqrt(length((!! rlang::sym(sigma)))))))
        )
    }

    # Make sure order of variables is identical for prior or posterior (facilitates processing of the
    # output of this function). For this we first add the group variable as a column if we're dealing
    # with the prior (which has no group variable since it's the same across groups). Then we sort the
    # columns.
    if (which == "prior") d.pars %<>% mutate((!! rlang::sym(group)) := factor("prior"))
    d.pars %<>% select(.chain, .iteration, .draw,
                       !! rlang::sym(group), !! rlang::sym(category),
                       kappa, nu, M, S, lapse_rate)

    if (wide) {
      if (which == "prior")
        d.pars %<>% gather(variable, value, c(M, S))
      else
        d.pars %<>% gather(variable, value, c(kappa, nu, M, S))

      d.pars %<>%
        unite(temp, !!! rlang::syms(pars.index), variable) %>%
        spread(temp, value)
    }

    return(d.pars)
  }
}






