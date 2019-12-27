#' @import assertthat purrr
#' @importFrom magrittr %<>%
#' @importFrom tidybayes spread_draws recover_types
#' @importFrom dplyr %>% mutate summarise
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

#' Add MCMC draws of IBBU parameters to a tibble.
#'
#' Add MCMC draws of all parameters from incremental Bayesian belief-updating (IBBU) to a tibble. Both wide
#' (wide=TRUE) or long format (wide=FALSE) can be chosen as output. By default all post-warmup draws are
#' returned, but if summarize=TRUE then just the mean of each parameter is returned instead. Users can
#' optionally provide a data frame with the unique values of the group, category, and cue variables (cue, cue2;
#' for example through \code{\link[tidyr]{crossing}}). If provided, this information will be used to recover
#' the types in the stanfit object, adding it to the resulting tibble with MCMC draws.
#'
#' @param fit mv-ibbu-stanfit object.
#' @param which Should parameters for the prior, posterior, or both be added? (default: posterior)
#' @param draws Vector with specific draw(s) to be returned, or NULL if all draws are to be returned. (default: NULL)
#' @param summarize Should the mean of the draws be returned instead of all of the draws? (default: FALSE)
#' @param wide Should all parameters be returned in one row? (default: FALSE)
#'
#' @return tibble with post-warmup (posterior) MCMC draws of the prior/posterior parameters of the IBBU model
#' (kappa, nu, M, S, lapse_rate). Kappa and nu are the pseudocounts that determine the strength of the beliefs
#' into the mean and covariance matrix, respectively. M is the mean of the multivariate normal distribution over category
#' means mu. S is the scatter matrix that determines both the covariance of the category means mu, and the
#' Inverse Wishart distribution over category covariance matrices Sigma.
#'
#' The expected value of the category mean mu is M. The expected value of the category covariance matrix Sigma
#' is S / (nu - D - 1), where D is the dimension of the multivariate Gaussian category. For details, see
#' Murphy (2012, p. 134).
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
#'
add_ibbu_draws = function(
  fit,
  which = c("prior", "posterior", "both")[2],
  draws = NULL,
  summarize = FALSE,
  wide = FALSE
) {
  assert_that(is.mv_ibbu_stanfit(fit))
  assert_that(which %in% c("prior", "posterior", "both"))
  assert_that(is.null(draws) | is.numeric(draws))
  assert_that(is.flag(summarize))
  assert_that(is.flag(wide))

  if (which == "both") {
    d.pars =
      rbind(
        add_ibbu_draws(fit = fit, which = "prior",
                          draws = if(!is.null(draws)) draws else NULL,
                          summarize = summarize, wide = wide),
        add_ibbu_draws(fit = fit, which = "posterior",
                          draws = if(!is.null(draws)) draws else NULL,
                          summarize = summarize, wide = wide)
      )
    return(d.pars)
  } else {
    # Parameters' names depend on whether prior or posterior is to be extracted.
    postfix = if (which == "prior") "_0" else "_n"
    kappa = paste0("kappa", postfix)
    nu = paste0("nu", postfix)
    mu = paste0("mu", postfix)
    sigma = paste0("sigma", postfix)

    # Variables by which parameters are indexed
    category = "category"
    group = "group"
    pars.index = if (which == "prior") category else c(category, group)

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

    # Make sure order of variables is identical for prior or posterior
    # (facilitates processing of the output of this function)
    if (which == "prior") d.pars %<>% mutate(subject = 0)
    d.pars %<>% select(.chain, .iteration, .draw, !! rlang::sym(group), category, kappa, nu, M, S, lapse_rate)

    if (wide) {
      if (which == "prior")
        d.pars %<>% gather(variable, value, c(mu, sigma))
      else
        d.pars %<>% gather(variable, value, c(kappa, nu, mu, sigma))

      d.pars %<>%
        unite(temp, !!! rlang::syms(pars.index), variable) %>%
        spread(temp, value)
    }

    return(d.pars)
  }
}






