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
#' Returns `ndraws` indices for random post-warmup MCMC draws (without replacement) from
#' `fit`.
#'
#' @param fit mv_ibbu_stanfit object.
#' @param ndraws Number of indices to be returned. Can't be larger than total number of
#' post-warmup samples across all MCMC chains in `fit`.
#'
#' @return A numeric vector.
get_random_draw_indices = function(fit, ndraws)
{
  n.all.draws = get_number_of_draws(fit)
  assert_that(ndraws <= n.all.draws,
              msg = paste0("Cannot return ", ndraws, " draws because there are only ", n.all.draws, " in the object."))

  draws = sample(1:n.all.draws, size = ndraws)
  return(draws)
}


#' Get the transform/untransform information from an NIW ideal adaptor stanfit.
#'
#' Returns the transform/untransform information handed to \code{stan} or \code{sampling} during the creation of the \code{stanfit}
#' object.
#'
#' @param fit \code{\link{NIW_ideal_adaptor_stanfit}} object.
#'
#' @return A function.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @rdname get_transform_information_from_stanfit
#' @export
get_transform_information_from_stanfit <- function(model) {
  assert_that(is.NIW_ideal_adaptor_stanfit(model))

  return(model@transform_information)
}

#' @rdname get_transform_information_from_stanfit
#' @export
get_transform_function_from_stanfit <- function(model) {
  return(get_transform_information_from_stanfit(model)$transform.function)
}

#' @rdname get_transform_information_from_stanfit
#' @export
get_untransform_function_from_stanfit <- function(model) {
  return(get_transform_information_from_stanfit(model)$untransform.function)
}



#' Get the input data from an NIW ideal adaptor stanfit.
#'
#' Returns the inputs handed to \code{stan} or \code{sampling} during the creation of the \code{stanfit}
#' object.
#'
#' @param fit \code{\link{NIW_ideal_adaptor_stanfit}} object.
#'
#' @return A list with element names and structures determined by the type of stanfit model.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @rdname get_ibbu_input
#' @export
get_input_from_stanfit = function(fit) {
  assert_that(is.NIW_ideal_adaptor_stanfit(fit))

  return(fit@input_data)
}


#' Get category sample mean or covariance matrix of exposure data from NIW IBBU stanfit.
#'
#' Returns the category means mu and/or category covariance matrix Sigma for the exposure data for an incremental
#' Bayesian belief-updating (IBBU) model from an NIW IBBU stanfit or NIW belief MCMC object.
#'
#' @param x \code{\link{NIW_ideal_adaptor_stanfit}} or NIW belief MCMC object.
#' @param category Character vector with categories (or category) for which category statistics are to be
#' returned.  If `NULL` then all categories are included. (default: `NULL`)
#' @param group Character vector with groups (or group) for which category statistics are to be
#' returned. If `NULL` then all groups are included. (default: `NULL`)
#' @param statistic Which exposure statistic should be returned? `n` for number of observations, `mean` for
#' category mean or `ss` for (uncentered) category sum-of-square matrix, or `c("mean", "ss")` for any combination
#' thereof. (default: all)
#'
#' @return If just one group and category was requested, a vector (for the mean) or matrix (for the covariance
#' matrix). If more than one group or category was requested, a tibble with one row for each unique combination
#' of group and category.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @rdname get_exposure_statistic_from_stanfit
#' @export
get_exposure_statistic_from_stanfit = function(x, category = NULL, group = NULL,
                                               statistic = c("n", "mean", "ss")) {
  assert_that(is.NIW_ideal_adaptor_input(x) | is.NIW_ideal_adaptor_stanfit(x))
  assert_that(all(statistic %in% c("n", "mean", "ss")),
              msg = "statistic must be one of 'mean' or 'ss'.")
  if (is.NIW_ideal_adaptor_stanfit(x)) x = get_input_from_stanfit(x)
  if (!is.null(category)) assert_that(all(category %in% unique(x$category)),
                                      msg = paste("Some categories were not found in the exposure data:",
                                                  paste(setdiff(category, unique(x$category)), collapse = ", ")))
  if (!is.null(group)) assert_that(all(group %in% unique(x$group)),
                                      msg = paste("Some groups were not found in the exposure data:",
                                                  paste(setdiff(group, unique(x$group)), collapse = ", ")))

  if ("n" %in% statistic) {
    stop("Not yet implementd for statistic = n.")
  }

  if ("mean" %in% statistic) {
    m <- x$x_mean
    d <- dim(m)
    dn <- dimnames(m)

    df.m <- tibble()
    for (c in 1:d[1]) { # category
      for (g in 1:d[2]) { # group/condition
        for (f in 1:d[3]) { # cue
          df.m <-
            rbind(
              df.m,
              tibble(
                group = dn[[2]][g],
                category = dn[[1]][c],
                cue = dn[[3]][f],
                value = m[c, g, f]))
        }
      }
    }

    df.m %<>%
      pivot_wider(names_from = "cue", values_from = "value") %>%
      make_vector_column(cols = dn[[3]], vector_col = "mean", .keep = "unused")

    df <- if (!is.null(df)) df %<>% left_join(df.m, by = c("group", "category")) else df.m
  }

  if ("ss" %in% statistic) {
    s <- x$x_ss
    d <- dim(s)
    dn <- dimnames(s)

    df.s <- tibble()
    for (c in 1:d[1]) { # category
      for (g in 1:d[2]) { # group/condition
        for (f1 in 1:d[3]) { # cue1
          for (f2 in 1:d[4]) { # cue2
            df.s <-
            rbind(
              df.s,
              tibble(
                group = dn[[2]][g],
                category = dn[[1]][c],
                cue = dn[[3]][f1],
                cue2 = dn[[4]][f2],
                value = s[c, g, f1, f2]))
          }
        }
      }
    }

    df.s %<>%
      group_by(category, group) %>%
      summarise(ss = list(matrix(value, nrow = sqrt(length(value)))))

    df <- if (!is.null(df)) df %<>% left_join(df.s, by = c("group", "category")) else df.s
  }

  df %<>%
    { if (!is.null(group)) filter(., group %in% group) else . } %>%
    { if (!is.null(category)) filter(., category %in% group) else . }

  return(df)
}

#' @rdname get_exposure_statistic_from_stanfit
#' @export
get_exposure_mean_from_stanfit = function(x, category, group) {
  return(get_exposure_statistic_from_stanfit(x, category, group, statistic = "mean"))
}

#' @rdname get_exposure_statistic_from_stanfit
#' @export
get_exposure_ss_from_stanfit = function(x, category, group) {
  return(get_exposure_statistic_from_stanfit(x, category, group, statistic = "ss"))
}


#' Get the test data from an NIW ideal adaptor stanfit.
#'
#' Returns the test data used during the creation of the \code{\link[rstan]{stanfit}}.
#' object.
#'
#' @param x \code{\link{NIW_ideal_adaptor_stanfit}} object.
#'
#' @return A \code{tibble} in which each row is a test token. Columns include the cues
#' and the response counts (one column per category) for all test tokens and all groups.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
get_test_data_from_stanfit = function(fit) {
  data <- get_input_from_stanfit(fit)
  data[["x_test"]] %>%
    cbind(data[["z_test_counts"]]) %>%
    mutate(
      group.id = data[["y_test"]],
      group = factor(attr(data[["y_test"]], "levels")[group.id],
                     levels = attr(data[["y_test"]], "levels")))
}




#' Get or restore the original group or category levels from an NIW ideal adaptor stanfit.
#'
#' Checks if information is available about the original values and order of the factor levels
#' for the category variable (for which beliefs about means and covariances are inferred) or
#' group variable (e.g., subject or exposure group), respectively. If available,
#' that information is returned. `get_category_levels_from_stanfit()` and `get_group_levels_from_stanfit()` are
#' convenience functions, calling `get_original_variable_levels_from_stanfit()`.
#'
#' @param fit \code{\link{NIW_ideal_adaptor_stanfit}} object.
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
#' @rdname get_original_variable_levels_from_stanfit
#' @export
get_original_variable_levels_from_stanfit = function(fit, variable = c("category", "group", "cue"), indices = NULL) {
  assert_that(is.null(indices) | all(indices > 0))
  f = get_constructor(fit, variable)

  if (is.null(indices)) return(levels(f(c()))) else return(f(indices))
}


#' @rdname get_original_variable_levels_from_stanfit
#' @export
get_category_levels_from_stanfit = function(fit, indices = NULL) {
  return(get_original_variable_levels_from_stanfit(fit, "category", indices))
}

#' @rdname get_original_variable_levels_from_stanfit
#' @export
get_group_levels_from_stanfit = function(fit, indices = NULL, include_prior = F) {
  groups = get_original_variable_levels_from_stanfit(fit, "group", indices)
  if (include_prior) groups = append(groups, "prior")

  return(groups)
}

#' @rdname get_original_variable_levels_from_stanfit
#' @export
get_cue_levels_from_stanfit = function(fit, indices = NULL) {
  cues = get_original_variable_levels_from_stanfit(fit, "cue", indices)

  return(cues)
}


#' Get tidybayes constructor from an NIW IBBU stanfit.
#'
#' Gets the tidybayes constructor function from the stanfit object. `get_category_constructor()` and
#' `get_group_constructor()` are convenience functions, calling `get_constructor()`. See \code{
#' \link[tidybayes]{recover_types}}. If variable is
#'
#' @param fit \code{\link{NIW_ideal_adaptor_stanfit}} object.
#' @param variable Either "category" or "group". If set to `NULL` then a list of all constructors is
#' returned. That list is `NULL` if not tidybayes constructors are found in fit. (default: c("category", "group"))
#'
#' @return A constructor function, a list of constructor functions, or `NULL`. If a specific constructor
#' function is requested but not found, a warning is shown.
#'
#' @seealso \code{\link[tidybayes]{recover_types}} from tidybayes, \code{\link{get_original_variable_levels_from_stanfit}}
#' @keywords TBD
#' @examples
#' TBD
#' @rdname get_constructor
#' @export
get_constructor = function(fit, variable = NULL) {
  available_constructors <- c("category", "group", "cue", "cue2")
  assert_NIW_ideal_adaptor_stanfit(fit)
  if (is.null(variable)) return(attr(fit, "tidybayes_constructors"))

  assert_that(variable %in% available_constructors,
              msg = paste0("Variable name must be one of ", paste(available_constructors, collapse = "or"), "."))

  if (is.null(attr(fit, "tidybayes_constructors")[[rlang::sym(variable)]])) {
    warning(paste0(class_name, " object does not contain type information about the variable ", variable,
                   ". Applying tidybayes::recover_types() to the object might fix this."))
    return(NULL)
  }

  f <- attr(fit, "tidybayes_constructors")[[rlang::sym(variable)]]

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

#' @rdname get_constructor
#' @export
get_cue_constructor = function(fit) {
  return(get_constructor(fit, "cue"))
}

#' @rdname get_constructor
#' @export
get_cue2_constructor = function(fit) {
  return(get_constructor(fit, "cue2"))
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
#' @param x An \code{\link[=is.NIW_ideal_adaptor_stanfit]{mv_ibbu_stanfit}} or \code{\link[=NIW_ideal_adaptor_MCMC]{NIW_ideal_adaptor_MCMC}} object.
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
#' @rdname get_expected_category_statistic_from_stanfit
#' @export
get_expected_category_statistic_from_stanfit = function(
  x,
  category = NULL,
  group = NULL,
  statistic = c("mu", "Sigma"),
  untransform_cues = T
) {
  assert_that(all(statistic %in% c("mu", "Sigma")))
  if (is.NIW_ideal_adaptor_stanfit(x)) {
    x = add_ibbu_stanfit_draws(x, which = "both", wide = F, nest = T, untransform_cues = untransform_cues)
  } else if (is.NIW_ideal_adaptor_MCMC(x)) {
    assert_that(is.NIW_ideal_adaptor_MCMC(x, is.nested = T, is.long = T),
                msg = "If x is an NIW_ideal_adaptor_MCMC object, it must be in nested long format.")
    assert_that(!untransform_cues,
                msg = "If untransform_cues = T, then x must be an NIW_ideal_adaptor_stanfit object.")
  }

  assert_that(any(is.null(category), is.character(category), is.numeric(category)))
  assert_that(any(is.null(group), is.character(group), is.numeric(group)))
  if (is.null(category)) category = unique(x$category)
  if (is.null(group)) group = unique(x$group)

  x %<>%
    filter(group %in% !! group, category %in% !! category) %>%
    mutate(Sigma = get_expected_Sigma_from_S(S, nu)) %>%
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

#' @rdname get_expected_category_statistic_from_stanfit
#' @export
get_expected_mu_from_stanfit = function(x, category, group, ...) {
  return(get_expected_category_statistic_from_stanfit(x, category, group, statistic = "mu", ...))
}

#' @rdname get_expected_category_statistic_from_stanfit
#' @export
get_expected_sigma_from_stanfit = function(x, category, group, ...) {
  return(get_expected_category_statistic_from_stanfit(x, category, group, statistic = "Sigma", ...))
}



#' @rdname get_NIW_categorization_function
#' @export
get_categorization_function_from_grouped_ibbu_stanfit_draws = function(fit, ...) {
  get_NIW_categorization_function(
    ms = fit$m,
    Ss = fit$S,
    kappas = fit$kappa,
    nus = fit$nu,
    lapse_rate = unique(unlist(fit$lapse_rate)),
    ...
  )
}


#' Add MCMC draws from an NIW IBBU stanfit to a tibble.
#'
#' Add MCMC draws of all parameters from incremental Bayesian belief-updating (IBBU) to a tibble. Both wide
#' (`wide=TRUE`) or long format (`wide=FALSE`) can be chosen as output. By default all post-warmup draws are
#' returned, but if `summarize=TRUE` then just the mean of each parameter is returned instead.
#'
#' By default, the category means and scatter matrices are nested, rather than each of their elements being
#' stored separately (`nest=TRUE`).
#'
#' @param fit \code{\link{NIW_ideal_adaptor_stanfit}} object.
#' @param which Should parameters for the prior, posterior, or both be added? (default: `"posterior"`)
#' @param ndraws Number of random draws or `NULL` if all draws are to be returned. Only `draws` or `ndraws` should be non-zero. (default: `NULL`)
#' @param draws Vector with specific draw(s) to be returned, or `NULL` if all draws are to be returned. (default: `NULL`)
#' @param untransform_cues Should m_0 and S_0 be transformed back into the original cue space? (default: `TRUE`)
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
  ndraws = NULL,
  draws = NULL,
  untransform_cues = TRUE,
  summarize = FALSE,
  wide = FALSE,
  nest = TRUE,
  category = "category",
  group = "group"
) {
  assert_NIW_ideal_adaptor_stanfit(fit)
  assert_that(which %in% c("prior", "posterior", "both"),
              msg = "which must be one of 'prior', 'posterior', or 'both'.")
  assert_that(any(is.null(ndraws), is.count(ndraws)),
              msg = "If not NULL, ndraw must be a count.")
  assert_that(any(all(is.null(draws), is.null(ndraws)), xor(!is.null(draws), !is.null(ndraws))),
              msg = "Only one of draws and ndraws can be non-NULL.")
  assert_that(any(is.null(draws), all(draws > 0)),
              msg = "If not NULL draws, must be a vector of positive integers.")
  assert_that(is.flag(summarize))
  assert_that(is.flag(wide))
  assert_that(!all(wide, !nest),
              msg = "Wide format is currently not implemented without nesting.")

  if (!is.null(ndraws)) draws = get_random_draw_indices(fit, ndraws)
  if (which == "both") {
    d.prior <-
      add_ibbu_stanfit_draws(
        fit = fit, which = "prior",
        ndraws = NULL, draws = draws,
        untransform_cues = untransform_cues,
        summarize = summarize, wide = wide, nest = nest)
    d.posterior <-
      add_ibbu_stanfit_draws(
        fit = fit, which = "posterior",
        ndraws = NULL, draws = draws,
        untransform_cues = untransform_cues,
        summarize = summarize, wide = wide, nest = nest)
    d.pars <-
      rbind(d.prior, d.posterior) %>%
      mutate(!! rlang::sym(group) :=
               factor(!! rlang::sym(group),
                      levels = c(with(d.prior, levels(!! rlang::sym(group))),
                                 with(d.posterior, levels(!! rlang::sym(group))))))
    return(d.pars)
  } else {
    assert_NIW_ideal_adaptor_stanfit(fit)

    # Parameters' names depend on whether prior or posterior is to be extracted.
    postfix <- if (which == "prior") "_0" else "_n"
    kappa <- paste0("kappa", postfix)
    nu <- paste0("nu", postfix)
    m <- paste0("m", postfix)
    S <- paste0("S", postfix)

    # Variables by which parameters are indexed
    pars.index <- if (which == "prior") category else c(category, group)

    # Get non-nested draws
    if (which == "prior") {
      d.pars <-
        fit %>%
        spread_draws(
          !! rlang::sym(kappa),
          !! rlang::sym(nu),
          (!! rlang::sym(m))[!!! rlang::syms(pars.index), cue],
          (!! rlang::sym(S))[!!! rlang::syms(pars.index), cue, cue2],
          lapse_rate,
          ndraws = ndraws) %>%
        { if (!is.null(draws)) filter(., .draw %in% draws) else . }
    } else {
      d.pars <-
        fit %>%
        spread_draws(
          (!! rlang::sym(kappa))[!!! rlang::syms(pars.index)],
          (!! rlang::sym(nu))[!!! rlang::syms(pars.index)],
          (!! rlang::sym(m))[!!! rlang::syms(pars.index), cue],
          (!! rlang::sym(S))[!!! rlang::syms(pars.index), cue, cue2],
          lapse_rate,
          ndraws = ndraws) %>%
        { if (!is.null(draws)) filter(., .draw %in% draws) else . }
    }

    d.pars %<>%
      rename_at(vars(ends_with(postfix)), ~ sub(postfix, "", .)) %>%
      { if (summarize) {
        group_by(., !!! syms(pars.index), cue, cue2) %>%
          summarise_at(., vars(kappa, nu, m, S, lapse_rate), mean) %>%
          mutate_at(., vars(.chain, .iteration, .draw), ~ "all")
      } else . } %>%
      # If group is prior, then add the group variable with value "prior" to d.pars first.
      { if (which == "prior") mutate(., (!! rlang::sym(group)) := "prior") else . } %>%
      # Make sure order of variables is identical for prior or posterior (facilitates processing
      # of the output of this function).
      select(.chain, .iteration, .draw,
             !! rlang::sym(group), !! rlang::sym(category),
             kappa, nu, m, S, lapse_rate)

    if (untransform_cues) {
      d.pars %<>%
        nest_cue_information_in_model() %>%
        untransform_model(transform = fit@transform_information)

      if (!nest) d.pars %<>% unnest_cue_information_in_model()
    } else if (nest) {
      d.pars %<>%
        nest_cue_information_in_model()
    }

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






