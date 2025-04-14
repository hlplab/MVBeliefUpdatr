#' DEPRECATED: get_transform_information_from_stanfit
#' @export
get_transform_information_from_stanfit <- function(...) get_transform_information.NIW_ideal_adaptor_stanfit(...)

#' DEPRECATED: get_transform_function_from_stanfit
#' @export
get_transform_function_from_stanfit <- function(...) get_transform_function.NIW_ideal_adaptor_stanfit(...)

#' DEPRECATED: get_untransform_function_from_stanfit
#' @export
get_untransform_function_from_stanfit <- function(...) get_untransform_function.NIW_ideal_adaptor_stanfit(...)

#' DEPRECATED: get_staninput_from_stanfit
#' @export
get_staninput_from_stanfit <- function(...) get_staninput.NIW_ideal_adaptor_stanfit(...)

# get_exposure_category_statistic_from_stanfit <- get_exposure_category_statistic.NIW_ideal_adaptor_stanfit
# get_exposure_mean_from_stanfit <- get_exposure_category_mean.NIW_ideal_adaptor_stanfit
# get_exposure_css_from_stanfit <- get_exposure_category_css.NIW_ideal_adaptor_stanfit
# get_exposure_uss_from_stanfit <- get_exposure_category_uss.NIW_ideal_adaptor_stanfit
# get_exposure_cov_from_stanfit <- get_exposure_category_cov.NIW_ideal_adaptor_stanfit

#' DEPRECATED: get_test_data_from_stanfit
#' @export
get_test_data_from_stanfit <- function(...) get_test_data(...)

# get_original_variable_levels_from_stanfit <- get_staninput_variable_levels
# get_category_levels_from_stanfit <- get_category_levels
# get_group_levels_from_stanfit <- get_group_levels
# get_cue_levels_from_stanfit <- get_cue_levels
#
# get_expected_category_statistic_from_stanfit <- get_expected_category_statistic
# get_expected_mu_from_stanfit <- get_expected_mu
# get_expected_sigma_from_stanfit <- get_expected_sigma

#' DEPRECATED: add_ibbu_stanfit_draw
#' @export
add_ibbu_stanfit_draw <- function(...) get_draws(...)

#' DEPRECATED: Infer prior beliefs
#'
#' Use \code{\link{infer_NIW_ideal_adaptor()}} instead, together with \code{\link{make_staninput_for_NIW_ideal_adaptor()}}.
#' @inheritParams make_staninput
#' @inheritParams fit_NIW_ideal_adaptor
#' @export
infer_prior_beliefs <- function(
  # arguments for make_staninput
  exposure, test,
  cues, category, response,
  group, group.unique = NULL,
  center.observations = TRUE, scale.observations = FALSE, pca.observations = FALSE, pca.cutoff = 1, transform_type = NULL,
  lapse_rate = NULL, mu_0 = NULL, Sigma_0 = NULL,
  tau_scale = NULL, L_omega_eta = 1,
  split_loglik_per_observation = 0,
  stanmodel = NULL,
  verbose = FALSE,
  silent = 1,
  ...
) {
  if (verbose) message("Entering verbose mode.")

  # Currently the make_staninput function is creating both the transforms *and* the data.
  # That's a bit confusing and should probably be split up in the future into separate
  # functions.
  if (any(!is.null(center.observations), !is.null(scale.observations), !is.null(pca.observations))) {
    staninput <-
      make_staninput_deprecated(
        exposure = exposure,
        test = test,
        cues = cues,
        category = category,
        response = response,
        group = group,
        group.unique = group.unique,
        center.observations = center.observations,
        scale.observations = scale.observations,
        pca.observations = pca.observations,
        pca.cutoff = pca.cutoff,
        lapse_rate = lapse_rate,
        mu_0 = mu_0,
        Sigma_0 = Sigma_0,
        tau_scale = tau_scale,
        L_omega_eta = L_omega_eta,
        split_loglik_per_observation = split_loglik_per_observation,
        use_univariate_updating = if (is.null(stanmodel)) { FALSE } else { stanmodel == 'uvg_conj_uninformative_priors_sufficient_stats_lapse'},
        verbose = verbose)
  } else  if (!is.null(transform_type)) {
    staninput <-
      make_staninput(
        exposure = exposure,
        test = test,
        cues = cues,
        category = category,
        response = response,
        group = group,
        group.unique = group.unique,
        transform_type = transform_type,
        lapse_rate = lapse_rate,
        mu_0 = mu_0,
        Sigma_0 = Sigma_0,
        tau_scale = if (is.null(tau_scale)) rep(5, length(cues)) else tau_scale,
        L_omega_eta = L_omega_eta,
        split_loglik_per_observation = split_loglik_per_observation,
        use_univariate_updating = if (is.null(stanmodel)) { FALSE } else { stanmodel == 'uvg_conj_uninformative_priors_sufficient_stats_lapse'},
        verbose = verbose,
        model_type = "NIW_ideal_adaptor")
  } else {
    stop2("Either transform or center.observations, scale.observations, or pca.observations must be specified.")
  }

  fit_NIW_ideal_adaptor(staninput = staninput, stanmodel = stanmodel, silent = silent, verbose = verbose, ...)
}

#' DEPRECATED: make_staninput
#'
#' @param center.observations Should the data be centered based on cues' means during exposure? Note that the cues' means
#' used for centering are calculated after aggregating the data to all unique combinations specified by \code{group.unique}.
#' These means are only expected to be the same as the standard deviations over the entire exposure data if the exposure data
#' are perfectly balanced with regard to \code{group.unique}. Centering will not affect the inferred correlation
#' or covariance matrices but it will affect the absolute position of the inferred means. The relative position of the inferred
#' means remains unaffected.
#' If \code{TRUE} and \code{mu_0} is specified, \code{mu_0} will also be centered (\code{Sigma_0} is not affected by centering and thus not changed).
#' (default: \code{TRUE})
#' @param scale.observations Should the data be standardized based on cues' standard deviation during exposure? Note that the
#' cues' standard deviations used for scaling are calculated after aggregating the data to all unique combinations specified
#' by \code{group.unique}. These standard deviations are only expected to be the same as the standard deviations over the entire
#' exposure data if the exposure data are perfectly balanced with regard to \code{group.unique}.
#' Scaling will not affect the inferred correlation matrix,
#' but it will affect the inferred covariance matrix because it affects the inferred standard deviations. It will also
#' affect the absolute position of the inferred means. The relative position of the inferred means remains unaffected.
#' If \code{TRUE} and \code{mu_0} and \code{Sigma_0} are specified, \code{mu_0} and \code{Sigma_0} will also be scaled.
#' (default: `FALSE`)
#' @param pca.observations Should the data be transformed into orthogonal principal components? (default: \code{FALSE})
#' @param pca.cutoff Determines which principal components are handed to the MVBeliefUpdatr Stan program: all
#' components necessary to explain at least the pca.cutoff of the total variance. (default: .95) Ignored if
#' \code{pca.observation = FALSE}. (default: 1)
make_staninput_deprecated <- function(
    exposure, test,
    cues, category = "category", response = "response",
    group = "group", group.unique = NULL,
    center.observations = T, scale.observations = F, pca.observations = F, pca.cutoff = 1,
    lapse_rate = NULL, mu_0 = NULL, Sigma_0 = NULL,
    tau_scale = NULL,
    L_omega_eta = 1,
    split_loglik_per_observation = 0,
    use_univariate_updating = FALSE,
    verbose = F
) {
  message("This variant of make_staninput() is DEPRECATED and is called internally because you used the DEPRECATED function infer_prior_beliefs().
          This function will be removed in a future version of MVBeliefUpdatr. Please use fit_NIW_ideal_adaptor() instead.")
  if (!center.observations)
    message("You did not center observations. Note that the prior of category means is symmetric around 0.")

  if (pca.observations)
    assert_that(between(pca.cutoff, 0, 1), msg = "pca.cutoff must be between 0 and 1.")
  if (!is.null(lapse_rate)) {
    assert_that(is.number(lapse_rate), msg = "If not NULL, lapse_rate must be a number.")
    assert_that(between(lapse_rate, 0, 1), msg = "If not NULL, lapse rate must be a number between 0 and 1.")
  }

  cues <- unique(cues)
  exposure <-
    check_exposure_test_data(
      data = exposure,
      cues = cues,
      category = category,
      response = NULL,
      group = group,
      which.data = "exposure",
      verbose = verbose)

  if (!is.null(group.unique)) {
    assert_that(group.unique %in% names(exposure),
                msg = paste("Column for group.unique ", group.unique, "not found in exposure data."))
    message(paste0("Collapsing *exposure* observations to unique values of group.unique (", group.unique, ") by
    discarding the data from all but the first group member. This means that each unique exposure condition will
    only be counted once. All test observations are still counted, but aggregated for each unique value of group.unique."))

    exposure %<>%
      mutate(across(all_of(group.unique), as.factor)) %>%
      group_by(!! sym(group.unique), !! sym(category), !!! syms(cues)) %>%
      filter(!! sym(group) == unique(!! sym(group))[1])

    group <- group.unique
  }

  # Make sure data is ungrouped so that transform_cues works correctly, and keep only the necessary columns
  exposure %<>%
    ungroup() %>%
    select(c(all_of(group), all_of(cues), all_of(category)))

  test <-
    check_exposure_test_data(
      data = test,
      cues = cues,
      category = NULL,
      response = response,
      group = group,
      which.data = "test",
      verbose = verbose) %>%
    select(c(!! group, !!! cues, !! response))

  assert_that(all(levels(exposure[[category]]) == levels(test[[response]])),
              msg = paste("category variable", category, "in exposure data and response variable", response, "in test data must be factors with the same levels in the same order. Either the levels do not match, or they are not in the same order."))
  assert_that(all(levels(exposure[[group]]) %in% levels(test[[group]])),
              msg = paste("All levels of the grouping variable", group, "found in exposure must also be present in test."))
  if (!all(levels(test[[group]]) %in% levels(exposure[[group]])))
    message(paste("Not all levels of the grouping variable", group, "that are present in test were found in exposure.
    This is expected if and only if the data contained a test prior to (or without any) exposure.
    Creating 0 exposure data for these groups."))
  exposure %<>%
    mutate(across(all_of(group), ~ factor(.x, levels = levels(test[[!! group]]))))

  if (!is.null(mu_0)) {
    if (nlevels(exposure[[category]]) == 1) {
      assert_that(is.vector(mu_0),
                  msg = "If mu_0 is not NULL and there is only one category, mu_0 must be a vector.")
    } else {
      assert_that(is.list(mu_0) & length(mu_0) == nlevels(exposure[[category]]),
                  msg = "If mu_0 is not NULL, mu_0 must be a list of vectors with as many elements as there are categories.")
    }
    assert_that(all(map_lgl(mu_0, is.numeric), map_lgl(mu_0, ~ is.null(dim(.x)) | length(dim(.x)) == 1)),
                msg = "If mu_0 is a list, each element must be a vector.")
    assert_that(all((map_int(mu_0, length)) == length(cues)),
                msg = paste(
                  "At least one element of mu_0 does not have the correct dimensionality. Observations have",
                  length(cues),
                  "dimensions. Dimensionality of mu_0 ranges from",
                  paste(map_int(mu_0, length) %>% range(), collapse = " to ")))
  }
  if (!is.null(Sigma_0)) {
    if (nlevels(exposure[[category]]) == 1) {
      assert_that(is.array(Sigma_0),
                  msg = "If Sigma_0 is not NULL and there is only one category, Sigma_0 must be a positive-definite  matrix.")
    } else {
      assert_that(is.list(Sigma_0) & length(Sigma_0) == nlevels(exposure[[category]]),
                  msg = "If Sigma_0 not NULL, Sigma_0 must be a list of positive-definite matrices with as many elements as there are categories.")
    }
    assert_that(all(map_lgl(Sigma_0, is.numeric), map_lgl(Sigma_0, ~ length(dim(.x)) == 2)),
                msg = "If Sigma_0 is a list, each element must be a k x k matrix.")
    assert_that(all(map_lgl(Sigma_0, ~ all(dim(.x) == length(cues)))),
                msg = paste(
                  "At least one element of Sigma_0 does not have the correct dimensionality. Observations have",
                  length(cues),
                  "dimensions. Sigma_0 includes matrices of dimension",
                  paste(paste(map(Sigma_0, ~ dim(.x) %>% paste(collapse = " x "))) %>% unique(), collapse = ", ")))
  }

  # Transform data
  transform <-
    transform_cues(
      data = exposure,
      cues = cues,
      center = center.observations,
      scale = scale.observations,
      pca = pca.observations,
      return.transformed.data = T,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T)

  if (pca.observations) {
    assert_that(all(is.null(mu_0), is.null(Sigma_0)),
                msg = "PCA is not yet implemented when mu_0 or Sigma_0 are specified.")
    s = summary(transform[["transform.parameters"]][["pca"]])$importance
    l = min(which(s["Cumulative Proportion",] >= pca.cutoff))
    assert_that(l >= 1, msg = "Specified pca.cutoff does not yield to any PCA component being included. Increase the
                pca.cutoff value.")
    if (length(cues) > l)
      message(paste("Given the specified pca.cutoff, only the first", l, "principal component(s) will be used as cues."))
    cues = colnames(s)[1:l]
  }

  # If data is not scaled, set tau_scale based on SD of cues in data
  if (is.null(tau_scale))
    tau_scale <- if (scale.observations) rep(5, length(cues)) else sapply(exposure[cues], sd) * 5

  # Stan doesn't recognize vectors of length 1
  if (length(tau_scale) == 1) tau_scale <- array(tau_scale, dim = 1)

  exposure <- transform[["data"]]
  test <- transform[["transform.function"]](test)
  if (!is.null(mu_0)) mu_0 <- map(mu_0, ~ transform_category_mean(m = .x, transform))
  if (!is.null(Sigma_0) & scale.observations) Sigma_0 <- map(Sigma_0, ~ transform_category_cov(S = .x, transform))

  test_counts <-
    get_test_counts(
      test = test,
      cues = cues,
      response = response,
      group = group,
      verbose = verbose)

  if (is.null(lapse_rate)) {
    lapse_rate <- numeric()
    lapse_rate_known <- 0
  } else {
    lapse_rate <- array(lapse_rate, dim = c(1))
    lapse_rate_known <- 1
  }

  n.cats <- nlevels(exposure[[category]])
  n.cues <- length(cues)
  if (is.null(mu_0)) {
    mu_0 <- array(numeric(), dim = c(0,0))
    mu_0_known <- 0
  } else {
    mu_0_known <- 1
    if (is.list(mu_0)) {
      temp <- array(dim = c(n.cats, n.cues))
      for (i in 1:length(mu_0)) temp[i,] <- mu_0[[i]]
      mu_0 <- temp
      rm(temp)
    }
  }

  if (is.null(Sigma_0)) {
    Sigma_0 <- array(numeric(), dim = c(0,0,0))
    Sigma_0_known <- 0
  } else {
    Sigma_0_known <- 1
    if (is.list(Sigma_0)) {
      temp <- array(dim = c(n.cats, n.cues, n.cues))
      for (i in 1:length(Sigma_0)) temp[i,,] <- Sigma_0[[i]]
      Sigma_0 <- temp
      rm(temp)
    }
  }

  if (length(cues) > 1 || !use_univariate_updating) {
    staninput <-
      exposure %>%
      get_sufficient_statistics_as_list_of_arrays(
        cues = cues,
        category = category,
        group = group,
        use_univariate_updating = F,
        verbose = verbose,
        # The part below currently is ignored by get_sufficient_statistics_as_list_of_arrays. If the same syntax as for univariate input could
        # also work for multivariate input to get_sufficient_statistics_as_list_of_arrays that would be more
        # elegant.
        x_mean = colMeans, N = length, x_ss = get_sum_of_uncentered_squares_from_df) %>%
      within({
        M <- dim(x_mean)[1]
        L <- dim(x_mean)[2]
        K <- length(cues)

        x_test <-
          test_counts %>%
          select(all_of(cues)) %>%
          as.matrix()
        y_test <-
          test_counts[[group]] %>%
          as.numeric() %T>%
          { attr(., which = "levels") <- levels(test[[group]]) }
        z_test_counts <-
          test_counts %>%
          mutate(rownames = paste0("group=", !! sym(group), "; ", paste(cues, collapse = "-"), "=", paste(!!! syms(cues), sep = ","))) %>%
          column_to_rownames("rownames") %>%
          select(levels(test[[response]])) %>%
          as.matrix()
        N_test <- nrow(x_test)

        lapse_rate_known <- lapse_rate_known
        lapse_rate_data <- lapse_rate
        mu_0_known <- mu_0_known
        mu_0_data <- mu_0
        Sigma_0_known <- Sigma_0_known
        Sigma_0_data <- Sigma_0

        tau_scale <- tau_scale
        L_omega_eta <- L_omega_eta

        split_loglik_per_observation <- split_loglik_per_observation
      })

    # Clean-up x_mean and x_ss for groups without exposure data. For reasons laid out in
    # get_sufficient_statistics_as_list_of_arrays, we had to set these means and sums of
    # squares to arbitrary values (since Stan doesn't accept typed NAs). But this can
    # create confusion when users try to retrieve the exposure statistics for those groups.
    # Here we're thus setting them to NAs.
    staninput$x_mean[staninput$N == 0] <- NA
    staninput$x_ss[staninput$N == 0] <- NA
  } else if (use_univariate_updating) {
    if (length(cues) > 1) stop2("Univariate updating is only implemented for univariate data.")

    staninput <-
      exposure %>%
      get_sufficient_statistics_as_list_of_arrays(
        cues = cues, category = category, group = group,
        use_univariate_updating = T,
        verbose = verbose,
        xbar = mean, n = length, xsd = sd) %>%
      within({
        m <- dim(xbar)[1]
        l <- dim(xbar)[2]

        x_test <-
          test_counts %>%
          select(all_of(cues)) %>%
          as.matrix()
        y_test <-
          test_counts[[group]] %>%
          as.numeric() %T>%
          { attr(., which = "levels") <- levels(test[[group]]) }
        z_test_counts <-
          test_counts %>%
          select(.dots = levels(test[[response]])) %>%
          as.matrix()
        n_test <- length(x_test)

        mu_0_sd <- 0
        sigma_0_sd <- 0
      })
  }

  data <-
    bind_rows(
      exposure[, c(group.unique, if (group == group.unique) NULL else group, category, cues)] %>%
        drop_na() %>%
        mutate(Phase = "exposure"),
      test[, c(group.unique, if (group == group.unique) NULL else group, response, cues)] %>%
        drop_na() %>%
        mutate(Phase = "test")) %>%
    relocate(Phase, all_of(c(group.unique, group, category, cues, response)))

  # Remove data from transform (for storage efficiency)
  transform$data <- NULL

  return(list(staninput = staninput, data = data, transform_information = transform))
}
