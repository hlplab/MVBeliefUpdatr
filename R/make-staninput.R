#' @export
check_exposure_test_data <- function(data, cues, category, response, group, which.data = "the", verbose = F) {
  assert_that(is_tibble(data) | is.data.frame(data))
  assert_cols_in_data(data, cues, which.data, scalar = F)
  assert_cols_in_data(data, group, which.data, scalar = T)

  data %<>%
    drop_na(.env$cues, .env$group) %>%
    mutate(across(.env$group, as.factor))

  if(!is.null(category)) {
    assert_cols_in_data(data, category, which.data, scalar = T)
    data %<>%
      mutate(across(.env$category, as.factor)) %>%
      drop_na(.env$category)
  }

  if (!is.null(response)) {
    assert_cols_in_data(data, response, which.data, scalar = T)
    data %<>%
      mutate(across(.env$response, as.factor)) %>%
      drop_na(.env$response)
  }
  assert_that(nrow(data) > 0,
              msg = paste("There must be at least one observation in", which.data, "data."))

  if (verbose){
    print("In check_exposure_test_data():")
    print (data)
  }

  return(data)
}


#' @export
get_test_counts <- function(test, cues, response, group, verbose = F) {
  test_counts <-
    test %>%
    as_tibble(.name_repair = "minimal") %>%
    group_by(
      !!! syms(cues),
      !! sym(response),
      !! sym(group)) %>%
    tally() %>%
    pivot_wider(
      names_from = !! response,
      values_from = n,
      values_fill = 0) %>%
    ungroup()

  if (verbose) {
    print("In get_test_counts():")
    print(test_counts)
  }

  return(test_counts)
}



#' Get sufficient statistics from a data set
#'
#' Get sufficient statistics from data. Calculates the means and sum of squares matrices for the specified cues for
#' any combination of groups (optional) and categories, and returns them as a list.
#'
#' @param data `tibble` or `data.frame` with the data. Each row should be an observation of a category,
#' and contain information about the category label, the cue values of the observation, and optionally grouping variables.
#' @param test `tibble` or `data.frame` with the test data. Each row should be an observation, and contain information
#' about the cue values of the test stimulus and the participant's response.
#' @param cues Names of columns with cue values.
#' @param category Name of column that contains the category label for the exposure data. Can be `NULL` for unsupervised updating
#' (not yet implemented). (default: "category")
#' @param group Name of column that contains information about which observations form a group. This could be individual
#' subjects or conditions in an experiment. The latter is more efficient, but should only be used if exposure is
#' identical for every individual within the group. Test does not have to be identical for every individual within
#' the same group. For example, one can group multiple groups of subjects that have received the same exposure
#' but were tested on different test tokens.
#'
#' @return A list.
#'
#' @keywords TBD
#' @export
get_sufficient_statistics_as_list_of_arrays <- function(
  data,
  cues,
  category,
  group,
  use_univariate_updating = F,
  verbose = F,
  ...
) {

  if (verbose) {
    print("In get_sufficient_statistics_as_list_of_arrays(), input data is:")
    print(data)
  }

  data_ss <-
    data %>%
    as_tibble(.name_repair = "minimal") %>%
    group_by(!! sym(category), !! sym(group))

  if (!use_univariate_updating) {
    # Multivariate observations
    data_ss %<>%
      summarise(
        N = length(!! sym(cues[1])),
        x_mean = list(colMeans(cbind(!!! syms(cues)))),
        x_ss = list(get_sum_of_uncentered_squares_from_df(cbind(!!! syms(cues)), verbose = verbose)))

    if (verbose) {
      print("In get_sufficient_statistics_as_list_of_arrays(), multivariate sum-of-uncentered-squares matrix:")
      print(data_ss)
    }

    ## -------------------------------------------------------------------------------
    # This is intended to map the elements of the tibble into arrays of the required
    # dimensionality. THERE PROBABLY IS A MUCH MORE CONCISE AND PERHAPS MORE EFFICIENT
    # WAY TO DO THIS. CHECK BACK.
    #
    # For helpful concise info on tibbles, see
    #   https://cran.r-project.org/web/packages/tibble/vignettes/tibble.html
    ## -------------------------------------------------------------------------------
    cats <- levels(data[[category]])
    groups <- levels(data[[group]])
    n_category <- length(cats)
    n_group <- length(groups)
    n_cues <- length(cues)

    N = array(dim = c(n_category,n_group))
    x_mean = array(dim = c(n_category,n_group,n_cues))
    x_ss = array(dim = c(n_category,n_group,n_cues,n_cues))

    for (i in 1:n_category) {
      for (j in 1:n_group) {
        temp.data_ss <-
          data_ss %>%
          ungroup() %>%
          filter(!! rlang::sym(category) == cats[i] &
                   !! rlang::sym(group) == groups[j])

        # Catch cases in which there is no data for a particular group/category combination
        # (this can happen for example, when the data contain a pre-exposure test, which has
        # no matching exposure statistics).
        if (nrow(temp.data_ss) > 0) {
          N[i,j] = temp.data_ss$N[[1]]
          x_mean[i,j,] = temp.data_ss$x_mean[[1]]
          x_ss[i,j,,] = temp.data_ss$x_ss[[1]]
        } else {
          # For groups without exposure data, we are setting the category means to 0 and the
          # sum-of-squares matrix to the identity matrix. This is a bit of a hack, that is
          # necessary because stan expects these variables to always be vectors/matrices of
          # the same type and dimensionality (even though they should really be NAs). So, in
          # order to avoid confusion, we're setting these quantities to NA *after* all required
          # input has been handed to stan.
          N[i,j] = 0
          x_mean[i,j,] = rep(0, length(cues))
          x_ss[i,j,,] = diag(length(cues))
        }
      }
    }

    dimnames(N) <- list(cats, groups)
    dimnames(x_mean) <- list(cats, groups, cues)
    dimnames(x_ss) <- list(cats, groups, cues, cues)
    data_ss <- list(N = N, x_mean = x_mean, x_ss = x_ss)
  } else {
    # Univariate observations
    data_ss %<>%
      summarise_at(cues, .funs = list(...))

    if (verbose) {
      print("In get_sufficient_statistics_as_list_of_arrays(), univariate sum-of-squares matrix (prior to map application):")
      print(data_ss)
    }

    stats <- names(list(...))
    data_ss <- map(stats, ~ acast(data_ss, as.list(c(category, group)),
                                  value.var = .x,
                                  drop = F,
                                  fill = 0)) %>%
      set_names(stats)
  }

  if (verbose) {
    print("In get_sufficient_statistics_as_list_of_arrays(), sum-of-squares matrix (uncentered for multivariate data, centered for univariate data):")
    print(data_ss)
  }

  return(data_ss)
}


#' Prepare data to fit NIW_ideal_adaptor_stanfit via rstan
#'
#' Take exposure and test data as input, and prepare the data for input into an MVBeliefUpdatr Stan program.
#'
#' It is important to use \code{group} to identify individuals that had a specific exposure (or no exposure at all)
#' and specific test trials. You should \emph{not} use \code{group} to identify exposure conditions. Setting \code{group} to an exposure condition
#' results in an exposure that concatenates the exposure observations from all subjects in that condition. Typically, this
#' is not what users intend, as it models exposure to the combination of exposure tokens across all subjects, rather than
#' exposure to one set of those exposure tokens. To achieve this intended outcome, use
#' \code{group.unique} to identify groups with identical exposure. This will correctly use only one unique instance of the
#' observations that any level of \code{group} receives during exposure.
#'
#' @param exposure `tibble` or `data.frame` with the exposure data. Each row should be an observation of a category,
#'   and contain information about the category label, the cue values of the observation, and optionally grouping variables.
#' @param test `tibble` or `data.frame` with the test data. Each row should be an observation, and contain information
#'   about the cue values of the test stimulus and the participant's response.
#' @param cues Names of columns with cue values. Must exist in both exposure and test data.
#' @param category Name of column in exposure data that contains the category label. Can be \code{NULL} for unsupervised updating
#'   (not yet implemented). (default: "category")
#' @param response Name of column in test data that contains participants' responses. (default: "response")
#' @param group Name of column that contains information about which observations form a group. Typically, this is
#'   a variable identifying subjects/participants. Must exist in both exposure and test data. (default: "group")
#' @param group.unique Name of column that uniquely identifies each group with identical exposure. This could be a
#'   variable indicating the different conditions in an experiment. Using group.unique is optional, but can be
#'   substantially more efficient if many groups share the same exposure. To ignore, set to \code{NULL}. (default: \code{NULL})
#' @param transform_type An affine transformation that can be applied to the data. See `type` in \code{\link{get_affine_transform}}
#'    for details. (default: "PCA whiten")
#' @param lapse_rate,mu_0,Sigma_0 Optionally, lapse rate, prior expected category means (\code{mu_0}) and/or prior expected
#'   category covariance matrices (\code{Sigma_0}) for all categories. Lapse rate should be a number between 0 and 1. For \code{mu_0}
#'   and \code{Sigma_0}, each should be a list, with each element being the expected mean/covariance matrix for a specific
#'   category prior to updating. Elements of \code{mu_0} and \code{Sigma_0} should be ordered in the same order as the levels of the
#'   category variable in \code{exposure} and \code{test}. These prior expected means and covariance matrices could be
#'   estimated, for example, from phonetically annotated speech recordings (see \code{\link{make_MVG_from_data}}
#'   for a convenient way to do so). Internally, m_0 is then set to \code{mu_0} (so that the expected value of the prior
#'   distribution of means is mu_0) and S_0 is set so that the expected value of the inverse-Wishart is \code{Sigma_0} given nu_0.
#'   Importantly, \strong{Sigma_0 should be convolved with perceptual noise (i.e., add perceptual noise covariance matrix to
#'   the category variability covariance matrices when you specify \code{Sigma_0})} since the stancode for the inference of the
#'   NIW ideal adaptor does \emph{not} infer category and noise variability separately.
#' @param tau_scale A vector of scales for the Cauchy priors for each cue's standard deviations. Used in
#'   both the prior for m_0 and the prior for S_0. (default: vector of 5s if scale.observations = TRUE, SD of each cue otherwise).
#' @param L_omega_eta A vector of etas of the LKJ prior for the correlations of the covariance matrix of \code{mu_0}. (default: 1,
#'   which corresponds to a uniform prior of correlation matrices)
#' @param split_loglik_per_observation Optionally, split the log likelihood per observation. This can be helpful of leave-one-out
#'   estimation in order to avoid high Pareto k, but it also makes the stored stanfit object much larger. (default: 0)
#'
#' @return A list consisting of:
#'   * `data`: A data.frame with the exposure and test data after exclusion of NAs and other checks.
#'   * `staninput`: A named list of variables and values to be handed to Stan.
#'   * `transform_information`: A list with information about the transformation that was applied to the data.
#'
#' @seealso \code{\link{is.NIW_ideal_adaptor_staninput}}
#' @keywords TBD
#'
#' @importFrom purrr map_lgl map_int
#' @importFrom reshape2 acast
#' @rdname make_staninput
#' @export
make_staninput <- function(..., model_type = "NIW_ideal_adaptor") {
  if (model_type == "NIW_ideal_adaptor") {
    make_staninput_for_NIW_ideal_adaptor(...)
  } else {
    message("Model type ", model_type, " not recognized. Supported model types are: NIW_ideal_adaptor.")
  }
}


#' @rdname make_staninput
#' @export
make_staninput_for_NIW_ideal_adaptor <- function(
    exposure, test,
    cues, category = "category", response = "response",
    group = "group", group.unique = NULL,
    lapse_rate = NULL, mu_0 = NULL, Sigma_0 = NULL,
    tau_scale = rep(5, length(cues)), L_omega_eta = 1,
    # THIS IS KEPT HERE JUST FOR NOW UNTIL I HAVE DETERMINED WHICH TRANSFORM IS BEST SUITED FOR FITTING.
    # (also remove the documentation for this once it's no longer needed AND remove transform_information
    # from the returned information AND change the documentation for the returned object above.)
    transform_type = c("identity", "center", "standardize", "PCA whiten", "ZCA whiten")[4],
    use_univariate_updating = F,
    split_loglik_per_observation = 0,
    verbose = F
) {
  if (!is.null(lapse_rate)) {
    assert_that(is.number(lapse_rate), msg = "If not NULL, lapse_rate must be a number.")
    assert_that(between(lapse_rate, 0, 1), msg = "If not NULL, lapse rate must be a number between 0 and 1.")
  }

  assert_that(
    length(tau_scale) == length(cues),
    msg = paste0("tau_scale must be a vector with the same number of elements as cues (", length(cues), ")."))
  assert_that(
    transform_type %in% c("identity", "center", "standardize", "PCA whiten", "ZCA whiten"),
    msg = paste("transform_type must be one of the following: identity, center, standardize, PCA whiten, ZCA whiten."))

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

  # Transform data and category representations (if provided)
  transform <- get_affine_transform(exposure, cues, transform_type)
  exposure_transformed <- transform[["transform.function"]](exposure)
  test_transformed <- transform[["transform.function"]](test)
  if (!is.null(mu_0)) mu_0_transformed <- map(mu_0, ~ transform_category_mean(m = .x, transform)) else mu_0_transformed <- mu_0
  if (!is.null(Sigma_0)) Sigma_0_transformed <- map(Sigma_0, ~ transform_category_cov(S = .x, transform)) else Sigma_0_transformed <- Sigma_0

  test_counts <-
    get_test_counts(
      test = test,
      cues = cues,
      response = response,
      group = group,
      verbose = verbose)

  test_counts_transformed <-
    get_test_counts(
      test = test_transformed,
      cues = cues,
      response = response,
      group = group,
      verbose = verbose)

  # Make sure that scales are converted into arrays for Stan
  # and check that any null values are set correctly.
  if (is.null(lapse_rate)) {
    lapse_rate <- numeric()
    lapse_rate_known <- 0
  } else {
    lapse_rate <- array(lapse_rate, dim = c(1))
    lapse_rate_known <- 1
  }

  n.cats <- nlevels(exposure[[category]])
  n.cues <- length(cues)
  if (is.null(mu_0_transformed)) {
    mu_0 <- array(numeric(), dim = c(0,0))
    mu_0_transformed <- array(numeric(), dim = c(0,0))
    mu_0_known <- 0
  } else {
    mu_0_known <- 1
    if (is.list(mu_0_transformed)) {
      temp <- array(dim = c(n.cats, n.cues))
      for (i in 1:length(mu_0)) temp[i,] <- mu_0[[i]]
      mu_0 <- temp

      temp <- array(dim = c(n.cats, n.cues))
      for (i in 1:length(mu_0_transformed)) temp[i,] <- mu_0[[i]]
      mu_0_transformed <- temp
      rm(temp)
    }
  }

  if (is.null(Sigma_0_transformed)) {
    Sigma_0 <- array(numeric(), dim = c(0,0,0))
    Sigma_0_transformed <- array(numeric(), dim = c(0,0,0))
    Sigma_0_known <- 0
  } else {
    Sigma_0_known <- 1
    if (is.list(Sigma_0_transformed)) {
      temp <- array(dim = c(n.cats, n.cues, n.cues))
      for (i in 1:length(Sigma_0)) temp[i,,] <- Sigma_0[[i]]
      Sigma_0 <- temp

      temp <- array(dim = c(n.cats, n.cues, n.cues))
      for (i in 1:length(Sigma_0_transformed)) temp[i,,] <- Sigma_0[[i]]
      Sigma_0_transformed <- temp
      rm(temp)
    }
  }

  if (length(tau_scale) == 1) tau_scale <- array(tau_scale, dim = 1)
  shift <- transform$transform.parameters[["shift"]]
  if (length(shift) == 1) shift <- array(shift, dim = 1)
  INV_SCALE <- transform$transform.parameters[["INV_SCALE"]]
  if (!is.matrix(INV_SCALE)) INV_SCALE <- as.matrix(INV_SCALE, ncol = length(cues))

  if (length(cues) > 1 | !use_univariate_updating) {
    # First get untransform input (to be stored in fit since it's helpful for plotting), and then get
    # transformed input below.
    staninput <-
      exposure %>%
      get_sufficient_statistics_as_list_of_arrays(
        cues = cues, category = category, group = group,
        use_univariate_updating = F,
        verbose = verbose,
        # The part below currently is ignored by get_sufficient_statistics_as_list_of_arrays. If the same syntax as for univariate input could
        # also work for multivariate input to get_sufficient_statistics_as_list_of_arrays that would be more
        # elegant.
        x_mean = colMeans, N = length, x_ss = get_sum_of_uncentered_squares_from_df) %>%
      within({
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

        mu_0_known <- mu_0_known
        mu_0_data <- mu_0
        Sigma_0_data <- Sigma_0
      })

    # Clean-up x_mean and x_ss for groups without exposure data. For reasons laid out in
    # get_sufficient_statistics_as_list_of_arrays, we had to set these means and sums of
    # squares to arbitrary values (since Stan doesn't accept typed NAs). But this can
    # create confusion when users try to retrieve the exposure statistics for those groups.
    # Here we're thus setting them to NAs.
    staninput$x_mean[staninput$N == 0] <- NA
    staninput$x_ss[staninput$N == 0] <- NA

    staninput_transformed <-
      exposure_transformed %>%
      get_sufficient_statistics_as_list_of_arrays(
        cues = cues, category = category, group = group,
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
          test_counts_transformed %>%
          select(all_of(cues)) %>%
          as.matrix()
        y_test <-
          test_counts_transformed[[group]] %>%
          as.numeric() %T>%
          { attr(., which = "levels") <- levels(test_transformed[[group]]) }
        z_test_counts <-
          test_counts_transformed %>%
          mutate(rownames = paste0("group=", !! sym(group), "; ", paste(cues, collapse = "-"), "=", paste(!!! syms(cues), sep = ","))) %>%
          column_to_rownames("rownames") %>%
          select(levels(test_transformed[[response]])) %>%
          as.matrix()
        N_test <- nrow(x_test)

        lapse_rate_known <- lapse_rate_known
        lapse_rate_data <- lapse_rate
        mu_0_known <- mu_0_known
        mu_0_data <- mu_0_transformed
        Sigma_0_known <- Sigma_0_known
        Sigma_0_data <- Sigma_0_transformed

        tau_scale <- tau_scale
        L_omega_eta <- L_omega_eta

        shift <- shift
        INV_SCALE <- INV_SCALE

        split_loglik_per_observation <- split_loglik_per_observation
      })

    # Clean-up x_mean and x_ss for groups without exposure data. For reasons laid out in
    # get_sufficient_statistics_as_list_of_arrays, we had to set these means and sums of
    # squares to arbitrary values (since Stan doesn't accept typed NAs). But this can
    # create confusion when users try to retrieve the exposure statistics for those groups.
    # Here we're thus setting them to NAs.
    staninput_transformed$x_mean[staninput_transformed$N == 0] <- NA
    staninput_transformed$x_ss[staninput_transformed$N == 0] <- NA
  } else if (use_univariate_updating) {
    if (length(cues) > 1) stop2("Univariate updating is only implemented for univariate data.")
    if (transform_type != "identity") stop2('Univariate updating is only implemented for transform_type = "identity".')

    staninput_transformed <- NULL
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
  } else {
    stop2("use_univariate_updating not specified.")
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

  return(
    list(
      staninput = list(
        transformed = staninput_transformed,
        untransformed = staninput),
      data = data,
      transform_information = transform))
}

#' Compose data to fit NIW_ideal_adaptor_stanfit via rstan
#'
#' DEPRECATED: Use \code{make_staninput(model_type = "NIW_ideal_adaptor")} instead.
#'
#' @export
compose_data_to_infer_NIW_ideal_adaptor <- function(model_type = "NIW_ideal_adaptor", ...)
  make_staninput(..., model_type = "NIW_ideal_adaptor")


