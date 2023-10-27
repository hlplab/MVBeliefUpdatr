#' @importFrom reshape2 acast

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
#' @examples
#' TBD
#' @export
get_sufficient_statistics_as_list_of_arrays <- function(
  data,
  cues,
  category,
  group,
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

  if (length(cues) > 1) {
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
        temp.data_ss <- data_ss %>%
          ungroup() %>%
          filter(!! rlang::sym(category) == cats[i] &
                   !! rlang::sym(group) == groups[j])

        if (nrow(temp.data_ss) > 0) {
          N[i,j] = temp.data_ss$N[[1]]
          x_mean[i,j,] = temp.data_ss$x_mean[[1]]
          x_ss[i,j,,] = temp.data_ss$x_ss[[1]]
        } else {
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


#' Compose data for input to RStan
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
#' and contain information about the category label, the cue values of the observation, and optionally grouping variables.
#' @param test `tibble` or `data.frame` with the test data. Each row should be an observation, and contain information
#' about the cue values of the test stimulus and the participant's response.
#' @param cues Names of columns with cue values. Must exist in both exposure and test data.
#' @param category Name of column in exposure data that contains the category label. Can be `NULL` for unsupervised updating
#' (not yet implemented). (default: "category")
#' @param response Name of column in test data that contains participants' responses. (default: "response")
#' @param group Name of column that contains information about which observations form a group. Typically, this is
#' a variable identifying subjects/participants. Must exist in both exposure and test data. (default: "group")
#' @param group.unique Name of column that uniquely identifies each group with identical exposure. This could be a
#' variable indicating the different conditions in an experiment. Using group.unique is optional, but can be
#' substantially more efficient if many groups share the same exposure. To ignore, set to `NULL`. (default: `NULL`)
#' @param center.observations Should the data be centered? Centering will not affect the inferred correlation or
#' covariance matrices but it will affect the absolute position of the inferred means. The relative position of
#' the inferred means remains unaffected. If `TRUE` and `mu_0` is specified, `mu_0` will also be centered (`Sigma_0` is not
#' affected by centering and thus not changed). (default: `TRUE`)
#' @param scale.observations Should the data be standardized? Scaling will not affect the inferred correlation matrix,
#' but it will affect the inferred covariance matrix because it affects the inferred standard deviations. It will also
#' affect the absolute position of the inferred means. The relative position of the inferred means remains unaffected.
#' If `TRUE` and `mu_0` and `Sigma_0` are specified, `mu_0` and `Sigma_0` will also be scaled.
#' (default: `TRUE`)
#' @param pca.observations Should the data be transformed into orthogonal principal components? (default: `FALSE`)
#' @param pca.cutoff Determines which principal components are handed to the MVBeliefUpdatr Stan program: all
#' components necessary to explain at least the pca.cutoff of the total variance. (default: .95) Ignored if
#' `pca.observation = FALSE`. (default: 1)
#' @param lapse_rate,mu_0,Sigma_0 Optionally, lapse rate, prior expected category means (mu_0) and/or prior expected
#' category covariance matrices (Sigma_0) for all categories. Lapse rate should be a number between 0 and 1. For mu_0
#' and Sigma_0, each should be a list, with each element being the expected mean/covariance matrix for a specific
#' category prior to updating. Elements of mu_0 and Sigma_0 should be ordered in the same order as the levels of the
#' category variable in \code{exposure} and \code{test}. These prior expected means and covariance matrices could be
#' estimated, for example, from phonetically annotated speech recordings (see \code{\link{make_MVG_from_data}}
#' for a convenient way to do so). Internally, m_0 is then set to mu_0 (so that the expected value of the prior
#' distribution of means is mu_0) and S_0 is set so that the expected value of the inverse-Wishart is Sigma_0 given nu_0.
#' @param tau_0_scales Optionally, a vector of scales for the Cauchy priors for each cue's standard deviations. Used in
#' both the prior for m_0 and the prior for S_0. (default: vector of 5s of length of cues, assumes scaled input)
#' @param omega_0_eta Optionally, etas the LKJ prior for the correlations of the covariance matrix of mu_0. Set to 0 to
#' ignore. (default: 0)
#'
#' @return A list that is an \code{NIW_ideal_adaptor_input}. In interpreting the inferred kappa_0 and nu_0, it should
#' be kept in mind that the \emph{inferred} scatter matrix S_0 includes variability from internal perceptual and/or
#' external environmental noise, \emph{in addition} to the motor noise that is reflected in production data. This also
#' implies that, if Sigma_0 is given, Sigma_0 and nu_0 mutually constrain each other, because the expected value of
#' Sigma_0 is determined by both S_0 and nu.
#'
#' @seealso \code{\link{is.NIW_ideal_adaptor_input}}
#' @keywords TBD
#' @examples
#' TBD
#' @rdname compose_data
#' @export
compose_data_to_infer_prior_via_conjugate_ibbu_w_sufficient_stats = function(
  exposure, test,
  cues, category = "category", response = "response", group = "group", group.unique = NULL,
  center.observations = T, scale.observations = T, pca.observations = F, pca.cutoff = 1,
  lapse_rate = NULL, mu_0 = NULL, Sigma_0 = NULL,
  tau_scale = 0, # rep(5, length(cues)),
  L_omega_scale = 0,
  verbose = F
) {
  if ((!center.observations | !scale.observations) & (tau_scale == 0 | L_omega_scale == 0))
    message("The tau_scale and L_omega_scale parameters are not specified (using defaults). Since you also did not ask for the input to be centered *and* scaled, this puts the priors assumed in the model on a scale that has no relation to the input. Unless you have manually centered and scaled the cues, this is strongly discouraged.")

  if (pca.observations)
    assert_that(between(pca.cutoff, 0, 1),
                msg = "pca.cutoff must be between 0 and 1.")

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

    group = group.unique
  }

  exposure %<>%
    # Make sure data is ungrouped so that transform_cues works correctly.
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
              msg = paste("category variable", category, "in exposure data and response variable", response, " in test data must be factors with the same levels in the same order. Either the levels do not match, or they are not in the same order."))
  assert_that(all(levels(exposure[[group]]) %in% levels(test[[group]])),
              msg = paste("All levels of the grouping variable", group, "found in exposure must also be present in test."))
  if (!all(levels(test[[group]]) %in% levels(exposure[[group]])))
    message(paste("Not all levels of the grouping variable", group, "that are present in test were found in exposure.
    Creating 0 exposure data for those groups."))
  exposure %<>%
    mutate(across(all_of(group), ~ factor(.x, levels = levels(test[[!! group]]))))

  if (!is.null(lapse_rate)) {
    assert_that(is.number(lapse_rate), msg = "If not NULL, lapse_rate must be a number.")
    assert_that(between(lapse_rate, 0, 1), msg = "If not NULL, lapse rate must be a number between 0 and 1.")
  }
  if (!is.null(mu_0)) {
    if (nlevels(exposure[[category]]) == 1) {
      assert_that(is.vector(mu_0),
                  msg = "If mu_0 is not NULL and there is only one category, mu_0 must be a vector.")
    } else {
      assert_that(is.list(mu_0) & length(mu_0) == nlevels(exposure[[category]]),
                msg = "If mu_0 is not NULL, mu_0 must be a list of vectors with as many elements as there are categories.")
    }
    assert_that(all(map(mu_0, is.numeric)  %>% unlist, map(mu_0, ~ is.null(dim(.x)) | length(dim(.x)) == 1) %>% unlist),
                msg = "If mu_0 is a list, each element must be a vector.")
    assert_that(all((map(mu_0, length) %>% unlist()) == length(cues)),
                msg = paste(
                  "At least one element of mu_0 does not have the correct dimensionality. Observations have",
                  length(cues),
                  "dimensions. Dimensionality of mu_0 ranges from",
                  paste(map(mu_0, length) %>% unlist() %>% range(), collapse = " to ")))
  }
  if (!is.null(Sigma_0)) {
    if (nlevels(exposure[[category]]) == 1) {
      assert_that(is.array(Sigma_0),
                msg = "If Sigma_0 is not NULL and there is only one category, Sigma_0 must be a positive-definite  matrix.")
    } else {
      assert_that(is.list(Sigma_0) & length(Sigma_0) == nlevels(exposure[[category]]),
                 msg = "If Sigma_0 not NULL, Sigma_0 must be a list of positive-definite matrices with as many elements as there are categories.")
    }
    assert_that(all(map(Sigma_0, is.numeric)  %>% unlist, map(Sigma_0, ~ length(dim(.x)) == 2) %>% unlist),
                msg = "If Sigma_0 is a list, each element must be a k x k matrix.")
    assert_that(all(map(Sigma_0, ~ dim(.x) == length(cues)) %>% unlist()),
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
      return.transform.function = T)

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

  if (length(cues) > 1) {
    data_list <-
      exposure %>%
      get_sufficient_statistics_as_list_of_arrays(
        cues = cues,
        category = category,
        group = group,
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
        L_omega_scale <- L_omega_scale
      })
  } else {
    message("For univariate input, beliefupdatr is run with legacy parameter names. This might change in the future.")
    data_list <- exposure %>%
      get_sufficient_statistics_as_list_of_arrays(
        cues = cues, category = category, group = group,
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

  if (verbose) {
    print("In compose_data_to_infer_prior_via_conjugate_ibbu_w_sufficient_stats():")
    print(data_list)
  }

  return(data_list)
}


attach_stanfit_input_data = function(stanfit, input) {
  assert_NIW_ideal_adaptor_stanfit(stanfit)
  assert_that(is.NIW_ideal_adaptor_input(input),
              msg = "input is not an acceptable input data.")
  slot(stanfit, "input_data", check = T) <- input

  return(stanfit)
}

attach_stanfit_transform = function(stanfit, transform_information) {
  assert_NIW_ideal_adaptor_stanfit(stanfit)
  slot(stanfit, "transform_information", check = T) <- transform_information

  return(stanfit)
}


