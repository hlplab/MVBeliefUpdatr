check_exposure_test_data <- function(data, cues, category, response, group, which.data = "the", verbose = F) {
  assert_that(is_tibble(data) | is.data.frame(data))
  assert_that(all(is_character(cues)),
              msg = "cues must be a column name or vector of column names.")
  assert_that(all(cues %in% names(data)),
              msg = paste("Cue column(s)", cues[which(cues %nin% names(data))], "not found in", which.data, "data." ))

  if(!is.null(category)) {
    assert_that(is_scalar_character(category),
                msg = "category must be a column name.")
    assert_that(category %in% names(data),
                msg = paste("Category column", category, "not found in", which.data, "data."))

    data %<>%
      mutate_at(category, as.factor)
  }

  if (!is.null(response)) {
    assert_that(is_scalar_character(response),
                msg = "response must be a column name.")
    assert_that(response %in% names(data),
                msg = paste("Response column", response, "not found in", which.data, "data."))

    data %<>%
      mutate_at(response, as.factor)
  }

  if (!is.null(group)) {
    assert_that(group %in% names(data),
                msg = paste("Group column", group, "not found in", which.data,"data."))

    data %<>%
      mutate_at(group, as.factor)
  }

  if (verbose){
    print("In check_exposure_test_data():")
    print (data)
  }


  return(data)
}



get_test_counts <- function(test, cues, category, response, group, verbose = F) {
  test_counts <- test %>%
    as_tibble() %>%
    group_by(!!! syms(cues),
             !! sym(response)) %>%
    { if (!is.null(group)) group_by(., !! sym(group), .add = T) } %>%
    tally() %>%
    pivot_wider(
      names_from = !! response,
      values_from = n,
      values_fill = 0
    ) %>%
    ungroup()

  if (verbose) {
    print("In test_counts():")
    print(test_counts)
  }

  return(test_counts)
}



#' Get sufficient statistics from a data set
#'
#' Get sufficient statistics from data. Calculates the means and covariance matrices for the specified cues for
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
#' @seealso
#' @keywords TBD
#' @examples
#' TBD
#' @export
get_sufficient_statistics_from_data <- function(data, cues, category, group, verbose = F, ...) {
  data = check_exposure_test_data(
    data = data,
    cues = cues,
    category = category,
    response = NULL,
    group = group,
    verbose = verbose)

  data_ss <- data %>%
    as_tibble() %>%
    group_by(!!! sym(category)) %>%
    { if (!is.null(group)) group_by(., !! sym(group), .add = T) }

  if (length(cues) > 1) {
    # Multivariate observations
    data_ss %<>%
      summarise(
        N = length(!! sym(cues[1])),
        x_mean = list(colMeans(cbind(!!! syms(cues)))),
        x_ss = list(get_sum_of_uncentered_squares(cbind(!!! syms(cues)), verbose = verbose))
      )

    if (verbose) {
      print("In get_sufficient_statistics_from_data(), multivariate sum-of-uncentered-squares matrix:")
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
    cats = levels(data[[category]])
    groups = levels(data[[group]])
    n_cat = length(cats)
    n_subj = length(groups)
    n_cues = length(cues)

    N = array(dim = c(n_cat,n_subj))
    x_mean = array(dim = c(n_cat,n_subj,n_cues))
    x_ss = array(dim = c(n_cat,n_subj,n_cues,n_cues))

    for (i in 1:n_cat) {
      for (j in 1:n_subj) {
        temp.data_ss = data_ss %>%
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

    dimnames(N) = list(cats, groups)
    dimnames(x_mean) = list(cats, groups, cues)
    dimnames(x_ss) = list(cats, groups, cues, cues)

    data_ss = list(N = N, x_mean = x_mean, x_ss = x_ss)
  } else {
    # Univariate observations
    data_ss %<>%
      summarise_at(cues, .funs=funs(...))

    if (verbose) {
      print("In get_sufficient_statistics_from_data(), univariate sum-of-squares matrix (prior to map application):")
      print(data_ss)
    }

    stats <- names(list(...))

    data_ss <- map(stats, ~ reshape2::acast(data_ss, as.list(c(category, group)), value.var=.x)) %>%
      set_names(stats)
  }


  if (verbose) {
    print("In get_sufficient_statistics_from_data(), sum-of-squares matrix (uncentered for multivariate data,
          centered for univariate data):")
    print(data_ss)
  }

  return(data_ss)
}




#' Compose data for input to RStan
#'
#' Take exposure and test data as input, and prepare the data for input into an MVBeliefUpdatr Stan program.
#'
#' @param exposure `tibble` or `data.frame` with the exposure data. Each row should be an observation of a category,
#' and contain information about the category label, the cue values of the observation, and optionally grouping variables.
#' @param test `tibble` or `data.frame` with the test data. Each row should be an observation, and contain information
#' about the cue values of the test stimulus and the participant's response.
#' @param cues Names of columns with cue values.
#' @param category Name of column that contains the category label for the exposure data. Can be `NULL` for unsupervised updating
#' (not yet implemented). (default: "category")
#' @param response Name of column that contains participants' responses for the test data. (default: "response")
#' @param group Name of column that contains information about which observations form a group. This could be individual
#' subjects or conditions in an experiment. The latter is more efficient, but should only be used if exposure is
#' identical for every individual within the group. Test does not have to be identical for every individual within
#' the same group. For example, one can group multiple groups of subjects that have received the same exposure
#' but were tested on different test tokens.
#' @param center.observations Should the data be centered? (default: `TRUE`)
#' @param scale.observations Should the data be standardized? (default: `TRUE`)
#' @param pca.observations Should the data be transformed into orthogonal principal components? (default: `FALSE`)
#' @param pca.cutoff Determines which principal components are handed to the MVBeliefUpdatr Stan program: all
#' components necessary to explain at least the pca.cutoff of the total variance. (default: .95) Ignored if
#' `pca.observation = FALSE`. (default: 1)
#' @param tau_scale,L_omega_scale Scale for the Cauchy prior for standard deviations of the covariance matrix of mu_0 and
#' scale for the LKJ prior for the correlations of the covariance matrix of mu_0. Set to 0 to ignore. (default: 0)
#'
#' @return A list that is an \code{mvg_ibbu_input}.
#'
#' @seealso \code{\link{is.mvg_ibbu_input}}
#' @keywords TBD
#' @examples
#' TBD
#' @rdname compose_data
#' @export
compose_data_to_infer_prior_via_conjugate_ibbu_w_sufficient_stats = function(
  exposure, test,
  cues, category = "category", response = "response", group = NULL,
  center.observations = T, scale.observations = T, pca.observations = F, pca.cutoff = 1,
  tau_scale = 0, L_omega_scale = 0,
  verbose = F
) {
  exposure <- check_exposure_test_data(
    data = exposure,
    cues = cues,
    category = category,
    response = NULL,
    group = group,
    which.data = "exposure",
    verbose = verbose)
  test <- check_exposure_test_data(
    data = test,
    cues = cues,
    category = NULL,
    response = response,
    group = group,
    which.data = "test",
    verbose = verbose)

  assert_that(all(levels(exposure[[category]]) == levels(test[[response]])),
              msg = paste("category variable", category, "in exposure and response colum", response, "must be factors
              with the same levels in the same order in. Either the levels do not match, or they are not in the same
              order."))
  if (!is.null(group))
    assert_that(all(levels(exposure[[group]]) %in% levels(test[[group]])),
                msg = paste("All levels of the group variable", group, "found in exposure must also be present in test."))
  if (!all(levels(test[[group]]) %in% levels(exposure[[group]])))
    message(paste("Not all levels of group variable", group, "that are present in test were found in exposure. Creating
                  0 exposure data for those groups."))
  exposure %<>%
    mutate_at(group, ~ factor(.x, levels = levels(test[[group]])))

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

  if (pca) {
    s = summary(transform[["transform.parameters"]][["pca"]])$importance
    cues = colnames(s)[1:min(which(s["Cumulative Proportion",] > pca.cutoff))]
  }

  exposure = transform[["data"]]
  test = transform[["transform.function"]](test)

  test_counts <- get_test_counts(
    test = test,
    cues = cues,
    category = category,
    response = response,
    group = group,
    verbose = verbose)

  print(exposure)
  print(test_counts)

  if (length(cues) > 1) {
    data_list <- exposure %>%
      get_sufficient_statistics_from_data(
        cues = cues, category = category, group = group,
        verbose = verbose,
        # The part below currently is ignored by get_sufficient_statistics_from_data. If the same syntax as for univariate input could
        # also work for multivariate input to get_sufficient_statistics_from_data that would be more
        # elegant.
        x_mean = colMeans, N = length, x_ss = get_sum_of_uncentered_squares) %>%
      within({
        M <- dim(x_mean)[1]
        L <- dim(x_mean)[2]
        K <- length(cues)

        x_test <- test_counts[[cues]]
        y_test <- as.numeric(test_counts[[group]])
        z_test_counts <-
          test_counts %>%
          select(.dots = levels(test[[response]])) %>%
          as.matrix()
        N_test <- nrow(x_test)

        tau_scale <- tau_scale
        L_omega_scale <- L_omega_scale
      })
  } else {
    message("For univariate input, beliefupdatr is run with legacy parameter names. This might change in
              the future.")
    data_list <- exposure %>%
      get_sufficient_statistics_from_data(
        cues = cues, category = category, group = group,
        verbose = verbose,
        xbar = mean, n = length, xsd = sd) %>%
      within({
        m <- dim(xbar)[1]
        l <- dim(xbar)[2]

        x_test <- test_counts[[cues]]
        y_test <- as.numeric(test_counts[[group]])
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
  assert_that(is.mvg_ibbu_stanfit(stanfit),
              msg = paste0("stanfit must be of class ", new_stanfit_class_name))
  assert_that(is.mvg_ibbu_input(input),
              msg = "input is not an acceptable input data.")

  message("Currently this function is only checking whether input is a list. Use at your own risk.")
  stanfit@input = input

  return(stanfit)
}



#' @rdname compose_data
#' @export
compose_data_to_infer_prior_kappanu_via_conjugate_ibbu_w_sufficient_stats = function() {
  # in composing the data and fitting the model make sure that the model inherits
  # variable names and values for e.g., the categories and cues, so that they can
  # can be used in spread_draws and alike.
  message("This function is not doing anything yet.")

  # Make sure to hand through for the the test data, too, for which group / condition
  # it was collected. SPECIFCIALLY, ANNOTATE Y_TEST WITH THE GROUP CHARACTER LABELS.

  # see also tidybayes::compose_data
}
