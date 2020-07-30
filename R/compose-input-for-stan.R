#' Transform and untransform cues by applying or undoing PCA, centering, and/or scaling.
#'
#' If the `transform.parameters` argument
#' is specified, the transforms in that object will be applied. This can be useful when the goal is to transform one
#' data set (e.g., test data) based on the statistics of the another data set (e.g.., training data). If no
#' `transform.parameters` are specified, then the transformations specified by the `center`, `scale`, and `pca`
#' flags will be applied. The transform and untransform functions can also return \emph{functions} that perform their
#' actions for the specific cues. This can be helpful if one wants to store those functions. For example, the
#' transform function can return both the transform and untransform functions necessary to perform the specified
#' centering, scaling, and/or pca \emph{and} to undo those transformations.
#'
#' @param data `tibble` or `data.frame`.
#' @param cues Vector of characters with names of cue variables.
#' @param center,uncenter Should the data be (un)centered? (default: `TRUE` unless `pca = TRUE`)
#' @param scale,unscale Should the data be (un)standardized? (default: `TRUE` unless `pca = TRUE`)
#' @param pca,unpca Should the data be transformed into/back from orthogonal principal components? If `TRUE` then \code{center} and
#' \code{scale} are default to `FALSE`. (default: `FALSE`)
#' @param transform.parameters List of transforms (default: `NULL`)
#' @param return.transformed.data,return.transformed.data Should the (un)transformed data be returned? (default: `TRUE`)
#' @param return.transform.parameters Should the list of transforms be returned? (default: `FALSE`)
#' @param return.transform.function,return.transform.function Should a function that applies the (un)transform be
#' returned? (default: `FALSE`)
#'
#' @return By default a \code{data.frame}. If `return.transform.parameters = TRUE`, a list of parameters. If
#' `return.transform.function = TRUE` a function. If more than one of these flags is `TRUE` then a list in which
#' the data element has name "data", the transform parameters have name "transform.parameters" and the transform
#' function has the name "transform.function".
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @rdname transform_cues
#' @export
#'
transform_cues = function(data, cues,
                          center =  if (pca) F else T, scale = if (pca) F else T, pca = F,
                          transform.parameters = NULL,
                          return.transformed.data = T, return.transform.parameters = F,
                          return.transform.function = F, return.untransform.function = F
) {
  assert_that(is.data.frame(data) | is_tibble(data))
  assert_that(is.null(transform) | is.list(transform))

  if (is.null(transform.parameters)) {
    transform.parameters = list()

    if (pca) {
      pca <- data %>%
        select(!!! rlang::syms(cues)) %>%
        prcomp(center = center, scale. = scale)

      print(summary(pca)$importance)
      transform.parameters[["pca"]] = pca
    }

    if (center) {
      transform.parameters[["center"]] = data %>%
        select(!!! rlang::syms(cues)) %>%
        summarise_all(mean)
    }

    if (scale) {
      transform.parameters[["scale"]] = data %>%
        select(!!! rlang::syms(cues)) %>%
        summarise_all(sd)
    }
  }

  if (pca) {
    data %<>%
      cbind(predict(transform.parameters[["pca"]], data))
    center = FALSE
    scale = FALSE
  }

  if (center) {
    newcues = data %>%
      select(!!! rlang::syms(cues)) %>%
      sweep(2, as.numeric(transform.parameters[["center"]]), FUN = "-")

    data %<>%
      select(-all_of(cues)) %>%
      cbind(newcues)
  }

  if (scale) {
    newcues = data %>%
      select(!!! rlang::syms(cues)) %>%
      sweep(2, as.numeric(transform.parameters[["scale"]]), FUN = "/")

    data %<>%
      select(-all_of(cues)) %>%
      cbind(newcues)
  }

  transform.function = if (!return.transform.function) NULL else {
    function(data,
             center = center, scale = scale, pca = pca) {
      cues = cues

      transform_cues(data, cues, center = center, scale = scale, pca = pca,
                     transform.parameters = transform.parameters,
                     return.transformed.data = T, return.transform.parameters = F, return.transform.function = F)

    }
  }

  untransform.function = if (!return.untransform.function) NULL else {
    untransform_cues(data, cues, uncenter = center, unscale = scale, unpca = pca,
                     transform.parameters = transform.parameters,
                     return.untransformed.data = F, return.untransform.function = T)
  }

  if (return.transformed.data & !return.transform.parameters & !return.transform.function & !return.untransform.function) return(data) else
    if (!return.transformed.data & return.transform.parameters & !return.transform.function & !return.untransform.function) return(transform.parameters) else
      if (!return.transformed.data & !return.transform.parameters & return.transform.function & !return.untransform.function) return(transform.function) else
        if (!return.transformed.data & !return.transform.parameters & !return.transform.function & return.untransform.function) return(untransform.function) else
          return(list(data = data, transform.parameters = transform.parameters, transform.function = transform.function, untransform.function = untransform.function))
}


#' @rdname transform_cues
#' @export
untransform_cues = function(data, cues,
                            uncenter = NULL, unscale = NULL, unpca = NULL,
                            transform.parameters = NULL,
                            return.untransformed.data = T, return.untransform.function = F
) {
  assert_that(is.data.frame(data) | is_tibble(data))
  assert_that(is.list(transform.parameters))

  # By default untransform all transformations available in transform object
  if (is.null(unpca)) pca = !is.null(transform.parameters[["pca"]])
  if (is.null(uncenter)) center = !is.null(transform.parameters[["center"]])
  if (is.null(unscale)) scale = !is.null(transform.parameters[["center"]])

  if (unpca) {
    stop("PCA untransform not yet implemented!")
    data %<>%
      cbind(predict(transform.parameters[["pca"]], data))
  }

  if (unscale) {
    newcues = data %>%
      select(!!! rlang::syms(cues)) %>%
      sweep(2, as.numeric(transform.parameters[["scale"]]), FUN = "*")

    data %<>%
      select(-all_of(cues)) %>%
      cbind(newcues)
  }

  if (uncenter) {
    newcues = data %>%
      select(!!! rlang::syms(cues)) %>%
      sweep(2, as.numeric(transform.parameters[["center"]]), FUN = "+")

    data %<>%
      select(-all_of(cues)) %>%
      cbind(newcues)
  }

  untransform.function = if (!return.transform.function) NULL else {
    function(data,
             uncenter = uncenter, unscale = unscale, unpca = unpca) {
      cues = cues

      untransform_cues(data, cues, uncenter = uncenter, unscale = unscale, unpca = unpca,
                     transform.parameters = transform.parameters,
                     return.untransformed.data = T, return.untransform.function = F)

    }
  }

  if (return.untransformed.data & !return.untransform.function) return(data) else
    if (!return.untransformed.data & return.untransform.function) return(untransform.function) else
        return(list(data = data, untransform.function = untransform.function))
}



check_exposure_test_data <- function(data, cues, category, response, group, which.data = "the") {
  assert_that(is_tibble(exposure) | is.data.frame(exposure))
  assert_that(all(is_character(cues)),
              msg = "cues must be a column name or vector of column names.")
  assert_that(cues %in% names(data),
              msg = paste("Cue column(s)", cues[which(cues %nin% names(data))], "not found in", which.data, "data." ))

  if(!null(category)) {
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
                msg = paste0("Response column", response, "not found in", which.data, "data."))

    data %<>%
      mutate_at(response, as.factor)
  }

  if (!is.null(group)) {
    assert_that(group %in% names(data),
                msg = paste0("Group column", group, "not found in", which.data,"data."))

    data %<>%
      mutate_at(group, as.factor)
  }

  if (verbose){
    print("In exposure_test_data():")
    print (data)
  }


  return(data)
}



get_test_counts <- function(training, test, cue, category, response, group, verbose) {
  test_counts <- test %>%
    as_tibble() %>%
    group_by(!!! rlang::syms(group),
             !!! rlang::syms(cue),
             !! rlang::sym(response)) %>%
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
get_sufficient_statistics_from_data <- function(exposure, cues, category, group, verbose = F, ...) {
  exposure = check_exposure_test_data(exposure, cues, category, NULL, group, verbose)

  if (length(cues) > 1) {
    # Multivariate observations
    exposure_ss <- exposure %>%
      as_tibble() %>%
      group_by(!!! rlang::syms(groupings)) %>%
      summarise(
        N = length(!! sym(cues[1])),
        x_mean = list(colMeans(cbind(!!! rlang::syms(cues)))),
        x_ss = list(get_sum_of_uncentered_squares(cbind(!!! rlang::syms(cues)), verbose = verbose))
      )

    if (verbose) {
      print("In get_sufficient_statistics_from_data(), multivariate sum-of-uncentered-squares matrix:")
      print(exposure_ss)
    }

    ## -------------------------------------------------------------------------------
    # This is intended to map the elements of the tibble into arrays of the required
    # dimensionality. THERE PROBABLY IS A MUCH MORE CONCISE AND PERHAPS MORE EFFICIENT
    # WAY TO DO THIS. CHECK BACK.
    #
    # For helpful concise info on tibbles, see
    #   https://cran.r-project.org/web/packages/tibble/vignettes/tibble.html
    ## -------------------------------------------------------------------------------
    cats = levels(exposure[[groupings[[1]]]])
    subjs = levels(exposure[[groupings[[2]]]])
    n_cat = length(cats)
    n_subj = length(subjs)
    n_cues = length(cues)

    N = array(dim = c(n_cat,n_subj))
    x_mean = array(dim = c(n_cat,n_subj,n_cues))
    x_ss = array(dim = c(n_cat,n_subj,n_cues,n_cues))

    for (i in 1:n_cat) {
      for (j in 1:n_subj) {
        temp.exposure_ss = exposure_ss %>%
          ungroup() %>%
          filter(!! rlang::sym(groupings[[1]]) == cats[i] &
                   !! rlang::sym(groupings[[2]]) == subjs[j])

        if (nrow(temp.exposure_ss) > 0) {
          N[i,j] = temp.exposure_ss$N[[1]]
          x_mean[i,j,] = temp.exposure_ss$x_mean[[1]]
          x_ss[i,j,,] = temp.exposure_ss$x_ss[[1]]
        } else {
          N[i,j] = 0
          x_mean[i,j,] = rep(0, length(cues))
          x_ss[i,j,,] = diag(length(cues))
        }
      }
    }

    exposure_ss = list(N = N, x_mean = x_mean, x_ss = x_ss)
  } else {
    # Univariate observations
    exposure_ss <-
      exposure %>%
      group_by(!!! rlang::syms(groupings)) %>%
      summarise_at(cues, .funs=funs(...))

    if (verbose) {
      print("In get_sufficient_statistics_from_data(), univariate sum-of-squares matrix (prior to map application):")
      print(exposure_ss)
    }

    stats <- names(list(...))

    exposure_ss <- map(stats, ~ reshape2::acast(exposure_ss, as.list(groupings), value.var=.x)) %>%
      set_names(stats)
  }


  if (debug) {
    print("In get_sufficient_statistics_from_data(), sum-of-squares matrix (uncentered for multivariate data, centered for univariate data):")
    print(exposure_ss)
  }

  return(exposure_ss)
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
#' `pca.observation = FALSE`
#'
#' @return A list that is an \code{mvg_ibbu_input}.
#'
#' @seealso \code{\link{is.mvg_ibbu_input}}
#' @keywords TBD
#' @examples
#' TBD
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


#' @rdname compose_data
#' @export
compose_data_to_infer_prior_via_conjugate_ibbu_w_sufficient_stats = function(
  exposure, test,
  cues, category = "category", response = "response", group,
  center.observations = T, scale.observations = T, pca.observations = F, pca.cutoff = .95,
  tau_scale, L_omega_scale,
  verbose = F
) {
  exposure <- check_exposure_test_data(exposure, cues, category, NULL, group, verbose)
  test <- check_exposure_test_data(test, cues, NULL, response, group, verbose)

  assert_that(all(levels(exposure[[category]]) == levels(test[[response]])),
              msg = paste("category column", category, "in exposure and response colum", response, "must be factors
              with the same levels in the same order in. Either the levels do not match, or they are not in the same
              order."))
  if (!is.null(group))
    assert_that(all(levels(exposure[[group]]) == levels(test[[group]])),
                msg = paste("group column", category, "must be a factor with the same levels in the same order in
                          the exposure and test data. Either the levels do not match, or they are not in the same
                          order."))

  test_counts <- get_test_counts(training, test, cue, category, response, group, verbose)

  if (length(cues) > 1) {
    data_list <- training %>%
      get_sufficient_statistics_from_data(
        list(category, group), cues = cues,
        verbose = verbose,
        x_mean = colMeans, N = length, x_ss = sum_uncentered_squares) %>%
      within({
        M <- dim(x_mean)[1]
        L <- dim(x_mean)[2]
        K <- length(cues)

        m <- dim(xbar)[1]
        l <- dim(xbar)[2]

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
    data_list <- training %>%
      get_sufficient_statistics_from_data(
        list(category, group), cues = cues,
        verbose = verbose,
        xbar = mean, n = length, xsd = sd) %>%
      within({
        m <- dim(xbar)[1]
        l <- dim(xbar)[2]

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

