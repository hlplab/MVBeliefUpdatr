#' Transform and untransform cues by applying or undoing PCA, centering, and/or scaling.
#'
#' If the `transform` argument
#' is specified, the transforms in that object will be applied. This can be useful when the goal is to transform one
#' data set (e.g., test data) based on the statistics of the another data set (e.g.., training data). If no `transform`
#' is specified, then the transformations specified by the `center`, `scale`, and `pca` flags will be applied.
#'
#' @param data `tibble` or `data.frame`.
#' @param cues Vector of characters with names of cue variables.
#' @param transform List of transforms (default: `NULL`)
#' @param return.transform Should the list of transforms be returned along with the data? (default: `FALSE`)
#' @param center Should the data be centered? (default: `TRUE`)
#' @param scale Should the data be standardized? (default: `TRUE`)
#' @param pca Should the data be transformed into orthogonal principal components? (default: `FALSE`)
#'
#' @return Data frame, unless `return.transform = T`. In that case, a list with two elements (`data` and `transform`).
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @rdname transform_cues
#' @export
#'
transform_cues = function(data, cues,
                          transform = NULL, return.transform = F,
                          center = T, scale = T, pca = F) {
  assert_that(is.data.frame(data) | is_tibble(data))
  assert_that(is.null(transform) | is.list(transform))

  if (is.null(transform)) {
    transform = list()

    if (pca) {
      pca <- data %>%
        select(!!! rlang::syms(cues)) %>%
        prcomp(center = center, scale. = scale)

      print(summary(pca)$importance)
      center = FALSE
      scale = FALSE

      transform[["pca"]] = pca
    }

    if (center) {
      transform[["center"]] = data %>%
        select(!!! rlang::syms(cues)) %>%
        summarise_all(mean)
    }

    if (scale) {
      transform[["scale"]] = data %>%
        select(!!! rlang::syms(cues)) %>%
        summarise_all(sd)
    }
  }

  if (pca) {
    data %<>%
      cbind(predict(transform[["pca"]], data))
    center = FALSE
    scale = FALSE
  }

  if (center) {
    newcues = data %>%
      select(!!! rlang::syms(cues)) %>%
      sweep(2, as.numeric(transform[["center"]]), FUN = "-")

    data %<>%
      select(-all_of(cues)) %>%
      cbind(newcues)
  }

  if (scale) {
    newcues = data %>%
      select(!!! rlang::syms(cues)) %>%
      sweep(2, as.numeric(transform[["scale"]]), FUN = "/")

    data %<>%
      select(-all_of(cues)) %>%
      cbind(newcues)
  }

  if (return.transform) return(list(data = data, transform = transform)) else return(data)
}


#' @rdname transform_cues
#' @export
untransform_cues = function(data, cues,
                            transform = NULL,
                            center = NULL, scale = NULL, pca = NULL) {
  assert_that(is.data.frame(data) | is_tibble(data))
  assert_that(is.list(transform))

  # By default untransform all transformations available in transform object
  if (is.null(pca)) pca = !is.null(transform[["pca"]])
  if (is.null(center)) center = !is.null(transform[["center"]])
  if (is.null(scale)) scale = !is.null(transform[["center"]])

  if (pca) {
    stop("PCA untransform not yet implemented!")
    data %<>%
      cbind(predict(transform[["pca"]], data))
    center = FALSE
    scale = FALSE
  }

  if (scale) {
    newcues = data %>%
      select(!!! rlang::syms(cues)) %>%
      sweep(2, as.numeric(transform[["scale"]]), FUN = "*")

    data %<>%
      select(-all_of(cues)) %>%
      cbind(newcues)
  }

  if (center) {
    newcues = data %>%
      select(!!! rlang::syms(cues)) %>%
      sweep(2, as.numeric(transform[["center"]]), FUN = "+")

    data %<>%
      select(-all_of(cues)) %>%
      cbind(newcues)
  }

  if (return.transform) return(list(data = data, transform = transform)) else return(data)
}


#' Compose data for input to RStan
#'
#' Take exposure and test data as input, and prepare the data for input into an MVBeliefUpdatr Stan program.
#'
#' @param exposure `tibble` or `data.frame` with the exposure data. Each row should be an observation of a category,
#' and contain information about the category label, the cue values of the observation, and optionally grouping variables.
#' @param test `tibble` or `data.frame` with the test data. Each row should be an observation, and contain information
#' about the cue values of the test stimulus and the participant's response.
#' @param cue.labels Vector of characters with names of cue variables.
#' @param group Column that contains information about which observations form a group. This could be individual
#' subjects or conditions in an experiment. The latter is more efficient, but should only be used if exposure is
#' identical for every individual within the group. Test does not have to be identical for every individual within
#' the same group. For example, one could group multiple groups of subjects that have received the same exposure
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
#' @seealso \code{link{is.mvg_ibbu_input}}
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
  training, test,
  cues, group, category, response,
  center.observations = T, scale.observations = T, pca.observations = F, pca.cutoff = .95,
  tau_scale, L_omega_scale,
  debug
) {
  # in composing the data and fitting the model make sure that the model inherits
  # variable names and values for e.g., the categories and cues, so that they can
  # can be used in spread_draws and alike.
  error("This function is not yet implemented.")

  # Make sure to hand through for the the test data, too, for which group / condition
  # it was collected. SPECIFCIALLY, ANNOTATE Y_TEST WITH THE GROUP CHARACTER LABELS.

  # see also tidybayes::compose_data
    training <- check_training_data(training, category, group, debug)
    test_counts <- test_counts(training, test, cue, category, response, group, debug)

    if (center.observations | scale.observations) {
      means.training = training %>% ungroup() %>% select(cue) %>% summarise_at(cue, .funs = mean, na.rm = T)
      sds.training = training %>% ungroup() %>% select(cue) %>% summarise_at(cue, .funs = sd, na.rm = T)

      ## -----------------------------------------------------------------------
      # There must be a more elegant way to subtract the training means
      ##  -----------------------------------------------------------------------
      for (i in cue) {
        # Note i here is an element of cue (i.e., a character)
        training[[i]] = training[[i]] - means.training[[i]]
        test_counts[[i]] = test_counts[[i]] - means.training[[i]]

        if (scale.observations) {
          training[[i]] = training[[i]] / sds.training[[i]]
          test_counts[[i]] = test_counts[[i]] / sds.training[[i]]
        }
      }
    }

    ## Category-by-subject/group matrices of sufficient stats
    if (useMultivariateUpdating) {
      # Multivariate input
      data_list <- training %>%
        training_ss_matrix(list(category, group), cue = cue,
                           useMultivariateUpdating = useMultivariateUpdating, debug = debug,
                           x_mean = colMeans, N = length, x_ss = sum_uncentered_squares) %>%
        within({
          M <- dim(x_mean)[1]
          L <- dim(x_mean)[2]
          K <- length(cue)

          x_test <- cbind(test_counts[,cue])
          y_test <- as.numeric(test_counts[[group]])
          z_test_counts <-
            test_counts %>%
            select(.dots=levels(training[[category]])) %>%
            as.matrix()

          N_test <- nrow(x_test)

          tau_scale <- tau_scale
          L_omega_scale <- L_omega_scale

        })
    } else {
      # Univariate input
      data_list <- training %>%
        training_ss_matrix(list(category, group), cue = cue,
                           useMultivariateUpdating = useMultivariateUpdating, debug = debug,
                           xbar = mean, n = length, xsd = sd) %>%
        within({
          m <- dim(xbar)[1]
          l <- dim(xbar)[2]

          x_test <- test_counts[[cue]]
          y_test <- as.numeric(test_counts[[group]])
          z_test_counts <-
            test_counts %>%
            select(.dots=levels(training[[category]])) %>%
            as.matrix()

          n_test <- length(x_test)
          mu_0_sd <- 0
          sigma_0_sd <- 0
        })
    }

    if (debug) {
      print("In prepare_data_conj_suff_stats_infer_prior():")
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

