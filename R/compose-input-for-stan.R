#' Transform and untransform cues by applying or undoing PCA, centering, and/or scaling.
#'
#' If the `transform` argument
#' is specified, the transforms in that object will be applied. This can be useful when the goal is to transform one
#' data set (e.g., test data) based on the statistics of the another data set (e.g.., training data). If no `transform`
#' is specified, then the transformations specified by the `center`, `scale`, and `pca` flags will be applied.
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

