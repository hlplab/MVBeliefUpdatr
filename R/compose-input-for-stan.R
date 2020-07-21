#' Transform and untransform cues by applying or undoing PCA, centering, and/or scaling.
#'
#' If the `transform` argument
#' is specified, the transforms in that object will be applied. This can be useful when the goal is to transform one
#' data set (e.g., test data) based on the statistics of the another data set (e.g.., training data). If no `transform`
#' is specified, then the transformations specified by the `center`, `scale`, and `PCA` flags will be applied.
#'
#' @param data `tibble` or `data.frame`.
#' @param cues Vector of characters with names of cue variables.
#' @param transform List of transforms (default: `NULL`)
#' @param return.transform Should the list of transforms be returned along with the data? (default: `FALSE`)
#' @param center Should the data be centered? (default: `TRUE`)
#' @param scale Should the data be standardized? (default: `TRUE`)
#' @param PCA Should the data be transformed into orthogonal principal components? (default: `FALSE`)
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
                          center = T, scale = T, PCA = F) {
  assert_that(is.data.frame(data) | is_tibble(data))
  assert_that(is.null(transform) | is.list(transform))

  if (is.null(transform)) {
    transform = list()

    if (PCA) {
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

  if (PCA) {
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
                            center = NULL, scale = NULL, PCA = NULL) {
  assert_that(is.data.frame(data) | is_tibble(data))
  assert_that(is.list(transform))

  # By default untransform all transformations available in transform object
  if (is.null(PCA)) PCA = !is.null(transform[["PCA"]])
  if (is.null(center)) center = !is.null(transform[["center"]])
  if (is.null(scale)) scale = !is.null(transform[["center"]])

  if (PCA) {
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

#' @export
compose_data = function() {
  # in composing the data and fitting the model make sure that the model inherits
  # variable names and values for e.g., the categories and cues, so that they can
  # can be used in spread_draws and alike.
  message("This function is not doing anything yet.")

  # Make sure to hand through for the the test data, too, for which group / condition
  # it was collected. SPECIFCIALLY, ANNOTATE Y_TEST WITH THE GROUP CHARACTER LABELS.

  # see also tidybayes::compose_data
}


make_standata = function() {
  message("This function is not doing anything yet.")
  # Check brms::make_standata
}


attach_stanfit_input_data = function(stanfit, input) {
  assert_that(is.mvg_ibbu_stanfit(stanfit),
              msg = paste0("stanfit must be of class ", new_stanfit_class_name))
  assert_that(is.list(input),
              msg = "input must be a list.")

    message("Currently this function is only checking whether input is a list.")
    stanfit@input = input

    return(stanfit)
}

