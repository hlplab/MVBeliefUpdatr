NULL

#' Get covariance matrix from correlation matrix and standard deviations
#'
#' Get covariance matrix from correlation matrix and vector of standard deviations. Based on code recommended
#' by Ben Bolder for efficiency.
#'
#' @param omega A correlation matrix.
#' @param tau A vector of standard deviations.
#'
#' @return A square matrix.
#'
#' @seealso \code{\link[stats]{cov2cor}}, \code{\link{cov2tau}}
#' @references \url{https://stats.stackexchange.com/questions/62850/obtaining-covariance-matrix-from-correlation-matrix}
#' @keywords TBD
#' @examples
#' TBD
#' @export
cor2cov = function(omega, tau) {
  assert_that(is.matrix(omega))
  assert_that(is.numeric(tau))
  assert_that(nrow(omega) == ncol(omega),
              msg = "omega must be a square matrix.")
  assert_that(length(tau) == nrow(omega),
              msg = "tau must have as many elements as omega has rows (and columns).")

  outer(tau,tau) * omega
}

#' Get vector of standard deviations from covariance matrix
#'
#' Get vector of standard deviations from covariance matrix.
#'
#' @param v A covariance matrix.
#'
#' @return A numeric vector.
#'
#' @seealso \code{\link[stats]{cov2cor}}, \code{\link{cor2cov}}
#' @keywords TBD
#' @examples
#' TBD
#' @export
cov2tau = function(v) {
  assert_that(is.matrix(v))
  assert_that(nrow(v) == ncol(v),
              msg = "omega must be a square matrix.")
  sqrt(diag(v))
}


make_named_vector = function(x, names) {
  x = as.vector(x)
  names(x) = names
  return(x)
}

#' Combine a number of columns into a new vector column
#'
#' Combine a number of columns into a new column in which each cell is the vector of values from the original columns.
#'
#' @param data `tibble` or `data.frame`.
#' @param cols Vector of characters with names of variables to combine.
#' @param vector_col Name of new column of vectors.
#'
#' @return Same as \code{data}.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
make_vector_column = function(data, cols, vector_col, transmute = F) {
  # CHECK: expand to also handle quo input. (each instance of cols then needs to change)
  # then make_NIW_prior_from... can use this function
  data %<>%
    mutate(!! sym(vector_col) := pmap(.l = list(!!! syms(cols)),
                                      .f = function(...) {
                                        x = c(...)
                                        names(x) = cols
                                        return(x)
                                      })) %>%
    { if (transmute) select(., vector_col) else . }

  return(data)
}


#' Get sum of uncentered squares
#'
#' Get sum of uncentered squares. This quantity is a sufficient statistic for, for example, multivariate Gaussian
#' belief-updating under an Normal-Inverse-Wishart prior.
#'
#' @param data A `tibble`, `data.frame`, or `matrix`. If data is a `tibble` or `data.frame`, the columns for
#' specified variables are extracted and (together) converted into a matrix with as many colums as there are
#' variables.
#' @param variables Only required if data is not already a `matrix`.
#'
#' @return A matrix.
#'
#' @seealso
#' @keywords TBD
#' @examples
#' TBD
#' @rdname get_sum_of_squares
#' @export
#'
get_sum_of_squares <- function(data, variables = NULL, center = T, verbose = F) {
  assert_that(is_tibble(data) | is.data.frame(data) | is.matrix(data))
  if (is_tibble(data) | is.data.frame(data))
    assert_that(variables %in% names(data),
                msg = paste("Variable column(s)", variables[which(variables %nin% names(data))], "not found in data."))

  data.matrix = if (is_tibble(data) | is.data.frame(data)) {
    # Assume that the variables are to be combined into a data.matrix
    data %>%
      mutate_at(variables, unlist) %>%
      pull(variables) %>%
      as.matrix()
  } else data

  if (center)
    data.matrix = data.matrix - colMeans(data.matrix)

  k = dim(data.matrix)[2]
  m = matrix(ncol = k, nrow = k)

  for (i in 1:k) {
    m[i,i] = sum(data.matrix[,i]**2)
    if (i < k) for (j in (i + 1):k) {
      m[j,i] = sum(data.matrix[,i] * data.matrix[,j])
      m[i,j] = m[j,i]
    }
  }

  return(m)
}

#' @rdname get_sum_of_squares
#' @export
get_sum_of_uncentered_squares <- function(data, variables = NULL, verbose = F) {
  get_sum_of_squares(data = data, variables = variables, center = F, verbose = verbose)
}

#' @rdname get_sum_of_squares
#' @export
get_sum_of_centered_squares <- function(data, variables = NULL, verbose = F) {
  get_sum_of_squares(data = data, variables = variables, center = T, verbose = verbose)
}


#' Get sufficient statistics from a data set
#'
#' Get sufficient statistics from data. Calculates functions for the specified cues for
#' any combination of groups (optional) and categories, and returns them as a tibble.
#'
#' @param data `tibble` or `data.frame` with the data. Each row should be an observation of a category,
#' and contain information about the category label, the cue values of the observation, and optionally grouping variables.
#' @param test `tibble` or `data.frame` with the test data. Each row should be an observation, and contain information
#' about the cue values of the test stimulus and the participant's response.
#' @param cues Names of columns with cue values.
#' @param category Name of column that contains the category label for the exposure data. Can be `NULL` for unsupervised updating
#' (not yet implemented). (default: "category")
#' @param groups Name of column(s) that contains information about which observations form a group. This could be individual
#' subjects or conditions in an experiment. The latter is more efficient, but should only be used if exposure is
#' identical for every individual within the group. Test does not have to be identical for every individual within
#' the same group. For example, one can group multiple groups of subjects that have received the same exposure
#' but were tested on different test tokens.
#'
#' @return A tibble of sufficient statistics for each combination of category and group.
#'
#' @seealso
#' @keywords TBD
#' @examples
#' TBD
#' @export
get_sufficient_category_statistics <- function(
  data,
  cues,
  category,
  groups,
  ...
) {
  message("rows with missing values for cues will be ignored in the calculation of the sufficient statistics.")
  data_ss <- data %>%
    as_tibble() %>%
    group_by(!! sym(category), !!! syms(groups))

  if (length(cues) > 1) {
    # Multivariate observations
    data_ss %<>%
      drop_na(!!! syms(cues)) %>%
      summarise(
        x_N = length(!! sym(cues[1])),
        x_mean = list(colMeans(cbind(!!! syms(cues)))),
        x_ss = list(get_sum_of_uncentered_squares(cbind(!!! syms(cues)), verbose = verbose)),
        x_cov = list(cov(cbind(!!! syms(cues)))))
  } else {
    # Univariate observations
    data_ss %<>%
      drop_na(!!! syms(cues)) %>%
      summarise(
        x_N = length(!! sym(cues)),
        x_mean = mean(!! sym(cues)),
        x_ss = as.numeric(get_sum_of_uncentered_squares(matrix(!! sym(cues)), verbose = verbose)),
        x_sd = sd(!!! syms(cues)))
  }

  return(data_ss)
}

#' Transform and untransform cues by applying or undoing PCA, centering, and/or scaling.
#'
#' If the `transform.parameters` argument
#' is specified, the transforms in that object will be applied. This can be useful when the goal is to transform one
#' data set (e.g., test data) based on the statistics of the another data set (e.g., exposure data). If no
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
#' @param pca,unpca Should the data be transformed into/back from orthogonal principal components? If `TRUE`
#' then \code{center} and \code{scale} are default to `FALSE`. If PCA is applied, it is applied after
#' centering and/or scaling. (default: `FALSE`)
#' @param attach Should the transformed cues be attached to \code{data} or should just the transformed cues be
#' returned? (default: `TRUE`)
#' @param transform.parameters List of transforms (default: `NULL`)
#' @param return.transformed.data,return.transformed.data Should the (un)transformed data be returned? (default: `TRUE`)
#' @param return.transform.parameters Should the list of transforms be returned? (default: `FALSE`)
#' @param return.transform.function,return.untransform.function Should a function that applies the (un)transform be
#' returned? (default: `FALSE`)
#'
#' @return By default a \code{data.frame}. If `return.transform.parameters = TRUE`, a list of parameters. If
#' `return.transform.function = TRUE` or `return.untransform.function = TRUE` a function. If more than one of
#' these flags is `TRUE` then a list in which
#' the data element has name "data", the transform parameters have name "transform.parameters" and the (un)transform
#' function(s) have the name "(un)transform.function".
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @rdname transform_cues
#' @export
transform_cues = function(data, cues,
                          center =  if (pca) T else F, scale = F, pca = F,
                          attach = T,
                          transform.parameters = NULL,
                          return.transformed.data = T, return.transform.parameters = F,
                          return.transform.function = F, return.untransform.function = F
) {
  assert_that(is.data.frame(data) | is_tibble(data))
  assert_that(is.null(transform.parameters) | is.list(transform.parameters))
  old_data = data
  groups = if (length(groups(data)) == 0) character() else groups(data) %>% as.character()

  if (is.null(transform.parameters)) {
    transform.parameters = list()
    transform.parameters[["cue.labels"]] = cues

    if (pca) {
      transform.parameters[["pca"]] = data %>%
        select(all_of(cues)) %>%
        prcomp(center = center, scale. = scale, retx = F)
    } else {
      if (center) {
        transform.parameters[["center"]] = data %>%
          select(all_of(cues)) %>%
          summarise_all(list(mean = mean))
      }

      if (scale) {
        transform.parameters[["scale"]] = data %>%
          select(all_of(cues)) %>%
          summarise_all(list(sd = sd))
      }
    }
  }

  if (return.transformed.data) {
    if (!is.null(transform.parameters[["pca"]])) {
      data %<>%
        predict(transform.parameters[["pca"]], .) %>%
        as_tibble()
    } else {
      if (!is.null(transform.parameters[["center"]])) {
        data %<>%
          ungroup() %>%
          select(cues) %>%
          { . - (data %>%
                   left_join(transform.parameters[["center"]], by = groups) %>%
                   ungroup() %>%
                   select(all_of(paste0(cues, "_mean"))))
          }
      }

      if (!is.null(transform.parameters[["scale"]])) {
        data %<>%
          ungroup() %>%
          select(cues) %>%
          { . / (data %>%
                   left_join(transform.parameters[["scale"]], by = groups) %>%
                   ungroup() %>%
                   select(all_of(paste0(cues, "_sd"))))
          }
      }
    }
  }

  transform.function = if (!return.transform.function) NULL else {
    function(data) {
      cues = cues
      center = center
      scale = scale
      pca = pca

      transform_cues(data, cues, center = center, scale = scale, pca = pca,
                     transform.parameters = transform.parameters,
                     return.transformed.data = T, return.transform.parameters = F,
                     return.transform.function = F, return.untransform.function = F)

    }
  }

  untransform.function = if (!return.untransform.function) NULL else {
    untransform_cues(data, cues, uncenter = center, unscale = scale, unpca = pca,
                     transform.parameters = transform.parameters,
                     return.untransformed.data = F, return.untransform.function = T)
  }

  if (!identical(groups, character())) data %<>% group_by(!! sym(groups))
  if (attach) {
    data %<>%
      cbind(
        old_data %>%
          select(-intersect(names(old_data), names(data))),
        .)
  }

  if (return.transformed.data & !return.transform.parameters & !return.transform.function & !return.untransform.function) return(data) else
    if (!return.transformed.data & return.transform.parameters & !return.transform.function & !return.untransform.function) return(transform.parameters) else
      if (!return.transformed.data & !return.transform.parameters & return.transform.function & !return.untransform.function) return(transform.function) else
        if (!return.transformed.data & !return.transform.parameters & !return.transform.function & return.untransform.function) return(untransform.function) else
          return(
            list(
              data = if (return.transformed.data) data else NULL,
              transform.parameters = if (return.transform.parameters) transform.parameters else NULL,
              transform.function = if (return.transform.function) transform.function else NULL,
              untransform.function = if (return.untransform.function) untransform.function else NULL))
}


#' @rdname transform_cues
#' @export
untransform_cues = function(data, cues,
                            uncenter = NULL, unscale = NULL, unpca = NULL,
                            transform.parameters = NULL,
                            return.untransformed.data = T, return.untransform.function = F
) {
  assert_that(is.data.frame(data) | is_tibble(data))
  assert_that(!is.null(transform.parameters) & is.list(transform.parameters),
              msg = "Must provide transform parameters.")

  # By default untransform all transformations available in transform object
  if (is.null(unpca)) unpca = !is.null(transform.parameters[["pca"]])
  if (is.null(uncenter)) uncenter = !is.null(transform.parameters[["center"]])
  if (is.null(unscale)) unscale = !is.null(transform.parameters[["scale"]])

  if (unpca) {
    # https://stackoverflow.com/questions/29783790/how-to-reverse-pca-in-prcomp-to-get-original-data
    newcues = data %>%
      select(!!! rlang::syms(cues)) %>%
      as.matrix() %>%
      { . %*% t(transform.parameters[["pca"]]$rotation) } %>%
      { if (unscale & uncenter)
        t(t(.) * transform.parameters[["pca"]]$scale + transform.parameters[["pca"]]$center) else
          if (unscale) t(t(.) * transform.parameters[["pca"]]$scale) else
            if (uncenter) t(t(.) + transform.parameters[["pca"]]$center) }

    data %<>%
      select(-all_of(cues)) %>%
      cbind(newcues)

    data %<>%
      rename_at(cues, ~ transform.parameters[["cue.labels"]])
    cues = transform.parameters[["cue.labels"]]
  } else {
    if (unscale) {
      data %<>%
        ungroup() %>%
        select(cues) %>%
        { . * (data %>%
                 left_join(transform.parameters[["scale"]], by = groups) %>%
                 ungroup() %>%
                 select(all_of(paste0(cues, "_sd"))))
        }
    }

    if (uncenter) {
      data %<>%
        ungroup() %>%
        select(cues) %>%
        { . + (data %>%
                 left_join(transform.parameters[["center"]], by = groups) %>%
                 ungroup() %>%
                 select(all_of(paste0(cues, "_mean"))))
        }
    }
  }

  untransform.function = if (!return.untransform.function) NULL else {
    function(data) {
      cues = cues
      uncenter = uncenter
      unscale = unscale
      unpca = unpca

      untransform_cues(data, cues, uncenter = uncenter, unscale = unscale, unpca = unpca,
                       transform.parameters = transform.parameters,
                       return.untransformed.data = T, return.untransform.function = F)

    }
  }

  if (return.untransformed.data & !return.untransform.function) return(data) else
    if (!return.untransformed.data & return.untransform.function) return(untransform.function) else
      return(list(data = data, untransform.function = untransform.function))
}




