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
#'
transform_cues = function(data, cues,
                          center =  T, scale = if (pca) F else T, pca = F,
                          transform.parameters = NULL,
                          return.transformed.data = T, return.transform.parameters = F,
                          return.transform.function = F, return.untransform.function = F
) {
  assert_that(is.data.frame(data) | is_tibble(data))
  assert_that(is.null(transform.parameters) | is.list(transform.parameters))

  if (is.null(transform.parameters)) {
    transform.parameters = list()

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

    if (pca) {
      PCA <- data %>%
        select(!!! rlang::syms(cues)) %>%
        prcomp(center = F, scale. = F)
      transform.parameters[["pca"]] = PCA
    }
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

  if (pca) {
    data %<>%
      cbind(predict(transform.parameters[["pca"]], data))
  }


  transform.function = if (!return.transform.function) NULL else {
    function(data) {
      cues = cues
      center = center
      scale = scale
      pca = pca

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
  assert_that(is.list(transform.parameters))

  # By default untransform all transformations available in transform object
  if (is.null(unpca)) pca = !is.null(transform.parameters[["pca"]])
  if (is.null(uncenter)) center = !is.null(transform.parameters[["center"]])
  if (is.null(unscale)) scale = !is.null(transform.parameters[["center"]])

  if (unpca) {
    # Since we're uncentering and unscaling separately, the following is not needed:
    # https://stackoverflow.com/questions/29783790/how-to-reverse-pca-in-prcomp-to-get-original-data
    # pca = transform.parameters[["pca"]]
    # t(t(pca$x %*% t(pca$rotation)) + pca$center)
    # If pca$scale is TRUE you will also need to re-scale
    # t(t(pca$x %*% t(pca$rotation)) * pca$scale + pca$center)
    newcues = data %>%
      select(!!! rlang::syms(cues)) %>%
      as.matrix() %>%
      { . %*% t(transform.parameters[["pca"]]$rotation) }

    data %<>%
      select(-all_of(cues)) %>%
      cbind(newcues)
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
#' @export
get_sum_of_uncentered_squares <- function(data, variables = NULL, verbose = F) {
  assert_that(is_tibble(data) | is.data.frame(data) | is.matrix(data))
  if (is_tibble(data) | is.data.frame(data))
    assert_that(variables %in% data,
                msg = paste("Variable column(s)", variables[which(variables %nin% names(data))], "not found in data."))

  data.matrix = if (is_tibble(data) | is.data.frame(data))
    data[[variables]] %>%
    as.matrix() else data

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

