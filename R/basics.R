NULL

#' Get sum-of-squares matrix
#'
#' Get sum-of-square matrix.
#'
#' @param x A matrix of observations, with each row being a vector observation.
#' @param centered Should the centered sum-of-squares be returned (`TRUE`) or the uncentered (`FALSE`)? (default: `TRUE`)
#'
#' @return A square matrix.
#'
#' @seealso \code{\link{css2cov}}, \code{\link{cov2css}}, \code{\link{uss2css}}
#' @keywords TBD
#' @examples
#' TBD
#' @export
ss <- function(x, center = TRUE) {
  if (center) {
    # Benchmarked to be more efficient than x - 1 %*% t(1) %*% x
    xm <- colMeans(x)
    xm <- matrix(xm, nrow = nrow(x), ncol = length(xm), byrow = T)
    x <- x - xm
  }

  # Benchmarked to be more efficient than t(x) %*% x
  k = dim(x)[2]
  m = matrix(ncol = k, nrow = k)
  for (i in 1:k) {
    m[i,i] = sum(x[,i]**2)
    if (i < k) for (j in (i + 1):k) {
      m[j,i] = sum(x[,i] * x[,j])
      m[i,j] = m[j,i]
    }
  }

  return(m)
}


#' Convert sum-of-square and covariance matrices
#'
#' Convert centered sum-of-square matrices (css), uncentered sum-of-square matrices (uss),
#' and covariance matrices (cov) into each other.
#'
#' @param uss,css,cov Uncentered sum-of-square, centered sum-of-square, and covariance matrix.
#' @param n Number of observations that have gone into the sum-of-square matrix.
#' @param mean Mean of the observations that have gone into the sum-of-square matrix.
#'
#' @return A square matrix.
#'
#' @seealso \code{\link{ss}}, \code{\link[stats]{cov2cor}}, \code{\link{cor2cov}}, \code{\link{cov2tau}}
#' @keywords TBD
#' @examples
#' TBD
#' @rdname uss2css
#' @export
uss2css <- function(uss, n, mean) {
  assert_that(is.matrix(uss))
  assert_that(length(dim(uss)) == 2)
  assert_that(dim(uss)[[1]] == dim(uss)[[2]],
              msg = "uss must be a square matrix.")
  assert_that(length(mean) == dim(uss)[[1]],
              msg = "uss and mean are not of compatible dimensions.")

  xm <- matrix(mean, nrow = n, ncol = length(mean), byrow = T)
  css <- uss - ss(xm, center = F)
  return(css)
}


#' @rdname uss2css
#' @export
css2uss <- function(css, n, mean) {
  assert_that(is.matrix(css))
  assert_that(length(dim(css)) == 2)
  assert_that(dim(css)[[1]] == dim(css)[[2]],
              msg = "uss must be a square matrix.")
  assert_that(length(mean) == dim(css)[[1]],
              msg = "uss and mean are not of compatible dimensions.")

  xm <- matrix(mean, nrow = n, ncol = length(mean), byrow = T)
  uss <- css + ss(xm, center = F)
  return(uss)
}

#' @rdname uss2css
#' @export
css2cov <- function(css, n) {
  assert_that(is.matrix(css))
  assert_that(length(dim(css)) == 2)
  assert_that(dim(css)[[1]] == dim(css)[[2]],
              msg = "css must be a square matrix.")

  return(css / (n - 1))
}

#' @rdname uss2css
#' @export
cov2css <- function(cov, n) {
  assert_that(is.matrix(cov))
  assert_that(length(dim(cov)) == 2)
  assert_that(dim(cov)[[1]] == dim(cov)[[2]],
              msg = "uss must be a square matrix.")

  return(cov * (n - 1))
}

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

make_named_square_matrix = function(x, names) {
  x = matrix(x, nrow = sqrt(length(x)), dimnames = list(names, names))
  return(x)
}

#' Combine a number of columns into a new vector column
#'
#' Combine a number of columns into a new column in which each cell is the vector of values from the original columns.
#'
#' @param data `tibble` or `data.frame`.
#' @param cols Vector of characters with names of variables to combine.
#' @param vector_col Name of new column of vectors.
#' @param .keep See \code{\link{dplyr::mutate}}.
#'
#' @return Same as \code{data}.
#'
#' @keywords TBD
#' @examples
#' TBD
#' @export
make_vector_column = function(data, cols, vector_col, .keep = "all") {
  # CHECK: expand to also handle quo input. (each instance of calls then needs to change)
  # then make_NIW_prior_from...  use this function
  data %<>%
    mutate(!! sym(vector_col) := pmap(
      .l = list(!!! syms(cols)),
      .f = function(...) {
        x = c(...)
        names(x) = cols
        return(x)
      }),
      .keep = .keep)

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
#' @keywords TBD
#' @examples
#' TBD
#' @rdname get_sum_of_squares_from_df
#' @export
get_sum_of_squares_from_df <- function(data, variables = NULL, center = T, verbose = F) {
  assert_that(is_tibble(data) | is.data.frame(data) | is.matrix(data))
  if (is_tibble(data) | is.data.frame(data))
    assert_that(all(variables %in% names(data)),
                msg = paste("Variable column(s)", variables[which(variables %nin% names(data))], "not found in data."))

  data.matrix <- if (is_tibble(data) | is.data.frame(data)) {
    # Assume that the variables are to be combined into a data.matrix
    data %>%
      mutate_at(variables, unlist) %>%
      select(all_of(variables)) %>%
      as.matrix()
  } else data

  ss(data.matrix, center = center)
}

#' @rdname get_sum_of_squares_from_df
#' @export
get_sum_of_uncentered_squares_from_df <- function(data, variables = NULL, verbose = F) {
  get_sum_of_squares_from_df(data = data, variables = variables, center = F, verbose = verbose)
}

#' @rdname get_sum_of_squares_from_df
#' @export
get_sum_of_centered_squares_from_df <- function(data, variables = NULL, verbose = F) {
  get_sum_of_squares_from_df(data = data, variables = variables, center = T, verbose = verbose)
}


#' Get sufficient statistics from a data set
#'
#' Get sufficient statistics from data. Calculates functions for the specified cues for
#' any combination of groups (optional) and categories, and returns them as a tibble.
#' Rows with missing values for cues will be ignored in the calculation of the sufficient statistics.
#'
#' @param data `tibble` or `data.frame` with the data. Each row should be an observation of a category,
#' and contain information about the category label, the cue values of the observation, and optionally grouping variables.
#' @param test `tibble` or `data.frame` with the test data. Each row should be an observation, and contain information
#' about the cue values of the test stimulus and the participant's response.
#' @param cues Names of columns with cue values.
#' @param category Name of column that contains the category label for the exposure data. Can be `NULL` for unsupervised updating
#' (not yet implemented). (default: "category")
#' @param group Name of column(s) that contains information about which observations form a group. This could be individual
#' subjects or conditions in an experiment. The latter is more efficient, but should only be used if exposure is
#' identical for every individual within the group. Test does not have to be identical for every individual within
#' the same group. For example, one can group multiple groups of subjects that have received the same exposure
#' but were tested on different test tokens. If `NULL` no grouping variable will be considered. (default: `NULL`)
#' @param categories,groups Character vector of categories/groups to be summarize. If `NULL`, all categories/groups will be
#' included. (default: `NULL`)
#'
#' @return A tibble of sufficient statistics for each combination of category and group. This includes the count, mean,
#' uncentered and centered sums-of-squares and the covariance matrix (or, for univariate, stimuli: the standard deviation).
#'
#' @keywords TBD
#' @examples
#' TBD
#' @export
get_sufficient_category_statistics <- function(
  data,
  cues,
  category = "category",
  group = NULL,
  categories = NULL,
  groups = NULL,
  ...
) {
  data_ss <- data %>%
    as_tibble(.name_repair = "minimal") %>%
    group_by(!! sym(category), !!! syms(group)) %>%
    { if(!is.null(categories)) filter(., !! sym(category) %in% categories) else . } %>%
    { if(!is.null(groups)) filter(., !! sym(group) %in% groups) else . } %>%
    droplevels()

  if (length(cues) > 1) {
    # Multivariate observations
    data_ss %<>%
      drop_na(!!! syms(cues)) %>%
      summarise(
        x_N = length(!! sym(cues[1])),
        x_mean = list(colMeans(cbind(!!! syms(cues)))),
        x_uss = list(get_sum_of_uncentered_squares_from_df(cbind(!!! syms(cues)), verbose = verbose)),
        x_css = list(get_sum_of_centered_squares_from_df(cbind(!!! syms(cues)), verbose = verbose)),
        x_cov = list(cov(cbind(!!! syms(cues)))))
  } else {
    # Univariate observations
    data_ss %<>%
      drop_na(!!! syms(cues)) %>%
      summarise(
        x_N = length(!! sym(cues)),
        x_mean = mean(!! sym(cues)),
        x_uss = as.numeric(get_sum_of_uncentered_squares_from_df(matrix(!! sym(cues)), verbose = verbose)),
        x_css = as.numeric(get_sum_of_centered_squares_from_df(matrix(!! sym(cues)), verbose = verbose)),
        x_sd = sd(!! sym(cues)))
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
#' @keywords TBD
#' @examples
#' TBD
#' @rdname transform_cues
#' @export
transform_cues <- function(
    data, cues,
    center =  if (pca) T else F, scale = F, pca = F,
    attach = T,
    transform.parameters = NULL,
    return.transformed.data = T, return.transform.parameters = F,
    return.transform.function = F, return.untransform.function = F
) {
  assert_that(is.data.frame(data) | is_tibble(data))
  assert_that(is.null(transform.parameters) | is.list(transform.parameters))
  old_data <- data
  groups <- if (length(groups(data)) == 0) character() else groups(data) %>% as.character()

  if (is.null(transform.parameters)) {
    transform.parameters = list()
    transform.parameters[["cue.labels"]] <- cues

    if (pca) {
      transform.parameters[["pca"]] <-
        data %>%
        select(all_of(cues)) %>%
        prcomp(center = center, scale. = scale, retx = F)
    } else {
      if (center) {
        transform.parameters[["center"]] <-
          data %>%
          select(all_of(cues)) %>%
          summarise_all(list(mean = mean))
      }

      if (scale) {
        transform.parameters[["scale"]] <-
          data %>%
          select(all_of(cues)) %>%
          summarise_all(list(sd = sd))
      }
    }
  }

  if (return.transformed.data) {
    if (!is.null(transform.parameters[["pca"]])) {
      data %<>%
        predict(transform.parameters[["pca"]], .) %>%
        as_tibble(.name_repair = "minimal")
    } else {
      if (!is.null(transform.parameters[["center"]])) {
        data %<>%
          ungroup() %>%
          select(all_of(cues)) %>%
          { . - (data %>%
                   { if (length(groups) > 0) left_join(., transform.parameters[["center"]], by = groups) else cross_join(., transform.parameters[["center"]]) } %>%
                   ungroup() %>%
                   select(all_of(paste0(cues, "_mean"))))
          }
      }

      if (!is.null(transform.parameters[["scale"]])) {
        data %<>%
          ungroup() %>%
          select(all_of(cues)) %>%
          { . / (data %>%
                   { if (length(groups) > 0) left_join(., transform.parameters[["scale"]], by = groups) else cross_join(., transform.parameters[["scale"]]) } %>%
                   ungroup() %>%
                   select(all_of(paste0(cues, "_sd"))))
          }
      }
    }
  }

  transform.function = if (!return.transform.function) NULL else {
    function(data) {
      cues <- cues
      center <- center
      scale <- scale
      pca <- pca
      attach <- attach

      transform_cues(data, cues, center = center, scale = scale, pca = pca,
                     attach = attach,
                     transform.parameters = transform.parameters,
                     return.transformed.data = T, return.transform.parameters = F,
                     return.transform.function = F, return.untransform.function = F)

    }
  }

  untransform.function = if (!return.untransform.function) NULL else {
    untransform_cues(data, cues, uncenter = center, unscale = scale, unpca = pca,
                     attach = attach,
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
untransform_cues <- function(
    data, cues,
    uncenter = NULL, unscale = NULL, unpca = NULL,
    attach = T,
    transform.parameters = NULL,
    return.untransformed.data = T, return.untransform.function = F
) {
  assert_that(is.data.frame(data) | is_tibble(data))
  assert_that(!is.null(transform.parameters) & is.list(transform.parameters),
              msg = "Must provide transform parameters.")
  old_data <- data
  groups <- if (length(groups(data)) == 0) character() else groups(data) %>% as.character()

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
        select(all_of(cues)) %>%
        { . * (data %>%
                 { if (length(groups) > 0) left_join(., transform.parameters[["scale"]], by = groups) else cross_join(., transform.parameters[["scale"]]) } %>%
                 ungroup() %>%
                 select(all_of(paste0(cues, "_sd"))))
        }
    }

    if (uncenter) {
      data %<>%
        ungroup() %>%
        select(all_of(cues)) %>%
        { . + (data %>%
                 { if (length(groups) > 0) left_join(., transform.parameters[["center"]], by = groups) else cross_join(., transform.parameters[["center"]]) } %>%
                 ungroup() %>%
                 select(all_of(paste0(cues, "_mean"))))
        }
    }
  }

  untransform.function = if (!return.untransform.function) NULL else {
    function(data) {
      cues <- cues
      uncenter <- uncenter
      unscale <- unscale
      unpca <- unpca
      attach <- attach

      untransform_cues(data, cues, uncenter = uncenter, unscale = unscale, unpca = unpca,
                       attach = attach,
                       transform.parameters = transform.parameters,
                       return.untransformed.data = T, return.untransform.function = F)

    }
  }

  if (attach) {
    data %<>%
      cbind(
        old_data %>%
          select(-intersect(names(old_data), names(data))),
        .)
  }

  if (return.untransformed.data & !return.untransform.function) return(data) else
    if (!return.untransformed.data & return.untransform.function) return(untransform.function) else
      return(list(data = data, untransform.function = untransform.function))
}

#' Transform and untransform category means and covariance matrices of a model
#'
#' Applies the transformation specified in \code{transform} object (e.g., centering, scaling, PCA) to the category
#' mean(s) and/or covariance(s) in the model.
#'
#' @param model A model with columns that specify category means (mu or m) and covariance information (Sigma or S).
#' @param m A single category mean or alike.
#' @param S A single covariance matrix or alike.
#' @param transform A transform or transform parameter object of the type returned by \code{\link{transform_cues}}.
#'
#' @return A model, category mean, or covariance matrix of the same type as the input.
#'
#' @keywords TBD
#' @examples
#' TBD
#' @rdname transform_model
#' @export
transform_model <- function(model, transform) {
  if (is.MVG(model) | is.MVG_ideal_observer(model)) {
    m = "mu"
    S = "Sigma"
  } else if (is.NIW_belief(model) | is.NIW_ideal_adaptor(model) | is.NIW_ideal_adaptor_stanfit(model) | is.NIW_ideal_adaptor_MCMC(model)) {
    m = "m"
    S = "S"
  } else {
    stop("Model type not recognized.")
  }

  model %>%
    mutate(
      !! sym(m) := map(!! sym(m), ~ transform_category_mean(.x, transform)),
      !! sym(S) := map(!! sym(S), ~ transform_category_cov(.x, transform)))
}

#' @rdname transform_model
#' @export
untransform_model <- function(model, transform) {
  if (is.MVG(model) | is.MVG_ideal_observer(model)) {
    m = "mu"
    S = "Sigma"
  } else if (is.NIW_belief(model) | is.NIW_ideal_adaptor(model) | is.NIW_ideal_adaptor_stanfit(model) | is.NIW_ideal_adaptor_MCMC(model)) {
    m = "m"
    S = "S"
  } else {
    stop("Model type not recognized.")
  }

  model %>%
    mutate(
      !! sym(m) := map(!! sym(m), ~ untransform_category_mean(.x, transform)),
      !! sym(S) := map(!! sym(S), ~ untransform_category_cov(.x, transform)))
}

#' @rdname transform_model
#' @export
transform_category_mean <- function(m, transform) {
  if (!is.null(transform[["transform.parameters"]])) transform <- transform[["transform.parameters"]]
  if (!is.null(transform$center)) {
    mean <- transform$center %>% as.numeric()
    m <- m - mean
  }

  if (!is.null(transform$scale)) {
    taus <- transform$scale %>% as.numeric()
    m <- m / taus
  }

  return(m)
}

#' @rdname transform_model
#' @export
untransform_category_mean <- function(m, transform) {
  if (!is.null(transform[["transform.parameters"]])) transform <- transform[["transform.parameters"]]
  if (!is.null(transform$scale)) {
    taus <- transform$scale %>% as.numeric()
    m <- m * taus
  }

  if (!is.null(transform$center)) {
    mean <- transform$center %>% as.numeric()
    m <- m + mean
  }

  return(m)
}

#' @rdname transform_model
#' @export
transform_category_cov <- function(S, transform) {
  if (!is.null(transform[["transform.parameters"]])) transform <- transform[["transform.parameters"]]
  if (!is.null(transform$scale)) {
    taus <- transform$scale %>% as.numeric()
    COVinv <- diag(taus) %>% solve()
    S <- COVinv %*% S %*% COVinv
  }
  return(S)
}

#' @rdname transform_model
#' @export
untransform_category_cov <- function(S, transform) {
  if (!is.null(transform[["transform.parameters"]])) transform <- transform[["transform.parameters"]]
  if (!is.null(transform$scale)) {
    taus <- transform$scale %>% as.numeric()
    COV <- diag(taus)

    S <- COV %*% S %*% COV
  }
  return(S)
}




