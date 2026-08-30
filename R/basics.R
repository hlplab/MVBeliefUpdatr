NULL

#' Infer default noise treatment from perceptual noise covariance
#'
#' @param Sigma_noise Perceptual noise matrix, scalar, or list of noise matrices.
#' @return A character scalar: `"marginalize"` if noise is present, `"no_noise"` otherwise.
#' @keywords internal
#' @noRd
.infer_noise_treatment <- function(Sigma_noise) {
  if (is.null(Sigma_noise)) {
    return("no_noise")
  }
  if (is.list(Sigma_noise)) {
    if (length(Sigma_noise) == 0L || any(vapply(Sigma_noise, is.null, logical(1)))) {
      return("no_noise")
    }
  }
  "marginalize"
}


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
#'
#' @export
ss <- function(x, center = TRUE) {
  if (center) {
    # Benchmarked to be more efficient than x - 1 %*% t(1) %*% x
    xm <- .colMeans(x)
    xm <- matrix(xm, nrow = nrow(x), ncol = length(xm), byrow = T)
    x <- x - xm
  }

  # Benchmarked to be more efficient than t(x) %*% x
  k <- dim(x)[2]
  m <- matrix(ncol = k, nrow = k)
  for (i in 1:k) {
    m[i, i] <- sum(x[, i]**2)
    if (i < k) {
      for (j in (i + 1):k) {
        m[j, i] <- sum(x[, i] * x[, j])
        m[i, j] <- m[j, i]
      }
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
#'
#' @rdname uss2css
#' @importFrom LaplacesDemon is.positive.semidefinite
#' @export
uss2css <- function(uss, n, mean) {
  .assert_numeric(uss, msg = "uss must be a numeric matrix.")
  if (.is_scalar(uss)) uss <- matrix(uss, nrow = 1, ncol = 1)
  .assert_true(is.positive.semidefinite(uss), msg = "uss must be positive definite.")
  .assert_that(length(mean) == dim(uss)[[1]],
    msg = "uss and mean are not of compatible dimensions."
  )

  css <- uss - n * mean %*% t(mean)
  return(css)
}


#' @rdname uss2css
#' @export
uss2cov <- function(uss, n, mean) {
  css <- uss2css(uss, n, mean)
  cov <- css2cov(css, n)
  return(css)
}

#' @rdname uss2css
#' @export
css2uss <- function(css, n, mean) {
  .assert_numeric(css, msg = "css must be a numeric matrix.")
  if (.is_scalar(css)) css <- matrix(css, nrow = 1, ncol = 1)
  .assert_true(is.positive.semidefinite(css), msg = "css must be positive definite.")
  .assert_that(length(mean) == dim(css)[[1]],
    msg = "uss and mean are not of compatible dimensions."
  )

  xm <- matrix(mean, nrow = n, ncol = length(mean), byrow = T)
  uss <- css + ss(xm, center = F)
  return(uss)
}

#' @rdname uss2css
#' @export
css2cov <- function(css, n) {
  .assert_numeric(css, msg = "css must be a numeric matrix.")
  if (.is_scalar(css)) css <- matrix(css, nrow = 1, ncol = 1)
  .assert_true(is.positive.semidefinite(css), msg = "css must be positive definite.")

  return(css / (n - 1))
}

#' @rdname uss2css
#' @export
cov2css <- function(cov, n) {
  .assert_numeric(cov, msg = "cov must be a numeric matrix.")
  if (.is_scalar(cov)) cov <- matrix(cov, nrow = 1, ncol = 1)
  .assert_true(is.positive.semidefinite(cov), msg = "cov must be positive definite.")

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
#' @seealso [stats::cov2cor], \code{\link{cov2tau}}
#' @references \url{https://stats.stackexchange.com/questions/62850/obtaining-covariance-matrix-from-correlation-matrix}
#' @keywords TBD
#'
#' @export
cor2cov <- function(omega, tau) {
  .assert_that(is.matrix(omega))
  .assert_that(is.numeric(tau))
  .assert_that(nrow(omega) == ncol(omega),
    msg = "omega must be a square matrix."
  )
  .assert_that(length(tau) == nrow(omega),
    msg = "tau must have as many elements as omega has rows (and columns)."
  )

  outer(tau, tau) * omega
}

#' Get vector of standard deviations from covariance matrix
#'
#' Get vector of standard deviations from covariance matrix.
#'
#' @param v A covariance matrix.
#'
#' @return A numeric vector.
#'
#' @seealso [stats::cov2cor], \code{\link{cor2cov}}
#' @keywords TBD
#'
#' @export
cov2tau <- function(v) {
  .assert_that(is.matrix(v))
  .assert_that(nrow(v) == ncol(v),
    msg = "omega must be a square matrix."
  )
  sqrt(diag(v))
}



#' Combine a number of columns into a new vector column
#'
#' Combine a number of columns into a new column in which each cell is the vector of values from the original columns.
#'
#' @param data A tibble or data.frame.
#' @param cols A character vector of variable names to combine.
#' @param vector_col Name of the new vector-valued column.
#' @param .keep A tidyselect option passed to dplyr::mutate.
#'
#' @return Same as \code{data}.
#'
#' @keywords TBD
#' @export
make_vector_column <- function(data, cols, vector_col, .keep = "all") {
  # CHECK: expand to also handle quo input. (each instance of calls then needs to change)
  # then make_NIW_prior_from...  use this function
  data %<>%
    mutate(
      !!sym(vector_col) := pmap(
        .l = list(!!!syms(cols)),
        .f = function(...) {
          x <- c(...)
          names(x) <- cols
          return(x)
        }
      ),
      .keep = .keep
    )

  return(data)
}


#' Get sum of uncentered squares
#'
#' Get sum of uncentered squares. This quantity is a sufficient statistic for, for example, multivariate Gaussian
#' belief-updating under an Normal-Inverse-Wishart prior.
#'
#' @param data A `tibble`, `data.frame`, or `matrix`. If data is a `tibble` or `data.frame`, the columns for
#' specified variables are extracted and (together) converted into a matrix with as many colums as there are
#' variables. If `NULL`, `NA` is returned.
#' @param variables Only required if data is not already a `matrix`.
#'
#' @return A matrix.
#'
#' @keywords TBD
#' @rdname get_sum_of_squares_from_df
#' @export
get_sum_of_squares_from_df <- function(data, variables = NULL, center = T, verbose = F) {
  if (is.null(data)) {
    return(NA)
  }

  .assert_that(is_tibble(data) | is.data.frame(data) | is.matrix(data))
  if (is_tibble(data) | is.data.frame(data)) {
    .assert_that(all(variables %in% names(data)),
      msg = paste("Variable column(s)", variables[which(variables %nin% names(data))], "not found in data.")
    )
  }

  data.matrix <- if (is_tibble(data) | is.data.frame(data)) {
    # Assume that the variables are to be combined into a data.matrix
    data %>%
      mutate(across(c(!!!syms(variables)), unlist)) %>%
      select(all_of(variables)) %>%
      as.matrix()
  } else {
    data
  }

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
#' Get sufficient statistics from data. Calculates functions for the specified
#' cues for any combination of groups (optional) and categories, and returns
#' them as a data.frame. Rows with missing values for cues will be ignored in the
#' calculation of the sufficient statistics.
#'
#' @param data A `data.frame` with the data. Each row should be an observation
#'   of a category, containing the category label, the cue values of the
#'   observation, and optionally grouping variables.
#' @param cues Names of columns with cue values.
#' @param category Name of column that contains the category label for the
#'   exposure data. (default: "category")
#' @param group Name of column(s) that contains information about which
#'   observations form a group. This could be individual subjects or conditions
#'   in an experiment. (default: `NULL`)
#' @param categories,groups Character vector of categories/groups to be
#'   summarized. If `NULL`, all categories/groups will be included.
#'   (default: `NULL`)
#' @param verbose Logical; if `TRUE`, produce verbose output. (default: `FALSE`)
#' @param ... Additional arguments (currently unused).
#'
#' @return A `data.frame` of sufficient statistics for each combination of
#'   category and group. This includes the count (`x_N`), mean vector
#'   (`x_mean`), uncentered and centered sums-of-squares (`x_uss`, `x_css`), and
#'   the covariance matrix (`x_cov`).
#'
#' @keywords TBD
#' @export
get_sufficient_category_statistics <- function(
  data,
  cues,
  category = "category",
  group = NULL,
  categories = NULL,
  groups = NULL,
  verbose = FALSE,
  ...
) {
  data <- as.data.frame(data)

  if (!is.null(categories)) {
    data <- data[data[[category]] %in% categories, , drop = FALSE]
  }
  if (!is.null(groups) && !is.null(group)) {
    if (length(group) == 1L) {
      data <- data[data[[group]] %in% groups, , drop = FALSE]
    } else {
      match_rows <- apply(
        data[, group, drop = FALSE],
        1L,
        function(r) any(r %in% groups)
      )
      data <- data[match_rows, , drop = FALSE]
    }
  }

  data <- droplevels(data)

  if (length(cues) >= 1L) {
    complete_idx <- stats::complete.cases(data[, cues, drop = FALSE])
    data <- data[complete_idx, , drop = FALSE]
  }

  group_cols <- c(category, group)

  if (nrow(data) == 0L) {
    res <- data[, group_cols, drop = FALSE]
    res$x_N <- integer(0)
    res$x_mean <- list()
    res$x_uss <- list()
    res$x_css <- list()
    res$x_cov <- list()
    return(res)
  }

  by_list <- lapply(group_cols, function(col) data[[col]])
  names(by_list) <- group_cols

  idx_split <- split(seq_len(nrow(data)), by_list, drop = TRUE, lex.order = TRUE)
  first_indices <- vapply(idx_split, function(idx) idx[1L], integer(1L))
  res <- data[first_indices, group_cols, drop = FALSE]
  rownames(res) <- NULL

  n_groups <- length(idx_split)
  x_N <- integer(n_groups)
  x_mean <- vector("list", n_groups)
  x_uss <- vector("list", n_groups)
  x_css <- vector("list", n_groups)
  x_cov <- vector("list", n_groups)

  cue_mat <- as.matrix(data[, cues, drop = FALSE])
  if (!is.numeric(cue_mat)) {
    storage.mode(cue_mat) <- "double"
  }

  for (i in seq_len(n_groups)) {
    idx <- idx_split[[i]]
    sub_mat <- cue_mat[idx, , drop = FALSE]
    n_obs <- length(idx)

    x_N[i] <- n_obs
    x_mean[[i]] <- colMeans(sub_mat)
    x_uss[[i]] <- ss(sub_mat, center = FALSE)
    x_css[[i]] <- ss(sub_mat, center = TRUE)
    cov_mat <- stats::cov(sub_mat)
    if (length(cues) == 1L) {
      cov_mat <- matrix(cov_mat, 1L, 1L, dimnames = list(cues, cues))
    }
    x_cov[[i]] <- cov_mat
  }

  res$x_N <- x_N
  res$x_mean <- x_mean
  res$x_uss <- x_uss
  res$x_css <- x_css
  res$x_cov <- x_cov

  return(res)
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
#'
#' @rdname transform_cues
#' @importFrom tidyselect all_of
#' @importFrom dplyr across cross_join
#' @importFrom rlang syms
#' @export
transform_cues <- function(
  data, cues,
  center = if (pca) T else F, scale = F, pca = F,
  attach = T,
  transform.parameters = NULL,
  return.transformed.data = T, return.transform.parameters = F,
  return.transform.function = F, return.untransform.function = F
) {
  .assert_data_frame_like(data)
  .assert_that(is.character(cues))
  .assert_that(all(cues %in% colnames(data)))
  .assert_that(all(
    is.logical(center), is.logical(scale), is.logical(pca), is.logical(attach),
    is.logical(return.transformed.data), is.logical(return.transform.parameters),
    is.logical(return.transform.function), is.logical(return.untransform.function)
  ))
  if (pca) center <- T
  .assert_optional_list(transform.parameters)
  old_data <- data
  groups <- if (length(dplyr::groups(data)) == 0) character() else dplyr::groups(data) %>% as.character()

  if (is.null(transform.parameters)) {
    transform.parameters <- list()
    transform.parameters[["cue.labels"]] <- cues

    if (pca) {
      transform.parameters[["pca"]] <-
        data %>%
        select(all_of(cues)) %>%
        # If centering and scaling is requested in addition to PCA, do it as part of
        # the PCA transformation
        prcomp(center = center, scale. = scale, retx = F)
    } else {
      if (center) {
        transform.parameters[["center"]] <-
          data %>%
          select(all_of(cues)) %>%
          summarise(across(everything(), list(mean = mean)))
      }

      if (scale) {
        transform.parameters[["scale"]] <-
          data %>%
          select(all_of(cues)) %>%
          summarise(across(everything(), list(sd = sd)))
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
          {
            . - (data %>%
              {
                if (length(groups) > 0) left_join(., transform.parameters[["center"]], by = groups) else cross_join(., transform.parameters[["center"]])
              } %>%
              ungroup() %>%
              select(all_of(paste0(cues, "_mean"))))
          }
      }

      if (!is.null(transform.parameters[["scale"]])) {
        data %<>%
          ungroup() %>%
          select(all_of(cues)) %>%
          {
            . / (data %>%
              {
                if (length(groups) > 0) left_join(., transform.parameters[["scale"]], by = groups) else cross_join(., transform.parameters[["scale"]])
              } %>%
              ungroup() %>%
              select(all_of(paste0(cues, "_sd"))))
          }
      }
    }
  }

  transform.function <- if (!return.transform.function) {
    NULL
  } else {
    function(data) {
      cues <- cues
      center <- center
      scale <- scale
      pca <- pca
      attach <- attach

      transform_cues(data, cues,
        center = center, scale = scale, pca = pca,
        attach = attach,
        transform.parameters = transform.parameters,
        return.transformed.data = T, return.transform.parameters = F,
        return.transform.function = F, return.untransform.function = F
      )
    }
  }

  untransform.function <- if (!return.untransform.function) {
    NULL
  } else {
    untransform_cues(
      data = data, cues = if (pca) colnames(transform.parameters[["pca"]]$rotation) else cues,
      uncenter = center, unscale = scale, unpca = pca,
      attach = attach,
      transform.parameters = transform.parameters,
      return.untransformed.data = F, return.untransform.function = T
    )
  }

  if (!identical(groups, character())) data %<>% group_by(!!sym(groups))
  if (attach) {
    data <-
      cbind(
        old_data %>%
          select(-intersect(names(old_data), names(data))),
        data
      )
  }

  if (!any(return.transformed.data, return.transform.parameters, return.transform.function, return.untransform.function)) {
    message("No return was requested.")
  } else if (return.transformed.data & !any(return.transform.parameters, return.transform.function, return.untransform.function)) {
    return(data)
  } else if (return.transform.parameters & !any(return.transformed.data, return.transform.function, return.untransform.function)) {
    return(transform.parameters)
  } else if (return.transform.function & !any(return.transformed.data, return.transform.parameters, return.untransform.function)) {
    return(transform.function)
  } else if (return.untransform.function & !any(return.transformed.data, return.transform.parameters, return.transform.function)) {
    return(untransform.function)
  } else {
    return(
      list(
        data = if (return.transformed.data) data else NULL,
        transform.parameters = if (return.transform.parameters) transform.parameters else NULL,
        transform.function = if (return.transform.function) transform.function else NULL,
        untransform.function = if (return.untransform.function) untransform.function else NULL
      )
    )
  }
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
  .assert_data_frame_like(data)
  .assert_that(!is.null(transform.parameters) & is.list(transform.parameters),
    msg = "Must provide transform parameters."
  )
  old_data <- data
  groups <- if (length(dplyr::groups(data)) == 0) character() else dplyr::groups(data) %>% as.character()

  # By default untransform all transformations available in transform object
  if (is.null(unpca)) unpca <- !is.null(transform.parameters[["pca"]])
  if (is.null(uncenter)) uncenter <- !is.null(transform.parameters[["center"]])
  if (is.null(unscale)) unscale <- !is.null(transform.parameters[["scale"]])

  if (return.untransformed.data) {
    if (unpca) {
      # https://stackoverflow.com/questions/29783790/how-to-reverse-pca-in-prcomp-to-get-original-data
      newcues <-
        data %>%
        select(all_of(cues)) %>%
        as.matrix() %>%
        {
          . %*% t(transform.parameters[["pca"]]$rotation)
        } %>%
        {
          if (unscale & uncenter) {
            t(t(.) * transform.parameters[["pca"]]$scale + transform.parameters[["pca"]]$center)
          } else if (unscale) t(t(.) * transform.parameters[["pca"]]$scale) else if (uncenter) t(t(.) + transform.parameters[["pca"]]$center)
        }

      data %<>%
        select(-all_of(cues)) %>%
        cbind(newcues)
    } else {
      if (unscale) {
        data %<>%
          ungroup() %>%
          select(all_of(cues)) %>%
          {
            . * (data %>%
              {
                if (length(groups) > 0) left_join(., transform.parameters[["scale"]], by = groups) else cross_join(., transform.parameters[["scale"]])
              } %>%
              ungroup() %>%
              select(all_of(paste0(cues, "_sd"))))
          }
      }

      if (uncenter) {
        data %<>%
          ungroup() %>%
          select(all_of(cues)) %>%
          {
            . + (data %>%
              {
                if (length(groups) > 0) left_join(., transform.parameters[["center"]], by = groups) else cross_join(., transform.parameters[["center"]])
              } %>%
              ungroup() %>%
              select(all_of(paste0(cues, "_mean"))))
          }
      }
    }

    if (attach) {
      data %<>%
        cbind(
          old_data %>%
            select(-intersect(names(old_data), names(data))),
          .
        )
    }
  }

  untransform.function <- if (!return.untransform.function) {
    NULL
  } else {
    function(data) {
      cues <- cues
      uncenter <- uncenter
      unscale <- unscale
      unpca <- unpca
      attach <- attach

      untransform_cues(data, cues,
        uncenter = uncenter, unscale = unscale, unpca = unpca,
        attach = attach,
        transform.parameters = transform.parameters,
        return.untransformed.data = T, return.untransform.function = F
      )
    }
  }

  if (!any(return.untransformed.data, return.untransform.function)) {
    message("No return was requested.")
  } else if (return.untransformed.data & !return.untransform.function) {
    return(data)
  } else if (!return.untransformed.data & return.untransform.function) {
    return(untransform.function)
  } else {
    return(list(data = data, untransform.function = untransform.function))
  }
}


#' Get affine transformation
#'
#' Get affine transformation for a set of `cues` in `data`. Returns a linear transformation of form
#' `f(x) = SCALE * (x + shift)`, the inverse of that transformation, and
#' all relevant transformation parameters. The transformation function and its inverse take data sets as input,
#' and return a new `data.frame` with only the (inverse) transformed `cue` columns. The transformation parameters
#' can be helpful when one wants to apply the affine transform to transform, for instance, covariance matrices.
#'
#' @param data A `tibble` or `data.frame`.
#' @param cues A character vector of column names in `data` that should be transformed.
#' @param type The type of transformation to apply. Can be one of "identity", "center", "standardize", "PCA whiten",
#'   or "ZCA whiten". Except for the identity transform, all transforms center. Standardize additionally divides by
#'   the standard deviations. PCA whitening rotates the data into the space of the principal components, followed
#'   by scaling along each principal axis to achieve unit variance. ZCA whitening also achieves unit variance along
#'   all dimensions, and---like PCA whitening---decorrelates the data but it aims to maintain the original orientation
#'   of the data as close as possible. If `cue` is a single cue, whitening reduces to standardization.
#'   (default : "identity")
#'
#' @return A list with the following elements:
#'  * `type`: A character vector of length 1, containing the type of transformation.
#'  * `transform.parameters`: A list with the following elements:
#'     * `cue.labels`: A character vector of cue labels for ease of recovery.
#'     * `shift` is a vector of `length(cues)`
#'     * `SCALE` is a `length(cues)` x `length(cues)` invertible square matrix. For types "identity" and "center",
#'        this is an identity matrix. For type "standardize", it is the identity matrix multiplied by the inverse
#'        of the vector of standard deviations of the cues.
#'     * `INV_SCALE` is the inverse of the `SCALE` matrix.
#'  * `transform.function`: A function f(data, return_type). The function applies the affine transformation to the
#'     columns of `data` specified in `cues`. The `return_type` argument determines what is returned. Can be one of
#'     "replace", "add", or "cues only". If "replace", the input data will be returned in full but with the cue
#'     columns replaced by the transformed cues. If "add", the input data will be returned with additional columns
#'     for the cues (labeled "*_transformed"). If "cues only", a new data frame with only the transformed cues will
#'     be returned. (default: "replace")
#'  * `untransform.function`: Same as the `transform.function` but return the inverse transformed data.
#'
#' @export
get_affine_transform <- function(
  data,
  cues,
  type = c("identity", "center", "scale", "PCA whiten", "ZCA whiten")[1]
) {
  .assert_that(is.character(cues))
  .assert_data_contains_cols(data, cues)

  # groups <- if (length(groups(data)) == 0) character() else groups(data) %>% as.character()

  transform.parameters <- list()
  transform.parameters[["cue.labels"]] <- cues

  n.cues <- length(cues)
  data <- as.matrix(data[, cues, drop = FALSE])

  if (type == "identity") {
    transform.parameters[["shift"]] <- rep(0, n.cues)
  } else {
    data <- scale(data, scale = F)
    transform.parameters[["shift"]] <- -(attr(data, "scaled:center"))
  }

  if (n.cues == 1) {
    if (type %in% c("identity", "center")) {
      transform.parameters[["SCALE"]] <- 1
    } else if (type %in% c("standardize", "PCA whiten", "ZCA whiten")) {
      transform.parameters[["SCALE"]] <- 1 / sd(data[, 1])
    } else {
      .stop("Unknown type.")
    }
  } else {
    if (type %in% c("identity", "center")) {
      transform.parameters[["SCALE"]] <- diag(n.cues)
    } else if (type == "standardize") {
      transform.parameters[["SCALE"]] <- diag(1 / apply(data, 2, sd))
    } else if (type %in% c("PCA whiten", "ZCA whiten")) {
      eig <- eigen(.cov(data))
      U <- eig$vectors
      D_inv_half <- diag(1 / sqrt(eig$values))

      if (type == "PCA whiten") {
        transform.parameters[["SCALE"]] <- D_inv_half %*% t(U)
      } else if (type == "ZCA whiten") {
        transform.parameters[["SCALE"]] <- U %*% D_inv_half %*% t(U)
      } else {
        .stop("Unknown whitening type.")
      }
    } else {
      .stop("Unknown type.")
    }
  }

  transform.parameters[["INV_SCALE"]] <- if (n.cues == 1) 1 / transform.parameters[["SCALE"]] else solve(transform.parameters[["SCALE"]])

  # Defining function that handles return of data that gets repeated below in the
  # different transform and untransform functions
  get_return_data <- function(data, newdata, return_type) {
    newdata <- as.data.frame(newdata)

    if (return_type == "replace") {
      data[, cues] <- newdata
    } else if (return_type == "add") {
      names(newdata) <- paste0(cues, "_transformed")
      data <- cbind(data, newdata)
    } else if (return_type %in% c("cues only", "only cues")) {
      names(newdata) <- cues
      data <- newdata
    } else {
      .stop("Unknown return type.")
    }

    return(data)
  }

  if (n.cues == 1) {
    transform.function <-
      function(data, return_type = c("replace", "add", "cues only")[1]) {
        get_return_data <- get_return_data
        cues <- cues
        shift <- transform.parameters[["shift"]]
        SCALE <- transform.parameters[["SCALE"]]

        newdata <- data[, cues]
        newdata <- (newdata + shift) * SCALE
        newdata <- get_return_data(data, newdata, return_type)

        return(newdata)
      }

    untransform.function <-
      function(data, return_type = c("replace", "add", "cues only")[1]) {
        get_return_data <- get_return_data
        cues <- cues
        shift <- transform.parameters[["shift"]]
        INV_SCALE <- transform.parameters[["INV_SCALE"]]

        newdata <- as.matrix(data[, cues])
        newdata <- (newdata * INV_SCALE) - shift
        newdata <- get_return_data(data, newdata, return_type)

        return(newdata)
      }
    # Handle cases with multiple cues
  } else {
    transform.function <-
      function(data, return_type = c("replace", "add", "cues only")[1]) {
        get_return_data <- get_return_data
        cues <- cues
        shift <- transform.parameters[["shift"]]
        SCALE <- transform.parameters[["SCALE"]]

        newdata <- as.matrix(data[, cues])
        newdata <- sweep(newdata, 2, -shift) %*% t(SCALE)
        newdata <- get_return_data(data, newdata, return_type)

        return(newdata)
      }

    untransform.function <-
      function(data, return_type = c("replace", "add", "cues only")[1]) {
        get_return_data <- get_return_data
        cues <- cues
        shift <- transform.parameters[["shift"]]
        INV_SCALE <- transform.parameters[["INV_SCALE"]]

        newdata <- as.matrix(data[, cues])
        newdata <- sweep(newdata %*% t(INV_SCALE), 2, +shift)
        newdata <- get_return_data(data, newdata, return_type)

        return(newdata)
      }
  }

  return(.nlist(type, transform.parameters, transform.function, untransform.function))
}


#' Transform and untransform category means and covariance matrices of a model
#'
#' Applies an affine (linear) transformation specified in a \code{transform} object
#' (e.g., centering, scaling, PCA, ZCA whitening from \code{\link{transform_cues}} or \code{\link{get_affine_transform}})
#' to the category mean(s) and covariance(s) within an S7 cognitive model.
#'
#' Only linear/affine transformations are supported because Gaussian and Normal-Inverse-Wishart
#' distributions transform analytically under affine maps \eqn{y = A x + b}:
#' \deqn{\mu' = A \mu + b, \quad \Sigma' = A \Sigma A^T}
#' Models containing \code{\link{MVG_CategoryRepresentation}}, \code{\link{UVG_CategoryRepresentation}},
#' \code{\link{NIW_CategoryRepresentation}}, \code{\link{NIX_CategoryRepresentation}}, or
#' \code{\link{MNIX_CategoryRepresentation}} are supported. Exemplar representations are transformed
#' by transforming their exemplar cues directly.
#'
#' Note that \code{(un)transform_category_cov} can also be used for centered sums-of-squares matrices, but not for
#' uncentered sums-of-squares matrices. The latter require additional adjustments that might be implemented in the
#' future (see https://chatgpt.com/c/683b3ed6-4924-800c-8f22-e61680f3360f).
#'
#'
#' @param model An \code{\link{MVBU_CognitiveModel}} object.
#' @param m A single category mean vector (or scalar for univariate representations).
#' @param S A single covariance matrix (or scalar variance for univariate representations).
#' @param transform A transform or transform parameter object returned by \code{\link{transform_cues}} or \code{\link{get_affine_transform}}.
#'
#' @return A model, category mean, or covariance matrix of the same class as the input.
#'
#' @keywords TBD
#' @rdname transform_model
#' @export
transform_model <- function(model, transform) {
  if (!S7::S7_inherits(model, MVBU_CognitiveModel)) {
    stop("model must be an MVBU_CognitiveModel object.")
  }

  reps <- model@category_template@representations
  for (cat_name in names(reps)) {
    r <- reps[[cat_name]]
    if (S7::S7_inherits(r, MVG_CategoryRepresentation)) {
      r@mu <- transform_category_mean(r@mu, transform)
      r@Sigma <- transform_category_cov(r@Sigma, transform)
    } else if (S7::S7_inherits(r, NIW_CategoryRepresentation)) {
      r@m <- transform_category_mean(r@m, transform)
      r@S <- transform_category_cov(r@S, transform)
    } else if (S7::S7_inherits(r, UVG_CategoryRepresentation)) {
      r@mu <- transform_category_mean(r@mu, transform)
      r@sigma2 <- as.numeric(transform_category_cov(matrix(r@sigma2, 1, 1), transform))
    } else if (S7::S7_inherits(r, NIX_CategoryRepresentation)) {
      r@m <- transform_category_mean(r@m, transform)
      r@S <- as.numeric(transform_category_cov(matrix(r@S, 1, 1), transform))
    } else if (S7::S7_inherits(r, MNIX_CategoryRepresentation)) {
      r@m <- transform_category_mean(r@m, transform)
      r@S <- transform_category_cov(r@S, transform)
    }
    reps[[cat_name]] <- r
  }
  model@category_template@representations <- reps
  return(model)
}

#' @rdname transform_model
#' @export
untransform_model <- function(model, transform) {
  if (!S7::S7_inherits(model, MVBU_CognitiveModel)) {
    stop("model must be an MVBU_CognitiveModel object.")
  }

  reps <- model@category_template@representations
  for (cat_name in names(reps)) {
    r <- reps[[cat_name]]
    if (S7::S7_inherits(r, MVG_CategoryRepresentation)) {
      r@mu <- untransform_category_mean(r@mu, transform)
      r@Sigma <- untransform_category_cov(r@Sigma, transform)
    } else if (S7::S7_inherits(r, NIW_CategoryRepresentation)) {
      r@m <- untransform_category_mean(r@m, transform)
      r@S <- untransform_category_cov(r@S, transform)
    } else if (S7::S7_inherits(r, UVG_CategoryRepresentation)) {
      r@mu <- untransform_category_mean(r@mu, transform)
      r@sigma2 <- as.numeric(untransform_category_cov(matrix(r@sigma2, 1, 1), transform))
    } else if (S7::S7_inherits(r, NIX_CategoryRepresentation)) {
      r@m <- untransform_category_mean(r@m, transform)
      r@S <- as.numeric(untransform_category_cov(matrix(r@S, 1, 1), transform))
    } else if (S7::S7_inherits(r, MNIX_CategoryRepresentation)) {
      r@m <- untransform_category_mean(r@m, transform)
      r@S <- untransform_category_cov(r@S, transform)
    }
    reps[[cat_name]] <- r
  }
  model@category_template@representations <- reps
  return(model)
}


.get_transform_params <- function(transform) {
  if (S7::S7_inherits(transform, MVBU_TransformInformation)) {
    return(transform@transform.parameters)
  }
  if (is.list(transform) && !is.null(transform[["transform.parameters"]])) {
    return(transform[["transform.parameters"]])
  }
  transform
}

#' @rdname transform_model
#' @export
transform_category_mean <- function(m, transform) {
  transform <- .get_transform_params(transform)
  m_names <- names(m)

  if (!is.null(transform$shift)) {
    m <- m + transform$shift
  }

  if (!is.null(transform$center)) {
    mean <- transform$center %>% as.numeric()
    m <- m - mean
  }

  if (!is.null(transform$scale)) {
    taus <- transform$scale %>% as.numeric()
    m <- m / taus
  }

  # For affine transforms
  if (!is.null(transform$SCALE)) {
    if (is.matrix(transform$SCALE) || (is.array(transform$SCALE) && length(dim(transform$SCALE)) > 1L)) {
      m <- m %*% t(transform$SCALE)
    } else {
      m <- as.numeric(m) * as.numeric(transform$SCALE)
    }
  }

  # Make sure m is a vector and that any potentially available cue names are maintained
  m <- as.vector(m)
  names(m) <- m_names

  return(m)
}

#' @rdname transform_model
#' @export
untransform_category_mean <- function(m, transform) {
  transform <- .get_transform_params(transform)
  m_names <- names(m)

  if (!is.null(transform$scale)) {
    taus <- transform$scale %>% as.numeric()
    m <- m * taus
  }

  if (!is.null(transform$center)) {
    mean <- transform$center %>% as.numeric()
    m <- m + mean
  }

  # For affine transforms
  if (!is.null(transform$SCALE)) {
    m <- m %*% t(transform$INV_SCALE)
  }

  if (!is.null(transform$shift)) {
    m <- m - transform$shift
  }

  # Make sure m is a vector and that any potentially available cue names are maintained
  m <- as.vector(m)
  names(m) <- m_names

  return(m)
}

#' @rdname transform_model
#' @export
transform_category_cov <- function(S, transform) {
  transform <- .get_transform_params(transform)
  S_names <- dimnames(S)

  if (!is.null(transform$scale)) {
    taus <- transform$scale %>% as.numeric()
    COVinv <- diag(taus, nrow = length(taus)) %>% solve()
    S <- COVinv %*% S %*% COVinv
  }

  # For affine transforms
  if (!is.null(transform$SCALE)) {
    if (is.matrix(transform$SCALE) || (is.array(transform$SCALE) && length(dim(transform$SCALE)) > 1L)) {
      S <- transform$SCALE %*% S %*% t(transform$SCALE)
    } else {
      S <- as.numeric(transform$SCALE)^2 * S
    }
  }

  dimnames(S) <- S_names

  return(S)
}

#' @rdname transform_model
#' @export
untransform_category_cov <- function(S, transform) {
  transform <- .get_transform_params(transform)
  S_names <- dimnames(S)

  if (!is.null(transform$scale)) {
    taus <- transform$scale %>% as.numeric()
    COV <- diag(taus, nrow = length(taus))

    S <- COV %*% S %*% COV
  }

  # For affine transforms
  if (!is.null(transform$SCALE)) {
    if (is.matrix(transform$INV_SCALE) || (is.array(transform$INV_SCALE) && length(dim(transform$INV_SCALE)) > 1L)) {
      S <- transform$INV_SCALE %*% S %*% t(transform$INV_SCALE)
    } else {
      S <- as.numeric(transform$INV_SCALE)^2 * S
    }
  }

  dimnames(S) <- S_names

  return(S)
}
