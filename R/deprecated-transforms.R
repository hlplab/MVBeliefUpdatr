#' @include basics.R
#' @importFrom lifecycle deprecate_warn
NULL

# =============================================================================
# deprecated
# =============================================================================

#' Deprecated: Transform and untransform cues by applying or undoing PCA, centering, and/or scaling
#'
#' @description `r lifecycle::badge("deprecated")`
#' `transform_cues()` and `untransform_cues()` are deprecated as of MVBeliefUpdatr 0.1.0
#' and will be removed in 0.2.0. Please use [get_affine_transform()] and [invert_transform()]
#' instead.
#'
#' @param data `tibble` or `data.frame`.
#' @param cues Vector of characters with names of cue variables.
#' @param center,uncenter Should the data be (un)centered? (default: `TRUE` unless `pca = TRUE`)
#' @param scale,unscale Should the data be (un)standardized? (default: `TRUE` unless `pca = TRUE`)
#' @param pca,unpca Should the data be transformed into/back from orthogonal principal components?
#' @param attach Should the transformed cues be attached to \code{data}?
#' @param transform.parameters List of transforms (default: `NULL`).
#' @param return.transformed.data,return.untransformed.data Should (un)transformed data be returned?
#' @param return.transform.parameters Should the list of transforms be returned?
#' @param return.transform.function,return.untransform.function Should a transform function be returned?
#'
#' @return By default a \code{data.frame}.
#'
#' @keywords internal
#' @rdname deprecated-transforms
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
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "transform_cues()",
    with = "get_affine_transform()"
  )
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
          summarise(across(all_of(cues), list(mean = mean)))
      }

      if (scale) {
        transform.parameters[["scale"]] <-
          data %>%
          select(all_of(cues)) %>%
          summarise(across(all_of(cues), list(sd = sd)))
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


#' @rdname deprecated-transforms
#' @export
untransform_cues <- function(
  data, cues,
  uncenter = NULL, unscale = NULL, unpca = NULL,
  attach = T,
  transform.parameters = NULL,
  return.untransformed.data = T, return.untransform.function = F
) {
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "untransform_cues()",
    with = "invert_transform()"
  )
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


#' Deprecated: Transform and untransform category means and covariances
#'
#' @description `r lifecycle::badge("deprecated")`
#' `transform_category_mean()`, `untransform_category_mean()`,
#' `transform_category_cov()`, and `untransform_category_cov()` are deprecated
#' as of MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. Please use
#' [transform_mean()], [untransform_mean()], [transform_cov()], and
#' [untransform_cov()] instead.
#'
#' @param m A mean vector.
#' @param S A covariance matrix.
#' @param transform An \code{\link{MVBU_TransformInformation}} object.
#'
#' @return A transformed or untransformed mean vector or covariance matrix.
#'
#' @keywords internal
#' @rdname deprecated-category-transforms
#' @export
transform_category_mean <- function(m, transform) {
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "transform_category_mean()",
    with = "transform_mean()"
  )
  transform_mean(m = m, transform = transform)
}

#' @rdname deprecated-category-transforms
#' @export
untransform_category_mean <- function(m, transform) {
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "untransform_category_mean()",
    with = "untransform_mean()"
  )
  untransform_mean(m = m, transform = transform)
}

#' @rdname deprecated-category-transforms
#' @export
transform_category_cov <- function(S, transform) {
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "transform_category_cov()",
    with = "transform_cov()"
  )
  transform_cov(S = S, transform = transform)
}

#' @rdname deprecated-category-transforms
#' @export
untransform_category_cov <- function(S, transform) {
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "untransform_category_cov()",
    with = "untransform_cov()"
  )
  untransform_cov(S = S, transform = transform)
}

