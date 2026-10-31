#' @include S7-generics.R S7-transform-information.R
NULL

#' Infer default noise treatment from perceptual noise covariance
#'
#' @param Sigma_noise Perceptual noise matrix, scalar, or list of noise
#'   matrices.
#' @return A character scalar: `"marginalize"` if noise is present, `"no_noise"`
#'   otherwise.
#' @keywords internal
#' @noRd
.infer_noise_treatment <- function(Sigma_noise) {
  if (is.null(Sigma_noise)) {
    return("no_noise")
  }
  if (is.list(Sigma_noise)) {
    if (length(Sigma_noise) == 0L ||
        any(vapply(Sigma_noise, is.null, logical(1)))) {
      return("no_noise")
    }
    if (all(vapply(
      Sigma_noise,
      function(m) is.numeric(m) && (length(m) == 0L || all(m == 0)),
      logical(1)
    ))) {
      return("no_noise")
    }
  }
  if (is.numeric(Sigma_noise) &&
      (length(Sigma_noise) == 0L || all(Sigma_noise == 0))) {
    return("no_noise")
  }
  "marginalize"
}

#' Infer default lapse treatment from lapse rate
#'
#' @param lapse_rate Lapse rate scalar or numeric vector.
#' @return A character scalar: `"marginalize"` if `lapse_rate > 0`,
#'   `"no_lapses"` otherwise.
#' @keywords internal
#' @noRd
.infer_lapse_treatment <- function(lapse_rate) {
  if (is.null(lapse_rate) || length(lapse_rate) == 0L ||
      !is.numeric(lapse_rate)) {
    return("no_lapses")
  }
  if (any(lapse_rate > 0)) {
    return("marginalize")
  }
  "no_lapses"
}


#' Get sum-of-squares matrix
#'
#' Calculate centered (\code{css()}) or uncentered (\code{uss()})
#' sum-of-squares matrix for a matrix of observations.
#'
#' @param x A matrix of observations, with each row being a vector observation.
#'
#' @return A square matrix.
#'
#' @seealso \code{\link{css2cov}}, \code{\link{cov2css}}, \code{\link{uss2css}}
#' @rdname css
#' @export
#' @examples
#' x <- matrix(rnorm(20), nrow = 10, ncol = 2)
#' # Centered sum-of-squares matrix
#' css(x)
#' # Uncentered sum-of-squares matrix
#' uss(x)
css <- function(x) {
  .ss(x, center = TRUE)
}

#' @rdname css
#' @export
uss <- function(x) {
  .ss(x, center = FALSE)
}

#' Internal sum-of-squares computation
#'
#' @param x A matrix of observations.
#' @param center Logical; whether to center observations (default: TRUE).
#' @return A square matrix.
#' @keywords internal
#' @noRd
.ss <- function(x, center = TRUE) {
  if (center) {
    # Benchmarked to be more efficient than x - 1 %*% t(1) %*% x
    xm <- .colMeans(x)
    xm <- matrix(xm, nrow = nrow(x), ncol = length(xm), byrow = TRUE)
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
#' Convert centered sum-of-square matrices (css), uncentered sum-of-square
#' matrices (uss), and covariance matrices (cov) into each other.
#'
#' @param uss,css,cov Uncentered sum-of-square, centered sum-of-square, and
#'   covariance matrix.
#' @param n Number of observations that have gone into the sum-of-square
#'   matrix.
#' @param mean Mean of the observations that have gone into the sum-of-square
#'   matrix.
#'
#' @return A square matrix.
#'
#' @seealso \code{\link{css}}, \code{\link{uss}}, \code{\link[stats]{cov2cor}},
#'   \code{\link{cor2cov}}, \code{\link{cov2tau}}
#' @keywords TBD
#'
#' @rdname uss2css
#' @importFrom LaplacesDemon is.positive.semidefinite
#' @export
#' @examples
#' x <- matrix(rnorm(20), nrow = 10, ncol = 2)
#' n <- nrow(x)
#' m <- colMeans(x)
#' u <- uss(x)
#' c_mat <- css(x)
#'
#' # Convert uncentered sum-of-squares to centered sum-of-squares
#' uss2css(u, n = n, mean = m)
#'
#' # Convert centered sum-of-squares to covariance matrix
#' css2cov(c_mat, n = n)
#'
#' # Convert covariance matrix back to centered sum-of-squares
#' cov2css(cov(x), n = n)
#'
#' # Convert centered sum-of-squares back to uncentered
#' css2uss(c_mat, n = n, mean = m)
#'
#' # Convert uncentered sum-of-squares directly to covariance
#' uss2cov(u, n = n, mean = m)
uss2css <- function(uss, n, mean) {
  .assert_numeric(uss, msg = "uss must be a numeric matrix.")
  if (.is_scalar(uss)) uss <- matrix(uss, nrow = 1, ncol = 1)
  .assert_that(n >= 0, msg = "n must be non-negative.")
  if (n == 0) {
    return(matrix(NA_real_,
      nrow = nrow(uss), ncol = ncol(uss),
      dimnames = dimnames(uss)
    ))
  }
  .assert_true(
    is.positive.semidefinite(uss),
    msg = "uss must be positive definite."
  )
  .assert_that(
    length(mean) == dim(uss)[[1]],
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
  return(cov)
}

#' @rdname uss2css
#' @export
css2uss <- function(css, n, mean) {
  .assert_numeric(css, msg = "css must be a numeric matrix.")
  if (.is_scalar(css)) css <- matrix(css, nrow = 1, ncol = 1)
  .assert_that(n >= 0, msg = "n must be non-negative.")
  if (n == 0) {
    return(matrix(NA_real_,
      nrow = nrow(css), ncol = ncol(css),
      dimnames = dimnames(css)
    ))
  }
  .assert_true(
    is.positive.semidefinite(css),
    msg = "css must be positive definite."
  )
  .assert_that(
    length(mean) == dim(css)[[1]],
    msg = "uss and mean are not of compatible dimensions."
  )

  xm <- matrix(mean, nrow = n, ncol = length(mean), byrow = TRUE)
  uss <- css + uss(xm)
  return(uss)
}

#' @rdname uss2css
#' @export
css2cov <- function(css, n) {
  .assert_numeric(css, msg = "css must be a numeric matrix.")
  if (.is_scalar(css)) css <- matrix(css, nrow = 1, ncol = 1)
  .assert_that(n >= 0, msg = "n must be non-negative.")
  if (n <= 1) {
    return(matrix(NA_real_,
      nrow = nrow(css), ncol = ncol(css),
      dimnames = dimnames(css)
    ))
  }
  .assert_true(
    is.positive.semidefinite(css),
    msg = "css must be positive definite."
  )

  return(css / (n - 1))
}

#' @rdname uss2css
#' @export
cov2css <- function(cov, n) {
  .assert_numeric(cov, msg = "cov must be a numeric matrix.")
  if (.is_scalar(cov)) cov <- matrix(cov, nrow = 1, ncol = 1)
  .assert_that(n >= 0, msg = "n must be non-negative.")
  if (n <= 1) {
    return(matrix(NA_real_,
      nrow = nrow(cov), ncol = ncol(cov),
      dimnames = dimnames(cov)
    ))
  }
  .assert_true(
    is.positive.semidefinite(cov),
    msg = "cov must be positive definite."
  )

  return(cov * (n - 1))
}

#' Get covariance matrix from correlation matrix and standard deviations
#'
#' Convert between a correlation matrix with standard deviations and a
#' covariance matrix.
#'
#' @param omega A correlation matrix.
#' @param tau A vector of standard deviations.
#' @param v A covariance matrix.
#'
#' @return For \code{cor2cov}, a covariance matrix. For \code{cov2tau}, a
#'   numeric vector of standard deviations.
#'
#' @seealso [stats::cov2cor], \code{\link{cor2cov}}, \code{\link{cov2tau}}
#' @references
#'   \url{https://stats.stackexchange.com/questions/62850/obtaining-covariance-matrix-from-correlation-matrix}
#' @keywords TBD
#' @rdname cor2cov
#' @export
#' @examples
#' omega <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' tau <- c(2, 3)
#' # Correlation matrix + standard deviations to covariance matrix
#' Sigma <- cor2cov(omega, tau)
#' Sigma
#'
#' # Extract standard deviations from covariance matrix
#' cov2tau(Sigma)
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

#' @rdname cor2cov
#' @export
cov2tau <- function(v) {
  .assert_that(is.matrix(v))
  .assert_that(nrow(v) == ncol(v),
    msg = "omega must be a square matrix."
  )
  sqrt(diag(v))
}

#' Get sufficient category statistics from a data frame
#'
#' Calculates sufficient statistics for specified cues for all categories and
#' optional groups in a data frame.
#'
#' @param data A \code{data.frame} or \code{tibble}.
#' @param cues Character vector of column names for cues. Required.
#' @param category Character string specifying category column name.
#'   Default: \code{"category"}.
#' @param group Optional character vector of group column name(s).
#'   Default: \code{NULL}.
#' @param stats Character vector of statistics to return. Allowed values are
#'   \code{"n"}, \code{"mean"}, \code{"cov"}, \code{"css"}, and \code{"uss"}.
#'   Default: \code{c("n", "mean", "cov", "css", "uss")}.
#' @param simplify Logical. If \code{TRUE} (default) and the data matches a
#'   single category and group with a single statistic, returns an atomic
#'   scalar, vector, or matrix. If \code{FALSE}, always returns a \code{tibble}.
#' @param ... Additional arguments (e.g. \code{model_family}).
#' @return If \code{simplify = TRUE} and a single statistic is selected for a
#'   single category and group, the atomic numeric value, vector, or matrix.
#'   Otherwise, a \code{tibble} with the grouping columns and requested
#' @seealso \code{\link{get_exposure_category_statistics}}
#' @export
#' @examples
#' # Example with synthetic phonetic data
#' df <- data.frame(
#'   category = rep(c("b", "p"), each = 20),
#'   vot = c(rnorm(20, 0, 10), rnorm(20, 50, 10)),
#'   f0 = c(rnorm(20, 200, 20), rnorm(20, 240, 20))
#' )
#' get_sufficient_category_statistics(df, cues = c("vot", "f0"))
get_sufficient_category_statistics <- function(
  data,
  cues,
  category = "category",
  group = NULL,
  stats = c("n", "mean", "cov", "css", "uss"),
  simplify = TRUE,
  ...
) {
  .assert_that(
    is.data.frame(data),
    msg = "data must be a data.frame or tibble."
  )
  .assert_non_NA_character(
    cues,
    msg = "cues must be a non-empty character vector."
  )
  .assert_true(
    length(cues) >= 1L,
    msg = "cues must contain at least one cue name."
  )
  .assert_true(
    all(cues %in% names(data)),
    msg = "All specified cues must be columns in data."
  )
  .assert_non_NA_character(
    category,
    msg = "category must be a character string specifying the category column."
  )
  .assert_true(
    length(category) == 1L,
    msg = "category must be a single column name."
  )
  .assert_true(
    category %in% names(data) || nrow(data) == 0L,
    msg = paste0("Column '", category, "' not found in data.")
  )
  if (!is.null(group)) {
    .assert_non_NA_character(
      group,
      msg = "group must be a character vector of group column names."
    )
    .assert_true(
      all(group %in% names(data)) || nrow(data) == 0L,
      msg = "All specified group columns must exist in data."
    )
  }

  stats_missing <- missing(stats)
  valid_stats <- c("n", "mean", "cov", "css", "uss")
  stats <- match.arg(stats, valid_stats, several.ok = TRUE)

  dots <- list(...)
  model_family <- if (!is.null(dots$model_family)) {
    dots$model_family
  } else {
    c("NIW", "NIX", "MNIX", "MUVG", "MVG", "UVG", "EXEMPLAR")
  }
  model_family <- match.arg(
    model_family,
    list_model_families(),
    several.ok = TRUE
  )

  compute_full_matrix <- any(c("NIW", "MVG") %in% model_family)
  compute_mnix_vector <- any(c("MNIX", "MUVG") %in% model_family)
  compute_nix_scalar <- any(c("NIX", "UVG") %in% model_family)

  if (!is.null(dots$model_family) && !compute_full_matrix && stats_missing) {
    stats <- setdiff(stats, c("cov", "uss"))
  }
  data <- as.data.frame(data)

  group_cols <- c(category, group)
  valid_group_cols <- intersect(group_cols, names(data))

  complete_idx <- stats::complete.cases(data[, cues, drop = FALSE])
  data <- data[complete_idx, , drop = FALSE]

  if (nrow(data) == 0L) {
    res <- data[, valid_group_cols, drop = FALSE]
    if ("n" %in% stats) res$n <- integer(0)
    if ("mean" %in% stats) res$mean <- list()
    if ("cov" %in% stats) res$cov <- list()
    if ("css" %in% stats) res$css <- list()
    if ("uss" %in% stats) res$uss <- list()
    if (isTRUE(simplify) && length(stats) == 1L) {
      return(res[[stats]])
    }
    return(tibble::as_tibble(res))
  }

  data <- droplevels(data)

  by_list <- lapply(group_cols, function(col) data[[col]])
  names(by_list) <- group_cols

  idx_split <- split(
    seq_len(nrow(data)),
    by_list,
    drop = TRUE,
    lex.order = TRUE
  )
  first_indices <- vapply(idx_split, function(idx) idx[1L], integer(1L))
  res <- data[first_indices, group_cols, drop = FALSE]
  rownames(res) <- NULL

  n_groups <- length(idx_split)

  need_cov <- "cov" %in% stats
  need_uss <- "uss" %in% stats
  need_css <- any(c("css", "cov", "uss") %in% stats)
  need_mean <- any(c("mean", "uss") %in% stats)
  need_n <- any("n" %in% stats, need_cov, need_uss)

  cue_mat <- as.matrix(data[, cues, drop = FALSE])
  if (!is.numeric(cue_mat)) {
    storage.mode(cue_mat) <- "double"
  }

  n_vals <- if (need_n) integer(n_groups) else NULL
  mean_list <- if (need_mean) vector("list", n_groups) else NULL
  css_list <- if (need_css) vector("list", n_groups) else NULL
  uss_list <- if (need_uss) vector("list", n_groups) else NULL
  cov_list <- if (need_cov) vector("list", n_groups) else NULL

  for (i in seq_len(n_groups)) {
    idx <- idx_split[[i]]
    sub_mat <- cue_mat[idx, , drop = FALSE]
    n_obs <- length(idx)

    if (need_n) n_vals[i] <- n_obs
    if (need_mean) mean_list[[i]] <- colMeans(sub_mat)

    if (compute_full_matrix) {
      if (need_uss) uss_list[[i]] <- uss(sub_mat)
      if (need_css || need_cov) {
        css_mat <- css(sub_mat)
        if (need_css) css_list[[i]] <- css_mat
      }
      if (need_cov) {
        cov_mat <- stats::cov(sub_mat)
        if (length(cues) == 1L) {
          cov_mat <- matrix(cov_mat, 1L, 1L, dimnames = list(cues, cues))
        }
        cov_list[[i]] <- cov_mat
      }
    } else if (compute_mnix_vector) {
      if (need_css || need_cov) {
        if (n_obs > 1L) {
          centered <- sweep(sub_mat, 2, colMeans(sub_mat), "-")
          vec_ss <- colSums(centered^2)
        } else {
          vec_ss <- rep(0, length(cues))
        }
        names(vec_ss) <- cues
        if (need_css) css_list[[i]] <- vec_ss
      }
      if (need_cov) {
        cov_mat <- if (n_obs > 1L) {
          diag(vec_ss / (n_obs - 1L), length(cues), length(cues))
        } else {
          matrix(NA_real_, length(cues), length(cues))
        }
        dimnames(cov_mat) <- list(cues, cues)
        cov_list[[i]] <- cov_mat
      }
      if (need_uss) {
        uss_vec <- colSums(sub_mat^2)
        names(uss_vec) <- cues
        uss_list[[i]] <- uss_vec
      }
    } else if (compute_nix_scalar) {
      if (need_css || need_cov) {
        scal_ss <- if (n_obs > 1L) {
          sum((sub_mat[, 1L] - mean(sub_mat[, 1L]))^2)
        } else {
          0
        }
        if (need_css) css_list[[i]] <- scal_ss
      }
      if (need_cov) {
        cov_mat <- if (n_obs > 1L) {
          matrix(
            scal_ss / (n_obs - 1L),
            1L,
            1L,
            dimnames = list(cues[1L], cues[1L])
          )
        } else {
          matrix(NA_real_, 1L, 1L, dimnames = list(cues[1L], cues[1L]))
        }
        cov_list[[i]] <- cov_mat
      }
      if (need_uss) {
        uss_list[[i]] <- sum(sub_mat[, 1L]^2)
      }
    } else {
      if (need_css) css_list[[i]] <- NULL
      if (need_uss) uss_list[[i]] <- NULL
      if (need_cov) cov_list[[i]] <- NULL
    }
  }

  if ("n" %in% stats) res$n <- n_vals
  if ("mean" %in% stats) res$mean <- mean_list
  if ("cov" %in% stats) res$cov <- cov_list
  if ("css" %in% stats) res$css <- css_list
  if ("uss" %in% stats) res$uss <- uss_list

  if (isTRUE(simplify) && nrow(res) == 1L && length(stats) == 1L) {
    return(res[[stats]][[1L]])
  }

  tibble::as_tibble(res)
}

#' Get affine transformation
#'
#' Computes an affine transformation for a set of \code{cues} in \code{data}.
#' Returns an \code{\link{MVBU_TransformInformation}} object representing a
#' linear map \eqn{f(x) = A(x + b)}, where \eqn{b} is the shift vector and
#' \eqn{A} is the scale/rotation matrix.
#'
#' @param data A \code{tibble} or \code{data.frame}.
#' @param cues A character vector of column names in \code{data} that should be
#'   transformed.
#' @param type The type of transformation to apply. Can be one of
#'   \code{"identity"}, \code{"center"}, \code{"standardize"},
#'   \code{"PCA whiten"}, or \code{"ZCA whiten"}. Except for the identity
#'   transform, all transforms center. Standardize additionally divides by the
#'   standard deviations. PCA whitening rotates the data into the space of
#'   principal components with unit variance. ZCA whitening achieves unit
#'   variance while maintaining original orientation as close as possible.
#'   (default: \code{"identity"})
#'
#' @return An \code{\link{MVBU_TransformInformation}} object that can be used to
#'   transform data or models.
#'
#' @export
#' @examples
#' df <- data.frame(
#'   vot = rnorm(50, 20, 10),
#'   f0 = rnorm(50, 220, 25)
#' )
#' # Standardize transformation
#' trans <- get_affine_transform(df, cues = c("vot", "f0"),
#'                               type = "standardize")
#' trans
get_affine_transform <- function(
  data,
  cues,
  type = c("identity", "center", "standardize", "PCA whiten", "ZCA whiten")[1]
) {
  .assert_that(is.character(cues))
  .assert_data_contains_cols(data, cues)

  if (identical(type, "scale")) {
    type <- "standardize"
  }

  n.cues <- length(cues)
  data_mat <- as.matrix(data[, cues, drop = FALSE])

  if (type == "identity") {
    shift <- rep(0, n.cues)
  } else {
    data_mat <- scale(data_mat, scale = FALSE)
    shift <- -(attr(data_mat, "scaled:center"))
  }

  if (n.cues == 1L) {
    if (type %in% c("identity", "center")) {
      SCALE <- 1
    } else if (type %in% c("standardize", "PCA whiten", "ZCA whiten")) {
      SCALE <- 1 / sd(data_mat[, 1])
    } else {
      .stop("Unknown type: ", type)
    }
  } else {
    if (type %in% c("identity", "center")) {
      SCALE <- diag(n.cues)
    } else if (type == "standardize") {
      SCALE <- diag(1 / apply(data_mat, 2, sd))
    } else if (type %in% c("PCA whiten", "ZCA whiten")) {
      eig <- eigen(.cov(data_mat))
      U <- eig$vectors
      D_inv_half <- diag(1 / sqrt(eig$values))

      if (type == "PCA whiten") {
        SCALE <- D_inv_half %*% t(U)
      } else if (type == "ZCA whiten") {
        SCALE <- U %*% D_inv_half %*% t(U)
      } else {
        .stop("Unknown whitening type: ", type)
      }
    } else {
      .stop("Unknown type: ", type)
    }
  }

  MVBU_TransformInformation(
    cue.labels = cues,
    type = type,
    shift = shift,
    SCALE = SCALE
  )
}


#' Transform and untransform means, covariances, and cognitive models
#'
#' Applies an affine (linear) transformation specified in an
#' \code{\link{MVBU_TransformInformation}} object (from
#' \code{\link{get_affine_transform}}) to a mean vector, covariance matrix, or
#' within an S7 cognitive model.
#'
#' Gaussian and Normal-Inverse-Wishart distributions transform analytically
#' under affine maps \eqn{y = A(x + b)}:
#' \deqn{\mu' = A(\mu + b), \quad \Sigma' = A \Sigma A^T}
#'
#' @param model An \code{\link{MVBU_CognitiveModel}} object.
#' @param m A mean vector (or scalar for univariate representations).
#' @param S A covariance matrix (or scalar variance for univariate
#'   representations).
#' @param transform An \code{\link{MVBU_TransformInformation}} object.
#'
#' @return A model, mean vector, or covariance matrix of the same class as the
#'   input.
#'
#' @rdname transform_model
#' @export
#' @examples
#' df <- data.frame(vot = rnorm(30, 20, 10), f0 = rnorm(30, 220, 25))
#' trans <- get_affine_transform(df, cues = c("vot", "f0"),
#'                               type = "standardize")
#'
#' # Transform mean vector and covariance matrix
#' m <- c("vot" = 20, "f0" = 220)
#' S <- matrix(c(100, 20, 20, 625), nrow = 2,
#'             dimnames = list(c("vot", "f0"), c("vot", "f0")))
#' m_trans <- transform_mean(m, trans)
#' S_trans <- transform_cov(S, trans)
#'
#' # Untransform back to original space
#' untransform_mean(m_trans, trans)
#' untransform_cov(S_trans, trans)
#'
#' # Transform an entire cognitive model
#' mod <- example_mvg_ideal_observer(n_cues = 2)
#' mod_trans <- transform_model(mod, trans)
#' untransform_model(mod_trans, trans)
transform_model <- function(model, transform) {
  if (!S7::S7_inherits(model, MVBU_CognitiveModel)) {
    stop("model must be an MVBU_CognitiveModel object.")
  }
  if (!S7::S7_inherits(transform, MVBU_TransformInformation)) {
    stop("transform must be an MVBU_TransformInformation object.")
  }

  reps <- model@category_template@representations
  for (cat_name in names(reps)) {
    r <- reps[[cat_name]]
    if (S7::S7_inherits(r, MVG_CategoryRepresentation)) {
      r@mu <- transform_mean(r@mu, transform)
      r@Sigma <- transform_cov(r@Sigma, transform)
    } else if (S7::S7_inherits(r, NIW_CategoryRepresentation)) {
      r@m <- transform_mean(r@m, transform)
      r@S <- transform_cov(r@S, transform)
    } else if (S7::S7_inherits(r, UVG_CategoryRepresentation)) {
      r@mu <- transform_mean(r@mu, transform)
      r@sigma2 <- as.numeric(transform_cov(matrix(r@sigma2, 1, 1), transform))
    } else if (S7::S7_inherits(r, MUVG_CategoryRepresentation)) {
      r@mu <- transform_mean(r@mu, transform)
      r@sigma2 <- diag(transform_cov(diag(r@sigma2), transform))
    } else if (S7::S7_inherits(r, NIX_CategoryRepresentation)) {
      r@m <- transform_mean(r@m, transform)
      r@sigma2 <- as.numeric(transform_cov(matrix(r@sigma2, 1, 1), transform))
    } else if (S7::S7_inherits(r, MNIX_CategoryRepresentation)) {
      r@m <- transform_mean(r@m, transform)
      r@sigma2 <- diag(transform_cov(diag(r@sigma2), transform))
    }
    reps[[cat_name]] <- r
  }
  model@category_template@representations <- reps
  return(model)
}

#' @rdname transform_model
#' @export
untransform_model <- function(model, transform) {
  transform_model(model, invert_transform(transform))
}

#' @rdname transform_model
#' @export
transform_mean <- function(m, transform) {
  if (!S7::S7_inherits(transform, MVBU_TransformInformation)) {
    stop("transform must be an MVBU_TransformInformation object.")
  }
  m_names <- names(m)
  shift_val <- transform@shift
  scale_val <- transform@SCALE
  k <- length(shift_val)

  if (k > 0L && length(shift_val) == length(m)) {
    m <- m + shift_val
  }

  if (k == 1L || (is.numeric(scale_val) && length(scale_val) == 1L)) {
    m <- as.numeric(m) * as.numeric(scale_val)
  } else if (is.matrix(scale_val)) {
    m <- as.vector(m %*% t(scale_val))
  }

  names(m) <- m_names
  return(m)
}

#' @rdname transform_model
#' @export
untransform_mean <- function(m, transform) {
  transform_mean(m, invert_transform(transform))
}

#' @rdname transform_model
#' @export
transform_cov <- function(S, transform) {
  if (!S7::S7_inherits(transform, MVBU_TransformInformation)) {
    stop("transform must be an MVBU_TransformInformation object.")
  }
  S_names <- dimnames(S)
  scale_val <- transform@SCALE
  k <- length(transform@shift)

  if (k == 1L || (is.numeric(scale_val) && length(scale_val) == 1L)) {
    S <- as.numeric(scale_val)^2 * S
  } else if (is.matrix(scale_val)) {
    S <- scale_val %*% S %*% t(scale_val)
  }

  dimnames(S) <- S_names
  return(S)
}

#' @rdname transform_model
#' @export
untransform_cov <- function(S, transform) {
  transform_cov(S, invert_transform(transform))
}
