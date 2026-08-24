#' Convert atomic values or lists into arrays with explicit inner/outer dimensions.
#'
#' This helper is used to shape prior information and other Stan inputs into the
#' dimensionality expected by the Stan programs.
#'
#' @param x Input value. Can be atomic, list, or NULL.
#' @param inner_dims Integer vector describing the inner dimensions of the array.
#' @param outer_dims Integer vector describing the outer dimensions of the array.
#' @param dimnames Optional dimnames for the resulting array.
#' @param simplify Logical controlling whether inner dimensions of length 1 are simplified.
#' @return An array.
#' @export
to_array <- function(
    x,
    inner_dims = NULL,
    outer_dims = NULL,
    dimnames = NULL,
    simplify = TRUE
) {
  .assert_true(
    is.null(inner_dims) || all(inner_dims == round(inner_dims)),
    msg = "inner_dims must be NULL or contain whole numbers."
  )
  .assert_true(
    is.null(outer_dims) || all(outer_dims == round(outer_dims)),
    msg = "outer_dims must be NULL or contain whole numbers."
  )

  if (is.null(x)) {
    if (simplify && !is.null(inner_dims) && all(inner_dims == 1)) {
      inner_dims <- NULL
    }

    total_dim_length <- length(inner_dims) + length(outer_dims)
    arr <- array(numeric(), dim = if (total_dim_length == 0) 0 else rep(0, total_dim_length))
    if (!is.null(dimnames)) dimnames(arr) <- dimnames
    return(arr)
  }

  if (is.atomic(x)) {
    found_inner_dims <- .dim(x)

    if (is.null(inner_dims)) inner_dims <- found_inner_dims
    if (any(found_inner_dims != inner_dims)) {
      .stop(paste0("Input's inner dimensions (", paste(found_inner_dims, collapse = ", "), ") do not match the provided inner dimension (", paste(inner_dims, collapse = ","), ")."))
    }

    if (length(found_inner_dims) != length(inner_dims)) {
      if (prod(inner_dims) == prod(found_inner_dims)) {
        if (length(inner_dims) == 1) {
          x <- as.vector(x)
        } else if (length(inner_dims) > 1) {
          x <- matrix(x, nrow = inner_dims[1], ncol = inner_dims[2])
        }
      } else {
        .stop(paste0("Input's inner dimensions (", paste(found_inner_dims, collapse = ", "), ") do not match the provided inner dimension (", paste(inner_dims, collapse = ","), ")."))
      }
    }

    if (simplify && all(inner_dims == 1)) {
      if (!is.null(outer_dims)) {
        inner_dims <- NULL
      } else {
        return(as.numeric(x))
      }
    }

    total_dim_length <- length(inner_dims) + length(outer_dims)
    arr <- rep(x, prod(outer_dims)) %>% array(dim = if (total_dim_length == 0) 0 else c(inner_dims, outer_dims))
    perm <- c(if (is.null(outer_dims)) NULL else seq(length(inner_dims) + 1, total_dim_length), if (is.null(inner_dims)) NULL else 1:length(inner_dims))
    arr <- aperm(a = arr, perm = perm)
    if (!is.null(dimnames)) dimnames(arr) <- dimnames
    return(arr)
  }

  if (is.list(x)) {
    found_inner_dims <- .dim(x[[1]])

    if (is.null(inner_dims)) inner_dims <- found_inner_dims
    if (any(found_inner_dims != inner_dims)) {
      .stop(paste0("Input's inner dimensions (", paste(found_inner_dims, collapse = ", "), ") do not match the provided inner dimension (", paste(inner_dims, collapse = ","), ")."))
    }

    expected_len <- if (length(outer_dims) == 0) 1 else prod(outer_dims)
    if (length(x) != expected_len) {
      .stop(paste0("Length of list (", length(x), ") input does not match product of provided outer dimensions (", expected_len, ")."))
    }

    arr <- simplify2array(x, except = NULL)
    if (!is.null(inner_dims)) {
      arr <- array(arr, dim = c(outer_dims, inner_dims))
    }
    if (!is.null(dimnames)) dimnames(arr) <- dimnames
    return(arr)
  }

  .stop("Unsupported input type for to_array.")
}
