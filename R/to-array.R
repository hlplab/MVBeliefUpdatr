#' Convert atomic values or lists into arrays with inner/outer dimensions
#'
#' This helper is used to shape prior information and other Stan inputs into
#' the dimensionality expected by the Stan programs.
#'
#' @param x Input value. Can be atomic, list, or NULL.
#' @param inner_dims Integer vector describing the inner dimensions of the
#'   array.
#' @param outer_dims Integer vector describing the outer dimensions of the
#'   array.
#' @param dimnames Optional dimnames for the resulting array.
#' @param simplify Logical controlling whether inner dimensions of length 1 are
#'   simplified.
#' @return An array.
#' @export
#' @examples
#' # Convert a vector into an array with outer dimension 3
#' to_array(c(1, 2), inner_dims = 2, outer_dims = 3)
#'
#' # Convert a list of matrices into an array
#' mat_list <- list(matrix(1:4, 2, 2), matrix(5:8, 2, 2))
#' to_array(mat_list, inner_dims = c(2, 2), outer_dims = 2)
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
    dim_vec <- if (total_dim_length == 0) 0 else rep(0, total_dim_length)
    arr <- array(numeric(), dim = dim_vec)
    if (!is.null(dimnames)) dimnames(arr) <- dimnames
    return(arr)
  }

  if (is.atomic(x)) {
    found_inner_dims <- .dim(x)

    if (is.null(inner_dims)) inner_dims <- found_inner_dims
    if (any(found_inner_dims != inner_dims)) {
      .stop(
        "Input's inner dimensions (",
        paste(found_inner_dims, collapse = ", "),
        ") do not match the provided inner dimension (",
        paste(inner_dims, collapse = ","),
        ")."
      )
    }

    if (length(found_inner_dims) != length(inner_dims)) {
      if (prod(inner_dims) == prod(found_inner_dims)) {
        if (length(inner_dims) == 1) {
          x <- as.vector(x)
        } else if (length(inner_dims) > 1) {
          x <- matrix(x, nrow = inner_dims[1], ncol = inner_dims[2])
        }
      } else {
        .stop(
          "Input's inner dimensions (",
          paste(found_inner_dims, collapse = ", "),
          ") do not match the provided inner dimension (",
          paste(inner_dims, collapse = ","),
          ")."
        )
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
    dim_vec <- if (total_dim_length == 0) 0 else c(inner_dims, outer_dims)
    arr <- array(rep(x, prod(outer_dims)), dim = dim_vec)
    perm <- c(
      if (is.null(outer_dims)) NULL else {
        seq(length(inner_dims) + 1, total_dim_length)
      },
      if (is.null(inner_dims)) NULL else 1:length(inner_dims)
    )
    arr <- aperm(a = arr, perm = perm)
    if (!is.null(dimnames)) dimnames(arr) <- dimnames
    return(arr)
  }

  if (is.list(x)) {
    found_inner_dims <- .dim(x[[1]])

    if (is.null(inner_dims)) inner_dims <- found_inner_dims
    if (any(found_inner_dims != inner_dims)) {
      .stop(
        "Input's inner dimensions (",
        paste(found_inner_dims, collapse = ", "),
        ") do not match the provided inner dimension (",
        paste(inner_dims, collapse = ","),
        ")."
      )
    }

    expected_len <- if (length(outer_dims) == 0) 1 else prod(outer_dims)
    if (length(x) != expected_len) {
      .stop(
        "Length of list (", length(x),
        ") input does not match product of provided outer dimensions (",
        expected_len, ")."
      )
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

