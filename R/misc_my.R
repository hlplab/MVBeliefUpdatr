get_current_versions <- function()
  list(
    MVBeliefUpdatr = utils::packageVersion("MVBeliefUpdatr"),
    rstan = utils::packageVersion("rstan"),
    stanHeaders = utils::packageVersion("StanHeaders"))


# dim that returns length of vector for vector (not from brms)
dim2 <- function(x) {
  if (is.null(dim(x))) return(length(x))
  return(dim(x))
}

replace_na_in_array <- function(x, fill = 0) {
  stopifnot(is.array(x))
  stopifnot(is.scalar(fill))

  x[is.na(x)] <- fill
  return(x)
}

#' Get aggregates from grouped data as list of lists or list of arrays
#'
#' Get aggregates from grouped data for any number aggregation functions specified in `...`. The columns specified in
#' `cols` jointly (essentially as a matrix) form the input to the aggregation functions. Aggregates are obtained
#' separately for each unique combinations values for the variables specified in the grouping variables `groups`
#' (which must be factors). For combinations of the grouping variables for which there is no data, the value in
#' `fill` will be substituted, repeated up to the necessary dimensionality (inferred from the other outputs of
#' the aggregation function).
#'
#' @param data `tibble` or `data.frame` with the data.
#' @param cols Character vector of column names for which aggregate values should be obtained. This can be more than
#'   one column, but only one aggregate will be returned \emph{across} columns, i.e., the aggregate will be calculated for
#'   each unique group in the data.
#' @param groups Character vector of column names that contain information about which observations form a group.
#'   These columns must factors.
#' @param fill What value should be substituted for `groups` for which there is no data (see `...`)? (default: `NA`)
#' @param ... A named list of aggregating functions of `cols` to calculate for each unique combination of grouping
#'   variables `groups` in `data`.
#'
#' @return A named list of length `...` (names are the names of the functions that have been computed). The list will
#'   be sorted by the `group` variables (in the order they are presented), in ascending order of the levels of those
#'   variables. The elements of the list will be lists.
#'   Missing values---resulting from combinations of grouping variables for which there is no data---will be filled
#'   with the value provided in `fill` and coerced into the same structure as all other outputs for that function
#'   (e.g., if the function f results in a 2-element vector for all combinations of `groups` variables for which
#'   there is data, then unobserved data results in a 2-element vector of `fill` values).
#'
#' @keywords TBD
#' @rdname get_aggregates_from_grouped_data_as_list_of_lists
#' @export
get_aggregates_from_grouped_data_as_list_of_lists <- function(
    data,
    groups,
    cols,
    fill = as.list(rep(NA, length(list(...)))),
    verbose = F,
    ...
) {
  stopifnot(all(groups %in% names(data)))
  stopifnot(all(cols %in% names(data)))
  if (!all(map_lgl(groups, ~ is.factor(data[[.x]])))) {
    stop2("Group variables must be factors.")
  }
  if (length(fill) != length(list(...))) stop2("fill must be a list of equal length as the number of functions provided in `...`.")

  fn_list <- list(...)

  data %<>%
    select(all_of(groups), all_of(cols)) %>%
    # Remove "." from group names to prevent it from leading to issues below
    # (since split uses "." as a separator)
    mutate(
      across(
        all_of(groups),
        ~ factor(
          gsub("\\.", "_", .x),
          levels = gsub("\\.", "_", levels(.x)))))

  # Get sorted levels
  group_levels <- map(groups, ~ levels(data[[.x]]))

  # Split the data by group and category
  split_data <- split(data, map(groups, ~ data[[.x]]), drop = FALSE)

  # Sort the split list by factor levels of group (first) and category (second)
  split_keys <-
    names(split_data)  %>%
    strsplit("\\.") %>%
    do.call(rbind, .) %>%
    as.data.frame(., stringsAsFactors = FALSE)

  split_order <-
    split_keys %>%
    Map(function(col, levs) factor(col, levels = levs),
        .,
        group_levels) %>%
    do.call(order, .)

  split_data <- split_data[split_order]

  # Initialize result list
  result <- lapply(fn_list, function(f) vector("list", length(split_data)))
  names(result) <- names(fn_list)

  # Apply functions
  for (i in seq_along(split_data)) {
    x <- do.call(rbind, lapply(seq_len(nrow(split_data[[i]])), function(j) {
      as.numeric(split_data[[i]][j, cols, drop = FALSE])
    }))
    for (fname in names(fn_list)) {
      # Call function and make sure that `NA` is returned instead if x is empty data frame
      result[[fname]][[i]] <- if (is.null(x)) NA else fn_list[[fname]](x)
    }
  }

  # For each function, get the highest dimensionality of any of its outputs, and then coerce
  # NA outputs into that same dimensionality.
  dim_target <- result %>% map(~ reduce((map(.x, dim2)), pmax))
  result <-
    pmap(
      .l = list(result, dim_target, fill),
      .f =
        ~ map(
          ..1,
          function(x) {
            if (any(is.na(x))) {
              if (verbose)
                message("Empty data found for a combination of grouping variables. Filling NA values with fill value (",
                        paste(..3, collapse = ", "),
                        "). If necessary and possible, an attempt will be made to coerce the fill value into the same ",
                        " dimensionality as the other outputs for this aggregation function (",
                        paste(..2, collapse = ", "), ").")
              if (all(dim2(..3) == ..2)) {
                return(..3)
              } else if (!is_scalar_double(..3)) {
                stop2("Attempt to coerce fill value into the same dimensionality as the other outputs for this aggregation function failed: fill value must either be a scalar or a vector/matrix/array with the same dimensionality as the other outputs for this aggregation function.")
              }

              x <-
                if (length(..2) == 1) {
                  # Coerce fill into vector
                  rep(..3, ..2)
                } else if (length(..2) == 2) {
                  # Coerce fill into matrix
                  matrix(..3, nrow = ..2[1], ncol = ..2[2])
                } else if (length(..2) > 2) {
                  # Coerce fill into array
                  array(rep(..3, prod(..2)), dim = ..2)
                } else ..3

              return(x)
            } else return(x)
          }))
  # Check that all elements have the same dimensionality as the dim_target
  if (!all(map_lgl(1:length(dim_target), ~ map_lgl(result[[.x]], function(x) all(dim2(x) == dim_target[[.x]])) %>% reduce(all))))
    stop2(
      paste0(
        "The results for some aggregate functions differ in their dimensionality (even after dealing with missing data): highest dimensionalities found (",
        paste(map_chr(1:length(dim_target), ~ paste0(names(dim_target)[.x], ": ", paste(dim_target[[.x]], collapse = ", "))), collapse = "; "), ")"))

  return(result)
}

#' @rdname get_aggregates_from_grouped_data_as_list_of_lists
#' @export
get_aggregates_from_grouped_data_as_list_of_arrays <- function(
    data,
    groups,
    cols,
    fill = as.list(rep(NA, length(list(...)))),
    simplify = as.list(rep(T, length(list(...)))),
    verbose = F,
    ...
) {
  if (!is.list(simplify) || !all(map_lgl(simplify, is.logical))) stop2("Argument simplify must be a list of logicals.")
  if (length(simplify) != length(list(...))) stop2("Simplify must be a list of equal length as the number of functions provided in `...`.")

  result <-
    get_aggregates_from_grouped_data_as_list_of_lists(
      data = data,
      groups = groups,
      cols = cols,
      fill = fill,
      verbose = verbose,
      ...)

  # (for now, only) name outer dimensions of array
  dimnames <- rev(map(groups, ~ levels(data[[.x]])))
  result %<>%
    map2(
      .y = simplify,
      .f =
        function(.x, .y) {
          to_array(
            .x,
            # reverse the outer dimensions since to_array iterates (first and thus fastest)
            # over the first group variable, followed by the second, etc., whereas the sorting
            # of get_aggregates_from_grouped_data_as_list_of_lists is the opposite (first sort
            # by first group variable, thus iterating *slowest* over it).
            outer_dims = rev(map_int(groups, ~ nlevels(data[[.x]]))),
            dimnames = dimnames,
            simplify = .y)
        })

  return(result)
}


#' Turn (lists of) numerics into arrays (for Stan input)
#'
#' Takes `NULL`, numeric atomics (single scalars, vectors, matrices) or lists of numeric elements as inputs and turns them into
#' numeric arrays. For lists, the \emph{outer} dimension of the array iterate over the list elements (for details, see argument
#' `outer_dims`).
#'
#' @param x The input to be turned into an array.
#' @param inner_dims Integer vector of intended dimensions of x.
#'   If `NULL`, `inner_dims` will be inferred from the input.
#'   If not `NULL`, `inner_dims` will be used for checks and to enforce the correct dimensions of `x`. (default: `NULL`)
#' @param outer_dims An integer vector of outer dimensions for the array. If `x` is not a list, an array of `x` will be repeated for each
#'   outer dimension. If `x` is list, combinations of the outer dimensions will each index one element of `x`.
#'   The total length of `x` must match the product of outer dimensions. The first dimension will be iterated
#'   over first (and thus alternating fastest), the second dimension will be iterated over second, etc. See details.
#' @param dimnames A list of character vectors with names for each dimension of the array. If not `NULL`, the length
#'   of this list must match the total number of dimensions, and each character vector must match the length of the
#'   corresponding dimension. Dimensions are named from the outside-in (i.e. first outer than inner dimensions).
#'   (default: `NULL`)
#' @param simplify Should the array be simplified by removing inner dimensions whenever the inner element(s) are just scalars?
#'   (e.g., if the inner elements are a one-element vector or a 1x1 matrix). Note that atomic values that can be simplified
#'   will be returned as scalars (and thus without dimnames) if there are no outer dimensions, since the returned object cannot
#'   be an array in that case. (default: `TRUE`)
#'
#' @return An array.
#'
#' @details
#' \itemize{
#'   \item{`NULL` will be turned in an `array(numeric(0), dim = c(outer_dims, inner_dims))`. This is sometimes required for
#'   optional Stan inputs.}
#'   \item{Atomics (scalars, vectors, and matrices) will be turned into an `array(x, dim = c(outer_dims, inner_dims))`.
#'   This essentially creates copies of `x` for each combination of `outer_dims`.}
#'   \item{Lists will be turned into an `array(x, dim = c(outer_dims, inner_dims))`. `prod(outer_dims)` must match the length
#'   of the list.}
#' }
#'
#' E.g., if the input is a list of six 2x2 covariance matrices, and `...` is `3, 2`, then output will
#' be an array of dimensionality `c(3, 2, 2, 2)`.
#'
#'
#' @keywords TBD
#' @export
to_array <- function(
    x,
    inner_dims = NULL,
    outer_dims = NULL,
    dimnames = NULL,
    simplify = T
) {
  stopifnot(is.null(inner_dims) || all(inner_dims == round(inner_dims)))
  stopifnot(is.null(outer_dims) || all(outer_dims == round(outer_dims)))

  if (is.null(x)) {
    # optionally simplify by removing inner_dims if all inner dimension have length 1 (i.e., if all elements are intended to be scalars)
    if (simplify) { inner_dims <- if (all(inner_dims == 1)) NULL else inner_dims }

    total_dim_length <- length(inner_dims) + length(outer_dims)
    arr <- array(numeric(), dim = if (total_dim_length == 0) 0 else rep(0, total_dim_length))
  } else
    # Handle atomic input
    if (is.atomic(x)) {
      # NOTE: no simplification is currently being applied for atomic elements.
      found_inner_dims <- dim2(x)

      if (is.null(inner_dims)) inner_dims <- found_inner_dims
      if (any(found_inner_dims != inner_dims))
        stop2(paste0("Input's inner dimensions (", paste(found_inner_dims, collapse = ", "), ") do not match the provided inner dimension (", paste(inner_dims, collapse = ","), ")."))

      # If desired inner_dims don't match found inner_dims, see whether the input can be coerced into the desired inner_dims format
      if (length(found_inner_dims) != length(inner_dims)) {
        if (prod(inner_dims) == prod(found_inner_dims)) {
          if (length(inner_dims) == 1) {
            # Coerce into scalar or vector format
            x <- as.vector(x)
          } else if (length(inner_dims) > 1) {
            # Coerce into matrix format
            x <- matrix(x, nrow = inner_dims[1], ncol = inner_dims[2])
          }
        } else {
          stop2(paste0("Input's inner dimensions (", paste(found_inner_dims, collapse = ", "), ") do not match the provided inner dimension (", paste(inner_dims, collapse = ","), ")."))
        }
      }

      # optionally simplify by removing inner_dims if all inner dimension have length 1 (i.e., if all elements are essentially scalars)
      # Since array with no dimensions are not allowed, special handling is required when there are no outer dimensions and we still
      # want to simplify. In that case, x as a scalar without dimnames.
      if (simplify) { inner_dims <- if (all(inner_dims == 1)) { if (!is.null(outer_dims)) NULL else return(as.numeric(x)) } else inner_dims }

      # First cast an array with outer_dims on the outside (necessary to get correct handling of matrices when outer_dims are requested)
      total_dim_length <- length(inner_dims) + length(outer_dims)
      arr <- rep(x, prod(outer_dims)) %>% array(dim = if (total_dim_length == 0) 0 else c(inner_dims, outer_dims))
      # Move outer dimensions to the front
      perm <- c(if (is.null(outer_dims)) NULL else seq(length(inner_dims) + 1, total_dim_length), if (is.null(inner_dims)) NULL else 1:length(inner_dims))
      arr <- aperm(a = arr, perm = perm)
    } else
      # Handle list input
      if (is.list(x)) {
        # Use first list element to infer detected dimensions
        found_inner_dims <- dim2(x[[1]])

        if (is.null(inner_dims)) inner_dims <- found_inner_dims
        if (any(found_inner_dims != inner_dims)) stop2(paste0("Input's inner dimensions (", paste(found_inner_dims, collapse = ","), ") do not match the provided inner dimension (", paste(inner_dims, collapse = ","), ")."))

        # Check that x has the right number of elements
        expected_len <- if (length(outer_dims) == 0) 1 else prod(outer_dims)
        if (length(x) != expected_len) {
          stop2(paste0("Length of list (", length(x), ") input does not match product of provided outer dimensions (", expected_len, ")."))
        }

        # Combine list elements into a matrix
        arr <- simplify2array(x, except = NULL)

        # Move last dimensions to the front
        total_dims <- dim(arr)
        perm <- c(length(total_dims), 1:(length(total_dims) - 1))
        arr <- aperm(a = arr, perm = perm)

        # optionally simplify by removing inner_dims if all inner dimension have length 1 (i.e., if all elements are essentially scalars)
        if (simplify) { inner_dims <- if (all(inner_dims == 1)) NULL else inner_dims }
        # split first dimension into all outer dimensions
        arr <- array(arr, dim = c(outer_dims, inner_dims))
      } else
        stop2("Unsupported input type.")

  # Name dimensions of array
  # (no special handling and checking necessary since dimnames already does that: it exhaustively names from the outer to
  # inner dimensions until it hits a mismatch)
  dimnames(arr) <- dimnames

  return(arr)
}
