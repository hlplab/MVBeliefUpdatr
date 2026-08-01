# Internal utility helpers shared across migration scaffolds.

#' @keywords internal
.is_numeric_vector <- function(x) {
  is.numeric(x) && is.null(dim(x))
}

is_scalar_character <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x[1])
}

is_character <- function(x) {
  is.character(x)
}

#' @keywords internal
.as_numeric_vector <- function(x) {
  if (!.is_numeric_vector(x)) {
    stop("Expected a numeric vector.", call. = FALSE)
  }
  as.numeric(x)
}


#' @export
is.Sigma <- function(x) {
  if (is.null(x)) stop2("Expected a covariance matrix, but got NULL.")
  if (is.matrix(x)) {
    if (all(x == 0) | is.positive.definite(x)) return(T) else return(F)
  } else {
    if (is_scalar_double(x)) return(T) else return(F)
  }
}


# dim that returns length of vector for vector (not from brms)
dim2 <- function(x) {
  if (is.null(dim(x))) return(length(x))
  return(dim(x))
}

replace_na_in_array <- function(x, fill = 0) {
  .assert_true(is.array(x), msg = "x must be an array.")
  .assert_true(is.scalar(fill), msg = "fill must be a scalar value.")

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
  .assert_true(all(groups %in% names(data)), msg = "All grouping columns must be present in data.")
  .assert_true(all(cols %in% names(data)), msg = "All aggregate columns must be present in data.")
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
