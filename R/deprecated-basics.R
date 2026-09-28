# =============================================================================
# deprecated
# =============================================================================

#' Deprecated: Combine columns into a vector column
#'
#' @description `r lifecycle::badge("deprecated")`
#' `make_vector_column()` is deprecated as of MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0.
#'
#' @param data A tibble or data.frame.
#' @param cols A character vector of variable names to combine.
#' @param vector_col Name of the new vector-valued column.
#' @param .keep A tidyselect option passed to dplyr::mutate.
#'
#' @return Same as \code{data}.
#'
#' @keywords internal
#' @export
make_vector_column <- function(data, cols, vector_col, .keep = "all") {
    lifecycle::deprecate_warn(
        when = "0.1.0",
        what = "make_vector_column()"
    )
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


#' Deprecated: Get sum of squares from a data frame or matrix
#'
#' @description `r lifecycle::badge("deprecated")`
#' `get_sum_of_squares_from_df()` and its aliases are deprecated as of MVBeliefUpdatr 0.1.0
#' and will be removed in 0.2.0. Please use [get_sufficient_category_statistics()] instead.
#'
#' @param data A `tibble`, `data.frame`, or `matrix`.
#' @param variables Only required if data is not already a `matrix`.
#' @param center Logical; if `TRUE`, compute centered sum of squares.
#' @param verbose Logical; produce verbose output.
#'
#' @return A matrix.
#'
#' @keywords internal
#' @rdname deprecated-sum-of-squares
#' @export
get_sum_of_squares_from_df <- function(data, variables = NULL, center = TRUE, verbose = FALSE) {
    lifecycle::deprecate_warn(
        when = "0.1.0",
        what = "get_sum_of_squares_from_df()",
        with = "get_sufficient_category_statistics()"
    )
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
        data %>%
            mutate(across(c(!!!syms(variables)), unlist)) %>%
            select(all_of(variables)) %>%
            as.matrix()
    } else {
        data
    }

    ss(data.matrix, center = center)
}

#' @rdname deprecated-sum-of-squares
#' @keywords internal
#' @export
get_sum_of_uncentered_squares_from_df <- function(data, variables = NULL, verbose = FALSE) {
    lifecycle::deprecate_warn(
        when = "0.1.0",
        what = "get_sum_of_uncentered_squares_from_df()",
        with = "get_sufficient_category_statistics()"
    )
    get_sum_of_squares_from_df(data = data, variables = variables, center = FALSE, verbose = verbose)
}

#' @rdname deprecated-sum-of-squares
#' @keywords internal
#' @export
get_sum_of_centered_squares_from_df <- function(data, variables = NULL, verbose = FALSE) {
    lifecycle::deprecate_warn(
        when = "0.1.0",
        what = "get_sum_of_centered_squares_from_df()",
        with = "get_sufficient_category_statistics()"
    )
    get_sum_of_squares_from_df(data = data, variables = variables, center = TRUE, verbose = verbose)
}
