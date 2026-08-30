#' @include S7-core-classes.R
#' @include S7-stanfit.R
#' @importFrom lifecycle deprecate_warn
NULL

# deprecated ------------------------------------------------------------------

.get_constructor <- function(x, variable = NULL) {
  available_constructors <- c("category", "group", "cue", "cue2")

  stanfit <- if (inherits(x, "stanfit")) {
    x
  } else if (S7::S7_inherits(x, MVBU_Stanfit)) {
    get_stanfit(x)
  } else {
    .assert_that(
      FALSE,
      msg = "x must be a stanfit or ideal adaptor fit object"
    )
  }

  constructors <- if (!is.null(stanfit) &&
    !is.null(attr(stanfit, "tidybayes_constructors"))) {
    attr(stanfit, "tidybayes_constructors")
  } else {
    NULL
  }

  if (is.null(variable)) {
    return(constructors)
  }

  .assert_that(
    variable %in% available_constructors,
    msg = paste0(
      "Variable name must be one of ",
      paste(available_constructors, collapse = ", "), "."
    )
  )

  if (is.null(constructors[[variable]])) {
    warning(
      paste0(
        class(x)[1], " object does not contain type information about ",
        variable, ". Applying recover_types() to the object might fix this."
      )
    )
    return(NULL)
  }

  constructors[[variable]]
}

#' Get tidybayes constructor from an ideal adaptor stanfit (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_constructor()`, `get_category_constructor()`, `get_group_constructor()`,
#' `get_cue_constructor()`, and `get_cue2_constructor()` are deprecated in favor of
#' \code{\link{get_category_labels}}, \code{\link{get_group_labels}}, and
#' \code{\link{get_cue_labels}}.
#'
#' @param x An \code{\link{IdealAdaptorStanfit}} or \code{\link{MVBU_Stanfit}} object.
#' @param variable Either "category", "group", "cue", or "cue2". If set to
#'   `NULL` then a list of all constructors is returned. (default: `NULL`)
#'
#' @return A constructor function, a list of constructor functions, or `NULL`.
#' @export
get_constructor <- function(x, variable = NULL) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_constructor()",
    details = "Use get_category_labels(), get_group_labels(), or get_cue_labels() instead."
  )
  .get_constructor(x, variable)
}

#' @rdname get_constructor
#' @export
get_category_constructor <- function(x) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_category_constructor()",
    "get_category_labels()"
  )
  .get_constructor(x, "category")
}

#' @rdname get_constructor
#' @export
get_group_constructor <- function(x) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_group_constructor()",
    "get_group_labels()"
  )
  .get_constructor(x, "group")
}

#' @rdname get_constructor
#' @export
get_cue_constructor <- function(x) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_cue_constructor()",
    "get_cue_labels()"
  )
  .get_constructor(x, "cue")
}

#' @rdname get_constructor
#' @export
get_cue2_constructor <- function(x) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_cue2_constructor()",
    "get_cue_labels()"
  )
  .get_constructor(x, "cue2")
}
