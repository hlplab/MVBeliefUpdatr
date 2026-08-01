#' Internal assertion helper.
#'
#' @param ... Conditions to evaluate.
#' @param msg Optional error message.
#' @return Invisibly TRUE if all conditions pass.
#' @keywords internal
.assert_that <- function(..., msg = NULL) {
  conditions <- list(...)
  if (length(conditions) == 0L) {
    return(invisible(TRUE))
  }

  cond_ok <- vapply(conditions, function(cond) {
    if (is.logical(cond) && length(cond) == 1L) {
      !is.na(cond) && isTRUE(cond)
    } else if (is.logical(cond)) {
      all(cond, na.rm = TRUE)
    } else {
      isTRUE(cond)
    }
  }, logical(1))

  if (!all(cond_ok)) {
    stop2(if (is.null(msg)) "Assertion failed." else msg)
  }

  invisible(TRUE)
}

#' Internal assertion helper.
#'
#' @param ... Conditions to evaluate.
#' @param msg Optional error message.
#' @return Invisibly TRUE if all conditions pass.
#' @keywords internal
.assert_all <- function(..., msg = NULL) {
  .assert_that(..., msg = msg)
}

#' Internal assertion helper.
#'
#' @param ... Conditions to evaluate.
#' @param msg Optional error message.
#' @return Invisibly TRUE if any condition passes.
#' @keywords internal
.assert_any <- function(..., msg = NULL) {
  conditions <- list(...)
  if (length(conditions) == 0L) {
    return(invisible(TRUE))
  }

  cond_ok <- vapply(conditions, function(cond) {
    if (is.logical(cond) && length(cond) == 1L) {
      !is.na(cond) && isTRUE(cond)
    } else if (is.logical(cond)) {
      all(cond, na.rm = TRUE)
    } else {
      isTRUE(cond)
    }
  }, logical(1))

  if (!any(cond_ok)) {
    stop2(if (is.null(msg)) "Assertion failed." else msg)
  }

  invisible(TRUE)
}

#' Internal assertion helper.
#'
#' @param cond Condition to evaluate.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the condition passes.
#' @keywords internal
.assert_true <- function(cond, msg = NULL) {
  .assert_that(cond, msg = if (is.null(msg)) "Assertion failed." else msg)
}

#' Internal assertion helper.
#'
#' @param cond Condition to evaluate.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the condition fails.
#' @keywords internal
.assert_false <- function(cond, msg = NULL) {
  .assert_true(!cond, msg = if (is.null(msg)) "Assertion failed." else msg)
}

#' Internal assertion helper.
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a logical flag.
#' @keywords internal
.assert_flag <- function(x, msg = NULL) {
  .assert_true(is.flag(x), msg = if (is.null(msg)) "Expected a logical flag." else msg)
}

#' Internal assertion helper.
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a data.frame-like object.
#' @keywords internal
.assert_data_frame_like <- function(x, msg = NULL) {
  .assert_true(is.data.frame(x) || is_tibble(x), msg = if (is.null(msg)) "Expected a data.frame or tibble." else msg)
}

#' Internal assertion helper.
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a stanfit object.
#' @keywords internal
.assert_stanfit <- function(x, msg = NULL) {
  .assert_true(inherits(x, "stanfit"), msg = if (is.null(msg)) "Expected a stanfit object." else msg)
}

#' Internal assertion helper.
#'
#' @param x Object to test.
#' @return Invisibly TRUE if the object contains draws.
#' @keywords internal
.assert_contains_draws <- function(x) {
  if (!.contains_draws(x)) stop2(paste("", deparse(substitute(x)), "does not contain any samples."))
}

#' Internal assertion helper.
#'
#' @param data Data object to inspect.
#' @param cols Column names to inspect.
#' @param which.data Description of the data source.
#' @param scalar Whether the columns should be scalar character names.
#' @return Invisibly TRUE if the columns are valid.
#' @keywords internal
.assert_cols_in_data <- function(data, cols, which.data = "the", scalar = TRUE) {
  if (scalar) {
    .assert_that(all(vapply(cols, is_scalar_character, logical(1))),
                 msg = paste0(paste(cols, collapse = ","), "must be a single column name."))
  } else {
    .assert_that(all(vapply(cols, is_character, logical(1))),
                 msg = paste0(paste(cols, collapse = ","), "must be column name or vector of column names."))
  }

  .assert_that(all(cols %in% names(data)),
               msg = paste("Column(s)", paste(cols[which(cols %nin% names(data))], collapse = ","), "not found in", which.data, "data."))

  if (length(cols) == 1L) {
    if (all(is.na(data[[cols[1]]]))) {
      warning(paste("The column(s)", paste(cols, collapse = ", "), "are present in", which.data, "data, but all values are NAs."))
    }
  } else {
    if (all(vapply(data[cols], function(x) all(is.na(x)), logical(1)))) {
      warning(paste("The column(s)", paste(cols, collapse = ", "), "are present in", which.data, "data, but all values are NAs."))
    }
  }

  invisible(TRUE)
}

# public-facing assert below this line

#' Assert that an object is an ideal adaptor Stanfit, Staninput, or Stanfit input container.
#'
#' @description These helpers validate whether an object is an ideal adaptor Stanfit,
#'   Staninput, or Stanfit input container.
#' @param x Object to test.
#' @param verbose Logical. Whether to emit additional diagnostics while testing the object.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the assertion passes.
#' @rdname assert_ideal_adaptor_stanfit
#' @export
assert_ideal_adaptor_stanfit <- function(x, verbose = FALSE) {
  .assert_true(
    is.ideal_adaptor_stanfit(x, verbose = verbose),
    msg = paste(deparse(substitute(x)), "must be of class ideal_adaptor_stanfit")
  )
}

#' @rdname assert_ideal_adaptor_stanfit
#' @export
assert_staninput <- function(x, msg = NULL) {
  .assert_true(S7::S7_inherits(x, IdealAdaptorStaninput), msg = if (is.null(msg)) "Expected an IdealAdaptorStaninput object." else msg)
}

#' @rdname assert_ideal_adaptor_stanfit
#' @export
assert_stanfit_input <- function(x, msg = NULL) {
  .assert_true(S7::S7_inherits(x, IdealAdaptorStanfitInput), msg = if (is.null(msg)) "Expected an IdealAdaptorStanfitInput object." else msg)
}

# Deprecated assertions below this line

#' @description Deprecated. Use \code{\link{assert_ideal_adaptor_stanfit}} for Stanfit objects, \code{\link{assert_staninput}} for Staninput objects, or \code{\link{assert_stanfit_input}} for Stanfit input containers instead.
#' @deprecated Use \code{\link{assert_ideal_adaptor_stanfit}}, \code{\link{assert_staninput}}, or \code{\link{assert_stanfit_input}} instead.
#' @keywords internal
#' @param x Object to test.
#' @param category Name of the category field to inspect.
#' @param verbose Logical. Whether to emit additional diagnostics while testing the object.
#' @return Invisibly TRUE if the assertion passes.
#' @export
assert_MVG_ideal_observer = function(x, category = "category", verbose = F) {
  .assert_that(is.MVG_ideal_observer(x, category = category, verbose = verbose),
               msg = paste(deparse(substitute(x)), "must be an MVG_ideal_observer object."))
}

#' @description Deprecated. Use \code{\link{assert_ideal_adaptor_stanfit}} for Stanfit objects, \code{\link{assert_staninput}} for Staninput objects, or \code{\link{assert_stanfit_input}} for Stanfit input containers instead.
#' @deprecated Use \code{\link{assert_ideal_adaptor_stanfit}}, \code{\link{assert_staninput}}, or \code{\link{assert_stanfit_input}} instead.
#' @keywords internal
#' @param x Object to test.
#' @param category Name of the category field to inspect.
#' @param verbose Logical. Whether to emit additional diagnostics while testing the object.
#' @param strict Logical. Whether to require the object to be an NIW belief strictly, rather than allowing an ideal adaptor.
#' @return Invisibly TRUE if the assertion passes.
#' @export
assert_NIW_belief = function(x, category = "category", verbose = F, strict = F) {
  if (strict) {
    .assert_that(is.NIW_belief(x, category = category, verbose = verbose),
                 msg = paste(deparse(substitute(x)), "must be an NIW_belief object."))
  } else {
    .assert_that(
      any(
        is.NIW_belief(x, category = category, verbose = verbose),
        is.NIW_ideal_adaptor(x, category = category, verbose = verbose)),
      msg = paste(deparse(substitute(x)), "must be an NIW_belief or NIW_ideal_adaptor object."))
  }
}

#' @description Deprecated. Use \code{\link{assert_ideal_adaptor_stanfit}} for Stanfit objects, \code{\link{assert_staninput}} for Staninput objects, or \code{\link{assert_stanfit_input}} for Stanfit input containers instead.
#' @deprecated Use \code{\link{assert_ideal_adaptor_stanfit}}, \code{\link{assert_staninput}}, or \code{\link{assert_stanfit_input}} instead.
#' @keywords internal
#' @param x Object to test.
#' @param category Name of the category field to inspect.
#' @param verbose Logical. Whether to emit additional diagnostics while testing the object.
#' @return Invisibly TRUE if the assertion passes.
#' @export
assert_NIW_ideal_adaptor = function(x, category = "category", verbose = F) {
  .assert_that(is.NIW_ideal_adaptor(x, category = category, verbose = verbose),
               msg = paste(deparse(substitute(x)), "must be an NIW_ideal_adaptor object."))
}

