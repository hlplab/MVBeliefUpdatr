#' @include internal-is.R internal-utils-imported.R
NULL

# Helper to format default assertion messages.
.default_assert_msg <- function(x, expectation) {
  var_name <- deparse(substitute(x))
  paste0("Expected ", var_name, " to be ", expectation, ".")
}

# Helper to interpret a single assertion condition.
.assert_condition_ok <- function(cond) {
  if (is.logical(cond) && length(cond) == 1L) {
    !is.na(cond) && isTRUE(cond)
  } else if (is.logical(cond)) {
    all(cond, na.rm = TRUE)
  } else {
    isTRUE(cond)
  }
}

# Helper to evaluate a quoted assertion condition safely.
.eval_assert_condition <- function(expr, env = parent.frame()) {
  tryCatch(
    .assert_condition_ok(eval(expr, envir = env)),
    error = function(e) FALSE
  )
}

#' Internal assertion helper
#'
#' @param ... Conditions to evaluate.
#' @param msg Optional error message.
#' @return Invisibly TRUE if all conditions pass.
#' @noRd
.assert_that <- function(..., msg = NULL) {
  conditions <- match.call(expand.dots = FALSE)$...
  if (length(conditions) == 0L) {
    return(invisible(TRUE))
  }

  env <- parent.frame()
  cond_ok <- vapply(conditions, function(cond) .eval_assert_condition(cond, env = env), logical(1))

  if (!all(cond_ok)) {
    .stop(if (is.null(msg)) "Assertion failed." else msg)
  }

  invisible(TRUE)
}

#' Internal assertion helper
#'
#' @param ... Conditions to evaluate.
#' @param msg Optional error message.
#' @return Invisibly TRUE if all conditions pass.
#' @noRd
.assert_all <- function(..., msg = NULL) {
  .assert_that(..., msg = msg)
}

#' Internal assertion helper
#'
#' @param ... Conditions to evaluate.
#' @param msg Optional error message.
#' @return Invisibly TRUE if any condition passes.
#' @noRd
.assert_any <- function(..., msg = NULL) {
  conditions <- match.call(expand.dots = FALSE)$...
  if (length(conditions) == 0L) {
    return(invisible(TRUE))
  }

  env <- parent.frame()
  cond_ok <- vapply(conditions, function(cond) .eval_assert_condition(cond, env = env), logical(1))

  if (!any(cond_ok)) {
    .stop(if (is.null(msg)) "Assertion failed." else msg)
  }

  invisible(TRUE)
}

#' Internal assertion helper
#'
#' @param cond Condition to evaluate.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the condition passes.
#' @noRd
.assert_true <- function(cond, msg = NULL) {
  .assert_that(cond, msg = if (is.null(msg)) "Assertion failed." else msg)
}

#' Internal assertion helper
#'
#' @param cond Condition to evaluate.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the condition fails.
#' @noRd
.assert_false <- function(cond, msg = NULL) {
  .assert_true(!cond, msg = if (is.null(msg)) "Assertion failed." else msg)
}

# Helper to describe a basic value of a given type.
.assert_type <- function(x, predicate, expectation, msg = NULL) {
  .assert_true(
    predicate(x),
    msg = if (is.null(msg)) .default_assert_msg(x, expectation) else msg
  )
}

# Helper to describe a basic scalar value of a given type.
.assert_scalar_type <- function(x, predicate, expectation, msg = NULL) {
  .assert_true(
    predicate(x),
    msg = if (is.null(msg)) .default_assert_msg(x, expectation) else msg
  )
}

# Helper to describe an optional value that may be missing, NULL, or a given type.
.assert_optional_type <- function(x, predicate, expectation, msg = NULL) {
  .assert_true(
    missing(x) || is.null(x) || predicate(x),
    msg = if (is.null(msg)) .default_assert_msg(x, expectation) else msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is missing, NULL, or a character value.
#' @noRd
.assert_optional_character <- function(x, msg = NULL) {
  .assert_optional_type(
    x,
    predicate = is.character,
    expectation = "missing, NULL, or a character value",
    msg = msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is missing, NULL, or a numeric value.
#' @noRd
.assert_optional_numeric <- function(x, msg = NULL) {
  .assert_optional_type(
    x,
    predicate = is.numeric,
    expectation = "missing, NULL, or a numeric value",
    msg = msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is missing, NULL, or a factor.
#' @noRd
.assert_optional_factor <- function(x, msg = NULL) {
  .assert_optional_type(
    x,
    predicate = is.factor,
    expectation = "missing, NULL, or a factor",
    msg = msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is missing, NULL, or a list.
#' @noRd
.assert_optional_list <- function(x, msg = NULL) {
  .assert_optional_type(
    x,
    predicate = is.list,
    expectation = "missing, NULL, or a list",
    msg = msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is missing, NULL, or a matrix.
#' @noRd
.assert_optional_matrix <- function(x, msg = NULL) {
  .assert_optional_type(
    x,
    predicate = is.matrix,
    expectation = "missing, NULL, or a matrix",
    msg = msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a character vector.
#' @noRd
.assert_character <- function(x, msg = NULL) {
  .assert_type(
    x,
    predicate = is.character,
    expectation = "a character vector",
    msg = msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a numeric vector.
#' @noRd
.assert_numeric <- function(x, msg = NULL) {
  .assert_type(
    x,
    predicate = is.numeric,
    expectation = "a numeric vector",
    msg = msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a logical vector.
#' @noRd
.assert_logical <- function(x, msg = NULL) {
  .assert_type(
    x,
    predicate = is.logical,
    expectation = "a logical vector",
    msg = msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a matrix.
#' @noRd
.assert_matrix <- function(x, msg = NULL) {
  .assert_type(
    x,
    predicate = is.matrix,
    expectation = "a matrix",
    msg = msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a factor.
#' @noRd
.assert_factor <- function(x, msg = NULL) {
  .assert_type(
    x,
    predicate = is.factor,
    expectation = "a factor",
    msg = msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a list.
#' @noRd
.assert_list <- function(x, msg = NULL) {
  .assert_type(
    x,
    predicate = is.list,
    expectation = "a list",
    msg = msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a character scalar.
#' @noRd
.assert_character_scalar <- function(x, msg = NULL) {
  .assert_scalar_type(
    x,
    predicate = .is_scalar_character,
    expectation = "a scalar character",
    msg = msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a numeric scalar.
#' @noRd
.assert_numeric_scalar <- function(x, msg = NULL) {
  .assert_scalar_type(
    x,
    predicate = .is_scalar_numeric,
    expectation = "a scalar numeric",
    msg = msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a logical scalar.
#' @noRd
.assert_logical_scalar <- function(x, msg = NULL) {
  .assert_scalar_type(
    x,
    predicate = .is_scalar_logical,
    expectation = "a scalar logical",
    msg = msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a scalar factor.
#' @noRd
.assert_factor_scalar <- function(x, msg = NULL) {
  .assert_scalar_type(
    x,
    predicate = .is_scalar_factor,
    expectation = "a scalar factor",
    msg = msg
  )
}

# Helper to describe a required value that must be non-NA and of a given type.
.assert_non_NA_type <- function(x, predicate, expectation, msg = NULL) {
  .assert_true(
    !missing(x) && !is.null(x) && predicate(x) && !anyNA(x) && length(x) > 0L,
    msg = if (is.null(msg)) {
      .default_assert_msg(x, paste0("a non-NA ", expectation))
    } else {
      msg
    }
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a non-NA character vector.
#' @noRd
.assert_non_NA_character <- function(x, msg = NULL) {
  .assert_non_NA_type(
    x,
    predicate = is.character,
    expectation = "character",
    msg = if (is.null(msg)) .default_assert_msg(x, "a non-NA character") else msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a non-NA numeric vector.
#' @noRd
.assert_non_NA_numeric <- function(x, msg = NULL) {
  .assert_non_NA_type(
    x,
    predicate = is.numeric,
    expectation = "numeric",
    msg = if (is.null(msg)) .default_assert_msg(x, "a non-NA numeric") else msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a non-NA logical vector.
#' @noRd
.assert_non_NA_logical <- function(x, msg = NULL) {
  .assert_non_NA_type(
    x,
    predicate = is.logical,
    expectation = "logical",
    msg = if (is.null(msg)) .default_assert_msg(x, "a non-NA logical") else msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a non-NA factor.
#' @noRd
.assert_non_NA_factor <- function(x, msg = NULL) {
  .assert_non_NA_type(
    x,
    predicate = is.factor,
    expectation = "factor",
    msg = if (is.null(msg)) .default_assert_msg(x, "a non-NA factor") else msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a non-NA list.
#' @noRd
.assert_non_NA_list <- function(x, msg = NULL) {
  .assert_non_NA_type(
    x,
    predicate = is.list,
    expectation = "list",
    msg = if (is.null(msg)) .default_assert_msg(x, "a non-NA list") else msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a non-NA matrix.
#' @noRd
.assert_non_NA_matrix <- function(x, msg = NULL) {
  .assert_non_NA_type(
    x,
    predicate = is.matrix,
    expectation = "matrix",
    msg = if (is.null(msg)) .default_assert_msg(x, "a non-NA matrix") else msg
  )
}


# Helper to describe a required non-NA scalar value of a given type.
.assert_non_NA_scalar_type <- function(x, predicate, expectation, msg = NULL) {
  .assert_true(
    predicate(x) && length(x) == 1L && !is.na(x[1]),
    msg = if (is.null(msg)) .default_assert_msg(x, expectation) else msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a non-NA, non-empty character scalar.
#' @noRd
.assert_non_NA_scalar_character <- function(x, msg = NULL) {
  .assert_true(
    .is_non_empty_scalar_character(x),
    msg = if (is.null(msg)) .default_assert_msg(x, "a non-NA, non-empty character scalar") else msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a non-NA numeric scalar.
#' @noRd
.assert_non_NA_scalar_numeric <- function(x, msg = NULL) {
  .assert_non_NA_scalar_type(
    x,
    predicate = .is_non_NA_scalar_numeric,
    expectation = "a non-NA numeric scalar",
    msg = msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a non-NA logical scalar.
#' @noRd
.assert_non_NA_scalar_logical <- function(x, msg = NULL) {
  .assert_non_NA_scalar_type(
    x,
    predicate = .is_non_NA_scalar_logical,
    expectation = "a non-NA logical scalar",
    msg = msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a non-NA scalar factor.
#' @noRd
.assert_non_NA_scalar_factor <- function(x, msg = NULL) {
  .assert_non_NA_scalar_type(
    x,
    predicate = .is_non_NA_scalar_factor,
    expectation = "a non-NA scalar factor",
    msg = msg
  )
}

# Helper to describe an optional value that may be missing, NULL, or a non-NA value of a given type.
.assert_missing_null_or_non_NA_type <- function(x, predicate, expectation, msg = NULL) {
  .assert_true(
    missing(x) || is.null(x) || (predicate(x) && !anyNA(x) && length(x) > 0L),
    msg = if (is.null(msg)) {
      .default_assert_msg(x, paste0("missing, NULL, or a non-empty ", expectation))
    } else {
      msg
    }
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is missing, NULL, or a non-empty character vector.
#' @noRd
.assert_missing_null_or_non_NA_character <- function(x, msg = NULL) {
  .assert_missing_null_or_non_NA_type(
    x,
    predicate = is.character,
    expectation = "character",
    msg = if (is.null(msg)) .default_assert_msg(x, "missing, NULL, or a non-empty character") else msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is missing, NULL, or a non-empty numeric vector.
#' @noRd
.assert_missing_null_or_non_NA_numeric <- function(x, msg = NULL) {
  .assert_missing_null_or_non_NA_type(
    x,
    predicate = is.numeric,
    expectation = "numeric",
    msg = if (is.null(msg)) .default_assert_msg(x, "missing, NULL, or a non-empty numeric") else msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is missing, NULL, or a non-empty factor.
#' @noRd
.assert_missing_null_or_non_NA_factor <- function(x, msg = NULL) {
  .assert_missing_null_or_non_NA_type(
    x,
    predicate = is.factor,
    expectation = "factor",
    msg = if (is.null(msg)) .default_assert_msg(x, "missing, NULL, or a non-empty factor") else msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is missing, NULL, or a non-empty matrix.
#' @noRd
.assert_missing_null_or_non_NA_matrix <- function(x, msg = NULL) {
  .assert_missing_null_or_non_NA_type(
    x,
    predicate = is.matrix,
    expectation = "matrix",
    msg = if (is.null(msg)) .default_assert_msg(x, "missing, NULL, or a non-empty matrix") else msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is non-negative.
#' @noRd
.assert_non_negative <- function(x, msg = NULL) {
  .assert_true(
    is.numeric(x) && length(x) >= 1L && all(!is.na(x)) && all(x >= 0),
    msg = if (is.null(msg)) .default_assert_msg(x, "non-negative") else msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param lower Lower bound.
#' @param upper Upper bound.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object lies between the bounds.
#' @noRd
.assert_between <- function(x, lower, upper, msg = NULL) {
  .assert_true(
    is.numeric(x) && length(x) >= 1L && all(!is.na(x)) && all(x >= lower) && all(x <= upper),
    msg = if (is.null(msg)) .default_assert_msg(x, paste0("between ", lower, " and ", upper)) else msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param choices Allowed values.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is one of the allowed choices.
#' @noRd
.assert_one_of <- function(x, choices, msg = NULL) {
  .assert_true(
    !missing(x) && !is.null(x) && length(x) >= 1L && !anyNA(x) && all(x %in% choices),
    msg = if (is.null(msg)) .default_assert_msg(x, paste0("one of ", paste(choices, collapse = ", "))) else msg
  )
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the object is a data.frame-like object.
#' @noRd
.assert_data_frame_like <- function(x, msg = NULL) {
  .assert_true(is.data.frame(x) || is_tibble(x), msg = if (is.null(msg)) "Expected a data.frame or tibble." else msg)
}

#' Internal assertion helper
#'
#' @param x Object to test.
#' @return Invisibly TRUE if the object contains draws.
#' @noRd
.assert_contains_draws <- function(x) {
  if (!.contains_draws(x)) .stop(paste("", deparse(substitute(x)), "does not contain any samples."))
}

#' Internal assertion helper
#'
#' @param data Data object to inspect.
#' @param cols Column names to inspect.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the columns are present.
#' @noRd
.assert_data_contains_cols <- function(data, cols, msg = NULL) {
  .assert_character(cols)
  .assert_data_frame_like(data)
  data_name <- deparse(substitute(data))
  if (is.null(msg)) {
    msg <- sprintf("Expected data %s to contain column(s): %s", data_name, paste(setdiff(cols, names(data)), collapse = ", "))
  }
  .assert_true(all(cols %in% names(data)), msg = msg)
  invisible(TRUE)
}