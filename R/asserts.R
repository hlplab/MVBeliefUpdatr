#' @export
assert_MVG_ideal_observer = function(x, category = "category", verbose = F) {
  assert_that(is.MVG_ideal_observer(x, category = category, verbose = verbose),
              msg = paste(deparse(substitute(x)), "must be an MVG_ideal_observer object."))
}

#' @export
assert_NIW_belief = function(x, category = "category", verbose = F, strict = F) {
  if (strict) {
    assert_that(is.NIW_belief(x, category = category, verbose = verbose),
                msg = paste(deparse(substitute(x)), "must be an NIW_belief object."))
  } else {
    assert_that(
      any(
        is.NIW_belief(x, category = category, verbose = verbose),
        is.NIW_ideal_adaptor(x, category = category, verbose = verbose)),
      msg = paste(deparse(substitute(x)), "must be an NIW_belief or NIW_ideal_adaptor object."))
  }
}

#' @export
assert_NIW_ideal_adaptor = function(x, category = "category", verbose = F) {
  assert_that(is.NIW_ideal_adaptor(x, category = category, verbose = verbose),
              msg = paste(deparse(substitute(x)), "must be an NIW_ideal_adaptor object."))
}

assert_ideal_adaptor_stanfit <- function(x, verbose = F) {
  assert_that(is.ideal_adaptor_stanfit(x, verbose = verbose),
              msg = paste(deparse(substitute(x)), "must be of class", "ideal_adaptor_stanfit"))
}


assert_cols_in_data <- function(data, cols, which.data = "the", scalar = T) {
  if (scalar) {
    assert_that(all(is_scalar_character(cols)),
                msg = paste0(paste(cols, collapse = ","), "must be a single column name."))
  } else {
    assert_that(all(is_character(cols)),
              msg = paste0(paste(cols, collapse = ","), "must be column name or vector of column names."))
  }

  assert_that(all(cols %in% names(data)),
              msg = paste("Column(s)", paste(cols[which(cols %nin% names(data))], collapse = ","), "not found in", which.data, "data." ))
  if (nrow(drop_na(data, all_of(cols))) == 0)
    warning(paste("The column(s)", paste(cols, collapse = ", "), "are present in", which.data, "data, but all values are NAs."))
}

assert_contains_draws <- function(x) {
  if (!contains_draws(x)) stop2(paste("", deparse(substitute(x)), "does not contain any samples."))
}
