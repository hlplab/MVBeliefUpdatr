assert_NIW_belief = function(x, category = "category") {
  assert_that(is.NIW_belief(x, category = category),
              msg = "x must be an NIW_belief object.")
}

assert_NIW_ideal_adaptor_stanfit = function(x) {
  assert_that(is.NIW_ideal_adaptor_stanfit(x),
              msg = paste0("x must be of class ", new_stanfit_class_name))
}

assert_cols_in_data = function(data, cols, which.data = "the", scalar = T) {
  if (scalar)
    assert_that(all(is_scalar_character(cols)),
                msg = paste0(paste(cols, collapse = ","), "must be a single column name.")) else
    assert_that(all(is_character(cols)),
              msg = paste0(paste(cols, collapse = ","), "must be column name or vector of column names."))
  assert_that(all(cols %in% names(data)),
              msg = paste("Column(s)", paste(cues[which(cols %nin% names(data))], collapse = ","), "not found in", which.data, "data." ))
  if (nrow(drop_na(data, cols)) == 0)
    warning(paste("The column(s)", paste(cols, collapse = ", "), "are present in", which.data, "data, but all values are NAs."))
}
