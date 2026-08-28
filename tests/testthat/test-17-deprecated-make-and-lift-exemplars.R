idahoans <- make_vowel_test_data()
legacy_exemplars <- tibble::tibble(
  category = c("A", "B"),
  exemplars = list(matrix(c(490, 500, 510), ncol = 1, dimnames = list(NULL, "F1")), matrix(c(690, 700, 710), ncol = 1, dimnames = list(NULL, "F1"))))

test_that("deprecated Exemplar make wrappers delegate to S7 constructors", {
  expect_warning(x <- make_exemplars_from_data(idahoans, category = "vowel", cues = c("F1", "F2")), "deprecated")
  expect_true(S7::S7_inherits(x, MVBU_CategoryRepresentationTemplate))
  expect_warning(x <- make_exemplar_model_from_data(idahoans, category = "vowel", cues = c("F1", "F2")), "deprecated")
  expect_true(S7::S7_inherits(x, Exemplar_Model))
  expect_error(make_exemplars_from_data(idahoans, group = "vowel", category = "vowel", cues = c("F1", "F2")))
  expect_warning(make_exemplars_from_data(idahoans, category = "vowel", cues = c("F1", "F2"), sim_function = function(x, y) 1), "no longer supported")
})

test_that("deprecated Exemplar lift wrapper delegates to S7 constructor", {
  expect_warning(model <- lift_exemplars_to_exemplar_model(legacy_exemplars), "deprecated")
  expect_true(S7::S7_inherits(model, Exemplar_Model))
  expect_error(lift_exemplars_to_exemplar_model(legacy_exemplars, group = "category"))
})
