idahoans <- make_vowel_test_data()
legacy_niw <- tibble::tibble(
  category = c("A", "B"),
  m = list(c(F1 = 500), c(F1 = 700)),
  kappa = c(10, 10),
  nu = c(10, 10),
  S = list(matrix(100, 1, 1, dimnames = list("F1", "F1")), matrix(120, 1, 1, dimnames = list("F1", "F1"))))

test_that("deprecated NIW make wrappers delegate to S7 constructors", {
  expect_warning(x <- make_NIW_belief_from_data(idahoans, category = "vowel", cues = c("F1", "F2")), "deprecated")
  expect_true(S7::S7_inherits(x, MVBU_CategoryRepresentationTemplate))
  expect_warning(x <- make_NIW_ideal_adaptor_from_data(idahoans, category = "vowel", cues = c("F1", "F2")), "deprecated")
  expect_true(S7::S7_inherits(x, NIW_IdealAdaptor))
  expect_error(make_NIW_belief_from_data(idahoans, group = "vowel", category = "vowel", cues = c("F1", "F2")))
})

test_that("deprecated NIW lift wrapper delegates to S7 constructor", {
  expect_warning(model <- lift_NIW_belief_to_NIW_ideal_adaptor(legacy_niw), "deprecated")
  expect_true(S7::S7_inherits(model, NIW_IdealAdaptor))
  expect_error(lift_NIW_belief_to_NIW_ideal_adaptor(legacy_niw, group = "category"))
})
