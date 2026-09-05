idahoans <- make_vowel_test_data()
niw_template <- new_niw_category_representation_template_from_data(idahoans, category = "vowel", cues = c("F1", "F2"))

test_that("deprecated NIW make wrappers delegate to S7 constructors", {
  expect_warning(x <- make_NIW_belief_from_data(idahoans, category = "vowel", cues = c("F1", "F2")), "deprecated")
  expect_true(S7::S7_inherits(x, MVBU_CategoryRepresentationTemplate))
  expect_warning(x <- make_NIW_ideal_adaptor_from_data(idahoans, category = "vowel", cues = c("F1", "F2")), "deprecated")
  expect_true(S7::S7_inherits(x, NIW_IdealAdaptor))
  expect_error(make_NIW_belief_from_data(idahoans, group = "vowel", category = "vowel", cues = c("F1", "F2")))
})

test_that("deprecated NIW lift wrapper delegates to S7 constructor", {
  expect_warning(model <- lift_NIW_belief_to_NIW_ideal_adaptor(niw_template), "deprecated")
  expect_true(S7::S7_inherits(model, NIW_IdealAdaptor))
  expect_error(lift_NIW_belief_to_NIW_ideal_adaptor(niw_template, group = "category"))
})
