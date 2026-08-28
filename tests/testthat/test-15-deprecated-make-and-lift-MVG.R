idahoans <- make_vowel_test_data()
legacy_mvg <- tibble::tibble(
  category = c("A", "B"),
  mu = list(c(F1 = 500), c(F1 = 700)),
  Sigma = list(matrix(100, 1, 1, dimnames = list("F1", "F1")), matrix(120, 1, 1, dimnames = list("F1", "F1"))))

test_that("deprecated MVG make wrappers delegate to S7 constructors", {
  expect_warning(x <- make_MVG_from_data(idahoans, category = "vowel", cues = c("F1", "F2")), "deprecated")
  expect_true(S7::S7_inherits(x, MVBU_CategoryRepresentationTemplate))
  expect_warning(x <- make_MVG_ideal_observer_from_data(idahoans, category = "vowel", cues = c("F1", "F2")), "deprecated")
  expect_true(S7::S7_inherits(x, MVG_IdealObserver))
  expect_error(make_MVG_from_data(idahoans, group = "vowel", category = "vowel", cues = c("F1", "F2")))
})

test_that("deprecated MVG lift wrapper delegates to S7 constructor", {
  expect_warning(model <- lift_MVG_to_MVG_ideal_observer(legacy_mvg), "deprecated")
  expect_true(S7::S7_inherits(model, MVG_IdealObserver))
  expect_error(lift_MVG_to_MVG_ideal_observer(legacy_mvg, group = "category"))
})
