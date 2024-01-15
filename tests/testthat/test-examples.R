context("make example objects")

test_that("creating examples", {
  expect_no_error(example_exemplar_model(1))
  expect_no_error(example_MVG_ideal_observer(1))
  expect_no_error(example_MVG_ideal_observer(2))
  expect_no_error(example_NIW_ideal_adaptor(1))
})
