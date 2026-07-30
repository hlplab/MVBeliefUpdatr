#' @import tidyverse
require(tidyverse)

context("class")

io <- example_MVG_ideal_observer()

test_that("is.MVBU_model - input", {
  expect_false(is.MVBU_model(NULL))
  expect_false(is.MVBU_model(1))
  expect_false(is.MVBU_model(list(1, 2)))
  expect_true(is.MVBU_model(example_MVG_ideal_observer(1)))
  expect_true(is.MVBU_model(example_MVG_ideal_observer(2)))
  expect_true(is.MVBU_model(example_exemplar_model(1)))
  expect_true(is.MVBU_model(example_exemplar_model(2)))
  expect_true(is.MVBU_model(example_exemplar_model(3)))
  expect_true(is.MVBU_model(example_NIW_ideal_adaptor(1)))
  expect_true(is.MVBU_model(example_NIW_ideal_adaptor(2)))
})


