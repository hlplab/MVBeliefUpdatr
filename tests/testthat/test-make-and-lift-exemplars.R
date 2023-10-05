context("make and lift exemplars and exemplar models")

library(curl)
if (has_internet()) remotes::install_github("joeystanley/joeysvowels")
library(joeysvowels)
data("idahoans")

test_that("make exemplars", {
  expect_true(is_tibble(make_exemplars_from_data(idahoans, category = "vowel", cues = c("F1"))))
  expect_true(is_tibble(make_exemplars_from_data(idahoans, category = "vowel", cues = c("F1", "F2"))))
  expect_true(is_tibble(make_exemplars_from_data(idahoans, category = "vowel", cues = c("F1", "F2", "F3"))))
})

x <- make_exemplars_from_data(idahoans, category = "vowel", cues = c("F1", "F2"))
test_that("recognize exemplars", {
  expect_true(is.exemplars(x))
  expect_false(is.exemplar_model(x))
  expect_false(is.MVG(x))
  expect_false(is.MVG_ideal_observer(x))
  expect_false(is.NIW_belief(x))
  expect_false(is.NIW_ideal_adaptor(x))
})

test_that("get category and cue information from exemplars", {
  expect_true(is.character(get_cue_labels_from_model(x)))
  expect_true(any(is.character(get_category_labels_from_model(x)),
                  is.factor(get_category_labels_from_model(x))))
  expect_true(is.numeric(get_nlevels_of_category_labels_from_model(x)))
})

test_that("lift exemplar model", {
  expect_true(is_tibble(lift_exemplars_to_exemplar_model(x)))
})

test_that("make exemplar model", {
  expect_true(is_tibble(make_exemplar_model_from_data(idahoans, category = "vowel", cues = c("F1"))))
  expect_true(is_tibble(make_exemplar_model_from_data(idahoans, category = "vowel", cues = c("F1", "F2"))))
  expect_true(is_tibble(make_exemplar_model_from_data(idahoans, category = "vowel", cues = c("F1", "F2", "F3"))))
  expect_true(is_tibble(make_exemplar_model_from_data(idahoans, category = "vowel", cues = c("F1"), prior = rep(1/11, 11), lapse_rate = .05, lapse_bias = rep(1/11, 11))))
  expect_message(make_exemplar_model_from_data(idahoans, category = "vowel", cues = c("F1")))
  expect_message(is_tibble(make_exemplar_model_from_data(idahoans, category = "vowel", cues = c("F1"), prior = rep(1/11, 11), lapse_rate = .05, lapse_bias = rep(1/11, 11))))
  expect_no_message(is_tibble(make_exemplar_model_from_data(idahoans %>% filter(vowel %in% c("AA", "OW")), category = "vowel", cues = c("F1"), prior = c("AA" = 1/2, "OW" = 1/2), lapse_rate = .05, lapse_bias = c("AA" = 1/2, "OW" = 1/2))))
  expect_error(make_exemplar_model_from_data(idahoans, category = "vowel", cues = c("F1", "F2"), prior = rep(0, 10)))
  expect_error(make_exemplar_model_from_data(idahoans, category = "vowel", cues = c("F1", "F2"), lapse_rate = rep(0, 10)))
  expect_error(make_exemplar_model_from_data(idahoans, category = "vowel", cues = c("F1", "F2"), lapse_bias = rep(0, 10)))
})

x <- make_exemplar_model_from_data(idahoans, category = "vowel", cues = c("F1", "F2"))
test_that("recognize exemplar model", {
  expect_true(is.exemplars(x))
  expect_true(is.exemplar_model(x))
  expect_false(is.MVG(x))
  expect_false(is.MVG_ideal_observer(x))
  expect_false(is.NIW_belief(x))
  expect_false(is.NIW_ideal_adaptor(x))
})

test_that("get category and cue information from exemplar model", {
  expect_true(is.character(get_cue_labels_from_model(x)))
  expect_true(any(is.character(get_category_labels_from_model(x)),
                  is.factor(get_category_labels_from_model(x))))
  expect_true(is.numeric(get_nlevels_of_category_labels_from_model(x)))
})

test_that("get lapse rate, bias, and category prior from exemplar model", {
  expect_true(is.numeric(get_priors_from_model(x)))
  expect_true(is.numeric(get_lapse_biases_from_model(x)))
  expect_true(is.numeric(get_lapse_rate_from_model(x)))
})


