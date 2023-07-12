context("get information from MVG - INITIAL TEST ONLY")

library(tidyverse)
library(curl)
if (has_internet()) remotes::install_github("joeystanley/joeysvowels")
library(joeysvowels)
data("idahoans")

model <- make_MVG_ideal_observer_from_data(idahoans, category = "vowel", cues = c("F1", "F2"))
x.1 <- idahoans %>% mutate(x = map(F1, ~ c(...))) %>% pull(x)
x.2 <- idahoans %>% mutate(x = map2(F1, F2, ~ c(...))) %>% pull(x)
x.3 <- idahoans %>% mutate(x = pmap(.l = list(F1, F2, F3), ~ c(...))) %>% pull(x)

test_that("Get likelihood - input check x", {
  expect_error(get_likelihood_from_MVG(x = x.1, model = model))
  expect_error(get_likelihood_from_MVG(x = x.3, model = model))
  expect_no_error(get_likelihood_from_MVG(x = x.2, model = model))
})


test_that("Get likelihood - input check x", {
  expect_no_error(get_categorization_from_MVG_ideal_observer(x = x.2, model = model,
                                                             noise_treatment = "no_noise",
                                                             lapse_treatment = "no_lapses",
                                                             decision_rule = "sampling"))
})
