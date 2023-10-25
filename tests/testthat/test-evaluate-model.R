context("get information from model - INITIAL TEST ONLY")

my_model <- example_MVG_ideal_observer(1)
my_data <-
  crossing(cue1 = seq(-2, 2, .25), cue2 = seq(-2, 2, .25)) %>%
  mutate(cues = map2(cue1, cue2, ~ c(...))) %>%
  mutate(response = map(cues, ~ get_categorization_from_MVG_ideal_observer(x = .x, model = my_model, decision_rule = "sampling", simplify = T)) %>% unlist())

test_that("evaluate_model - out check", {
  expect_true(is_scalar_atomic(evaluate_model(my_model, x = my_data$cues, correct_category = my_data$response, decision_rule = "proportional")))
  expect_true(is_scalar_atomic(evaluate_model(my_model, x = my_data$cues, correct_category = my_data$response, decision_rule = "proportional", method = "accuracy")))
  expect_true(is.infinite(evaluate_model(my_model, x = my_data$cues, correct_category = my_data$response, decision_rule = "criterion")))
  expect_true(is_scalar_atomic(evaluate_model(my_model, x = my_data$cues, correct_category = my_data$response, decision_rule = "criterion", method = "accuracy")))
})






