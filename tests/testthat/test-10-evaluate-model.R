
my_model <- example_MVG_ideal_observer(5)
my_data <-
  crossing(cue1 = seq(-2, 2, .25), cue2 = seq(-2, 2, .25)) %>%
  mutate(cues = map2(cue1, cue2, ~ c(...))) %>%
  mutate(response = map(cues, ~ expect_warning(
    get_categorization_from_MVG_ideal_observer(
      x = .x,
      model = my_model,
      decision_rule = "sampling",
      simplify = TRUE
    ),
    "deprecated"
  )) %>% unlist())

test_that("evaluate_model - output check (method accuracy)", {
  accuracy_proportional <- expect_warning(
    evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "proportional", method = "accuracy"),
    "deprecated"
  )
  accuracy_criterion <- expect_warning(
    evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "criterion", method = "accuracy"),
    "deprecated"
  )

  expect_true(is.atomic(accuracy_proportional))
  expect_true(is.atomic(accuracy_criterion))
})

test_that("evaluate_model - output check (method likelihood-up-to-constant)", {
  likelihood_constant_proportional <- expect_warning(
    evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "proportional", method = "likelihood-up-to-constant"),
    "deprecated"
  )
  likelihood_constant_criterion <- expect_warning(
    evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "criterion", method = "likelihood-up-to-constant"),
    "deprecated"
  )

  expect_true(is.atomic(likelihood_constant_proportional))
  expect_true(is.infinite(likelihood_constant_criterion))
})

test_that("evaluate_model - output check (method likelihood)", {
  likelihood_proportional <- expect_warning(
    evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "proportional", method = "likelihood"),
    "deprecated"
  )
  likelihood_criterion <- expect_warning(
    evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "criterion", method = "likelihood"),
    "deprecated"
  )

  expect_true(is.atomic(likelihood_proportional))
  expect_true(is.infinite(likelihood_criterion))
})

test_that("evaluate_model - output check (method default)", {
  default_proportional <- expect_warning(
    evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "proportional"),
    "deprecated"
  )
  default_criterion <- expect_warning(
    evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "criterion"),
    "deprecated"
  )

  expect_true(is.atomic(default_proportional))
  expect_true(is.infinite(default_criterion))
})

