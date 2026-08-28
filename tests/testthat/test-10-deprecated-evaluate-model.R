
make_evaluate_model_test_data <- function() {
  my_model <- suppressMessages(suppressWarnings(example_mvg_ideal_observer(n_cues = 2)))
  my_data <- tidyr::crossing(cue1 = seq(-2, 2, .25), cue2 = seq(-2, 2, .25))
  my_data$cues <- purrr::map2(my_data$cue1, my_data$cue2, ~ c(...))
  my_data$response <- vapply(my_data$cues, function(x) {
    as.character(
      suppressWarnings(
        get_categorization_from_MVG_ideal_observer(
          x = x,
          model = my_model,
          decision_rule = "sampling",
          simplify = TRUE
        )))
  }, character(1))

  list(model = my_model, data = my_data)
}

test_that("evaluate_model - output check (method accuracy)", {
  setup <- make_evaluate_model_test_data()
  my_model <- setup$model
  my_data <- setup$data
  accuracy_proportional <- NULL
  expect_warning(
    {
      accuracy_proportional <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "proportional", method = "accuracy")
    },
    "deprecated"
  )
  accuracy_criterion <- NULL
  expect_warning(
    {
      accuracy_criterion <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "criterion", method = "accuracy")
    },
    "deprecated"
  )

  expect_true(is.atomic(accuracy_proportional))
  expect_true(is.atomic(accuracy_criterion))
})

test_that("evaluate_model - output check (method likelihood-up-to-constant)", {
  setup <- make_evaluate_model_test_data()
  my_model <- setup$model
  my_data <- setup$data
  likelihood_constant_proportional <- NULL
  expect_warning(
    {
      likelihood_constant_proportional <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "proportional", method = "likelihood-up-to-constant")
    },
    "deprecated"
  )
  likelihood_constant_criterion <- NULL
  expect_warning(
    {
      likelihood_constant_criterion <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "criterion", method = "likelihood-up-to-constant")
    },
    "deprecated"
  )

  expect_true(is.atomic(likelihood_constant_proportional))
  expect_true(is.infinite(likelihood_constant_criterion))
})

test_that("evaluate_model - output check (method likelihood)", {
  setup <- make_evaluate_model_test_data()
  my_model <- setup$model
  my_data <- setup$data
  likelihood_proportional <- NULL
  expect_warning(
    {
      likelihood_proportional <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "proportional", method = "likelihood")
    },
    "deprecated"
  )
  likelihood_criterion <- NULL
  expect_warning(
    {
      likelihood_criterion <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "criterion", method = "likelihood")
    },
    "deprecated"
  )

  expect_true(is.atomic(likelihood_proportional))
  expect_true(is.infinite(likelihood_criterion))
})

test_that("evaluate_model - output check (method default)", {
  setup <- make_evaluate_model_test_data()
  my_model <- setup$model
  my_data <- setup$data
  default_proportional <- NULL
  expect_warning(
    {
      default_proportional <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "proportional")
    },
    "deprecated"
  )
  default_criterion <- NULL
  expect_warning(
    {
      default_criterion <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "criterion")
    },
    "deprecated"
  )

  expect_true(is.atomic(default_proportional))
  expect_true(is.infinite(default_criterion))
})
