
make_evaluate_model_test_data <- function() {
  my_model <- example_mvg_ideal_observer(n_cues = 2)
  my_data <- tidyr::crossing(cue1 = seq(-2, 2, 0.5), cue2 = seq(-2, 2, 0.5))
  my_data$cues <- purrr::map2(my_data$cue1, my_data$cue2, ~ c(.x, .y))
  cat_res <- categorize(my_model, as.matrix(my_data[, c("cue1", "cue2")]), decision_rule = "sampling")
  my_data$response <- cat_res$category

  list(model = my_model, data = my_data)
}



test_that("evaluate_model - output check (method accuracy)", {
  setup <- make_evaluate_model_test_data()
  my_model <- setup$model
  my_data <- setup$data

  accuracy_proportional <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "proportional", method = "accuracy")
  accuracy_criterion <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "criterion", method = "accuracy")

  expect_true(is.atomic(accuracy_proportional))
  expect_true(accuracy_proportional >= 0 && accuracy_proportional <= 1)
  expect_true(is.atomic(accuracy_criterion))
  expect_true(accuracy_criterion >= 0 && accuracy_criterion <= 1)

  # Check return_by_x
  acc_by_x <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, method = "accuracy", return_by_x = TRUE)
  expect_s3_class(acc_by_x, "data.frame")
  expect_true("accuracy" %in% names(acc_by_x))
})

test_that("evaluate_model - output check (method log_lik)", {
  setup <- make_evaluate_model_test_data()
  my_model <- setup$model
  my_data <- setup$data

  log_lik_proportional <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "proportional", method = "log_lik")
  log_lik_criterion <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "criterion", method = "log_lik")

  expect_true(is.atomic(log_lik_proportional))
  expect_true(is.finite(log_lik_proportional))
  expect_true(is.atomic(log_lik_criterion))

  # Check return_by_x
  ll_by_x <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, method = "log_lik", return_by_x = TRUE)
  expect_s3_class(ll_by_x, "data.frame")
  expect_true("log_lik" %in% names(ll_by_x))
})

test_that("evaluate_model - output check (method log_lik_permutation_constant)", {
  setup <- make_evaluate_model_test_data()
  my_model <- setup$model
  my_data <- setup$data

  const_val <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, method = "log_lik_permutation_constant")
  expect_true(is.atomic(const_val))
  expect_true(is.finite(const_val))
  expect_true(const_val >= 0)

  # Check return_by_x
  const_by_x <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, method = "log_lik_permutation_constant", return_by_x = TRUE)
  expect_s3_class(const_by_x, "data.frame")
  expect_true("log_lik_permutation_constant" %in% names(const_by_x))
})

test_that("evaluate_model - output check (method likelihood-up-to-constant deprecated alias)", {
  setup <- make_evaluate_model_test_data()
  my_model <- setup$model
  my_data <- setup$data

  expect_warning(
    likelihood_constant_proportional <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "proportional", method = "likelihood-up-to-constant"),
    class = "lifecycle_warning_deprecated"
  )
  expect_true(is.atomic(likelihood_constant_proportional))
  expect_true(is.finite(likelihood_constant_proportional))

  # Check return_by_x
  ll_by_x <- suppressWarnings(evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, method = "likelihood-up-to-constant", return_by_x = TRUE))
  expect_s3_class(ll_by_x, "data.frame")
  expect_true("log_likelihood" %in% names(ll_by_x))
})

test_that("evaluate_model - output check (method default)", {
  setup <- make_evaluate_model_test_data()
  my_model <- setup$model
  my_data <- setup$data

  default_proportional <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "proportional")
  default_criterion <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "criterion")

  expect_true(is.atomic(default_proportional))
  expect_true(is.atomic(default_criterion))

  # Default method should equal method = "log_lik"
  ll_explicit <- evaluate_model(model = my_model, x = my_data$cues, response_category = my_data$response, decision_rule = "proportional", method = "log_lik")
  expect_equal(default_proportional, ll_explicit)
})

test_that("evaluate_model on MVBU_Stanfit works", {
  skip_if_not_installed("rstan")
  model_name <- "minimal-nix_ideal_adaptor-NIX"
  fit <- tryCatch(load_ideal_adaptor_fit_model(model_name), error = function(e) NULL)
  if (!is.null(fit)) {
    # If no test data embedded, expecting informative error
    expect_error(evaluate_model(fit), "No test data found")

    # With explicit x and response_category
    x_test <- matrix(c(0.2, 0.8), ncol = 1)
    resp_test <- c("A", "B")
    eval_res <- evaluate_model(fit, x = x_test, response_category = resp_test)
    expect_true(is.atomic(eval_res))
    expect_true(is.finite(eval_res))
    expect_true(eval_res < 0)

    eval_acc <- evaluate_model(fit, x = x_test, response_category = resp_test, method = "accuracy")
    expect_true(is.atomic(eval_acc))
    expect_true(eval_acc >= 0 && eval_acc <= 1)
  }
})

