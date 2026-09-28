make_evaluate_model_test_data <- function() {
  my_model <- example_mvg_ideal_observer(
    n_cues = 2,
    categories = c("/d/", "/t/")
  )
  exp_data <- example_exposure_test_data(
    model_family = "NIW",
    n_cues = 2L,
    categories = c("/d/", "/t/"),
    seed = 42L
  )
  test_data <- exp_data$test
  cues_list <- purrr::map2(
    test_data$VOT,
    test_data$f0_semitones,
    ~ c(.x, .y)
  )
  my_data <- list(
    cues = cues_list,
    response = as.character(test_data$response_category)
  )
  list(model = my_model, data = my_data)
}

test_that("evaluate_model - output check (method accuracy)", {
  setup <- make_evaluate_model_test_data()
  my_model <- setup$model
  my_data <- setup$data

  accuracy_proportional <- evaluate_model(
    model = my_model,
    x = my_data$cues,
    response_category = my_data$response,
    decision_rule = "proportional",
    method = "accuracy"
  )
  accuracy_criterion <- evaluate_model(
    model = my_model,
    x = my_data$cues,
    response_category = my_data$response,
    decision_rule = "criterion",
    method = "accuracy"
  )

  expect_true(is.atomic(accuracy_proportional))
  expect_true(accuracy_proportional >= 0 && accuracy_proportional <= 1)
  expect_true(is.atomic(accuracy_criterion))
  expect_true(accuracy_criterion >= 0 && accuracy_criterion <= 1)

  acc_by_x <- evaluate_model(
    model = my_model,
    x = my_data$cues,
    response_category = my_data$response,
    method = "accuracy",
    return_by_x = TRUE
  )
  expect_s3_class(acc_by_x, "data.frame")
  expect_true("accuracy" %in% names(acc_by_x))
})

test_that("evaluate_model - output check (method log_lik)", {
  setup <- make_evaluate_model_test_data()
  my_model <- setup$model
  my_data <- setup$data

  log_lik_proportional <- evaluate_model(
    model = my_model,
    x = my_data$cues,
    response_category = my_data$response,
    decision_rule = "proportional",
    method = "log_lik"
  )
  log_lik_criterion <- evaluate_model(
    model = my_model,
    x = my_data$cues,
    response_category = my_data$response,
    decision_rule = "criterion",
    method = "log_lik"
  )

  expect_true(is.atomic(log_lik_proportional))
  expect_true(is.finite(log_lik_proportional))
  expect_true(is.atomic(log_lik_criterion))

  ll_by_x <- evaluate_model(
    model = my_model,
    x = my_data$cues,
    response_category = my_data$response,
    method = "log_lik",
    return_by_x = TRUE
  )
  expect_s3_class(ll_by_x, "data.frame")
  expect_true("log_lik" %in% names(ll_by_x))
})

test_that(
  "evaluate_model - output check (method log_lik_permutation_constant)",
  {
    setup <- make_evaluate_model_test_data()
    my_model <- setup$model
    my_data <- setup$data

    const_val <- evaluate_model(
      model = my_model,
      x = my_data$cues,
      response_category = my_data$response,
      method = "log_lik_permutation_constant"
    )
    expect_true(is.atomic(const_val))
    expect_true(is.finite(const_val))
    expect_true(const_val >= 0)

    const_by_x <- evaluate_model(
      model = my_model,
      x = my_data$cues,
      response_category = my_data$response,
      method = "log_lik_permutation_constant",
      return_by_x = TRUE
    )
    expect_s3_class(const_by_x, "data.frame")
    expect_true("log_lik_permutation_constant" %in% names(const_by_x))
  }
)

test_that(
  "evaluate_model - output check (method likelihood-up-to-constant alias)",
  {
    setup <- make_evaluate_model_test_data()
    my_model <- setup$model
    my_data <- setup$data

    expect_warning(
      likelihood_constant_proportional <- evaluate_model(
        model = my_model,
        x = my_data$cues,
        response_category = my_data$response,
        decision_rule = "proportional",
        method = "likelihood-up-to-constant"
      ),
      class = "lifecycle_warning_deprecated"
    )
    expect_true(is.atomic(likelihood_constant_proportional))
    expect_true(is.finite(likelihood_constant_proportional))

    ll_by_x <- suppressWarnings(
      evaluate_model(
        model = my_model,
        x = my_data$cues,
        response_category = my_data$response,
        method = "likelihood-up-to-constant",
        return_by_x = TRUE
      )
    )
    expect_s3_class(ll_by_x, "data.frame")
    expect_true("log_likelihood" %in% names(ll_by_x))
  }
)

test_that("evaluate_model - output check (method default)", {
  setup <- make_evaluate_model_test_data()
  my_model <- setup$model
  my_data <- setup$data

  default_proportional <- evaluate_model(
    model = my_model,
    x = my_data$cues,
    response_category = my_data$response,
    decision_rule = "proportional"
  )
  default_criterion <- evaluate_model(
    model = my_model,
    x = my_data$cues,
    response_category = my_data$response,
    decision_rule = "criterion"
  )

  expect_true(is.atomic(default_proportional))
  expect_true(is.atomic(default_criterion))

  ll_explicit <- evaluate_model(
    model = my_model,
    x = my_data$cues,
    response_category = my_data$response,
    decision_rule = "proportional",
    method = "log_lik"
  )
  expect_equal(default_proportional, ll_explicit)
})

test_that("evaluate_model on MVBU_Stanfit works", {
  skip_if_not_installed("rstan")
  model_path <- testthat::test_path(
    "models",
    "minimal-nix_ideal_adaptor-NIX.rds"
  )
  fit <- example_ideal_adaptor_stanfit(
    model_family = "NIX",
    n_cues = 1L,
    file = model_path,
    file_refit = "on_change"
  )

  fit_no_test <- fit
  fit_no_test@data <- fit@data[fit@data$Phase != "test", , drop = FALSE]
  expect_error(evaluate_model(fit_no_test), "No test data found")

  eval_default <- evaluate_model(fit)
  expect_true(is.atomic(eval_default))
  expect_true(is.finite(eval_default))
  expect_true(eval_default < 0)

  x_test <- matrix(c(0.2, 0.8), ncol = 1)
  resp_test <- c("/d/", "/t/")
  eval_res <- evaluate_model(
    fit,
    x = x_test,
    response_category = resp_test
  )
  expect_true(is.atomic(eval_res))
  expect_true(is.finite(eval_res))
  expect_true(eval_res < 0)

  eval_acc <- evaluate_model(
    fit,
    x = x_test,
    response_category = resp_test,
    method = "accuracy"
  )
  expect_true(is.atomic(eval_acc))
  expect_true(eval_acc >= 0 && eval_acc <= 1)
})
