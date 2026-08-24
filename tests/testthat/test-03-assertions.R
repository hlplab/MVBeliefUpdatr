test_that(".assert_that handles simple and multi-condition checks", {
  expect_true(MVBeliefUpdatr:::.assert_that(TRUE))
  expect_true(MVBeliefUpdatr:::.assert_that(1 > 0, 2 > 1))
  expect_error(
    MVBeliefUpdatr:::.assert_that(FALSE, msg = "custom message"),
    "custom message"
  )
})

test_that(".assert_true and .assert_false provide clear wrappers", {
  expect_true(MVBeliefUpdatr:::.assert_true(1 == 1, "must hold"))
  expect_true(MVBeliefUpdatr:::.assert_false(FALSE, "must be false"))
  expect_error(
    MVBeliefUpdatr:::.assert_true(FALSE, "must fail"),
    "must fail"
  )
  expect_error(
    MVBeliefUpdatr:::.assert_false(TRUE, "must fail"),
    "must fail"
  )
})

test_that("internal predicates distinguish basic scalar types", {
  expect_true(MVBeliefUpdatr:::.is_scalar(1))
  expect_true(MVBeliefUpdatr:::.is_scalar(NA))
  expect_true(MVBeliefUpdatr:::.is_scalar(list(1)))
  expect_true(MVBeliefUpdatr:::.is_scalar(matrix(1, 1, 1)))
  expect_false(MVBeliefUpdatr:::.is_scalar(NULL))
  expect_false(MVBeliefUpdatr:::.is_scalar(c(1, 2)))

  expect_true(MVBeliefUpdatr:::.is_scalar_numeric(1))
  expect_true(MVBeliefUpdatr:::.is_scalar_numeric(1L))
  expect_true(MVBeliefUpdatr:::.is_scalar_integer(1L))
  expect_false(MVBeliefUpdatr:::.is_scalar_integer(1))
  expect_true(MVBeliefUpdatr:::.is_scalar_double(1))
  expect_false(MVBeliefUpdatr:::.is_scalar_double(1L))
  expect_true(MVBeliefUpdatr:::.is_scalar_character("x"))
  expect_true(MVBeliefUpdatr:::.is_scalar_factor(factor("x")))
  expect_true(MVBeliefUpdatr:::.is_scalar_logical(TRUE))
  expect_true(MVBeliefUpdatr:::.is_scalar_count(0))
  expect_true(MVBeliefUpdatr:::.is_scalar_count(2L))
  expect_true(MVBeliefUpdatr:::.is_scalar_numeric(NA_real_))
  expect_true(MVBeliefUpdatr:::.is_scalar_integer(NA_integer_))
  expect_true(MVBeliefUpdatr:::.is_scalar_double(NA_real_))
  expect_true(MVBeliefUpdatr:::.is_scalar_character(NA_character_))
  expect_true(MVBeliefUpdatr:::.is_scalar_factor(factor(NA)))
  expect_true(MVBeliefUpdatr:::.is_scalar_logical(NA))
  expect_true(MVBeliefUpdatr:::.is_scalar_count(NA_real_))

  expect_false(MVBeliefUpdatr:::.is_scalar_character(c("x", "y")))
  expect_false(MVBeliefUpdatr:::.is_scalar_count(-1))
  expect_false(MVBeliefUpdatr:::.is_scalar_count(1.5))

  expect_false(MVBeliefUpdatr:::.is_non_NA_scalar_numeric(NA_real_))
  expect_false(MVBeliefUpdatr:::.is_non_NA_scalar_integer(NA_integer_))
  expect_false(MVBeliefUpdatr:::.is_non_NA_scalar_double(NA_real_))
  expect_false(MVBeliefUpdatr:::.is_non_NA_scalar_character(NA_character_))
  expect_false(MVBeliefUpdatr:::.is_non_NA_scalar_factor(factor(NA)))
  expect_false(MVBeliefUpdatr:::.is_non_NA_scalar_logical(NA))
  expect_false(MVBeliefUpdatr:::.is_non_NA_scalar_count(NA_real_))

  expect_true(MVBeliefUpdatr:::.is_equal(1, 1))
  expect_true(MVBeliefUpdatr:::.is_like_factor("x"))
  expect_true(MVBeliefUpdatr:::.is_try_error(try(stop("error"), silent = TRUE)))
  expect_true(MVBeliefUpdatr:::.is_sigma(diag(2)))
  expect_false(MVBeliefUpdatr:::.is_sigma(matrix(c(1, 2, 2, 1), 2)))
})

test_that(".assert_* helpers validate basic variable types", {
  expect_true(MVBeliefUpdatr:::.assert_character("x"))
  expect_true(MVBeliefUpdatr:::.assert_numeric(1:3))
  expect_true(MVBeliefUpdatr:::.assert_logical(c(TRUE, FALSE)))
  expect_true(MVBeliefUpdatr:::.assert_factor(factor("x")))
  expect_true(MVBeliefUpdatr:::.assert_list(list(a = 1)))

  expect_true(MVBeliefUpdatr:::.assert_character_scalar("x"))
  expect_true(MVBeliefUpdatr:::.assert_numeric_scalar(1))
  expect_true(MVBeliefUpdatr:::.assert_logical_scalar(NA))
  expect_true(MVBeliefUpdatr:::.assert_factor_scalar(factor("x")))

  expect_true(MVBeliefUpdatr:::.assert_optional_list(NULL))

  expect_error(MVBeliefUpdatr:::.assert_character(1), "Expected x to be a character vector")
  expect_error(MVBeliefUpdatr:::.assert_numeric("x"), "Expected x to be a numeric vector")
  expect_error(MVBeliefUpdatr:::.assert_logical(1), "Expected x to be a logical vector")
  expect_error(MVBeliefUpdatr:::.assert_factor(1), "Expected x to be a factor")
  expect_error(MVBeliefUpdatr:::.assert_list("x"), "Expected x to be a list")

  expect_error(MVBeliefUpdatr:::.assert_character_scalar(c("x", "y")), "Expected x to be a scalar character")
  expect_error(MVBeliefUpdatr:::.assert_numeric_scalar(c(1, 2)), "Expected x to be a scalar numeric")
  expect_error(MVBeliefUpdatr:::.assert_logical_scalar(c(TRUE, FALSE)), "Expected x to be a scalar logical")
  expect_error(MVBeliefUpdatr:::.assert_factor_scalar(factor(c("x", "y"))), "Expected x to be a scalar factor")

  expect_error(MVBeliefUpdatr:::.assert_optional_list("x"), "Expected x to be missing, NULL, or a list")
})

test_that("public S7 class assertions validate the expected S7 classes", {
  representation <- MVBU_CategoryRepresentation(
    category_likelihood_function = function(x) x,
    metadata = list(label_information = list(category = c("A"), cue = c("x")))
  )

  template <- MVBU_CategoryRepresentationTemplate(
    representations = list(representation),
    metadata = list()
  )

  model <- MVBU_CognitiveModel(
    category_template = template,
    category_posterior_functions = list(),
    decision_rule = "max",
    category_prior = c(A = 1),
    lapse_behavior = list(lapse_rate = 0, lapse_bias = c(A = 1), lapse_treatment = "no_lapses"),
    noise_behavior = list(Sigma_noise = NULL, noise_treatment = "no_noise"),
    metadata = list()
  )

  distribution <- MVBU_ModelDistribution(
    model_family = "MVG",
    cache = list(),
    metadata = list(),
    group_label = "group"
  )

  stanfit <- IdealAdaptorStanfit()

  expect_true(MVBeliefUpdatr::assert_MVBU_CategoryRepresentation(representation))
  expect_true(MVBeliefUpdatr::assert_MVBU_CategoryRepresentationTemplate(template))
  expect_true(MVBeliefUpdatr::assert_MVBU_CognitiveModel(model))
  expect_true(MVBeliefUpdatr::assert_MVBU_ModelDistribution(distribution))
  expect_true(MVBeliefUpdatr::assert_IdealAdaptorStanfit(stanfit))

  expect_error(MVBeliefUpdatr::assert_MVBU_CategoryRepresentation(list()), "must inherit from MVBU_CategoryRepresentation")
  expect_error(MVBeliefUpdatr::assert_IdealAdaptorStanfit(list()), "must inherit from IdealAdaptorStanfit")
})

test_that("required non-NA helpers reject missing and NA values", {
  expect_true(MVBeliefUpdatr:::.assert_non_NA_character("x"))
  expect_true(MVBeliefUpdatr:::.assert_non_NA_numeric(1))
  expect_true(MVBeliefUpdatr:::.assert_non_NA_logical(TRUE))
  expect_true(MVBeliefUpdatr:::.assert_non_NA_factor(factor("x")))
  expect_true(MVBeliefUpdatr:::.assert_non_NA_list(list(a = 1)))
  expect_true(MVBeliefUpdatr:::.assert_non_NA_scalar_character("x"))
  expect_true(MVBeliefUpdatr:::.assert_non_NA_scalar_numeric(1))
  expect_true(MVBeliefUpdatr:::.assert_non_NA_scalar_logical(TRUE))
  expect_true(MVBeliefUpdatr:::.assert_non_NA_scalar_factor(factor("x")))

  expect_true(MVBeliefUpdatr:::.assert_character("x"))
  expect_true(MVBeliefUpdatr:::.assert_character(NA_character_))
  expect_true(MVBeliefUpdatr:::.assert_optional_character(NULL))
  expect_true(MVBeliefUpdatr:::.assert_matrix(matrix(1, 1, 1)))
  expect_true(MVBeliefUpdatr:::.assert_non_negative(0))
  expect_true(MVBeliefUpdatr:::.assert_between(0.5, 0, 1))
  expect_true(MVBeliefUpdatr:::.assert_one_of("sample", c("sample", "marginalize")))

  expect_error(MVBeliefUpdatr:::.assert_non_NA_character(), "Expected x to be a non-NA character")
  expect_error(MVBeliefUpdatr:::.assert_non_NA_character(NA_character_), "Expected x to be a non-NA character")
  expect_error(MVBeliefUpdatr:::.assert_non_NA_numeric(NA_real_), "Expected x to be a non-NA numeric")
  expect_error(MVBeliefUpdatr:::.assert_non_NA_logical(NA), "Expected x to be a non-NA logical")
  expect_error(MVBeliefUpdatr:::.assert_non_NA_factor(NA), "Expected x to be a non-NA factor")
  expect_error(MVBeliefUpdatr:::.assert_non_NA_scalar_character(NA_character_), "Expected x to be a non-NA, non-empty character scalar")
  expect_error(MVBeliefUpdatr:::.assert_non_NA_scalar_numeric(NA_real_), "Expected x to be a non-NA numeric scalar")
  expect_error(MVBeliefUpdatr:::.assert_non_NA_scalar_factor(factor(NA)), "Expected x to be a non-NA scalar factor")
  expect_error(MVBeliefUpdatr:::.assert_non_NA_scalar_logical(NA), "Expected x to be a non-NA logical scalar")
  expect_error(MVBeliefUpdatr:::.assert_matrix(1), "Expected x to be a matrix")
  expect_error(MVBeliefUpdatr:::.assert_non_negative(-1), "Expected x to be non-negative")
  expect_error(MVBeliefUpdatr:::.assert_between(2, 0, 1), "Expected x to be between 0 and 1")
  expect_error(MVBeliefUpdatr:::.assert_one_of("other", c("sample", "marginalize")), "Expected x to be one of")
})

test_that(".assert_data_contains_cols handles column presence checks", {
  data <- data.frame(a = 1, b = 2)

  expect_true(MVBeliefUpdatr:::.assert_data_contains_cols(data, c("a", "b")))
  expect_error(
    MVBeliefUpdatr:::.assert_data_contains_cols(data, c("a", "missing")),
    "Expected data data to contain column\\(s\\): missing"
  )
})
