idahoans <- make_vowel_test_data()

my_model <- suppressMessages(suppressWarnings(
  make_MVG_ideal_observer_from_data(idahoans, category = "vowel", cues = c("F1", "F2"), verbose = FALSE)))
x.1 <- idahoans %>% mutate(x = map(F1, ~ c(...))) %>% pull(x)
x.2 <- idahoans %>% mutate(x = map2(F1, F2, ~ c(...))) %>% pull(x)
x.3 <- idahoans %>% mutate(x = pmap(.l = list(F1, F2, F3), ~ c(...))) %>% pull(x)

test_that("Test is.MVG_ideal_observer", {
  expect_false(suppressWarnings(is.MVG_ideal_observer(NULL)))
  expect_false(suppressWarnings(is.MVG_ideal_observer(NA)))
  expect_false(suppressWarnings(is.MVG_ideal_observer(1)))
  expect_false(suppressWarnings(is.MVG_ideal_observer("1")))
  expect_false(suppressWarnings(is.MVG_ideal_observer(TRUE)))
  expect_false(suppressWarnings(is.MVG_ideal_observer(list(1))))
  expect_false(suppressWarnings(is.MVG_ideal_observer(suppressMessages(example_exemplar_model(n_cues = 1)))))
  expect_true(suppressWarnings(is.MVG_ideal_observer(suppressMessages(example_mvg_ideal_observer(n_cues = 1)))))
  expect_false(suppressWarnings(is.MVG_ideal_observer(suppressMessages(example_niw_ideal_adaptor(n_cues = 1)))))
})

test_that("Get MVG likelihood - input check x (single element, non-list)", {
  expect_error(get_MVG_likelihood(
    x = 1,
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
  expect_error(get_MVG_likelihood(
    x = matrix(c(1,0,1,2,1,2), nrow = 6),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
  expect_error(get_MVG_likelihood(
    x = matrix(c(1,0,1,2,1,2), nrow = 2),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
  expect_no_error(get_MVG_likelihood(
    x = matrix(c(1,0,1,2,1,2), nrow = 3),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
})

test_that("Get MVG likelihood - input check x (single-element list)", {
  expect_error(get_MVG_likelihood(
    x = list(1),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
  expect_error(get_MVG_likelihood(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 6)),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
  expect_error(get_MVG_likelihood(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 2)),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
  expect_no_error(get_MVG_likelihood(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 3)),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
})

test_that("Get MVG likelihood - input check x (multi-element list)", {
  expect_error(get_MVG_likelihood(
    x = map(rep(1, 3), ~ .x),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
  expect_error(get_MVG_likelihood(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 6), matrix(c(1,0,1,2,1,2), nrow = 6), matrix(c(1,0,1,2,1,2), nrow = 6)),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
  expect_error(get_MVG_likelihood(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 2), matrix(c(1,0,1,2,1,2), nrow = 2), matrix(c(1,0,1,2,1,2), nrow = 2)),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
  expect_no_error(get_MVG_likelihood(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 3), matrix(c(1,0,1,2,1,2), nrow = 3), matrix(c(1,0,1,2,1,2), nrow = 3)),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
})

test_that("Get likelihood from MVG - input check x (single element, non-list)", {
  expect_error(suppressWarnings(get_likelihood_from_MVG(
    x = 1,
    model = my_model)))
  expect_error(suppressWarnings(get_likelihood_from_MVG(
    x = matrix(c(1,0,1,2,1,2), nrow = 6),
    model = my_model)))
  expect_error(suppressWarnings(get_likelihood_from_MVG(
    x = matrix(c(1,0,1,2,1,2), nrow = 2),
    model = my_model)))
  expect_no_error(suppressWarnings(get_likelihood_from_MVG(
    x = matrix(c(1,0,1,2,1,2), nrow = 3),
    model = my_model)))
})

test_that("Get likelihood from MVG - input check x (single-element list)", {
  expect_error(suppressWarnings(get_likelihood_from_MVG(
    x = list(1),
    model = my_model)))
  expect_error(suppressWarnings(get_likelihood_from_MVG(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 6)),
    model = my_model)))
  expect_error(suppressWarnings(get_likelihood_from_MVG(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 2)),
    model = my_model)))
  expect_no_error(suppressWarnings(get_likelihood_from_MVG(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 3)),
    model = my_model)))
})

test_that("Get likelihood from MVG - input check x (multi-element list)", {
  expect_error(suppressWarnings(get_likelihood_from_MVG(
    x = map(rep(1, 3), ~ .x),
    model = my_model)))
  expect_error(suppressWarnings(get_likelihood_from_MVG(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 6), matrix(c(1,0,1,2,1,2), nrow = 6), matrix(c(1,0,1,2,1,2), nrow = 6)),
    model = my_model)))
  expect_error(suppressWarnings(get_likelihood_from_MVG(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 2), matrix(c(1,0,1,2,1,2), nrow = 2), matrix(c(1,0,1,2,1,2), nrow = 2)),
    model = my_model)))
  expect_no_error(suppressWarnings(get_likelihood_from_MVG(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 3), matrix(c(1,0,1,2,1,2), nrow = 3), matrix(c(1,0,1,2,1,2), nrow = 3)),
    model = my_model)))
})

test_that("Get likelihood from MVG - input check x", {
  expect_error(suppressWarnings(get_likelihood_from_MVG(x = x.1, model = my_model)))
  expect_error(suppressWarnings(get_likelihood_from_MVG(x = x.3, model = my_model)))
  expect_no_error(suppressWarnings(get_likelihood_from_MVG(x = x.2, model = my_model)))
})

test_that("Get categorization from MVG ideal observer - input check x", {
  result <- expect_warning(
    get_categorization_from_MVG_ideal_observer(
      x = x.2,
      model = my_model,
      noise_treatment = "no_noise",
      lapse_treatment = "no_lapses",
      decision_rule = "sampling"
    ),
    "deprecated"
  )

  expect_true(is.list(result))
})

test_that("MVG categorization aligns categories with observations", {
  x <- x.2[1:3]
  category_labels <- as.character(get_category_labels(my_model))

  result <- suppressWarnings(get_categorization_from_MVG_ideal_observer(
    x = x,
    model = my_model,
    decision_rule = "proportional"
  ))

  expect_equal(result$observationID, rep(seq_along(x), each = length(category_labels)))
  expect_equal(as.character(result$category), rep(category_labels, times = length(x)))
  expect_equal(result$x, rep(x, each = length(category_labels)))
})
