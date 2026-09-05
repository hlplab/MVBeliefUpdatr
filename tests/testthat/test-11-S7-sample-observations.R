test_that("sample_observations works on MVBU_CategoryRepresentation", {
  rep <- new_uvg_category_representation("A", "F1", mu = 100, sigma2 = 25)
  draws <- sample_observations(rep, n = 20L)
  expect_s3_class(draws, "data.frame")
  expect_equal(nrow(draws), 20L)
  expect_true("F1" %in% names(draws))
  expect_true("category" %in% names(draws))
  expect_true(all(draws$category == "A"))
})

test_that("sample_observations works uniformly on MVBU_CategoryRepresentationTemplate", {
  repA <- new_uvg_category_representation("A", "F1", mu = 100, sigma2 = 25)
  repB <- new_uvg_category_representation("B", "F1", mu = 200, sigma2 = 25)
  tpl <- new_category_representation_template(list(A = repA, B = repB))

  set.seed(123)
  draws <- sample_observations(tpl, n = 100L)
  expect_s3_class(draws, "data.frame")
  expect_equal(nrow(draws), 100L)
  expect_true(all(c("A", "B") %in% unique(draws$category)))

  # Vector Ns
  draws_vec <- sample_observations(tpl, n = c(A = 30L, B = 70L), randomize_order = FALSE)
  expect_equal(sum(draws_vec$category == "A"), 30L)
  expect_equal(sum(draws_vec$category == "B"), 70L)
})

test_that("sample_observations samples proportionally to category prior on MVBU_CognitiveModel", {
  repA <- new_uvg_category_representation("A", "F1", mu = 100, sigma2 = 25)
  repB <- new_uvg_category_representation("B", "F1", mu = 200, sigma2 = 25)
  tpl <- new_category_representation_template(list(A = repA, B = repB))
  mod <- new_uvg_ideal_observer(tpl, category_prior = c(A = 0.9, B = 0.1))

  set.seed(42)
  draws <- sample_observations(mod, n = 500L)
  expect_s3_class(draws, "data.frame")
  expect_equal(nrow(draws), 500L)
  # High probability category A should dominate
  prop_A <- mean(draws$category == "A")
  expect_gt(prop_A, 0.80)
})

test_that("sample_observations handles with_replacement FALSE and error guards on exemplar models", {
  # Create exemplar representation with small exemplar set
  exA <- data.frame(F1 = c(10, 11, 12, 13, 14), category = "A")
  exB <- data.frame(F1 = c(20, 21, 22), category = "B")
  repA <- new_exemplar_category_representation_from_data(exA, cues = "F1", category = "category")
  repB <- new_exemplar_category_representation_from_data(exB, cues = "F1", category = "category")
  tpl <- new_category_representation_template(list(A = repA, B = repB))

  # Sampling within limits without replacement succeeds
  draws_no_rep <- sample_observations(tpl, n = c(A = 3L, B = 2L), with_replacement = FALSE)
  expect_equal(nrow(draws_no_rep), 5L)

  # Requesting more than total available fails with informative error
  expect_error(
    sample_observations(tpl, n = 10L, with_replacement = FALSE),
    "Requested 10 samples without replacement, but template only contains 8 total exemplars"
  )

  # Requesting single representation beyond size without replacement fails
  expect_error(
    sample_observations(repB, n = 5L, with_replacement = FALSE),
    "Cannot sample 5 observations without replacement"
  )
})
