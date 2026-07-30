context("S7 Stan wrapper classes")

test_that("ideal_adaptor_staninput creates an S7 object with list slots", {
  obj <- IdealAdaptorStaninput(
    transformed = list(M = 2, D = 1),
    untransformed = list(raw = TRUE)
  )

  expect_true(S7::S7_inherits(obj, IdealAdaptorStaninput))
  expect_equal(obj@transformed$M, 2)
  expect_equal(obj@untransformed$raw, TRUE)
})

test_that("transform_information accepts list and function properties", {
  identity_fn <- function(x) x
  obj <- transform_information(
    transform.parameters = list(scale = 1),
    transform.function = identity_fn,
    untransform.function = identity_fn
  )

  expect_true(S7::S7_inherits(obj, transform_information_class))
  expect_equal(obj@transform.parameters$scale, 1)
  expect_equal(obj@transform.function(3), 3)
})

test_that("ideal_adaptor_stanfit stores staninput and transform information", {
  staninput_obj <- IdealAdaptorStaninput(transformed = list(M = 2), untransformed = list())
  transform_obj <- transform_information(
    transform.parameters = list(),
    transform.function = function(x) x,
    untransform.function = function(x) x
  )

  obj <- IdealAdaptorStanfit(
    data = data.frame(x = 1:3),
    staninput = staninput_obj,
    transform_information = transform_obj,
    backend = "rstan"
  )

  expect_true(S7::S7_inherits(obj, IdealAdaptorStanfit))
  expect_equal(obj@backend, "rstan")
  expect_true(S7::S7_inherits(obj@staninput, IdealAdaptorStaninput))
  expect_true(S7::S7_inherits(obj@transform_information, transform_information_class))
})

test_that("S7 accessors expose staninput and transform data from legacy fit objects", {
  staninput_obj <- IdealAdaptorStaninput(transformed = list(M = 2), untransformed = list(raw = TRUE))
  obj <- IdealAdaptorStanfit(staninput = staninput_obj)

  expect_equal(get_staninput(obj, which = "transformed")$M, 2)
  expect_equal(get_staninput(obj, which = "untransformed")$raw, TRUE)
  expect_equal(get_transform_function(obj)(4), 4)
  expect_equal(get_untransform_function(obj)(4), 4)
})

test_that("family-style Stan wrapper classes inherit from the shared S7 base classes", {
  staninput_obj <- IdealAdaptorStaninput(
    transformed = list(M = 2),
    untransformed = list(raw = TRUE)
  )
  transform_obj <- MVBU_TransformInformation(
    transform.parameters = list(scale = 1),
    transform.function = function(x) x,
    untransform.function = function(x) x
  )
  fit_obj <- IdealAdaptorStanfit(
    data = data.frame(x = 1:3),
    staninput = staninput_obj,
    transform_information = transform_obj,
    backend = "rstan"
  )

  expect_true(S7::S7_inherits(staninput_obj, IdealAdaptorStaninput))
  expect_true(S7::S7_inherits(staninput_obj, MVBU_Staninput))
  expect_true(S7::S7_inherits(transform_obj, MVBU_TransformInformation))
  expect_true(S7::S7_inherits(fit_obj, IdealAdaptorStanfit))
  expect_true(S7::S7_inherits(fit_obj, MVBU_Stanfit))
})
