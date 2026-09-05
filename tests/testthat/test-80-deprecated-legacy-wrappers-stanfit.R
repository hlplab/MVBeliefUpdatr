
test_that("ideal_adaptor_staninput creates an S7 object with list slots", {
  obj <- IdealAdaptorStaninput(values = list(
    transformed = list(M = 2, D = 1),
    untransformed = list(raw = TRUE)
  ))

  expect_true(S7::S7_inherits(obj, IdealAdaptorStaninput))
  expect_equal(obj@values$transformed$M, 2)
  expect_equal(obj@values$untransformed$raw, TRUE)
})

test_that("transform_information accepts list and function properties", {
  identity_fn <- function(x) x
  obj <- MVBU_TransformInformation(
    transform.parameters = list(scale = 1),
    transform.function = identity_fn,
    untransform.function = identity_fn
  )

  expect_true(S7::S7_inherits(obj, MVBU_TransformInformation))
  expect_equal(obj@transform.parameters$scale, 1)
  expect_equal(obj@transform.function(3), 3)
})

test_that("ideal_adaptor_stanfit stores staninput and transform information", {
  staninput_obj <- IdealAdaptorStaninput(values = list(transformed = list(M = 2), untransformed = list()))
  transform_obj <- MVBU_TransformInformation(
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
  expect_true(S7::S7_inherits(obj@transform_information, MVBU_TransformInformation))
})

test_that("S7 accessors expose staninput and transform data", {
  staninput_obj <- IdealAdaptorStaninput(values = list(transformed = list(M = 2), untransformed = list(raw = TRUE)))
  obj <- IdealAdaptorStanfit(staninput = staninput_obj)

  expect_equal(get_staninput(obj)@values$transformed$M, 2)
  expect_equal(get_staninput(obj)@values$untransformed$raw, TRUE)
  expect_equal(get_transform_function(obj)(4), 4)
  expect_equal(get_untransform_function(obj)(4), 4)
})

test_that("MVBU_Stanfit creates a placeholder staninput before fitting", {
  obj <- MVBU_Stanfit(data = data.frame(x = 1:3))

  expect_true(S7::S7_inherits(obj, MVBU_Stanfit))
  expect_true(S7::S7_inherits(obj@staninput, MVBU_Staninput))
  expect_true(is.list(obj@staninput@values))
})

test_that("family-specific fit constructors use the appropriate S7 subclass", {
  staninput_obj <- NIX_IdealAdaptorStaninput(values = list(
    transformed = list(M = 2),
    untransformed = list(raw = TRUE)
  ))

  fit_obj <- ideal_adaptor_stanfit(staninput = staninput_obj)

  expect_true(S7::S7_inherits(fit_obj, NIX_IdealAdaptorStanfit))
  expect_true(S7::S7_inherits(fit_obj, IdealAdaptorStanfit))
})
