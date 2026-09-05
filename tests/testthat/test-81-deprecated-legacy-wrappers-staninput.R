
test_that("deprecated wrappers forward to the S7 fit-input constructor", {
  exposure <- data.frame(
    category = factor(c("A", "B")),
    group = factor(c("g1", "g1")),
    cue1 = c(0.1, 0.4),
    cue2 = c(0.2, 0.3)
  )
  test <- data.frame(
    response = factor(c("A", "B")),
    group = factor(c("g1", "g1")),
    cue1 = c(0.15, 0.45),
    cue2 = c(0.25, 0.35)
  )

  expect_warning(
    obj <- make_staninput(
      exposure = exposure,
      test = test,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "response",
      group = "group",
      stanmodel = "NIW_ideal_adaptor"
    ),
    "new_ideal_adaptor_stanfit_input"
  )

  expect_true(S7::S7_inherits(obj, IdealAdaptorStanfitInput))
  expect_true(S7::S7_inherits(obj@staninput, IdealAdaptorStaninput))
})

test_that("legacy wrappers reject invalid fixed parameters", {
  exposure <- data.frame(
    category = factor(c("A", "B")),
    group = factor(c("g1", "g1")),
    cue1 = c(0, 1),
    cue2 = c(0, 1)
  )
  test <- data.frame(
    response = factor(c("A", "B")),
    group = factor(c("g1", "g1")),
    cue1 = c(0.1, 0.9),
    cue2 = c(0.1, 0.9)
  )

  expect_warning(
    expect_error(
      make_ideal_adaptor_stanfit_input(
        exposure = exposure,
        test = test,
        cues = c("cue1", "cue2"),
        category = "category",
        response = "response",
        group = "group",
        lapse_rate = 1.5,
        control = control_staninput(transform_type = "identity"),
        stanmodel = "NIW_ideal_adaptor"
      ),
      "between 0 and 1"
    ),
    "new_ideal_adaptor_stanfit_input"
  )
})
