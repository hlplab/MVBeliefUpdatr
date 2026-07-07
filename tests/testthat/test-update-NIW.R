context("update NIW models")

test_that("incremental updating runs for all methods (1 cue)", {
  model <- example_NIW_ideal_adaptor(example = 1, verbose = FALSE)  # 1 cue: VOT
  cats  <- get_category_labels_from_model(model)

  exposure <- tibble::tibble(
    category = rep(cats, each = 2),
    VOT      = c(10, 20, 40, 50)
  )

  for (m in c("label-certain", "nolabel-criterion", "nolabel-sampling",
              "nolabel-proportional", "nolabel-uniform")) {
    expect_no_error(
      update_NIW_ideal_adaptor_incrementally(
        model, exposure = exposure,
        exposure.category = "category", exposure.cues = "VOT",
        method = m, verbose = FALSE)
    )
  }

  # a bad method should still error clearly
  expect_error(
    update_NIW_ideal_adaptor_incrementally(
      model, exposure = exposure,
      exposure.category = "category", exposure.cues = "VOT",
      method = "not-a-method"),
    regexp = "not an acceptable updating method"
  )
})
