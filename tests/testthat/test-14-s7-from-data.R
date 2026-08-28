idahoans <- make_vowel_test_data()

representation_types <- c("UVG", "NIX", "MUVG", "MNIX", "MVG", "NIW", "EXEMPLAR")
template_types <- representation_types
model_types <- representation_types

representation_cues <- function(type) if (type %in% c("UVG", "NIX")) "F1" else c("F1", "F2")

representation_args <- function(type) {
  args <- list(data = idahoans, category = "vowel", cues = representation_cues(type))
  if (type %in% c("NIX", "MNIX", "NIW")) args[c("kappa", "nu")] <- list(kappa = 10, nu = 30)
  args
}

test_that("category representation from-data constructors cover all families", {
  for (type in representation_types) {
    args <- representation_args(type)
    args$data <- args$data[args$data$vowel == "AA", ]
    args$category <- "vowel"
    result <- do.call(new_category_representation_from_data, c(args, type = type))
    expect_true(S7::S7_inherits(result, MVBU_CategoryRepresentation))
  }
  expect_error(new_category_representation_from_data(idahoans, type = "unknown", cues = "F1"))
})

test_that("category representation from-data constructors require one category", {
  expect_true(S7::S7_inherits(
    new_mvg_category_representation_from_data(idahoans[idahoans$vowel == "AA", ], category = "vowel", cues = "F1"),
    MVG_CategoryRepresentation))
  expect_error(new_mvg_category_representation_from_data(idahoans, category = "vowel", cues = "F1"))
})

test_that("category representation template from-data constructors cover all families", {
  for (type in template_types) {
    result <- do.call(
      new_category_representation_template_from_data,
      c(list(data = idahoans, category = "vowel", cues = representation_cues(type)), type = type,
        if (type %in% c("NIX", "MNIX", "NIW")) list(kappa = 10, nu = 30) else list())
    )
    expect_true(S7::S7_inherits(result, MVBU_CategoryRepresentationTemplate))
  }
  expect_error(new_category_representation_template_from_data(idahoans, type = "unknown", cues = "F1"))
})

test_that("model from-data constructors cover all families", {
  for (type in model_types) {
    result <- do.call(
      new_model_from_data,
      c(list(data = idahoans, category = "vowel", cues = representation_cues(type)), type = type,
        if (type %in% c("NIX", "MNIX", "NIW")) list(kappa = 10, nu = 30) else list())
    )
    expect_true(S7::S7_inherits(result, MVBU_CognitiveModel))
  }
  expect_error(new_model_from_data(idahoans, type = "unknown", cues = "F1"))
})

test_that("one-cue family constraints are enforced", {
  expect_error(new_category_representation_from_data(idahoans, type = "UVG", cues = c("F1", "F2")))
  expect_error(new_category_representation_from_data(idahoans, type = "NIX", cues = c("F1", "F2")))
  expect_true(S7::S7_inherits(
    new_category_representation_from_data(idahoans[idahoans$vowel == "AA", ], type = "MUVG", category = "vowel", cues = c("F1", "F2")),
    MUVG_CategoryRepresentation))
  expect_true(S7::S7_inherits(
    new_category_representation_from_data(idahoans[idahoans$vowel == "AA", ], type = "MNIX", category = "vowel", cues = c("F1", "F2"), kappa = 10, nu = 30),
    MNIX_CategoryRepresentation))
})
