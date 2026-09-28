utils::data("pb52", package = "MVBeliefUpdatr")
pb52_clean <- pb52[stats::complete.cases(pb52), ]
ex_template <- new_exemplar_category_representation_template_from_data(
  pb52_clean,
  category = "vowel",
  cues = c("f1", "f2")
)

test_that("deprecated Exemplar make wrappers delegate to S7 constructors", {
  expect_warning(
    x <- make_exemplars_from_data(
      pb52_clean,
      category = "vowel",
      cues = c("f1", "f2")
    ),
    "deprecated"
  )
  expect_true(S7::S7_inherits(x, MVBU_CategoryRepresentationTemplate))
  expect_warning(
    x <- make_exemplar_model_from_data(
      pb52_clean,
      category = "vowel",
      cues = c("f1", "f2")
    ),
    "deprecated"
  )
  expect_true(S7::S7_inherits(x, Exemplar_Model))
  expect_error(
    make_exemplars_from_data(
      pb52_clean,
      group = "vowel",
      category = "vowel",
      cues = c("f1", "f2")
    )
  )
  expect_warning(
    make_exemplars_from_data(
      pb52_clean,
      category = "vowel",
      cues = c("f1", "f2"),
      sim_function = function(x, y) 1
    ),
    "no longer supported"
  )
})

test_that("deprecated Exemplar lift wrapper delegates to S7 constructor", {
  expect_warning(
    model <- lift_exemplars_to_exemplar_model(ex_template),
    "deprecated"
  )
  expect_true(S7::S7_inherits(model, Exemplar_Model))
  expect_error(
    lift_exemplars_to_exemplar_model(ex_template, group = "category")
  )
})

