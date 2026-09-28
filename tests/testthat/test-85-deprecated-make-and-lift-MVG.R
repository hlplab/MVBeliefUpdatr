utils::data("pb52", package = "MVBeliefUpdatr")
pb52_clean <- pb52[stats::complete.cases(pb52), ]
mvg_template <- new_mvg_category_representation_template_from_data(
  pb52_clean,
  category = "vowel",
  cues = c("f1", "f2")
)

test_that("deprecated MVG make wrappers delegate to S7 constructors", {
  expect_warning(
    x <- make_MVG_from_data(
      pb52_clean,
      category = "vowel",
      cues = c("f1", "f2")
    ),
    "deprecated"
  )
  expect_true(S7::S7_inherits(x, MVBU_CategoryRepresentationTemplate))
  expect_warning(
    x <- make_MVG_ideal_observer_from_data(
      pb52_clean,
      category = "vowel",
      cues = c("f1", "f2")
    ),
    "deprecated"
  )
  expect_true(S7::S7_inherits(x, MVG_IdealObserver))
  expect_error(
    make_MVG_from_data(
      pb52_clean,
      group = "vowel",
      category = "vowel",
      cues = c("f1", "f2")
    )
  )
})

test_that("deprecated MVG lift wrapper delegates to S7 constructor", {
  expect_warning(
    model <- lift_MVG_to_MVG_ideal_observer(mvg_template),
    "deprecated"
  )
  expect_true(S7::S7_inherits(model, MVG_IdealObserver))
  expect_error(lift_MVG_to_MVG_ideal_observer(mvg_template, group = "category"))
})

