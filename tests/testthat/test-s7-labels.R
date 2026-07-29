context("S7 category-label access across representations, templates, and models")

test_that("models delegate category and cue labels through their category likelihood", {
  rep_a <- new_category_representation(category_labels = "A", cue_labels = c("F1", "F2"))
  rep_b <- new_category_representation(category_labels = "B", cue_labels = c("F1", "F2"))
  template <- new_category_representation_template(representations = list(A = rep_a, B = rep_b))
  model <- new_cognitive_model(category_template = template)

  expect_equal(get_category_labels(rep_a), "A")
  expect_equal(get_category_labels(template), c("A", "B"))
  expect_equal(get_category_labels(model), c("A", "B"))
  expect_equal(get_cue_labels(rep_a), c("F1", "F2"))
  expect_equal(get_cue_labels(model), c("F1", "F2"))
})

test_that("label information is nested under metadata$label_information", {
  rep <- new_category_representation(category_labels = "A", cue_labels = c("F1", "F2"))
  template <- new_category_representation_template(representations = list(A = rep))

  expect_equal(rep@metadata$label_information$category, "A")
  expect_equal(rep@metadata$label_information$cue, c("F1", "F2"))
  expect_equal(template@metadata$label_information$category, "A")
  expect_equal(template@metadata$label_information$cue, c("F1", "F2"))
})
