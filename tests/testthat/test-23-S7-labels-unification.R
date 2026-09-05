test_that("label handling is unified and consistent across core S7 objects", {
  # 1. Category representation
  rep1 <- MVG_CategoryRepresentation(
    mu = c(VOT = 10, F0 = 200),
    Sigma = diag(c(25, 400)),
    metadata = list(label_information = list(cue = c("VOT", "F0"), category = "b"))
  )
  expect_equal(get_cue_labels(rep1), c("VOT", "F0"))
  expect_equal(get_category_labels(rep1), "b")
  expect_equal(get_labels(rep1), list(cue = c("VOT", "F0"), category = "b", group = character(0)))

  # 2. Category representation template
  rep2 <- MVG_CategoryRepresentation(
    mu = c(VOT = 50, F0 = 250),
    Sigma = diag(c(30, 450)),
    metadata = list(label_information = list(cue = c("VOT", "F0"), category = "p"))
  )
  tmpl <- new_category_representation_template(list(b = rep1, p = rep2))
  expect_equal(get_cue_labels(tmpl), c("VOT", "F0"))
  expect_equal(get_category_labels(tmpl), c("b", "p"))
  expect_equal(get_labels(tmpl), list(cue = c("VOT", "F0"), category = c("b", "p"), group = character(0)))

  # 3. Cognitive model
  model <- new_cognitive_model(category_template = tmpl)
  expect_equal(get_cue_labels(model), c("VOT", "F0"))
  expect_equal(get_category_labels(model), c("b", "p"))
  expect_equal(get_labels(model), list(cue = c("VOT", "F0"), category = c("b", "p"), group = character(0)))
})

test_that("label handling is unified and consistent for IdealAdaptorStanfitInput", {
  labels_expected <- list(
    cue = c("c1", "c2"),
    category = c("cat1", "cat2"),
    group = c("grp1", "grp2")
  )

  # Construction with unified labels list
  input_obj <- IdealAdaptorStanfitInput(
    data = data.frame(c1 = 1:4, c2 = 5:8),
    metadata = list(label_information = labels_expected)
  )

  expect_equal(input_obj@metadata$label_information, labels_expected)
  expect_false("labels" %in% names(S7::props(input_obj)))
  expect_false("cues" %in% names(S7::props(input_obj)))
  expect_false("category_levels" %in% names(S7::props(input_obj)))
  expect_false("group_levels" %in% names(S7::props(input_obj)))
  expect_equal(get_cue_labels(input_obj), c("c1", "c2"))
  expect_equal(get_category_labels(input_obj), c("cat1", "cat2"))
  expect_equal(get_group_labels(input_obj), c("grp1", "grp2"))
  expect_equal(get_group_labels(input_obj, include_prior = TRUE), c("prior", "grp1", "grp2"))
  expect_equal(get_labels(input_obj), labels_expected)
})

test_that("label handling is unified and consistent for MVBU_Stanfit", {
  labels_expected <- list(
    cue = c("c1", "c2"),
    category = c("cat1", "cat2"),
    group = c("grp1", "grp2")
  )

  fit_obj <- MVBU_Stanfit(
    data = data.frame(c1 = 1:4),
    metadata = list(label_information = labels_expected)
  )

  expect_equal(fit_obj@metadata$label_information, labels_expected)
  expect_false("labels" %in% names(S7::props(fit_obj)))
  expect_equal(get_cue_labels(fit_obj), c("c1", "c2"))
  expect_equal(get_category_labels(fit_obj), c("cat1", "cat2"))
  expect_equal(get_group_labels(fit_obj), c("grp1", "grp2"))
  expect_equal(get_group_labels(fit_obj, include_prior = TRUE), c("prior", "grp1", "grp2"))
  expect_equal(get_labels(fit_obj), labels_expected)
})

