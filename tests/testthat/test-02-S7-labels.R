test_that("models delegate category and cue labels through their category likelihood", {
  rep_a <- new_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2")
  )
  rep_b <- new_category_representation(
    category_labels = "B",
    cue_labels = c("F1", "F2")
  )
  template <- new_category_representation_template(
    representations = list(A = rep_a, B = rep_b)
  )
  model <- new_cognitive_model(category_template = template)

  expect_equal(get_category_labels(rep_a), "A")
  expect_equal(get_category_labels(template), c("A", "B"))
  expect_equal(get_category_labels(model), c("A", "B"))
  expect_equal(get_cue_labels(rep_a), c("F1", "F2"))
  expect_equal(get_cue_labels(model), c("F1", "F2"))
})

test_that("label information is nested under metadata$label_information", {
  rep <- new_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2")
  )
  template <- new_category_representation_template(
    representations = list(A = rep)
  )

  expect_equal(rep@metadata$label_information$category, "A")
  expect_equal(rep@metadata$label_information$cue, c("F1", "F2"))
  expect_equal(template@metadata$label_information$category, "A")
  expect_equal(template@metadata$label_information$cue, c("F1", "F2"))
})

test_that("label handling is unified and consistent across core S7 objects", {
  # 1. Category representation
  rep1 <- MVG_CategoryRepresentation(
    mu = c(VOT = 10, F0 = 200),
    Sigma = diag(c(25, 400)),
    metadata = list(label_information = list(
      cue = c("VOT", "F0"),
      category = "b",
      response_category = "b",
      group = character(0)
    ))
  )
  expect_equal(get_cue_labels(rep1), c("VOT", "F0"))
  expect_equal(get_category_labels(rep1), "b")
  expect_equal(get_response_category_labels(rep1), "b")
  expect_equal(
    get_labels(rep1),
    list(cue = c("VOT", "F0"), category = "b", response_category = "b", group = character(0))
  )

  # 2. Category representation template
  rep2 <- MVG_CategoryRepresentation(
    mu = c(VOT = 50, F0 = 250),
    Sigma = diag(c(30, 450)),
    metadata = list(label_information = list(
      cue = c("VOT", "F0"),
      category = "p",
      response_category = "p",
      group = character(0)
    ))
  )
  tmpl <- new_category_representation_template(list(b = rep1, p = rep2))
  expect_equal(get_cue_labels(tmpl), c("VOT", "F0"))
  expect_equal(get_category_labels(tmpl), c("b", "p"))
  expect_equal(get_response_category_labels(tmpl), c("b", "p"))
  expect_equal(
    get_labels(tmpl),
    list(cue = c("VOT", "F0"), category = c("b", "p"), response_category = c("b", "p"), group = character(0))
  )

  # 3. Cognitive model
  model <- new_cognitive_model(category_template = tmpl)
  expect_equal(get_cue_labels(model), c("VOT", "F0"))
  expect_equal(get_category_labels(model), c("b", "p"))
  expect_equal(get_response_category_labels(model), c("b", "p"))
  expect_equal(
    get_labels(model),
    list(cue = c("VOT", "F0"), category = c("b", "p"), response_category = c("b", "p"), group = character(0))
  )
})

test_that("label handling is unified and consistent for IdealAdaptorStanfitInput", {
  labels_expected <- list(
    cue = c("c1", "c2"),
    category = c("cat1", "cat2"),
    response_category = c("rcat1", "rcat2"),
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
  expect_equal(get_response_category_labels(input_obj), c("rcat1", "rcat2"))
  expect_equal(get_group_labels(input_obj), c("grp1", "grp2"))
  expect_equal(
    get_group_labels(input_obj, include_prior = TRUE),
    c("prior", "grp1", "grp2")
  )
  expect_equal(get_labels(input_obj), labels_expected)
})

test_that("label handling is unified and consistent for MVBU_Stanfit", {
  labels_expected <- list(
    cue = c("c1", "c2"),
    category = c("cat1", "cat2"),
    response_category = c("rcat1", "rcat2"),
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
  expect_equal(get_response_category_labels(fit_obj), c("rcat1", "rcat2"))
  expect_equal(get_group_labels(fit_obj), c("grp1", "grp2"))
  expect_equal(
    get_group_labels(fit_obj, include_prior = TRUE),
    c("prior", "grp1", "grp2")
  )
  expect_equal(get_labels(fit_obj), labels_expected)
})

test_that("set_labels updates label information correctly", {
  rep <- new_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2")
  )
  rep_updated <- set_labels(rep, category = "X", cue = c("cue1", "cue2"))
  expect_equal(get_category_labels(rep_updated), "X")
  expect_equal(get_cue_labels(rep_updated), c("cue1", "cue2"))

  meta <- list(other = 123)
  meta_updated <- set_labels(meta, cue = "F1", category = "A")
  expect_equal(meta_updated$other, 123)
  expect_equal(meta_updated$label_information$cue, "F1")
  expect_equal(meta_updated$label_information$category, "A")
})

test_that(".validate_requested_labels validates and errors informatively", {
  avail <- c("control", "exposure")
  expect_equal(.validate_requested_labels(NULL, avail), avail)
  expect_equal(.validate_requested_labels("control", avail), "control")
  expect_equal(.validate_requested_labels(c("control", "exposure"), avail), c("control", "exposure"))

  expect_error(
    .validate_requested_labels("c99", avail, label_type = "group"),
    "Requested group.*not found.*c99.*Available group.*control, exposure"
  )
})

