
test_that("Phase 2 adapters migrate legacy MVG/NIW/Exemplar rows to S7 templates", {
  legacy_mvg <- data.frame(category = factor(c("A", "B")))
  legacy_mvg$mu <- list(c(F1 = 0, F2 = 1), c(F1 = 2, F2 = 3))
  legacy_mvg$Sigma <- list(diag(c(1, 2)), diag(c(1, 1)))

  mvg_template <- as_s7_category_representation_template(legacy_mvg, family = "MVG")
  expect_true(S7::S7_inherits(mvg_template, MVBU_CategoryRepresentationTemplate))
  expect_equal(length(mvg_template@representations), 2)
  expect_true(all(vapply(mvg_template@representations, function(r) S7::S7_inherits(r, MVG_CategoryRepresentation), logical(1))))

  legacy_niw <- data.frame(category = factor(c("A", "B")))
  legacy_niw$m <- list(c(F1 = 0, F2 = 1), c(F1 = 2, F2 = 3))
  legacy_niw$kappa <- list(1, 1)
  legacy_niw$nu <- list(4, 4)
  legacy_niw$S <- list(diag(c(1, 1)), diag(c(2, 2)))

  niw_template <- as_s7_category_representation_template(legacy_niw, family = "NIW")
  expect_true(all(vapply(niw_template@representations, function(r) S7::S7_inherits(r, NIW_CategoryRepresentation), logical(1))))

  legacy_ex <- data.frame(category = factor(c("A", "B")))
  legacy_ex$exemplars <- list(
    matrix(c(0, 1, 1, 2), nrow = 2, byrow = TRUE, dimnames = list(NULL, c("F1", "F2"))),
    matrix(c(2, 3, 3, 4), nrow = 2, byrow = TRUE, dimnames = list(NULL, c("F1", "F2")))
  )

  ex_template <- as_s7_category_representation_template(legacy_ex, family = "EXEMPLAR")
  expect_true(all(vapply(ex_template@representations, function(r) S7::S7_inherits(r, Exemplar_CategoryRepresentation), logical(1))))
})

test_that("Phase 2 adapters migrate legacy model rows to S7 models", {
  legacy_mvg_model <- data.frame(category = factor(c("A", "B")))
  legacy_mvg_model$mu <- list(c(F1 = 0, F2 = 1), c(F1 = 2, F2 = 3))
  legacy_mvg_model$Sigma <- list(diag(c(1, 2)), diag(c(1, 1)))
  legacy_mvg_model$prior <- c(0.4, 0.6)
  legacy_mvg_model$lapse_rate <- c(0.1, 0.1)
  legacy_mvg_model$lapse_bias <- c(0.5, 0.5)

  mvg_model <- as_s7_mvg_ideal_observer(legacy_mvg_model)
  expect_true(S7::S7_inherits(mvg_model, MVG_IdealObserver))
  expect_equal(unname(get_category_prior(mvg_model)), c(0.4, 0.6), tolerance = MVBU_PROB_TOL)

  legacy_niw_model <- data.frame(category = factor(c("A", "B")))
  legacy_niw_model$m <- list(c(F1 = 0, F2 = 1), c(F1 = 2, F2 = 3))
  legacy_niw_model$kappa <- list(1, 1)
  legacy_niw_model$nu <- list(4, 4)
  legacy_niw_model$S <- list(diag(c(1, 1)), diag(c(2, 2)))
  legacy_niw_model$prior <- c(0.5, 0.5)
  legacy_niw_model$lapse_rate <- c(0, 0)
  legacy_niw_model$lapse_bias <- c(0.5, 0.5)

  niw_model <- as_s7_niw_ideal_adaptor(legacy_niw_model)
  expect_true(S7::S7_inherits(niw_model, NIW_IdealAdaptor))

  legacy_ex_model <- data.frame(category = factor(c("A", "B")))
  legacy_ex_model$exemplars <- list(
    matrix(c(0, 1, 1, 2), nrow = 2, byrow = TRUE, dimnames = list(NULL, c("F1", "F2"))),
    matrix(c(2, 3, 3, 4), nrow = 2, byrow = TRUE, dimnames = list(NULL, c("F1", "F2")))
  )

  ex_model <- as_s7_exemplar_model(legacy_ex_model)
  expect_true(S7::S7_inherits(ex_model, Exemplar_Model))
})



test_that("Phase 2 adapters normalize and validate constructor arguments consistently", {
  legacy_mvg <- data.frame(category = factor(c("A", "B")))
  legacy_mvg$mu <- list(c(F1 = 0, F2 = 1), c(F1 = 2, F2 = 3))
  legacy_mvg$Sigma <- list(diag(c(1, 2)), diag(c(1, 1)))

  # Lower-case family names are normalized.
  mvg_template <- as_s7_category_representation_template(legacy_mvg, family = "mvg")
  expect_true(S7::S7_inherits(mvg_template, MVBU_CategoryRepresentationTemplate))

  expect_error(
    as_s7_category_representation_template(legacy_mvg, family = c("MVG", "NIW")),
    "family must be a non-empty scalar character value"
  )

  expect_error(
    as_s7_mvg_ideal_observer(legacy_mvg, decision_rule = c("sampling", "map")),
    "decision_rule must be a non-empty scalar character value"
  )



  expect_error(
    as_s7_mvg_representations(legacy_mvg, category = ""),
    "category must be a non-empty scalar character value"
  )
})

test_that("Phase 2 migration pattern is reusable for MUVG/MNIX prototype families", {
  legacy_muvg <- data.frame(category = factor(c("A", "B")))
  legacy_muvg$component_mu <- list(c(F1 = 0, F2 = 1), c(F1 = 2, F2 = 3))
  legacy_muvg$component_sigma2 <- list(c(F1 = 1, F2 = 2), c(F1 = 1, F2 = 1))

  muvg_template <- as_s7_category_representation_template(legacy_muvg, family = "MUVG")
  expect_true(S7::S7_inherits(muvg_template, MVBU_CategoryRepresentationTemplate))
  expect_true(all(vapply(muvg_template@representations, function(r) S7::S7_inherits(r, MUVG_CategoryRepresentation), logical(1))))

  muvg_model <- as_s7_muvg_ideal_observer(legacy_muvg)
  expect_true(S7::S7_inherits(muvg_model, MUVG_IdealObserver))

  legacy_mnix <- data.frame(category = factor(c("A", "B")))
  legacy_mnix$component_m <- list(c(0, 1), c(2, 3))
  legacy_mnix$component_kappa <- list(c(1, 2), c(1, 2))
  legacy_mnix$component_nu <- list(c(3, 4), c(3, 4))
  legacy_mnix$component_sigma2 <- list(c(1, 2), c(1, 1))

  mnix_template <- as_s7_category_representation_template(legacy_mnix, family = "MNIX")
  expect_true(S7::S7_inherits(mnix_template, MVBU_CategoryRepresentationTemplate))
  expect_true(all(vapply(mnix_template@representations, function(r) S7::S7_inherits(r, MNIX_CategoryRepresentation), logical(1))))

  mnix_model <- as_s7_mnix_ideal_adaptor(legacy_mnix)
  expect_true(S7::S7_inherits(mnix_model, MNIX_IdealAdaptor))


})

test_that("Phase 2 migrated adapter outputs are S7-only (no S4 construction path)", {
  legacy_mvg <- data.frame(category = factor(c("A", "B")))
  legacy_mvg$mu <- list(c(F1 = 0, F2 = 1), c(F1 = 2, F2 = 3))
  legacy_mvg$Sigma <- list(diag(c(1, 2)), diag(c(1, 1)))

  legacy_niw <- data.frame(category = factor(c("A", "B")))
  legacy_niw$m <- list(c(F1 = 0, F2 = 1), c(F1 = 2, F2 = 3))
  legacy_niw$kappa <- list(1, 1)
  legacy_niw$nu <- list(4, 4)
  legacy_niw$S <- list(diag(c(1, 1)), diag(c(2, 2)))

  legacy_ex <- data.frame(category = factor(c("A", "B")))
  legacy_ex$exemplars <- list(
    matrix(c(0, 1, 1, 2), nrow = 2, byrow = TRUE, dimnames = list(NULL, c("F1", "F2"))),
    matrix(c(2, 3, 3, 4), nrow = 2, byrow = TRUE, dimnames = list(NULL, c("F1", "F2")))
  )

  objs <- list(
    as_s7_category_representation_template(legacy_mvg, family = "MVG"),
    as_s7_mvg_ideal_observer(legacy_mvg),
    as_s7_category_representation_template(legacy_niw, family = "NIW"),
    as_s7_niw_ideal_adaptor(legacy_niw),
    as_s7_category_representation_template(legacy_ex, family = "EXEMPLAR"),
    as_s7_exemplar_model(legacy_ex)
  )

  expect_true(all(vapply(objs, function(x) S7::S7_inherits(x, MVBU_Object), logical(1))))
  expect_false(any(vapply(objs, methods::is, logical(1), class2 = "S4")))
  expect_false(any(vapply(objs, isS4, logical(1))))
})

test_that("S7 accessors expose category priors and lapse biases directly", {
  rep_a <- new_category_representation(category_labels = "A", cue_labels = c("F1", "F2"))
  rep_b <- new_category_representation(category_labels = "B", cue_labels = c("F1", "F2"))
  template <- new_category_representation_template(representations = list(A = rep_a, B = rep_b))
  model <- new_cognitive_model(
    category_template = template,
    category_prior = c(A = 0.7, B = 0.3),
    lapse_rate = 0.1,
    lapse_bias = c(A = 0.8, B = 0.2)
  )

  expect_equal(get_category_prior(model, categories = c("A", "B")), c(0.7, 0.3))
  expect_equal(get_lapse_rate(model), 0.1)
  expect_equal(get_lapse_bias(model, categories = c("A", "B")), c(0.8, 0.2))
  expect_equal(get_cue_labels(model), c("F1", "F2"))
  expect_equal(get_category_labels(model), c("A", "B"))
  expect_equal(length(get_category_labels(model)), 2L)
  expect_equal(get_cue_labels(model, indices = 1), "F1")
  expect_equal(get_category_labels(model, indices = 1), "A")
})

test_that("legacy info helpers work with S7 objects through the new accessors", {
  rep_a <- new_category_representation(category_labels = "A", cue_labels = c("F1", "F2"))
  rep_b <- new_category_representation(category_labels = "B", cue_labels = c("F1", "F2"))
  template <- new_category_representation_template(representations = list(A = rep_a, B = rep_b))
  model <- new_cognitive_model(
    category_template = template,
    category_prior = c(A = 0.7, B = 0.3),
    lapse_rate = 0.1,
    lapse_bias = c(A = 0.8, B = 0.2)
  )

  expect_warning(expect_equal(get_priors_from_model(model, categories = c("A", "B")), c(0.7, 0.3)), class = "lifecycle_warning_deprecated")
  expect_warning(expect_equal(get_lapse_rate_from_model(model), 0.1), class = "lifecycle_warning_deprecated")
  expect_warning(expect_equal(get_lapse_biases_from_model(model, categories = c("A", "B")), c(0.8, 0.2)), class = "lifecycle_warning_deprecated")
  expect_warning(expect_equal(get_cue_labels_from_model(model), c("F1", "F2")), class = "lifecycle_warning_deprecated")
  expect_warning(expect_equal(get_category_labels_from_model(model), c("A", "B")), class = "lifecycle_warning_deprecated")
  expect_warning(expect_equal(get_nlevels_of_category_labels_from_model(model), 2L), class = "lifecycle_warning_deprecated")
  expect_warning(expect_equal(get_cue_labels_from_model(model, indices = 1), "F1"), class = "lifecycle_warning_deprecated")
  expect_warning(expect_equal(get_category_labels_from_model(model, indices = 1), "A"), class = "lifecycle_warning_deprecated")
  expect_warning(expect_equal(get_cue_dimensionality_from_model(model), 2L), class = "lifecycle_warning_deprecated")
  expect_warning(expect_null(get_perceptual_noise_from_model(model)), class = "lifecycle_warning_deprecated")
  expect_null(get_noise(model))

  model_with_noise <- new_cognitive_model(
    category_template = template,
    Sigma_noise = diag(2)
  )
  expect_equal(get_noise(model_with_noise), diag(2))
  expect_warning(expect_equal(get_perceptual_noise_from_model(model_with_noise), diag(2)), class = "lifecycle_warning_deprecated")
})
