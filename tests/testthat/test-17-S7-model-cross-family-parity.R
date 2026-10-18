# Cross-family method parity test matrix
# Verifies that all concrete cognitive model families (UVG, NIX, MUVG, MNIX,
# MVG, NIW, Exemplar) adhere to the unified S7 method interface contract.

test_that("all cognitive model families satisfy the core generic method parity matrix", {
  families <- c("UVG", "NIX", "MUVG", "MNIX", "MVG", "NIW", "EXEMPLAR")

  for (fam in families) {
    nc <- if (fam %in% c("UVG", "NIX")) 1L else 2L
    m <- example_model(fam, n_cues = nc)

    # 1. Type and structure introspection
    expect_true(S7::S7_inherits(m, MVBU_CognitiveModel), info = fam)
    expect_true(is.character(get_family_name(m)), info = fam)
    expect_true(S7::S7_inherits(get_family(m), MVBU_Family), info = fam)

    cues <- get_cue_labels(m)
    expect_equal(length(cues), nc, info = fam)
    cats <- get_category_labels(m)
    expect_true(length(cats) >= 2L, info = fam)

    # 2. Parameters extraction
    pars <- get_parameters(m)
    expect_true(is.list(pars) || is.data.frame(pars), info = fam)

    # Generate test observation points
    set.seed(42)
    dat <- as.data.frame(matrix(rnorm(nc * 6L), nrow = 6L))
    colnames(dat) <- cues

    # 3. Likelihood computation
    lik <- likelihood(m, dat)
    expect_true(is.data.frame(lik) || is.matrix(lik), info = fam)
    expect_equal(nrow(lik), 6L, info = fam)
    expect_true(all(lik >= 0, na.rm = TRUE), info = fam)

    # 4. Posterior computation
    post <- posterior(m, dat)
    expect_true(is.data.frame(post) || is.matrix(post), info = fam)
    expect_equal(nrow(post), 6L, info = fam)
    # Probabilities should sum to 1 across categories for each observation
    row_sums <- rowSums(post[, cats, drop = FALSE])
    expect_equal(row_sums, rep(1, 6L), tolerance = 1e-6, info = fam)

    # 5. Categorization under decision rules
    # 5a. Deterministic MAP criterion
    cat_crit <- categorize(m, dat, decision_rule = "criterion")
    expect_true(is.character(cat_crit), info = fam)
    expect_equal(length(cat_crit), 6L, info = fam)
    expect_true(all(cat_crit %in% cats), info = fam)

    # 5b. Probability matching (proportional) with detailed return format
    cat_prop <- categorize(m, dat, decision_rule = "proportional", simplify = FALSE)
    expect_s3_class(cat_prop, "data.frame")
    expect_equal(nrow(cat_prop), 6L, info = fam)
    expect_true(all(c("response_category", "response_probability") %in% names(cat_prop)), info = fam)
    expect_true(all(cat_prop$response_probability >= 0 & cat_prop$response_probability <= 1), info = fam)

    # 5c. Posterior sampling
    set.seed(123)
    cat_samp <- categorize(m, dat, decision_rule = "sampling")
    expect_true(is.character(cat_samp), info = fam)
    expect_equal(length(cat_samp), 6L, info = fam)
    expect_true(all(cat_samp %in% cats), info = fam)

    # 6. Sample observations (generative synthesis)
    samp <- sample_observations(m, n = 8L)
    expect_s3_class(samp, "data.frame")
    expect_equal(nrow(samp), 8L, info = fam)
    expect_true("category" %in% names(samp), info = fam)
    expect_true(all(cues %in% names(samp)), info = fam)

    # 7. Model evaluation
    resp_vec <- sample(cats, size = 6L, replace = TRUE)
    ev_loglik <- evaluate_model(m, dat, response_category = resp_vec, method = "log_lik")
    expect_true(is.numeric(ev_loglik), info = fam)
    expect_true(is.finite(ev_loglik), info = fam)

    ev_acc <- evaluate_model(m, dat, response_category = resp_vec, method = "accuracy")
    expect_true(is.numeric(ev_acc), info = fam)
    expect_true(ev_acc >= 0 && ev_acc <= 1, info = fam)

    # 8. Printing and summary
    out_print <- capture.output(print(m))
    expect_true(length(out_print) > 0L, info = fam)
    out_sum <- capture.output(summary(m))
    expect_true(length(out_sum) > 0L, info = fam)

    # 9. Plotting methods
    p_cat <- plot_categories(m, interactive = FALSE)
    expect_s3_class(p_cat, "ggplot")

    p_dec <- plot_categorization_functions(m, interactive = FALSE)
    expect_s3_class(p_dec, "ggplot")
  }
})

test_that("cross-family template parity holds across representation families", {
  families <- c("UVG", "NIX", "MUVG", "MNIX", "MVG", "NIW", "EXEMPLAR")

  for (fam in families) {
    nc <- if (fam %in% c("UVG", "NIX")) 1L else 2L
    tpl <- example_category_representation_template(fam, n_cues = nc)

    expect_true(S7::S7_inherits(tpl, MVBU_CategoryRepresentationTemplate), info = fam)
    expect_true(is.character(get_family_name(tpl)), info = fam)
    expect_true(S7::S7_inherits(get_family(tpl), MVBU_Family), info = fam)
    expect_equal(length(get_cue_labels(tpl)), nc, info = fam)
    expect_true(length(get_category_labels(tpl)) >= 2L, info = fam)

    # Check that individual representations inherit from MVBU_CategoryRepresentation
    for (cat_name in get_category_labels(tpl)) {
      rep_obj <- tpl@representations[[cat_name]]
      expect_true(S7::S7_inherits(rep_obj, MVBU_CategoryRepresentation), info = paste(fam, cat_name))
      expect_equal(get_cue_labels(rep_obj), get_cue_labels(tpl), info = paste(fam, cat_name))
    }
  }
})

test_that("early and late multi-cue integration models satisfy core generic method parity", {
  set.seed(42)
  cues <- c("cue1", "cue2")
  d_train <- data.frame(
    category = factor(rep(c("A", "B"), each = 20)),
    cue1 = c(rnorm(20, -1, 1), rnorm(20, 1, 1)),
    cue2 = c(rnorm(20, -0.5, 1), rnorm(20, 0.5, 1))
  )

  mnix_tpl <- new_mnix_category_representation_template_from_data(
    d_train,
    category = "category", cues = cues
  )
  d_proj <- d_train
  d_proj$.comp <- 0.6 * d_train$cue1 + 0.4 * d_train$cue2
  nix_tpl <- new_nix_category_representation_template_from_data(
    d_proj,
    category = "category", cues = ".comp"
  )

  models <- list(
    MNIX_IdealIntegrator = new_mnix_ideal_adaptor_ideal_integrator(category_template = mnix_tpl),
    MNIX_WeightedIntegrator = new_mnix_ideal_adaptor_weighted_integrator(
      category_template = mnix_tpl, weights = c(0.6, 0.4)
    ),
    MNIX_CategoryWeightedIntegrator = new_mnix_ideal_adaptor_category_weighted_integrator(
      category_template = mnix_tpl, weights = rbind(A = c(0.6, 0.4), B = c(0.3, 0.7))
    ),
    IdealIntegrator_MNIX = new_ideal_integrator_mnix_ideal_adaptor(
      category_template = mnix_tpl, cues = cues
    ),
    WeightedIntegrator_MNIX = new_weighted_integrator_mnix_ideal_adaptor(
      category_template = mnix_tpl, weights = c(0.6, 0.4), cues = cues
    ),
    CategoryWeightedIntegrator_MNIX = new_category_weighted_integrator_mnix_ideal_adaptor(
      category_template = mnix_tpl, weights = rbind(A = c(0.6, 0.4), B = c(0.3, 0.7)), cues = cues
    )
  )

  dat <- data.frame(
    cue1 = c(-1.2, 0.0, 1.1, -0.5),
    cue2 = c(-0.8, 0.2, 0.9, -0.1)
  )
  cats <- c("A", "B")

  for (mod_name in names(models)) {
    m <- models[[mod_name]]
    expect_true(S7::S7_inherits(m, MVBU_CognitiveModel), info = mod_name)

    # 1. Likelihood
    lik <- likelihood(m, dat)
    expect_true(is.data.frame(lik) || is.matrix(lik), info = mod_name)
    expect_equal(nrow(lik), 4L, info = mod_name)
    expect_true(all(lik >= 0, na.rm = TRUE), info = mod_name)

    # 2. Posterior (sums to 1)
    post <- posterior(m, dat)
    expect_true(is.data.frame(post) || is.matrix(post), info = mod_name)
    expect_equal(nrow(post), 4L, info = mod_name)
    expect_equal(rowSums(post[, cats, drop = FALSE]), rep(1, 4L), tolerance = 1e-6, info = mod_name)

    # 3. Decision rules
    cat_crit <- categorize(m, dat, decision_rule = "criterion")
    expect_equal(length(cat_crit), 4L, info = mod_name)
    expect_true(all(cat_crit %in% cats), info = mod_name)

    cat_prop <- categorize(m, dat, decision_rule = "proportional", simplify = FALSE)
    expect_s3_class(cat_prop, "data.frame")
    expect_equal(nrow(cat_prop), 4L, info = mod_name)
    expect_true(all(cat_prop$response_probability >= 0 & cat_prop$response_probability <= 1), info = mod_name)

    cat_samp <- categorize(m, dat, decision_rule = "sampling")
    expect_equal(length(cat_samp), 4L, info = mod_name)
    expect_true(all(cat_samp %in% cats), info = mod_name)

    # 4. Generative sampling
    samp <- sample_observations(m, n = 6L)
    expect_s3_class(samp, "data.frame")
    expect_equal(nrow(samp), 6L, info = mod_name)
    expect_true("category" %in% names(samp), info = mod_name)

    # 5. Model evaluation
    resp_vec <- c("A", "A", "B", "B")
    ev_loglik <- evaluate_model(m, dat, response_category = resp_vec, method = "log_lik")
    expect_true(is.numeric(ev_loglik) && is.finite(ev_loglik), info = mod_name)

    ev_acc <- evaluate_model(m, dat, response_category = resp_vec, method = "accuracy")
    expect_true(is.numeric(ev_acc) && ev_acc >= 0 && ev_acc <= 1, info = mod_name)
  }
})

test_that("early and late integration models support analytical updating via update_template", {
  set.seed(123)
  cues <- c("cue1", "cue2")
  d_train <- data.frame(
    category = factor(rep(c("A", "B"), each = 20)),
    cue1 = c(rnorm(20, -1, 1), rnorm(20, 1, 1)),
    cue2 = c(rnorm(20, -0.5, 1), rnorm(20, 0.5, 1))
  )

  mnix_tpl <- new_mnix_category_representation_template_from_data(
    d_train,
    category = "category", cues = cues
  )
  d_proj <- d_train
  d_proj$.comp <- 0.5 * d_train$cue1 + 0.5 * d_train$cue2
  nix_tpl <- new_nix_category_representation_template_from_data(
    d_proj,
    category = "category", cues = ".comp"
  )

  d_new <- data.frame(
    category = factor(c("A", "A", "B", "B")),
    cue1 = c(-0.8, -1.2, 1.2, 0.8),
    cue2 = c(-0.4, -0.6, 0.6, 0.4)
  )

  models <- list(
    MNIX_IdealIntegrator = new_mnix_ideal_adaptor_ideal_integrator(category_template = mnix_tpl),
    MNIX_WeightedIntegrator = new_mnix_ideal_adaptor_weighted_integrator(
      category_template = mnix_tpl, weights = c(0.5, 0.5)
    ),
    MNIX_CategoryWeightedIntegrator = new_mnix_ideal_adaptor_category_weighted_integrator(
      category_template = mnix_tpl, weights = rbind(A = c(0.5, 0.5), B = c(0.5, 0.5))
    ),
    IdealIntegrator_MNIX = new_ideal_integrator_mnix_ideal_adaptor(
      category_template = mnix_tpl, cues = cues
    ),
    WeightedIntegrator_MNIX = new_weighted_integrator_mnix_ideal_adaptor(
      category_template = mnix_tpl, weights = c(0.5, 0.5), cues = cues
    ),
    CategoryWeightedIntegrator_MNIX = new_category_weighted_integrator_mnix_ideal_adaptor(
      category_template = mnix_tpl, weights = rbind(A = c(0.5, 0.5), B = c(0.5, 0.5)), cues = cues
    )
  )

  for (mod_name in names(models)) {
    m <- models[[mod_name]]
    m_up <- update_template(m, d_new)
    expect_true(S7::S7_inherits(m_up, S7::S7_class(m)), info = mod_name)
    # Ensure kappa increased by exposure observations
    reps_orig <- get_category_representations(m)
    reps_up <- get_category_representations(m_up)
    expect_true(all(reps_up$A@kappa > reps_orig$A@kappa), label = paste(mod_name, "A kappa updated"))
    expect_true(all(reps_up$B@kappa > reps_orig$B@kappa), label = paste(mod_name, "B kappa updated"))
  }
})

test_that("mathematical parity holds across integration variants under matched weights", {
  set.seed(789)
  cues <- c("cue1", "cue2")
  d_train <- data.frame(
    category = factor(rep(c("A", "B"), each = 25)),
    cue1 = c(rnorm(25, -1, 1), rnorm(25, 1, 1)),
    cue2 = c(rnorm(25, -0.5, 1), rnorm(25, 0.5, 1))
  )

  mnix_tpl <- new_mnix_category_representation_template_from_data(
    d_train,
    category = "category", cues = cues
  )
  w <- c(0.7, 0.3)
  d_proj <- d_train
  d_proj$.comp <- w[1] * d_train$cue1 + w[2] * d_train$cue2
  nix_tpl <- new_nix_category_representation_template_from_data(
    d_proj,
    category = "category", cues = ".comp"
  )

  dat <- data.frame(
    cue1 = c(-1.5, 0.0, 1.2, 0.5),
    cue2 = c(-0.5, 0.1, 0.8, -0.3)
  )

  # Early integration parity: Ideal vs Weighted vs CategoryWeighted (with matching weights)
  m_early_ideal <- new_ideal_integrator_mnix_ideal_adaptor(
    category_template = mnix_tpl, cues = cues
  )
  m_early_weighted_equal <- new_weighted_integrator_mnix_ideal_adaptor(
    category_template = mnix_tpl, weights = c(0.5, 0.5), cues = cues
  )
  m_early_weighted <- new_weighted_integrator_mnix_ideal_adaptor(
    category_template = mnix_tpl, weights = w, cues = cues
  )
  m_early_cat <- new_category_weighted_integrator_mnix_ideal_adaptor(
    category_template = mnix_tpl, weights = rbind(A = w, B = w), cues = cues
  )

  lik_ideal <- likelihood(m_early_ideal, dat)
  lik_weighted_equal <- likelihood(m_early_weighted_equal, dat)
  lik_weighted <- likelihood(m_early_weighted, dat)
  lik_cat <- likelihood(m_early_cat, dat)

  expect_equal(lik_ideal, lik_weighted_equal, tolerance = 1e-10)
  expect_equal(lik_weighted, lik_cat, tolerance = 1e-10)

  post_ideal <- posterior(m_early_ideal, dat)
  post_weighted_equal <- posterior(m_early_weighted_equal, dat)
  post_weighted <- posterior(m_early_weighted, dat)
  post_cat <- posterior(m_early_cat, dat)

  expect_equal(post_ideal, post_weighted_equal, tolerance = 1e-10)
  expect_equal(post_weighted, post_cat, tolerance = 1e-10)

  # Late integration parity: Weighted vs CategoryWeighted (with matching rows)
  m_late_weighted <- new_mnix_ideal_adaptor_weighted_integrator(
    category_template = mnix_tpl, weights = w
  )
  m_late_cat <- new_mnix_ideal_adaptor_category_weighted_integrator(
    category_template = mnix_tpl, weights = rbind(A = w, B = w)
  )

  lik_lw <- likelihood(m_late_weighted, dat)
  lik_lc <- likelihood(m_late_cat, dat)
  expect_equal(lik_lw, lik_lc, tolerance = 1e-10)

  post_lw <- posterior(m_late_weighted, dat)
  post_lc <- posterior(m_late_cat, dat)
  expect_equal(post_lw, post_lc, tolerance = 1e-10)
})

test_that(
  paste0(
    ".new_staninput constructs typed Staninput for all integration models"
  ),
  {
    set.seed(999)
    cues <- c("cue1", "cue2")
    exposure <- data.frame(
      category = factor(rep(c("A", "B"), each = 15)),
      group = factor(rep("G1", 30)),
      cue1 = rnorm(30),
      cue2 = rnorm(30)
    )
    test_data <- data.frame(
      group = factor(rep("G1", 10)),
      cue1 = rnorm(10),
      cue2 = rnorm(10),
      response_category = factor(sample(c("A", "B"), 10, replace = TRUE))
    )

    sfi <- new_ideal_adaptor_stanfit_input(
      exposure = exposure,
      test = test_data,
      cues = cues
    )

    # 1. MNIX Ideal Integrator
    si_mnix_ideal <- MVBeliefUpdatr:::.new_staninput(
      sfi,
      stanmodel = "MNIX_ideal_adaptor_ideal_integrator"
    )
    expect_true(S7::S7_inherits(si_mnix_ideal, IdealAdaptorStaninput))
    expect_equal(
      get_family_name(si_mnix_ideal),
      "MNIX_ideal_integrator"
    )
    expect_equal(si_mnix_ideal@values$K, 2L)

    # 2. MNIX Weighted Integrator
    si_mnix_weighted <- MVBeliefUpdatr:::.new_staninput(
      sfi,
      stanmodel = "MNIX_ideal_adaptor_weighted_integrator"
    )
    expect_true(S7::S7_inherits(si_mnix_weighted, IdealAdaptorStaninput))
    expect_equal(
      get_family_name(si_mnix_weighted),
      "MNIX_weighted_integrator"
    )

    # 3. MNIX Category Weighted Integrator
    si_mnix_cat <- MVBeliefUpdatr:::.new_staninput(
      sfi,
      stanmodel = "MNIX_ideal_adaptor_category_weighted_integrator"
    )
    expect_true(S7::S7_inherits(si_mnix_cat, IdealAdaptorStaninput))
    expect_equal(
      get_family_name(si_mnix_cat),
      "MNIX_category_weighted_integrator"
    )

    # 4. Early MNIX Ideal Integrator
    si_early_ideal <- MVBeliefUpdatr:::.new_staninput(
      sfi,
      stanmodel = "ideal_integrator_MNIX_ideal_adaptor"
    )
    expect_true(S7::S7_inherits(si_early_ideal, IdealAdaptorStaninput))
    expect_equal(
      get_family_name(si_early_ideal),
      "ideal_integrator_MNIX"
    )
    # Should have 2x2 scatter matrix for NIW summary
    expect_equal(dim(si_early_ideal@values$css_exposure)[3:4], c(2L, 2L))

    # Retired stanmodel name is rejected
    expect_error(
      MVBeliefUpdatr:::.new_staninput(
        sfi,
        stanmodel = "MNIX_ideal_adaptor"
      )
    )

    # 5. Early MNIX Weighted Integrator
    si_early_weighted <- MVBeliefUpdatr:::.new_staninput(
      sfi,
      stanmodel = "weighted_integrator_MNIX_ideal_adaptor",
      fixed_parameters_transformed = list(weights = c(0.7, 0.3))
    )
    expect_true(S7::S7_inherits(si_early_weighted, IdealAdaptorStaninput))
    expect_equal(
      get_family_name(si_early_weighted),
      "weighted_integrator_MNIX"
    )
    expect_equal(
      as.numeric(si_early_weighted@values$weights),
      c(0.7, 0.3)
    )
    expect_equal(dim(si_early_weighted@values$css_exposure)[3:4], c(2L, 2L))

    # 6. Early MNIX Category Weighted Integrator
    si_early_cat <- MVBeliefUpdatr:::.new_staninput(
      sfi,
      stanmodel = "category_weighted_integrator_MNIX_ideal_adaptor"
    )
    expect_true(S7::S7_inherits(si_early_cat, IdealAdaptorStaninput))
    expect_equal(
      get_family_name(si_early_cat),
      "category_weighted_integrator_MNIX"
    )
    expect_equal(dim(si_early_cat@values$css_exposure)[3:4], c(2L, 2L))
  }
)
