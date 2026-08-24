
test_that("foundation base classes instantiate", {
  base_obj <- new_mvbu_object()
  expect_true(S7::S7_inherits(base_obj, MVBU_Object))
  expect_equal(class(base_obj)[1], "MVBU_Object")

  rep_obj <- new_category_representation(
    category_labels = c("A", "B"),
    cue_labels = c("F1", "F2")
  )
  expect_true(S7::S7_inherits(rep_obj, MVBU_CategoryRepresentation))
  expect_equal(class(rep_obj)[1], "MVBU_CategoryRepresentation")

  dist_obj <- new_model_distribution(model_family = "NIW")
  expect_true(S7::S7_inherits(dist_obj, MVBU_ModelDistribution))
  expect_equal(get_model_family(dist_obj), "NIW")
})

test_that("validation helpers behave as expected", {
  rep_obj <- new_category_representation(
    category_labels = c("A"),
    cue_labels = c("F1")
  )

  expect_true(validate_object(rep_obj))
  expect_true(.is_valid(rep_obj))
  expect_false(.is_valid(1L))
})

test_that("generic stubs are available and family-specific constructor policy is enforced", {
  rep_obj <- new_category_representation(
    category_labels = c("A"),
    cue_labels = c("F1")
  )

  expect_true(is.function(get_model_family))
  expect_true(is.function(get_cue_labels))
  expect_equal(get_category_labels(rep_obj), "A")
  expect_equal(get_cue_labels(rep_obj), "F1")
})

test_that("core generic aliases dispatch on base classes", {
  rep_obj <- new_category_representation(
    category_labels = c("A"),
    cue_labels = c("F1")
  )

  rep_set <- new_category_representation_template(representations = list(A = rep_obj))
  model <- new_cognitive_model(category_template = rep_set)

  expect_identical(construct_mvbu(rep_obj), rep_obj)
  expect_true(validate_mvbu(rep_obj))

  summary_stub <- summarize_mvbu(rep_obj)
  expect_equal(summary_stub$class, "MVBU_CategoryRepresentation")
  expect_equal(summary_stub$model_family, "MVBU_CategoryRepresentation")

  expect_invisible(print_mvbu(rep_obj))

  expect_error(
    posterior(model, data.frame(x = 1), categories = NULL),
    "category_likelihood not implemented|not yet implemented|Can't find method"
  )
  expect_error(
    categorize(model, data.frame(x = 1), decision_rule = NULL),
    "category_likelihood not implemented|not yet implemented|Can't find method"
  )
  expect_error(
    plot_prep_mvbu(rep_obj),
    "not yet implemented"
  )
})

test_that("cognitive model constructor resolves template arguments without ambiguity", {
  rep_a <- new_category_representation(category_labels = "A", cue_labels = "F1")
  rep_b <- new_category_representation(category_labels = "B", cue_labels = "F1")
  rep_set <- new_category_representation_template(representations = list(A = rep_a, B = rep_b))

  expect_error(
    new_cognitive_model(category_template = NULL),
    "category_template must be supplied"
  )
})

test_that("cognitive model defaults and distribution group labels are consistent", {
  rep_a <- new_category_representation(category_labels = "A", cue_labels = "F1")
  rep_b <- new_category_representation(category_labels = "B", cue_labels = "F1")
  rep_set <- new_category_representation_template(representations = list(A = rep_a, B = rep_b))

  model <- new_cognitive_model(
    category_template = rep_set,
    decision_rule = "sampling",
    lapse_rate = 0.1
  )

  expect_true(S7::S7_inherits(model, MVBU_CognitiveModel))
  expect_true(S7::S7_inherits(model@category_template, MVBU_CategoryRepresentationTemplate))
  expect_equal(length(model@category_template@representations), 2)
  expect_identical(model@category_template, get_category_template(model))
  expect_equal(length(get_category_prior(model)), 2)
  expect_equal(unname(get_category_prior(model)), c(0.5, 0.5), tolerance = MVBU_PROB_TOL)
  expect_equal(unname(model@lapse_behavior$lapse_bias), unname(get_category_prior(model)), tolerance = MVBU_PROB_TOL)
  expect_equal(get_group_labels(model), character(0))

  distribution <- new_model_distribution(model_family = "NIW", group_label = "fit-1")
  expect_equal(get_group_labels(distribution), "fit-1")
})

test_that("accessor generics provide consistent introspection for representations, templates, models, and distributions", {
  rep <- new_category_representation(category_labels = "A", cue_labels = "F1")
  template <- new_category_representation_template(representations = list(A = rep))
  model <- new_cognitive_model(
    category_template = template,
    category_prior = c(A = 1),
    lapse_rate = 0.1,
    lapse_bias = c(A = 1)
  )
  distribution <- new_model_distribution(model_family = "NIW", group_label = "fit-1")

  expect_identical(get_metadata(rep), rep@metadata)
  expect_identical(get_category_representations(template), template@representations)
  expect_equal(get_category_template(model), model@category_template)
  expect_identical(get_category_template(model), model@category_template)
  expect_identical(get_category_representations(model), get_category_representations(get_category_template(model)))
  expect_equal(get_category_prior(model), c(A = 1))
  expect_equal(get_lapse_rate(model), 0.1)
  expect_equal(get_lapse_bias(model), c(A = 1))
  expect_equal(get_group_labels(distribution), "fit-1")
  expect_equal(get_cue_labels(rep, indices = 1), "F1")
  expect_equal(get_category_labels(template, indices = 1), "A")
})

test_that("cognitive model validators enforce category_prior and lapse_bias constraints", {
  rep_a <- new_category_representation(category_labels = "A", cue_labels = "F1")
  rep_b <- new_category_representation(category_labels = "B", cue_labels = "F1")
  rep_set <- new_category_representation_template(representations = list(A = rep_a, B = rep_b))

  expect_error(
    new_cognitive_model(
      category_template = rep_set,
      category_prior = c(1),
      lapse_bias = c(0.5, 0.5)
    ),
    "category_prior length must match"
  )

  expect_error(
    new_cognitive_model(
      category_template = rep_set,
      category_prior = c(0.5, 0.5),
      lapse_bias = c(0.8, 0.3)
    ),
    "lapse_bias entries must sum to 1"
  )
})

test_that("models store noise and lapse behavior and expose posterior functions by treatment", {
  rep_a <- new_category_representation(category_labels = "A", cue_labels = c("F1", "F2"))
  rep_b <- new_category_representation(category_labels = "B", cue_labels = c("F1", "F2"))
  rep_set <- new_category_representation_template(representations = list(A = rep_a, B = rep_b))

  model <- new_cognitive_model(
    category_template = rep_set,
    category_prior = c(A = 0.7, B = 0.3),
    lapse_rate = 0.1,
    lapse_bias = c(A = 0.8, B = 0.2),
    Sigma_noise = c(0.1, 0.2),
    noise_treatment = "marginalize",
    lapse_treatment = "sample"
  )

  expect_equal(model@noise_behavior$Sigma_noise, diag(c(0.1, 0.2)))
  expect_equal(model@noise_behavior$noise_treatment, "marginalize")
  expect_equal(model@lapse_behavior$lapse_rate, 0.1)
  expect_equal(model@lapse_behavior$lapse_bias, c(A = 0.8, B = 0.2))
  expect_equal(model@lapse_behavior$lapse_treatment, "sample")
  expect_true(is.function(get_category_posterior_function(model)))
  expect_true(is.function(get_category_posterior_function(model, noise_treatment = "no_noise", lapse_treatment = "no_lapses")))
  expect_true(is.function(get_category_posterior_function(model, noise_treatment = "marginalize", lapse_treatment = "sample")))
  expect_true(is.list(model@category_posterior_functions))
  expect_true(all(c("marginalize__sample", "no_noise__no_lapses") %in% names(model@category_posterior_functions)))

  expect_error(
    new_cognitive_model(category_template = rep_set, Sigma_noise = c(0.1)),
    "Sigma_noise length must match"
  )
})

test_that("sample-based noise treatment is available in posterior computation", {
  rep_a <- new_mvg_category_representation(
    category_labels = "A",
    cue_labels = "F1",
    mu = 0,
    Sigma = matrix(1, nrow = 1, ncol = 1)
  )
  rep_b <- new_mvg_category_representation(
    category_labels = "B",
    cue_labels = "F1",
    mu = 1,
    Sigma = matrix(1, nrow = 1, ncol = 1)
  )
  rep_set <- new_category_representation_template(representations = list(A = rep_a, B = rep_b))

  model <- new_mvg_ideal_observer(
    category_template = rep_set,
    category_prior = c(A = 0.5, B = 0.5),
    Sigma_noise = c(0.5),
    noise_treatment = "sample"
  )

  posterior_fn <- get_category_posterior_function(model, noise_treatment = "sample", lapse_treatment = "no_lapses")
  set.seed(123)
  posterior_matrix <- posterior_fn(data.frame(F1 = 0))

  expect_true(is.matrix(posterior_matrix))
  expect_equal(ncol(posterior_matrix), 2)
  expect_equal(colnames(posterior_matrix), c("A", "B"))
  expect_true(all(abs(rowSums(posterior_matrix) - 1) < MVBU_PROB_TOL))
  expect_true(all(posterior_matrix >= 0))
})

test_that("sample-based noise treatment inflates the effective variance", {
  rep <- new_uvg_category_representation(
    category_labels = "A",
    cue_labels = "F1",
    mu = 0,
    sigma2 = 1
  )
  likelihood_fn <- get_category_likelihood_function(rep)

  set.seed(123)
  sample_lik <- likelihood_fn(matrix(0, nrow = 1, ncol = 1), log = TRUE, noise_treatment = "sample", Sigma_noise = matrix(1, 1, 1))

  set.seed(123)
  noise_draw <- mvtnorm::rmvnorm(n = 1, mean = rep(0, 1), sigma = matrix(1, 1, 1))
  expected_lik <- stats::dnorm(as.numeric(noise_draw), mean = 0, sd = sqrt(2), log = TRUE)

  expect_equal(sample_lik, expected_lik)
})

test_that("MUVG representation schema validators and constructor defaults work", {
  muvg_rep <- new_muvg_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    component_mu = c(0, 1),
    component_sigma2 = c(1, 1)
  )

  expect_true(S7::S7_inherits(muvg_rep, MUVG_CategoryRepresentation))
  params <- get_parameters(muvg_rep)
  expect_equal(params$component_weights, c(0.5, 0.5), tolerance = MVBU_PROB_TOL)

  expect_error(
    new_muvg_category_representation(
      category_labels = "A",
      cue_labels = c("F1", "F2"),
      component_mu = c(0, 1),
      component_sigma2 = c(1, 1),
      component_weights = c(0.9, 0.2)
    ),
    "component_weights entries must sum to 1"
  )

  expect_error(
    new_muvg_category_representation(
      category_labels = "A",
      cue_labels = "F1",
      component_mu = c(0, 1),
      component_sigma2 = c(1, 1)
    ),
    "cue_labels length must match"
  )
})

test_that("MNIX representation schema validators and typed model constructor work", {
  mnix_rep <- new_mnix_category_representation(
    category_labels = "A",
    cue_labels = "F1",
    component_m = c(0, 1),
    component_kappa = c(1, 2),
    component_nu = c(3, 4),
    component_sigma2 = c(1, 2)
  )

  expect_true(S7::S7_inherits(mnix_rep, MNIX_CategoryRepresentation))
  mnix_params <- get_parameters(mnix_rep)
  expect_equal(mnix_params$component_weights, c(0.5, 0.5), tolerance = MVBU_PROB_TOL)

  mnix_rep_set <- new_category_representation_template(representations = list(A = mnix_rep))
  mnix_model <- new_mnix_ideal_adaptor(category_template = mnix_rep_set)
  expect_true(S7::S7_inherits(mnix_model, MNIX_IdealAdaptor))
  expect_equal(get_category_prior(mnix_model), c(A = 1), tolerance = MVBU_PROB_TOL)
  muvg_rep <- new_muvg_category_representation(
    category_labels = "B",
    cue_labels = c("F1", "F2"),
    component_mu = c(0, 1),
    component_sigma2 = c(1, 1)
  )

  expect_error(
    new_mnix_ideal_adaptor(
      category_template = new_category_representation_template(representations = list(A = muvg_rep))
    ),
    "must inherit from MNIX category representation class"
  )

  expect_error(
    new_mnix_category_representation(
      category_labels = "A",
      cue_labels = "F1",
      component_m = c(0, 1),
      component_kappa = c(1),
      component_nu = c(3, 4),
      component_sigma2 = c(1, 2),
      component_weights = c(0.5, 0.5)
    ),
    "MNIX component parameter vectors must all have equal length"
  )

  mnix_rep_b <- new_mnix_category_representation(
    category_labels = "B",
    cue_labels = "F1",
    component_m = c(2, 3),
    component_kappa = c(1, 2),
    component_nu = c(3, 4),
    component_sigma2 = c(1, 1)
  )
  mnix_model_bi <- new_mnix_ideal_adaptor(
    category_template = new_category_representation_template(
      representations = list(A = mnix_rep, B = mnix_rep_b)
    )
  )
  expect_true(S7::S7_inherits(mnix_model_bi, MNIX_IdealAdaptor))
})

test_that("MUVG typed ideal observer constructor enforces representation class", {
  muvg_rep_a <- new_muvg_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    component_mu = c(0, 1),
    component_sigma2 = c(1, 1)
  )
  muvg_rep_b <- new_muvg_category_representation(
    category_labels = "B",
    cue_labels = c("F1", "F2"),
    component_mu = c(2, 3),
    component_sigma2 = c(1, 1)
  )

  muvg_rep_set <- new_category_representation_template(representations = list(A = muvg_rep_a, B = muvg_rep_b))
  muvg_model <- new_muvg_ideal_observer(category_template = muvg_rep_set)
  expect_true(S7::S7_inherits(muvg_model, MUVG_IdealObserver))
  expect_equal(get_category_prior(muvg_model), c(A = 0.5, B = 0.5), tolerance = MVBU_PROB_TOL)

  plain_rep <- new_category_representation(category_labels = "X", cue_labels = "F1")
  expect_error(
    new_muvg_ideal_observer(
      category_template = new_category_representation_template(representations = list(X = plain_rep))
    ),
    "must inherit from MUVG category representation class"
  )
})

test_that("MVG/NIW/UVG/NIX/Exemplar typed constructors and validators work", {
  mvg_rep <- new_mvg_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    mu = c(0, 1),
    Sigma = matrix(c(1, 0, 0, 2), nrow = 2)
  )
  expect_true(S7::S7_inherits(mvg_rep, MVG_CategoryRepresentation))
  mvg_params <- get_parameters(mvg_rep)
  expect_equal(mvg_params$mu, c(0, 1))

  niw_rep <- new_niw_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    m = c(0, 1),
    kappa = 1,
    nu = 3,
    S = matrix(c(1, 0, 0, 1), nrow = 2)
  )
  expect_true(S7::S7_inherits(niw_rep, NIW_CategoryRepresentation))
  niw_params <- get_parameters(niw_rep)
  expect_equal(niw_params$kappa, 1)

  uvg_rep <- new_uvg_category_representation(
    category_labels = "A",
    cue_labels = "F1",
    mu = 0,
    sigma2 = 1
  )
  expect_true(S7::S7_inherits(uvg_rep, UVG_CategoryRepresentation))

  nix_rep <- new_nix_category_representation(
    category_labels = "A",
    cue_labels = "F1",
    m = 0,
    kappa = 1,
    nu = 2,
    sigma2 = 1
  )
  expect_true(S7::S7_inherits(nix_rep, NIX_CategoryRepresentation))

  mnix_rep_lik <- new_mnix_category_representation(
    category_labels = "A",
    cue_labels = "F1",
    component_m = c(0, 1),
    component_kappa = c(1, 2),
    component_nu = c(3, 4),
    component_sigma2 = c(1, 2)
  )
  expect_true(S7::S7_inherits(mnix_rep_lik, MNIX_CategoryRepresentation))

  ex_rep <- new_exemplar_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    exemplars = matrix(c(0, 1, 1, 2), nrow = 2, byrow = TRUE)
  )
  expect_true(S7::S7_inherits(ex_rep, Exemplar_CategoryRepresentation))
  ex_params <- get_parameters(ex_rep)
  expect_equal(ex_params$exemplar_weights, c(0.5, 0.5), tolerance = MVBU_PROB_TOL)

  mvg_model <- new_mvg_ideal_observer(
    category_template = new_category_representation_template(representations = list(A = mvg_rep))
  )
  expect_true(S7::S7_inherits(mvg_model, MVG_IdealObserver))

  niw_model <- new_niw_ideal_adaptor(
    category_template = new_category_representation_template(representations = list(A = niw_rep))
  )
  expect_true(S7::S7_inherits(niw_model, NIW_IdealAdaptor))

  uvg_model <- new_uvg_ideal_observer(
    category_template = new_category_representation_template(representations = list(A = uvg_rep))
  )
  expect_true(S7::S7_inherits(uvg_model, UVG_IdealObserver))

  nix_model <- new_nix_ideal_adaptor(
    category_template = new_category_representation_template(representations = list(A = nix_rep))
  )
  expect_true(S7::S7_inherits(nix_model, NIX_IdealAdaptor))

  ex_model <- new_exemplar_model(
    category_template = new_category_representation_template(representations = list(A = ex_rep))
  )
  expect_true(S7::S7_inherits(ex_model, Exemplar_Model))

  expect_error(
    new_mvg_category_representation(
      category_labels = "A",
      cue_labels = c("F1", "F2"),
      mu = c(0, 1),
      Sigma = matrix(c(1, 2, 0, 1), nrow = 2)
    ),
    "must be symmetric"
  )

  expect_error(
    new_uvg_category_representation(
      category_labels = "A",
      cue_labels = c("F1", "F2"),
      mu = 0,
      sigma2 = 1
    ),
    "single cue dimension"
  )

  expect_error(
    new_exemplar_category_representation(
      category_labels = "A",
      cue_labels = c("F1", "F2"),
      exemplars = matrix(c(0, 1, 1, 2), nrow = 2, byrow = TRUE),
      exemplar_weights = c(0.8, 0.3)
    ),
    "entries must sum to 1"
  )

  expect_error(
    new_mvg_ideal_observer(
      category_template = new_category_representation_template(representations = list(A = niw_rep))
    ),
    "must inherit from MVG category representation class"
  )

  uvg_d <- uvg_rep@category_likelihood_function(c(0, 1))
  uvg_ld <- uvg_rep@category_likelihood_function(c(0, 1), log = TRUE)
  expect_length(uvg_d, 2)
  expect_true(all(is.finite(uvg_d)))
  expect_equal(uvg_ld, log(uvg_d), tolerance = 1e-10)

  nix_d <- nix_rep@category_likelihood_function(c(0, 1))
  nix_ld <- nix_rep@category_likelihood_function(c(0, 1), log = TRUE)
  expect_length(nix_d, 2)
  expect_true(all(is.finite(nix_d)))
  expect_equal(nix_ld, log(nix_d), tolerance = 1e-10)

  niw_d <- niw_rep@category_likelihood_function(matrix(c(0, 1, 1, 2), ncol = 2, byrow = TRUE))
  niw_ld <- niw_rep@category_likelihood_function(matrix(c(0, 1, 1, 2), ncol = 2, byrow = TRUE), log = TRUE)
  expect_length(niw_d, 2)
  expect_true(all(is.finite(niw_d)))
  expect_equal(niw_ld, log(niw_d), tolerance = 1e-10)

  mvg_d <- mvg_rep@category_likelihood_function(matrix(c(0, 1, 1, 2), ncol = 2, byrow = TRUE))
  mvg_ld <- mvg_rep@category_likelihood_function(matrix(c(0, 1, 1, 2), ncol = 2, byrow = TRUE), log = TRUE)
  expect_length(mvg_d, 2)
  expect_true(all(is.finite(mvg_d)))
  expect_equal(mvg_ld, log(mvg_d), tolerance = 1e-10)

  ex_d <- ex_rep@category_likelihood_function(matrix(c(0, 1, 1, 2), ncol = 2, byrow = TRUE))
  ex_ld <- ex_rep@category_likelihood_function(matrix(c(0, 1, 1, 2), ncol = 2, byrow = TRUE), log = TRUE)
  expect_length(ex_d, 2)
  expect_true(all(is.finite(ex_d)))
  expect_equal(ex_ld, log(ex_d), tolerance = 1e-10)

  mnix_d <- mnix_rep_lik@category_likelihood_function(c(0, 1))
  mnix_ld <- mnix_rep_lik@category_likelihood_function(c(0, 1), log = TRUE)
  expect_length(mnix_d, 2)
  expect_true(all(is.finite(mnix_d)))
  expect_equal(mnix_ld, log(mnix_d), tolerance = 1e-10)

  muvg_rep_lik <- new_muvg_category_representation(
    category_labels = "C",
    cue_labels = c("F1", "F2"),
    component_mu = c(0, 1),
    component_sigma2 = c(1, 2)
  )
  muvg_d <- muvg_rep_lik@category_likelihood_function(matrix(c(0, 1, 1, 2), ncol = 2, byrow = TRUE))
  muvg_ld <- muvg_rep_lik@category_likelihood_function(matrix(c(0, 1, 1, 2), ncol = 2, byrow = TRUE), log = TRUE)
  expect_length(muvg_d, 2)
  expect_true(all(is.finite(muvg_d)))
  expect_equal(muvg_ld, log(muvg_d), tolerance = 1e-10)
})

test_that("unified S7 categorization and prediction methods support single and batch inputs", {
  mvg_rep_a <- new_mvg_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    mu = c(0, 0),
    Sigma = diag(c(1, 1))
  )
  mvg_rep_b <- new_mvg_category_representation(
    category_labels = "B",
    cue_labels = c("F1", "F2"),
    mu = c(3, 3),
    Sigma = diag(c(1, 1))
  )

  model <- new_mvg_ideal_observer(
    category_template = new_category_representation_template(representations = list(A = mvg_rep_a, B = mvg_rep_b)),
    category_prior = c(A = 0.5, B = 0.5)
  )

  x_single <- matrix(c(0, 0, 3, 3), ncol = 2, byrow = TRUE)
  post_single <- posterior(model, x_single, categories = NULL)
  expect_true(is.matrix(post_single))
  expect_equal(dim(post_single), c(2, 2))
  expect_equal(rowSums(post_single), c(1, 1), tolerance = 1e-10)
  expect_equal(colnames(post_single), c("A", "B"))

  pred_single <- categorize(model, x_single, decision_rule = NULL)
  expect_true(is.data.frame(pred_single))
  expect_equal(names(pred_single), c("category", "probability"))
  expect_equal(nrow(pred_single), 2)
  expect_true(all(pred_single$probability >= 0 & pred_single$probability <= 1))

  x_batch <- list(
    matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE),
    matrix(c(3, 3, 2, 2), ncol = 2, byrow = TRUE)
  )

  post_batch <- posterior(model, x_batch, categories = NULL)
  pred_batch <- categorize(model, x_batch, decision_rule = NULL)
  posterior_batch <- list(category_posterior = post_batch, category = pred_batch)

  expect_true(is.list(post_batch))
  expect_true(is.list(pred_batch))
  expect_true(is.list(posterior_batch))
  expect_length(post_batch, 2)
  expect_length(pred_batch, 2)
  expect_true(all(vapply(post_batch, is.matrix, logical(1))))
  expect_true(all(vapply(pred_batch, is.data.frame, logical(1))))
  expect_true(all(vapply(posterior_batch, is.list, logical(1))))
})

test_that("family-typed model distributions and registry hooks work", {
  niw_dist <- new_niw_model_distribution(group_label = "fit-niw")
  mvg_dist <- new_mvg_model_distribution(group_label = "fit-mvg")
  muvg_dist <- new_muvg_model_distribution(group_label = "fit-muvg")
  mnix_dist <- new_mnix_model_distribution(group_label = "fit-mnix")
  uvg_dist <- new_uvg_model_distribution(group_label = "fit-uvg")
  nix_dist <- new_nix_model_distribution(group_label = "fit-nix")
  ex_dist <- new_exemplar_model_distribution(group_label = "fit-ex")

  expect_true(S7::S7_inherits(niw_dist, NIW_IdealAdaptorDistribution))
  expect_true(S7::S7_inherits(mvg_dist, MVG_IdealObserverDistribution))
  expect_true(S7::S7_inherits(muvg_dist, MUVG_IdealObserverDistribution))
  expect_true(S7::S7_inherits(mnix_dist, MNIX_IdealAdaptorDistribution))
  expect_true(S7::S7_inherits(uvg_dist, UVG_IdealObserverDistribution))
  expect_true(S7::S7_inherits(nix_dist, NIX_IdealAdaptorDistribution))
  expect_true(S7::S7_inherits(ex_dist, Exemplar_ModelDistribution))

  expect_equal(get_model_family(niw_dist), "NIW")
  expect_equal(get_model_family(mvg_dist), "MVG")
  expect_equal(get_model_family(muvg_dist), "MUVG")
  expect_equal(get_model_family(mnix_dist), "MNIX")
  expect_equal(get_model_family(uvg_dist), "UVG")
  expect_equal(get_model_family(nix_dist), "NIX")
  expect_equal(get_model_family(ex_dist), "EXEMPLAR")

  baseline_families <- get_registered_model_families()
  expect_true(all(c("EXEMPLAR", "MNIX", "MUVG", "MVG", "NIW", "NIX", "UVG") %in% baseline_families))

  register_model_family(
    family = "custommix",
    category_representation_class = "CustomMix_CategoryRepresentation",
    cognitive_model_class = "CustomMix_Model",
    model_distribution_class = "CustomMix_ModelDistribution"
  )

  expect_true("CUSTOMMIX" %in% get_registered_model_families())
  custom_registration <- get_model_family_registration("CUSTOMMIX")
  expect_equal(custom_registration$category_representation, "CustomMix_CategoryRepresentation")
  expect_equal(custom_registration$cognitive_model, "CustomMix_Model")
  expect_equal(custom_registration$model_distribution, "CustomMix_ModelDistribution")
})

test_that("Stan-family extension hooks can be registered and retrieved", {
  niw_hooks <- get_stan_family_hooks("NIW")
  expect_true("get_stanfit" %in% niw_hooks$bridge_methods)
  expect_true("posterior::as_draws_df" %in% niw_hooks$bridge_methods)

  register_stan_family_hooks(
    family = "uvg",
    stanfit_class = "UVG_Stanfit",
    bridge_methods = c("get_stanfit", "summary"),
    dependency_rationale = "Minimal bridge for prototype workflows."
  )

  uvg_hooks <- get_stan_family_hooks("UVG")
  expect_equal(uvg_hooks$stanfit_class, "UVG_Stanfit")
  expect_equal(uvg_hooks$bridge_methods, c("get_stanfit", "summary"))
})
