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
})

test_that("validation helpers behave as expected", {
  rep_obj <- new_category_representation(
    category_labels = c("A"),
    cue_labels = c("F1")
  )

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

  expect_invisible(summary(rep_obj))
  expect_invisible(print(rep_obj))

  expect_error(
    posterior(model, data.frame(x = 1), categories = NULL),
    "category_likelihood not implemented|not yet implemented|Can't find method"
  )
  expect_error(
    categorize(model, data.frame(x = 1), decision_rule = NULL),
    "category_likelihood not implemented|not yet implemented|Can't find method"
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
  expect_identical(get_metadata(rep), rep@metadata)
  expect_identical(get_category_representations(template), template@representations)
  expect_equal(get_category_template(model), model@category_template)
  expect_identical(get_category_template(model), model@category_template)
  expect_identical(get_category_representations(model), get_category_representations(get_category_template(model)))
  expect_equal(get_category_prior(model), c(A = 1))
  expect_equal(get_lapse_rate(model), 0.1)
  expect_equal(get_lapse_bias(model), c(A = 1))
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
    mu = c(0, 1),
    sigma2 = c(1, 1)
  )

  expect_true(S7::S7_inherits(muvg_rep, MUVG_CategoryRepresentation))
  params <- get_parameters(muvg_rep)
  expect_equal(params$weights, c(0.5, 0.5), tolerance = MVBU_PROB_TOL)

  expect_error(
    new_muvg_category_representation(
      category_labels = "A",
      cue_labels = c("F1", "F2"),
      mu = c(0, 1),
      sigma2 = c(1, 1),
      weights = c(0.9, 0.2)
    ),
    "weights entries must sum to 1"
  )

  expect_error(
    new_muvg_category_representation(
      category_labels = "A",
      cue_labels = "F1",
      mu = c(0, 1),
      sigma2 = c(1, 1)
    ),
    "cue_labels length must match"
  )
})

test_that("MNIX representation schema validators and typed model constructor work", {
  mnix_rep <- new_mnix_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    m = c(0, 1),
    kappa = c(1, 2),
    nu = c(3, 4),
    sigma2 = c(1, 2)
  )

  expect_true(S7::S7_inherits(mnix_rep, MNIX_CategoryRepresentation))
  mnix_params <- get_parameters(mnix_rep)
  expect_equal(mnix_params$weights, c(0.5, 0.5), tolerance = MVBU_PROB_TOL)

  mnix_rep_set <- new_category_representation_template(representations = list(A = mnix_rep))
  mnix_model <- new_mnix_ideal_adaptor(category_template = mnix_rep_set)
  expect_true(S7::S7_inherits(mnix_model, MNIX_IdealAdaptor))
  expect_equal(get_category_prior(mnix_model), c(A = 1), tolerance = MVBU_PROB_TOL)
  muvg_rep <- new_muvg_category_representation(
    category_labels = "B",
    cue_labels = c("F1", "F2"),
    mu = c(0, 1),
    sigma2 = c(1, 1)
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
      m = c(0, 1),
      kappa = c(1),
      nu = c(3, 4),
      sigma2 = c(1, 2),
      weights = c(0.5, 0.5)
    ),
    "MNIX parameter vectors must all have equal length"
  )

  mnix_rep_b <- new_mnix_category_representation(
    category_labels = "B",
    cue_labels = c("F1", "F2"),
    m = c(2, 3),
    kappa = c(1, 2),
    nu = c(3, 4),
    sigma2 = c(1, 1)
  )
  mnix_model_bi <- new_mnix_ideal_adaptor(
    category_template = new_category_representation_template(
      representations = list(A = mnix_rep, B = mnix_rep_b)
    )
  )
  expect_true(S7::S7_inherits(mnix_model_bi, MNIX_IdealAdaptor))
})

test_that("new_mnix_*_from_data supports scalar, vector, and dots fallback", {
  df <- data.frame(
    category = rep(c("A", "B"), each = 20),
    F1 = rnorm(40),
    F2 = rnorm(40)
  )
  # Scalar kappa and nu
  rep_scalar <- new_mnix_category_representation_from_data(
    df[df$category == "A", ],
    category = "category",
    cues = c("F1", "F2"),
    kappa = 5,
    nu = 4
  )
  expect_equal(rep_scalar@kappa, c(5, 5))
  expect_equal(rep_scalar@nu, c(4, 4))

  # Vector kappa and nu
  rep_vector <- new_mnix_category_representation_from_data(
    df[df$category == "A", ],
    category = "category",
    cues = c("F1", "F2"),
    kappa = c(10, 2),
    nu = c(6, 4)
  )
  expect_equal(rep_vector@kappa, c(10, 2))
  expect_equal(rep_vector@nu, c(6, 4))

  # Generic factory fallback via kappa / nu in dots
  mod_dots <- new_model_from_data(
    data = df,
    type = "MNIX",
    category = "category",
    cues = c("F1", "F2"),
    kappa = 8,
    nu = 5
  )
  expect_true(S7::S7_inherits(mod_dots, MNIX_IdealAdaptor))
  rep_a <- mod_dots@category_template@representations$A
  expect_equal(rep_a@kappa, c(8, 8))
  expect_equal(rep_a@nu, c(5, 5))
})

test_that("MUVG typed ideal observer constructor enforces representation class", {
  muvg_rep_a <- new_muvg_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    mu = c(0, 1),
    sigma2 = c(1, 1)
  )
  muvg_rep_b <- new_muvg_category_representation(
    category_labels = "B",
    cue_labels = c("F1", "F2"),
    mu = c(2, 3),
    sigma2 = c(1, 1)
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
    cue_labels = c("F1", "F2"),
    m = c(0, 1),
    kappa = c(1, 2),
    nu = c(3, 4),
    sigma2 = c(1, 2)
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

  # Test get_parameter_names for representation objects
  expect_equal(get_parameter_names(uvg_rep), c("mu", "sigma2"))
  expect_equal(get_parameter_names(nix_rep), c("m", "kappa", "nu", "sigma2"))
  expect_equal(get_parameter_names(mvg_rep), c("mu", "Sigma"))
  expect_equal(get_parameter_names(niw_rep), c("m", "kappa", "nu", "S"))
  expect_equal(
    get_parameter_names(ex_rep),
    c("exemplars", "exemplar_weights", "c")
  )
  expect_true(ex_rep@c > 0)

  # Test custom c vs default Silverman c
  ex_custom <- new_exemplar_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    exemplars = matrix(c(0, 1, 1, 2), nrow = 2, byrow = TRUE),
    c = 2.5
  )
  expect_equal(ex_custom@c, 2.5)

  # Test get_parameter_names for cognitive model objects
  expect_equal(
    get_parameter_names(uvg_model),
    c("mu", "sigma2", "category_prior", "lapse_rate", "lapse_bias", "Sigma_noise")
  )
  expect_equal(
    get_parameter_names(nix_model),
    c("m", "kappa", "nu", "sigma2", "category_prior", "lapse_rate", "lapse_bias", "Sigma_noise")
  )
  expect_equal(
    get_parameter_names(mvg_model),
    c("mu", "Sigma", "category_prior", "lapse_rate", "lapse_bias", "Sigma_noise")
  )
  expect_equal(
    get_parameter_names(niw_model),
    c("m", "kappa", "nu", "S", "category_prior", "lapse_rate", "lapse_bias", "Sigma_noise")
  )
  expect_equal(
    get_parameter_names(ex_model),
    c("exemplars", "exemplar_weights", "c", "category_prior", "lapse_rate", "lapse_bias", "Sigma_noise")
  )

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
    new_exemplar_category_representation(
      category_labels = "A",
      cue_labels = c("F1", "F2"),
      exemplars = matrix(c(0, 1, 1, 2), nrow = 2, byrow = TRUE),
      c = -1
    ),
    "c must be a single positive numeric scalar"
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

  mnix_d <- mnix_rep_lik@category_likelihood_function(matrix(c(0, 1, 1, 2), ncol = 2, byrow = TRUE))
  mnix_ld <- mnix_rep_lik@category_likelihood_function(matrix(c(0, 1, 1, 2), ncol = 2, byrow = TRUE), log = TRUE)
  expect_length(mnix_d, 2)
  expect_true(all(is.finite(mnix_d)))
  expect_equal(mnix_ld, log(mnix_d), tolerance = 1e-10)

  muvg_rep_lik <- new_muvg_category_representation(
    category_labels = "C",
    cue_labels = c("F1", "F2"),
    mu = c(0, 1),
    sigma2 = c(1, 2)
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
  expect_type(pred_single, "character")
  expect_length(pred_single, 2)
  expect_equal(pred_single, c("A", "B"))

  pred_single_df <- categorize(model, x_single, decision_rule = NULL, simplify = FALSE)
  expect_true(is.data.frame(pred_single_df))
  expect_equal(names(pred_single_df), c("response_category", "response_probability"))
  expect_equal(nrow(pred_single_df), 2)
  expect_equal(pred_single_df$response_probability, c(1, 1))

  pred_prop <- categorize(model, x_single, decision_rule = "proportional")
  expect_true(is.data.frame(pred_prop))
  expect_equal(names(pred_prop), c("response_category", "response_probability"))
  expect_true(all(pred_prop$response_probability >= 0 & pred_prop$response_probability <= 1))

  x_batch <- list(
    matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE),
    matrix(c(3, 3, 2, 2), ncol = 2, byrow = TRUE)
  )

  post_batch <- posterior(model, x_batch, categories = NULL)
  pred_batch <- categorize(model, x_batch, decision_rule = NULL, simplify = FALSE)
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

test_that("model family registry and hooks work", {
  baseline_families <- list_model_families()
  expect_true(
    all(
      c("EXEMPLAR", "MNIX", "MUVG", "MVG", "NIW", "NIX", "UVG") %in%
        baseline_families
    )
  )

  register_model_family(
    family = "custommix",
    category_representation_class = "CustomMix_CategoryRepresentation",
    cognitive_model_class = "CustomMix_Model"
  )

  expect_true("CUSTOMMIX" %in% list_model_families())
  custom_registration <- get_model_family_registration("CUSTOMMIX")
  expect_equal(
    custom_registration$category_representation,
    "CustomMix_CategoryRepresentation"
  )
  expect_equal(custom_registration$cognitive_model, "CustomMix_Model")
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

test_that("get_model_type and get_representation_type work across core objects", {
  rep_uvg <- example_category_representation("UVG")
  rep_mvg <- example_category_representation("MVG")
  rep_niw <- example_category_representation("NIW")
  rep_nix <- example_category_representation("NIX")
  rep_ex <- example_category_representation("EXEMPLAR")

  expect_equal(get_representation_type(rep_uvg), "UVG")
  expect_equal(get_representation_type(rep_mvg), "MVG")
  expect_equal(get_representation_type(rep_niw), "NIW")
  expect_equal(get_representation_type(rep_nix), "NIX")
  expect_equal(get_representation_type(rep_ex), "EXEMPLAR")

  tpl_mvg <- example_category_representation_template("MVG")
  tpl_niw <- example_category_representation_template("NIW")
  expect_equal(get_representation_type(tpl_mvg), "MVG")
  expect_equal(get_representation_type(tpl_niw), "NIW")

  mod_mvg <- example_mvg_ideal_observer()
  mod_niw <- example_niw_ideal_adaptor()
  mod_ex <- example_exemplar_model()
  expect_equal(get_model_type(mod_mvg), "MVG")
  expect_equal(get_model_type(mod_niw), "NIW")
  expect_equal(get_model_type(mod_ex), "EXEMPLAR")
})


test_that("add_category_representation generic works on templates, cognitive models, and model lists", {
  repA <- new_niw_category_representation(
    category_labels = "A",
    cue_labels = c("c1", "c2"),
    m = c(0, 0),
    S = diag(2, 2L),
    kappa = 2,
    nu = 6
  )
  repB <- new_niw_category_representation(
    category_labels = "B",
    cue_labels = c("c1", "c2"),
    m = c(2, 2),
    S = diag(2, 2L),
    kappa = 2,
    nu = 6
  )
  repC <- new_niw_category_representation(
    category_labels = "C",
    cue_labels = c("c1", "c2"),
    m = c(4, 4),
    S = diag(2, 2L),
    kappa = 2,
    nu = 6
  )

  tpl <- new_category_representation_template(list(repA, repB))
  expect_equal(get_category_labels(tpl), c("A", "B"))

  # Add to template
  tpl_with_c <- add_category_representation(tpl, repC)
  expect_equal(get_category_labels(tpl_with_c), c("A", "B", "C"))

  # Add to cognitive model with default prior and bias (default 0)
  mod <- new_niw_ideal_adaptor(
    category_template = tpl,
    category_prior = c(A = 0.6, B = 0.4),
    lapse_bias = c(A = 0.7, B = 0.3)
  )
  expect_equal(get_category_labels(mod), c("A", "B"))
  mod_with_c <- add_category_representation(mod, repC)
  expect_equal(get_category_labels(mod_with_c), c("A", "B", "C"))
  expect_equal(
    get_category_prior(mod_with_c),
    c(A = 0.6, B = 0.4, C = 0),
    tolerance = MVBU_PROB_TOL
  )
  expect_equal(
    get_lapse_bias(mod_with_c),
    c(A = 0.7, B = 0.3, C = 0),
    tolerance = MVBU_PROB_TOL
  )
  expect_true(S7::S7_inherits(mod_with_c, NIW_IdealAdaptor))

  # Add with scalar category_prior and lapse_bias (rescaling old proportionally)
  mod_scalar <- add_category_representation(
    mod,
    repC,
    category_prior = 0.2,
    lapse_bias = 0.5
  )
  expect_equal(
    get_category_prior(mod_scalar),
    c(A = 0.6 * 0.8, B = 0.4 * 0.8, C = 0.2),
    tolerance = MVBU_PROB_TOL
  )
  expect_equal(
    get_lapse_bias(mod_scalar),
    c(A = 0.7 * 0.5, B = 0.3 * 0.5, C = 0.5),
    tolerance = MVBU_PROB_TOL
  )

  # Add with full vector for all categories
  mod_vector <- add_category_representation(
    mod,
    repC,
    category_prior = c(A = 0.2, B = 0.5, C = 0.3),
    lapse_bias = c(A = 0.1, B = 0.1, C = 0.8)
  )
  expect_equal(
    get_category_prior(mod_vector),
    c(A = 0.2, B = 0.5, C = 0.3),
    tolerance = MVBU_PROB_TOL
  )
  expect_equal(
    get_lapse_bias(mod_vector),
    c(A = 0.1, B = 0.1, C = 0.8),
    tolerance = MVBU_PROB_TOL
  )

  # Validation errors for invalid prior / bias
  expect_error(
    add_category_representation(mod, repC, category_prior = 1.5),
    "category_prior for the new category must be between 0 and 1"
  )
  expect_error(
    add_category_representation(mod, repC, category_prior = c(0.1, 0.2)),
    "category_prior must be a scalar.*or a numeric vector of length 3"
  )
  expect_error(
    add_category_representation(
      mod,
      repC,
      category_prior = c(A = 0.2, B = 0.2, C = 0.2)
    ),
    "category_prior entries must sum to 1"
  )

  # Add to model list
  mlist <- as_model_list(list(M1 = mod, M2 = mod))
  mlist_with_c <- add_category_representation(
    mlist,
    repC,
    category_prior = 0.25
  )
  expect_equal(get_category_labels(mlist_with_c[[1]]), c("A", "B", "C"))
  expect_equal(get_category_labels(mlist_with_c[[2]]), c("A", "B", "C"))
  expect_equal(
    get_category_prior(mlist_with_c[[1]]),
    c(A = 0.6 * 0.75, B = 0.4 * 0.75, C = 0.25),
    tolerance = MVBU_PROB_TOL
  )
})
