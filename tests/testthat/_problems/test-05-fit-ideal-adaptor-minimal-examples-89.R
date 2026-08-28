# Extracted from test-05-fit-ideal-adaptor-minimal-examples.R:89

# test -------------------------------------------------------------------------
skip_if_not_installed("rstan")
model_name <- "minimal-nix_ideal_adaptor-nix"
fit <- load_ideal_adaptor_fit_model(model_name)
tmp_file <- tempfile(fileext = ".rds")
saveRDS(fit, tmp_file)
data <- build_shifted_prior_fit_example(
    cues = "VOT",
    n_exposure_per_condition = 20L,
    n_test_per_condition = 10L,
    seed = 321L
  )
input <- new_ideal_adaptor_stanfit_input(
    exposure = data$exposure,
    test = data$test,
    cues = "VOT",
    category = "category",
    response = "response",
    group = "group",
    group.unique = "Condition",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "NIX_ideal_adaptor"
  )
reloaded_fit <- fit_ideal_adaptor(
    stanfit_input = input,
    file = tmp_file,
    file_refit = "never",
    chains = 0,
    iter = 0,
    warmup = 0,
    refresh = 0,
    stanmodel = "NIX_ideal_adaptor"
  )
