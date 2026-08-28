# Extracted from test-04-stanfit-input-compatibility-with-stanprogram.R:26

# test -------------------------------------------------------------------------
for (n_obs_exposure in c(0L, 1L, 3L)) {
    data <- make_minimal_staninput_data(
      cues = "cue1",
      n_obs_exposure = n_obs_exposure,
      n_obs_test = 2L,
      n_group = 1L,
      n_category = 2L
    )

    input <- new_ideal_adaptor_stanfit_input(
      exposure = data$exposure,
      test = data$test,
      cues = "cue1",
      category = "category",
      response = "response",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIX_ideal_adaptor"
    )

    run_fixed_param_stan_program(
      stan_file = system.file("stan", "NIX_ideal_adaptor.stan", package = "MVBeliefUpdatr"),
      input = input,
      model_name = "NIX_compat"
    )
  }
