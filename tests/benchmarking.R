source("functions-to-make-or-load-models.R")

library(microbenchmark)

example <- 1
seed <- 42
.data <- make_data_for_stanfit(example, seed = seed)

make_staninput_for_NIW_ideal_adaptor_with_transform <-
  function(transform_type)
    make_staninput_for_NIW_ideal_adaptor(
      exposure = .data %>% filter(Phase == "exposure"),
      test = .data %>% filter(Phase == "test"),
      cues =
        if (example %in% 1:3)
        {
          c("VOT", "f0_semitones", "vowel_duration")[1:(example)]
        } else if (example %in% 4:6) {
          c("cue1", "cue2", "cue3")[1:(example %% 4 + 1)]
        },
      category = "category",
      response = "Response",
      group = "Subject",
      group.unique = "Condition",
      transform_type = transform_type)

mb <-
  microbenchmark(
    # Cholesky with transformation and back-transformation
    "Cholesky, center" =
      fit_NIW_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "center"),
        stanmodel = "mvg_conj_sufficient_stats_lapse_cholesky",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.85)),
    "Cholesky, standardize" =
      fit_NIW_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "standardize"),
        stanmodel = "mvg_conj_sufficient_stats_lapse_cholesky",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.85)),
    "Cholesky, PCA whiten" =
      fit_NIW_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "PCA whiten"),
        stanmodel = "mvg_conj_sufficient_stats_lapse_cholesky",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.85)),
    # No Cholesky with transformation and back-transformation
    "no cholesky, no transform" =
      fit_NIW_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "identity"),
        stanmodel = "mvg_conj_sufficient_stats_lapse",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.85)),
    "no cholesky, center" =
      fit_NIW_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "center"),
        stanmodel = "mvg_conj_sufficient_stats_lapse",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.85)),
    "no cholesky, standardize" =
      fit_NIW_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "standardize"),
        stanmodel = "mvg_conj_sufficient_stats_lapse",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.85)),
    "no cholesky, PCA whiten" =
      fit_NIW_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "PCA whiten"),
        stanmodel = "mvg_conj_sufficient_stats_lapse",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.85)),
    "no cholesky, ZCA whiten" =
      fit_NIW_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "ZCA whiten"),
        stanmodel = "mvg_conj_sufficient_stats_lapse",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.85)),
    times = 25,
    control = list(order = "random"))

plot(mb)
autoplot(mb)
print(mb)
