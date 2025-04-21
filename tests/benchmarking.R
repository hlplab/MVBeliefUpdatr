source("functions-to-make-or-load-models.R")

library(microbenchmark)

example <- 1
seed <- 42
.data <- make_data_for_stanfit(example, seed = seed)

make_staninput_for_NIW_ideal_adaptor_with_transform <-
  function(transform_type, example = 1)
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

mb1 <-
  microbenchmark(
    # Cholesky with transformation and back-transformation
    "Cholesky, center" =
      fit_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "center"),
        stanmodel = "NIW_ideal_adaptor_cholesky",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.92)),
    "Cholesky, standardize" =
      fit_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "standardize"),
        stanmodel = "NIW_ideal_adaptor_cholesky",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.92)),
    "Cholesky, PCA whiten" =
      fit_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "PCA whiten"),
        stanmodel = "NIW_ideal_adaptor_cholesky",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.92)),
    # No Cholesky with transformation and back-transformation
    "no cholesky, no transform" =
      fit_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "identity"),
        stanmodel = "NIW_ideal_adaptor",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.92)),
    "no cholesky, center" =
      fit_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "center"),
        stanmodel = "NIW_ideal_adaptor",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.92)),
    "no cholesky, standardize" =
      fit_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "standardize"),
        stanmodel = "NIW_ideal_adaptor",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.92)),
    "no cholesky, PCA whiten" =
      fit_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "PCA whiten"),
        stanmodel = "NIW_ideal_adaptor",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.92)),
    "no cholesky, ZCA whiten" =
      fit_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "ZCA whiten"),
        stanmodel = "NIW_ideal_adaptor",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.92)),
    times = 25,
    control = list(order = "random"))

plot(mb1)
autoplot(mb1)
print(mb1)

saveRDS(mb1, "benchmark-1cue.rds", compress = T)


example <- 3
mb3 <-
  microbenchmark(
    # Cholesky with transformation and back-transformation
    "Cholesky, center" =
      fit_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "center", example = 3),
        stanmodel = "NIW_ideal_adaptor_cholesky",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.999)),
    "Cholesky, standardize" =
      fit_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "standardize", example = 3),
        stanmodel = "NIW_ideal_adaptor_cholesky",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.999)),
    # No Cholesky with transformation and back-transformation
    "no cholesky, standardize" =
      fit_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "standardize", example = 3),
        stanmodel = "NIW_ideal_adaptor",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.999)),
    "no cholesky, PCA whiten" =
      fit_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "PCA whiten", example = 3),
        stanmodel = "NIW_ideal_adaptor",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.999)),
    "no cholesky, ZCA whiten" =
      fit_ideal_adaptor(
        staninput = make_staninput_for_NIW_ideal_adaptor_with_transform(transform_type = "ZCA whiten", example = 3),
        stanmodel = "NIW_ideal_adaptor",
        file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.999)),
    times = 25,
    control = list(order = "random"))

plot(mb3)
autoplot(mb3)
print(mb3)

saveRDS(mb3, "benchmark-3cues.rds", compress = T)
