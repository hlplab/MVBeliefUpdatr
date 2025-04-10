source("functions-to-make-or-load-models.R")

library(microbenchmark)

example <- 1
seed <- 42
.data <- make_data_for_stanfit(example, seed = seed)

mb <- microbenchmark(
  # Cholesky with transformation but without back-transformation
  "cholesky, no transform" = infer_prior_beliefs(
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
    center.observations = F, scale.observations = F, pca.observations = F, transform_type = NULL, stanmodel = "mvg_conj_sufficient_stats_lapse_cholesky",
    file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.85)),
  "cholesky, center" = infer_prior_beliefs(
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
    center.observations = T, scale.observations = F, pca.observations = F, transform_type = NULL, stanmodel = "mvg_conj_sufficient_stats_lapse_cholesky",
    file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.85)),
  "cholesky, standardize" = infer_prior_beliefs(
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
    center.observations = T, scale.observations = T, pca.observations = F, transform_type = NULL, stanmodel = "mvg_conj_sufficient_stats_lapse_cholesky",
    file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.85)),
  # No Cholesky with transformation and back-transformation
  "no cholesky, no transform" = infer_prior_beliefs(
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
    center.observations = NULL, scale.observations = NULL, pca.observations = NULL, transform_type = "identity", stanmodel = "mvg_conj_sufficient_stats_lapse",
    file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.85)),
  "no cholesky, center" = infer_prior_beliefs(
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
    center.observations = NULL, scale.observations = NULL, pca.observations = NULL, transform_type = "center", stanmodel = "mvg_conj_sufficient_stats_lapse",
    file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.85)),
  "no cholesky, standardize" = infer_prior_beliefs(
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
    center.observations = NULL, scale.observations = NULL, pca.observations = NULL, transform_type = "standardize", stanmodel = "mvg_conj_sufficient_stats_lapse",
    file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.85)),
  "no cholesky, PCA whiten" = infer_prior_beliefs(
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
    center.observations = NULL, scale.observations = NULL, pca.observations = NULL, transform_type = "PCA whiten", stanmodel = "mvg_conj_sufficient_stats_lapse",
    file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.85)),
  "no cholesky, ZCA whiten" = infer_prior_beliefs(
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
    center.observations = NULL, scale.observations = NULL, pca.observations = NULL, transform_type = "ZCA whiten", stanmodel = "mvg_conj_sufficient_stats_lapse",
    file = NULL, refresh = 400, iter = 2000, cores = 4, control = list(adapt_delta = 0.85)),
  times = 20,
  control = list(order = "random")
)

plot(mb)
autoplot(mb)
print(mb)
