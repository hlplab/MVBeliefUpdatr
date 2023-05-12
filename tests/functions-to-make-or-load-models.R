get_example_stanfit <- function() {
  filename <- "../example-stanfit.rds"
  if (file.exists(filename)) {
    fit <- readRDS(filename)
  } else {
    require(tidyverse)
    require(MVBeliefUpdatr)

    n_subject <- 30
    # number of trials in exposure per category per subject
    # (and there will be 2 * n_trials trials in test per subject)
    n_trial <- 50
    .cues <- c("cue1", "cue2")

    # Make 5 ideal observers
    .io <-
      example_MVG_ideal_observer(1) %>%
      mutate(Sigma = map(Sigma, ~ .x * 5))
    .io.2.2 <-
      .io %>%
      mutate(mu = map(mu, ~ .x + c(2, 2)))
    .io.4.4 <-
      .io %>%
      mutate(mu = map(mu, ~ .x + c(4, 4)))
    .io.2.4 <-
      .io %>%
      mutate(mu = map(mu, ~ .x + c(2, 4)))
    .io.4.2 <-
      .io %>%
      mutate(mu = map(mu, ~ .x + c(4, 2)))

    # Sample responses for subjects that have converged against those five states
    .data <-
      bind_rows(
        sample_MVG_data_from_model(
          model = .io,
          Ns = n_trial,
          keep.input_parameters = F) %>%
          make_vector_column(cols = c("cue1", "cue2"), vector_col = "cue") %>%
          crossing(Subject = 1:n_subject) %>%
          mutate(Condition = "baseline",
                 Response = get_categorization_from_MVG_ideal_observer(
                   x = cue,
                   model = .io,
                   decision_rule = "sampling",
                   simplify = T,
                   noise_treatment = "no_noise",
                   lapse_treatment = "no_lapses")),
        sample_MVG_data_from_model(
          model = .io.2.2,
          Ns = n_trial,
          keep.input_parameters = F) %>%
          make_vector_column(cols = c("cue1", "cue2"), vector_col = "cue") %>%
          crossing(Subject = 1:n_subject) %>%
          mutate(Condition = "plus2.2",
                 Response = get_categorization_from_MVG_ideal_observer(
                   x = cue,
                   model = .io.2.2,
                   decision_rule = "sampling",
                   simplify = T,
                   noise_treatment = "no_noise",
                   lapse_treatment = "no_lapses")),
        sample_MVG_data_from_model(
          model = .io.4.4,
          Ns = n_trial,
          keep.input_parameters = F) %>%
          make_vector_column(cols = c("cue1", "cue2"), vector_col = "cue") %>%
          crossing(Subject = 1:n_subject) %>%
          mutate(Condition = "plus4.4",
                 Response = get_categorization_from_MVG_ideal_observer(
                   x = cue,
                   model = .io.4.4,
                   decision_rule = "sampling",
                   simplify = T,
                   noise_treatment = "no_noise",
                   lapse_treatment = "no_lapses")),
        sample_MVG_data_from_model(
          model = .io.2.4,
          Ns = n_trial,
          keep.input_parameters = F) %>%
          make_vector_column(cols = c("cue1", "cue2"), vector_col = "cue") %>%
          crossing(Subject = 1:n_subject) %>%
          mutate(Condition = "plus2.4",
                 Response = get_categorization_from_MVG_ideal_observer(
                   x = cue,
                   model = .io.2.4,
                   decision_rule = "sampling",
                   simplify = T,
                   noise_treatment = "no_noise",
                   lapse_treatment = "no_lapses")),
        sample_MVG_data_from_model(
          model = .io.4.2,
          Ns = n_trial,
          keep.input_parameters = F) %>%
          make_vector_column(cols = c("cue1", "cue2"), vector_col = "cue") %>%
          crossing(Subject = 1:n_subject) %>%
          mutate(Condition = "plus4.2",
                 Response = get_categorization_from_MVG_ideal_observer(
                   x = cue,
                   model = .io.4.2,
                   decision_rule = "sampling",
                   simplify = T,
                   noise_treatment = "no_noise",
                   lapse_treatment = "no_lapses")))

    fit <- infer_prior_beliefs(
      exposure = .data,
      test = .data,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "Response",
      group = "Subject",
      group.unique = "Condition",
      file = filename,
      cores = 4)
  }

  return(fit)
}
