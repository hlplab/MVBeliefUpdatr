make_data_for_stanfit <- function(example = 1, seed = NULL) {
  require(tidyverse)
  require(magrittr)
  require(MVBeliefUpdatr)

  if (!is.null(seed)) set.seed(seed)

  if (example == 1) {
    return(make_data_for_1Dstanfit_with_exposure())
  } else if (example == 2) {
    return(make_data_for_2Dstanfit_with_exposure())
  } else if (example == 3) {
    return(make_data_for_3Dstanfit_with_exposure())
  } else if (example == 4) {
    return(make_data_for_1Dstanfit_without_exposure())
  } else if (example == 5) {
    return(make_data_for_2Dstanfit_without_exposure())
  } else if (example == 6) {
    return(make_data_for_3Dstanfit_without_exposure())
  }
}

make_data_for_1Dstanfit_with_exposure <- function() {
  n_subject <- 30
  # number of trials in exposure per category per subject
  # (and there will be 2 * n_trials trials in test per subject)
  n_trial <- 50
  .cues <- c("VOT")

  # Make 5 ideal observers
  .io <-
    example_MVG_ideal_observer(1) %>%
    mutate(Sigma = map(Sigma, ~ .x * 5))
  .io.p20 <-
    .io %>%
    mutate(mu = map(mu, ~ .x + c(20)))
  .io.p40 <-
    .io %>%
    mutate(mu = map(mu, ~ .x + c(40)))
  .io.m20 <-
    .io %>%
    mutate(mu = map(mu, ~ .x + c(-20)))
  .io.m40 <-
    .io %>%
    mutate(mu = map(mu, ~ .x + c(-40)))

  # Sample responses for subjects that have converged against those five states
  .data <-
    bind_rows(
      sample_MVG_data_from_model(
        model = .io,
        Ns = n_trial,
        keep.input_parameters = F) %>%
        make_vector_column(cols = c("VOT"), vector_col = "cue") %>%
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
        model = .io.p20,
        Ns = n_trial,
        keep.input_parameters = F) %>%
        make_vector_column(cols = c("VOT"), vector_col = "cue") %>%
        crossing(Subject = 1:n_subject) %>%
        mutate(Condition = "plus20",
               Response = get_categorization_from_MVG_ideal_observer(
                 x = cue,
                 model = .io.p20,
                 decision_rule = "sampling",
                 simplify = T,
                 noise_treatment = "no_noise",
                 lapse_treatment = "no_lapses")),
      sample_MVG_data_from_model(
        model = .io.p40,
        Ns = n_trial,
        keep.input_parameters = F) %>%
        make_vector_column(cols = c("VOT"), vector_col = "cue") %>%
        crossing(Subject = 1:n_subject) %>%
        mutate(Condition = "plus40",
               Response = get_categorization_from_MVG_ideal_observer(
                 x = cue,
                 model = .io.p40,
                 decision_rule = "sampling",
                 simplify = T,
                 noise_treatment = "no_noise",
                 lapse_treatment = "no_lapses")),
      sample_MVG_data_from_model(
        model = .io.m20,
        Ns = n_trial,
        keep.input_parameters = F) %>%
        make_vector_column(cols = c("VOT"), vector_col = "cue") %>%
        crossing(Subject = 1:n_subject) %>%
        mutate(Condition = "minus20",
               Response = get_categorization_from_MVG_ideal_observer(
                 x = cue,
                 model = .io.m20,
                 decision_rule = "sampling",
                 simplify = T,
                 noise_treatment = "no_noise",
                 lapse_treatment = "no_lapses")),
      sample_MVG_data_from_model(
        model = .io.m40,
        Ns = n_trial,
        keep.input_parameters = F) %>%
        make_vector_column(cols = c("VOT"), vector_col = "cue") %>%
        crossing(Subject = 1:n_subject) %>%
        mutate(Condition = "minus40",
               Response = get_categorization_from_MVG_ideal_observer(
                 x = cue,
                 model = .io.m40,
                 decision_rule = "sampling",
                 simplify = T,
                 noise_treatment = "no_noise",
                 lapse_treatment = "no_lapses")))

  return(.data %>%
           crossing(Phase = c("exposure", "test")))
}


make_data_for_2Dstanfit_with_exposure <- function() {
  n_subject <- 30
  # number of trials in exposure per category per subject
  # (and there will be 2 * n_trials trials in test per subject)
  n_trial <- 50
  .cues <- c("VOT", "f0_semitones")

  # Make 5 ideal observers
  .io <-
    example_MVG_ideal_observer(2) %>%
    mutate(Sigma = map(Sigma, ~ .x * 5))
  .io.20.20 <-
    .io %>%
    mutate(mu = map(mu, ~ .x + c(20, 20)))
  .io.40.40 <-
    .io %>%
    mutate(mu = map(mu, ~ .x + c(40, 40)))
  .io.20.40 <-
    .io %>%
    mutate(mu = map(mu, ~ .x + c(20, 40)))
  .io.40.20 <-
    .io %>%
    mutate(mu = map(mu, ~ .x + c(40, 20)))

  # Sample responses for subjects that have converged against those five states
  .data <-
    bind_rows(
      sample_MVG_data_from_model(
        model = .io,
        Ns = n_trial,
        keep.input_parameters = F) %>%
        make_vector_column(cols = .cues, vector_col = "cue") %>%
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
        model = .io.20.20,
        Ns = n_trial,
        keep.input_parameters = F) %>%
        make_vector_column(cols = .cues, vector_col = "cue") %>%
        crossing(Subject = 1:n_subject) %>%
        mutate(Condition = "plus20.20",
               Response = get_categorization_from_MVG_ideal_observer(
                 x = cue,
                 model = .io.20.20,
                 decision_rule = "sampling",
                 simplify = T,
                 noise_treatment = "no_noise",
                 lapse_treatment = "no_lapses")),
      sample_MVG_data_from_model(
        model = .io.40.40,
        Ns = n_trial,
        keep.input_parameters = F) %>%
        make_vector_column(cols = .cues, vector_col = "cue") %>%
        crossing(Subject = 1:n_subject) %>%
        mutate(Condition = "plus40.40",
               Response = get_categorization_from_MVG_ideal_observer(
                 x = cue,
                 model = .io.40.40,
                 decision_rule = "sampling",
                 simplify = T,
                 noise_treatment = "no_noise",
                 lapse_treatment = "no_lapses")),
      sample_MVG_data_from_model(
        model = .io.20.40,
        Ns = n_trial,
        keep.input_parameters = F) %>%
        make_vector_column(cols = .cues, vector_col = "cue") %>%
        crossing(Subject = 1:n_subject) %>%
        mutate(Condition = "plus20.40",
               Response = get_categorization_from_MVG_ideal_observer(
                 x = cue,
                 model = .io.20.40,
                 decision_rule = "sampling",
                 simplify = T,
                 noise_treatment = "no_noise",
                 lapse_treatment = "no_lapses")),
      sample_MVG_data_from_model(
        model = .io.40.20,
        Ns = n_trial,
        keep.input_parameters = F) %>%
        make_vector_column(cols = .cues, vector_col = "cue") %>%
        crossing(Subject = 1:n_subject) %>%
        mutate(Condition = "plus40.20",
               Response = get_categorization_from_MVG_ideal_observer(
                 x = cue,
                 model = .io.40.20,
                 decision_rule = "sampling",
                 simplify = T,
                 noise_treatment = "no_noise",
                 lapse_treatment = "no_lapses")))

  return(.data %>%
           crossing(Phase = c("exposure", "test")))
}

make_data_for_3Dstanfit_with_exposure <- function() {
  n_subject <- 30
  # number of trials in exposure per category per subject
  # (and there will be 2 * n_trials trials in test per subject)
  n_trial <- 50
  .cues <- c("VOT", "f0_semitones", "vowel_duration")

  # Make 5 ideal observers
  .io <-
    example_MVG_ideal_observer(3) %>%
    mutate(Sigma = map(Sigma, ~ .x * 5))
  .io.20.20.20 <-
    .io %>%
    mutate(mu = map(mu, ~ .x + c(20, 20, 20)))
  .io.40.40.40 <-
    .io %>%
    mutate(mu = map(mu, ~ .x + c(40, 40, 40)))
  .io.20.40.60 <-
    .io %>%
    mutate(mu = map(mu, ~ .x + c(20, 40, 60)))
  .io.40.20.0 <-
    .io %>%
    mutate(mu = map(mu, ~ .x + c(40, 20, 0)))

  # Sample responses for subjects that have converged against those five states
  .data <-
    bind_rows(
      sample_MVG_data_from_model(
        model = .io,
        Ns = n_trial,
        keep.input_parameters = F) %>%
        make_vector_column(cols = .cues, vector_col = "cue") %>%
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
        model = .io.20.20.20,
        Ns = n_trial,
        keep.input_parameters = F) %>%
        make_vector_column(cols = .cues, vector_col = "cue") %>%
        crossing(Subject = 1:n_subject) %>%
        mutate(Condition = "plus20.20.20",
               Response = get_categorization_from_MVG_ideal_observer(
                 x = cue,
                 model = .io.20.20.20,
                 decision_rule = "sampling",
                 simplify = T,
                 noise_treatment = "no_noise",
                 lapse_treatment = "no_lapses")),
      sample_MVG_data_from_model(
        model = .io.40.40.40,
        Ns = n_trial,
        keep.input_parameters = F) %>%
        make_vector_column(cols = .cues, vector_col = "cue") %>%
        crossing(Subject = 1:n_subject) %>%
        mutate(Condition = "plus40.40.40",
               Response = get_categorization_from_MVG_ideal_observer(
                 x = cue,
                 model = .io.40.40.40,
                 decision_rule = "sampling",
                 simplify = T,
                 noise_treatment = "no_noise",
                 lapse_treatment = "no_lapses")),
      sample_MVG_data_from_model(
        model = .io.20.40.60,
        Ns = n_trial,
        keep.input_parameters = F) %>%
        make_vector_column(cols = .cues, vector_col = "cue") %>%
        crossing(Subject = 1:n_subject) %>%
        mutate(Condition = "plus20.40.60",
               Response = get_categorization_from_MVG_ideal_observer(
                 x = cue,
                 model = .io.20.40.60,
                 decision_rule = "sampling",
                 simplify = T,
                 noise_treatment = "no_noise",
                 lapse_treatment = "no_lapses")),
      sample_MVG_data_from_model(
        model = .io.40.20.0,
        Ns = n_trial,
        keep.input_parameters = F) %>%
        make_vector_column(cols = .cues, vector_col = "cue") %>%
        crossing(Subject = 1:n_subject) %>%
        mutate(Condition = "plus40.20.0",
               Response = get_categorization_from_MVG_ideal_observer(
                 x = cue,
                 model = .io.40.20.0,
                 decision_rule = "sampling",
                 simplify = T,
                 noise_treatment = "no_noise",
                 lapse_treatment = "no_lapses")))

  return(.data %>%
           crossing(Phase = c("exposure", "test")))
}


make_data_for_2Dstanfit_without_exposure <- function() {
  n_subject <- 60
  # number of trials in exposure per category per subject
  # (and there will be 2 * n_trials trials in test per subject)
  n_trial.exposure <- 90
  .cues <- c("cue1", "cue2")

  # Make ideal adaptor for prior
  .ia_0 <-
    example_MVG_ideal_observer(5) %>%
    lift_MVG_ideal_observer_to_NIW_ideal_adaptor(kappa = 10, nu = 100)
  # Update that ideal adaptor with shifted exposure
  # Shift 1
  .exposure_1 <-
    sample_MVG_data_from_model(
      model =
        example_MVG_ideal_observer(5) %>%
        mutate(
          mu = map(mu, ~ .x + c(-1, 3)),
          Sigma = ifelse(category == "B", map(Sigma, ~ .x * 2), Sigma)),
      Ns = n_trial.exposure,
      keep.input_parameters = F) %>%
    make_vector_column(cols = .cues, vector_col = "cue")
  .ia_1 <-
    .ia_0 %>%
    update_NIW_ideal_adaptor_batch(
      prior = .,
      exposure = .exposure_1,
      noise_treatment = "no_noise")
  # Shift 2
  .exposure_2 <-
    sample_MVG_data_from_model(
      model =
        example_MVG_ideal_observer(5) %>%
        mutate(mu = map(mu, ~ .x + c(4, -1))),
      Ns = n_trial.exposure,
      keep.input_parameters = F) %>%
    make_vector_column(cols = .cues, vector_col = "cue")
  .ia_2 <-
    .ia_0 %>%
    update_NIW_ideal_adaptor_batch(
      prior = .,
      exposure = .exposure_2,
      noise_treatment = "no_noise")
  # Shift 3
  .exposure_3 <-
    sample_MVG_data_from_model(
      model =
        example_MVG_ideal_observer(5) %>%
        mutate(
          mu = map(mu, ~ .x + c(-4, -2)),
          Sigma = ifelse(category == "B", map(Sigma, ~ .x * 2), Sigma)),
      Ns = n_trial.exposure,
      keep.input_parameters = F) %>%
    make_vector_column(cols = .cues, vector_col = "cue")
  .ia_3 <-
    .ia_0 %>%
    update_NIW_ideal_adaptor_batch(
      prior = .,
      exposure = .exposure_3,
      noise_treatment = "no_noise")

  # store exposure data
  df.exposure <-
    bind_rows(
      .exposure_1 %>% mutate(Condition = "shift_1"),
      .exposure_2 %>% mutate(Condition = "shift_2"),
      .exposure_3 %>% mutate(Condition = "shift_3")) %>%
    mutate(Phase = "exposure")

  # define a test grid
  df.test <-
    crossing(
      cue1 = seq(-5, 5, length.out = 10),
      cue2 = seq(-5, 5, length.out = 10),
      category = NA) %>%
    make_vector_column(cols = c("cue1", "cue2"), vector_col = "cue") %>%
    mutate(Phase = "test")

  # plot_expected_categories_contour2D(.ia_0) +
  #   geom_point(
  #     data = df.exposure,
  #     aes(x = cue1, y = cue2, shape = category, color = Condition)) +
  #   geom_point(
  #     data = df.test,
  #     aes(x = cue1, y = cue2), shape = 3, color = "black") +
  #   scale_color_manual(values = c("red", "blue", "green")) +
  #   theme_bw()

  # Sample tests responses for subjects after the four exposure conditions
  # (one of which is no_exposure)
  df.exposure %<>% crossing(Subject = 1:n_subject)
  df.test %<>% crossing(Subject = 1:n_subject)
  .data <-
    df.exposure %>%
    bind_rows(
      df.test %>%
        mutate(
          Condition = "no_exposure",
          Response =
            get_categorization_from_NIW_ideal_adaptor(
              x = cue,
              model = .ia_0,
              decision_rule = "sampling",
              simplify = T,
              noise_treatment = "no_noise",
              lapse_treatment = "no_lapses")),
      df.test %>%
        mutate(
          Condition = "shift_1",
          Response =
            get_categorization_from_NIW_ideal_adaptor(
              x = cue,
              model = .ia_1,
              decision_rule = "sampling",
              simplify = T,
              noise_treatment = "no_noise",
              lapse_treatment = "no_lapses")),
      df.test %>%
        mutate(
          Condition = "shift_2",
          Response =
            get_categorization_from_NIW_ideal_adaptor(
              x = cue,
              model = .ia_2,
              decision_rule = "sampling",
              simplify = T,
              noise_treatment = "no_noise",
              lapse_treatment = "no_lapses")),
      df.test %>%
        mutate(
          Condition = "shift_3",
          Response =
            get_categorization_from_NIW_ideal_adaptor(
              x = cue,
              model = .ia_3,
              decision_rule = "sampling",
              simplify = T,
              noise_treatment = "no_noise",
              lapse_treatment = "no_lapses")))

  .data %<>%
    arrange(Condition, Subject, Phase)

  return(.data)
}

make_data_for_3Dstanfit_without_exposure <- function() {
  n_subject <- 60
  # number of trials in exposure per category per subject
  # (and there will be 2 * n_trials trials in test per subject)
  n_trial.exposure <- 90
  .cues <- c("cue1", "cue2", "cue3")

  # Make ideal adaptor for prior
  .ia_0 <-
    example_MVG_ideal_observer(3) %>%
    lift_MVG_ideal_observer_to_NIW_ideal_adaptor(kappa = 10, nu = 100)
  # Update that ideal adaptor with shifted exposure
  # Shift 1
  .exposure_1 <-
    sample_MVG_data_from_model(
      model =
        example_MVG_ideal_observer(3) %>%
        mutate(
          mu = map(mu, ~ .x + c(-1, 3, 2)),
          Sigma = ifelse(category == "B", map(Sigma, ~ .x * 2), Sigma)),
      Ns = n_trial.exposure,
      keep.input_parameters = F) %>%
    make_vector_column(cols = .cues, vector_col = "cue")
  .ia_1 <-
    .ia_0 %>%
    update_NIW_ideal_adaptor_batch(
      prior = .,
      exposure = .exposure_1,
      noise_treatment = "no_noise")
  # Shift 2
  .exposure_2 <-
    sample_MVG_data_from_model(
      model =
        example_MVG_ideal_observer(3) %>%
        mutate(mu = map(mu, ~ .x + c(4, -1))),
      Ns = n_trial.exposure,
      keep.input_parameters = F) %>%
    make_vector_column(cols = .cues, vector_col = "cue")
  .ia_2 <-
    .ia_0 %>%
    update_NIW_ideal_adaptor_batch(
      prior = .,
      exposure = .exposure_2,
      noise_treatment = "no_noise")
  # Shift 3
  .exposure_3 <-
    sample_MVG_data_from_model(
      model =
        example_MVG_ideal_observer(3) %>%
        mutate(
          mu = map(mu, ~ .x + c(-4, -2)),
          Sigma = ifelse(category == "B", map(Sigma, ~ .x * 2), Sigma)),
      Ns = n_trial.exposure,
      keep.input_parameters = F) %>%
    make_vector_column(cols = .cues, vector_col = "cue")
  .ia_3 <-
    .ia_0 %>%
    update_NIW_ideal_adaptor_batch(
      prior = .,
      exposure = .exposure_3,
      noise_treatment = "no_noise")

  # store exposure data
  df.exposure <-
    bind_rows(
      .exposure_1 %>% mutate(Condition = "shift_1"),
      .exposure_2 %>% mutate(Condition = "shift_2"),
      .exposure_3 %>% mutate(Condition = "shift_3")) %>%
    mutate(Phase = "exposure")

  # define a test grid
  df.test <-
    crossing(
      cue1 = seq(-5, 5, length.out = 10),
      cue2 = seq(-5, 5, length.out = 10),
      category = NA) %>%
    make_vector_column(cols = c("cue1", "cue2"), vector_col = "cue") %>%
    mutate(Phase = "test")

  plot_expected_categories_contour2D(.ia_0) +
    geom_point(
      data = df.exposure,
      aes(x = cue1, y = cue2, shape = category, color = Condition)) +
    geom_point(
      data = df.test,
      aes(x = cue1, y = cue2), shape = 3, color = "black") +
    scale_color_manual(values = c("red", "blue", "green")) +
    theme_bw()

  # Sample tests responses for subjects after the four exposure conditions
  # (one of which is no_exposure)
  df.exposure %<>% crossing(Subject = 1:n_subject)
  df.test %<>% crossing(Subject = 1:n_subject)
  .data <-
    df.exposure %>%
    bind_rows(
      df.test %>%
        mutate(
          Condition = "no_exposure",
          Response =
            get_categorization_from_NIW_ideal_adaptor(
              x = cue,
              model = .ia_0,
              decision_rule = "sampling",
              simplify = T,
              noise_treatment = "no_noise",
              lapse_treatment = "no_lapses")),
      df.test %>%
        mutate(
          Condition = "shift_1",
          Response =
            get_categorization_from_NIW_ideal_adaptor(
              x = cue,
              model = .ia_1,
              decision_rule = "sampling",
              simplify = T,
              noise_treatment = "no_noise",
              lapse_treatment = "no_lapses")),
      df.test %>%
        mutate(
          Condition = "shift_2",
          Response =
            get_categorization_from_NIW_ideal_adaptor(
              x = cue,
              model = .ia_2,
              decision_rule = "sampling",
              simplify = T,
              noise_treatment = "no_noise",
              lapse_treatment = "no_lapses")),
      df.test %>%
        mutate(
          Condition = "shift_3",
          Response =
            get_categorization_from_NIW_ideal_adaptor(
              x = cue,
              model = .ia_3,
              decision_rule = "sampling",
              simplify = T,
              noise_treatment = "no_noise",
              lapse_treatment = "no_lapses")))

  .data %<>%
    arrange(Condition, Subject, Phase)

  return(.data)
}

get_example_stanfit <- function(
    example = 1,
    center.observations = T, scale.observations = F,
    stanmodel = "mvg_conj_sufficient_stats_lapse_cholesky",
    silent = 2, refresh = 0,
    seed = 42,
    ...
) {
  filename <-
    paste0(
      "../example-stanfit", example, "-",
      paste(c(if (!is.null(stanmodel)) stanmodel else "", seed, center.observations, scale.observations), collapse = "-"),
      ".rds")

  if (file.exists(filename)) {
    fit <- readRDS(filename)
  } else {
    .data <- make_data_for_stanfit(example, seed = seed)
    fit <-
      infer_prior_beliefs(
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
        center.observations = center.observations, scale.observations = scale.observations,
        stanmodel = stanmodel,
        control = list(adapt_delta = 0.9),
        file = filename,
        refresh = refresh,
        silent = silent,
        cores = 4,
        ...)
  }

  return(fit)
}
