make_data_for_stanfit <- function(example = 1) {
  require(tidyverse)
  require(MVBeliefUpdatr)

  if (example == 1) {
    return(make_data_for_2Dstanfit_with_exposure())
  } else if (example == 2) {
    return(make_data_for_2Dstanfit_without_exposure())
  } else if (example == 3) {
    return(make_data_for_1Dstanfit_with_exposure())
  } else if (example == 5) {
    return(make_data_for_3Dstanfit_with_exposure())
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
    example_MVG_ideal_observer(2) %>%
    mutate(Sigma = map(Sigma, ~ .x * 5))
  .io.p2 <-
    .io %>%
    mutate(mu = map(mu, ~ .x + c(2)))
  .io.p4 <-
    .io %>%
    mutate(mu = map(mu, ~ .x + c(4)))
  .io.m2 <-
    .io %>%
    mutate(mu = map(mu, ~ .x + c(-2)))
  .io.m4 <-
    .io %>%
    mutate(mu = map(mu, ~ .x + c(-4)))

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
        model = .io.p2,
        Ns = n_trial,
        keep.input_parameters = F) %>%
        make_vector_column(cols = c("VOT"), vector_col = "cue") %>%
        crossing(Subject = 1:n_subject) %>%
        mutate(Condition = "plus2",
               Response = get_categorization_from_MVG_ideal_observer(
                 x = cue,
                 model = .io.p2,
                 decision_rule = "sampling",
                 simplify = T,
                 noise_treatment = "no_noise",
                 lapse_treatment = "no_lapses")),
      sample_MVG_data_from_model(
        model = .io.p4,
        Ns = n_trial,
        keep.input_parameters = F) %>%
        make_vector_column(cols = c("VOT"), vector_col = "cue") %>%
        crossing(Subject = 1:n_subject) %>%
        mutate(Condition = "plus4",
               Response = get_categorization_from_MVG_ideal_observer(
                 x = cue,
                 model = .io.p4,
                 decision_rule = "sampling",
                 simplify = T,
                 noise_treatment = "no_noise",
                 lapse_treatment = "no_lapses")),
      sample_MVG_data_from_model(
        model = .io.m2,
        Ns = n_trial,
        keep.input_parameters = F) %>%
        make_vector_column(cols = c("VOT"), vector_col = "cue") %>%
        crossing(Subject = 1:n_subject) %>%
        mutate(Condition = "minus2",
               Response = get_categorization_from_MVG_ideal_observer(
                 x = cue,
                 model = .io.m2,
                 decision_rule = "sampling",
                 simplify = T,
                 noise_treatment = "no_noise",
                 lapse_treatment = "no_lapses")),
      sample_MVG_data_from_model(
        model = .io.m4,
        Ns = n_trial,
        keep.input_parameters = F) %>%
        make_vector_column(cols = c("VOT"), vector_col = "cue") %>%
        crossing(Subject = 1:n_subject) %>%
        mutate(Condition = "minus4",
               Response = get_categorization_from_MVG_ideal_observer(
                 x = cue,
                 model = .io.m4,
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
    example_MVG_ideal_observer(1) %>%
    lift_MVG_ideal_observer_to_NIW_ideal_adaptor(kappa = 10, nu = 100)
  # Update that ideal adaptor with shifted exposure
  # Shift 1
  .exposure_1 <-
    sample_MVG_data_from_model(
      model =
        example_MVG_ideal_observer(1) %>%
        mutate(
          mu = map(mu, ~ .x + c(-1, 3)),
          Sigma = ifelse(category == "B", map(Sigma, ~ .x * 2), Sigma)),
      Ns = n_trial.exposure,
      keep.input_parameters = F) %>%
    make_vector_column(cols = c("cue1", "cue2"), vector_col = "cue")
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
        example_MVG_ideal_observer(1) %>%
        mutate(mu = map(mu, ~ .x + c(4, -1))),
      Ns = n_trial.exposure,
      keep.input_parameters = F) %>%
    make_vector_column(cols = c("cue1", "cue2"), vector_col = "cue")
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
        example_MVG_ideal_observer(1) %>%
        mutate(
          mu = map(mu, ~ .x + c(-4, -2)),
          Sigma = ifelse(category == "B", map(Sigma, ~ .x * 2), Sigma)),
      Ns = n_trial.exposure,
      keep.input_parameters = F) %>%
    make_vector_column(cols = c("cue1", "cue2"), vector_col = "cue")
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

make_data_for_3Dstanfit_without_exposure <- function() {
  n_subject <- 60
  # number of trials in exposure per category per subject
  # (and there will be 2 * n_trials trials in test per subject)
  n_trial.exposure <- 90
  .cues <- c("cue1", "cue2", "cue3")

  # Make ideal adaptor for prior
  .ia_0 <-
    example_MVG_ideal_observer(1) %>%
    lift_MVG_ideal_observer_to_NIW_ideal_adaptor(kappa = 10, nu = 100)
  # Update that ideal adaptor with shifted exposure
  # Shift 1
  .exposure_1 <-
    sample_MVG_data_from_model(
      model =
        example_MVG_ideal_observer(1) %>%
        mutate(
          mu = map(mu, ~ .x + c(-1, 3)),
          Sigma = ifelse(category == "B", map(Sigma, ~ .x * 2), Sigma)),
      Ns = n_trial.exposure,
      keep.input_parameters = F) %>%
    make_vector_column(cols = c("cue1", "cue2"), vector_col = "cue")
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
        example_MVG_ideal_observer(1) %>%
        mutate(mu = map(mu, ~ .x + c(4, -1))),
      Ns = n_trial.exposure,
      keep.input_parameters = F) %>%
    make_vector_column(cols = c("cue1", "cue2"), vector_col = "cue")
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
        example_MVG_ideal_observer(1) %>%
        mutate(
          mu = map(mu, ~ .x + c(-4, -2)),
          Sigma = ifelse(category == "B", map(Sigma, ~ .x * 2), Sigma)),
      Ns = n_trial.exposure,
      keep.input_parameters = F) %>%
    make_vector_column(cols = c("cue1", "cue2"), vector_col = "cue")
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

get_example_stanfit <- function(example = 1) {
  filename <- paste0("../example-stanfit", example, ".rds")
  if (file.exists(filename)) {
    fit <- readRDS(filename)
  } else {
    .data <- make_data_for_stanfit(example)
    fit <-
      infer_prior_beliefs(
        exposure = .data %>% filter(Phase == "exposure"),
        test = .data %>% filter(Phase == "test"),
        cues = unique(case_when(
          example %in% 1:2 ~  c("cue1", "cue2"),
          example %in% 3 ~ c("VOT"),
          example %in% 4 ~ c("cue1", "cue2", "cue3"),
          T ~ "")),
        category = "category",
        response = "Response",
        group = "Subject",
        group.unique = "Condition",
        file = filename,
        cores = 4)
  }

  return(fit)
}
