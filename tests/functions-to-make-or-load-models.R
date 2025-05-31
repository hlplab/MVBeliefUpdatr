n_subject <- 30
# number of trials in exposure per category per subject
n_exposure_trial <- 50
n_test_trial <- 15

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

get_test_responses_after_updating_based_on_exposure <- function(.io, .exposure, .test, .cues) {
  .kappa <- 5
  .nu <- 100

  # Define NIW and update it based on exposures, and then sample responses over test grid
  .ia <-
    .io %>%
    lift_MVG_ideal_observer_to_NIW_ideal_adaptor(kappa = .kappa, nu = .nu) %>%
    crossing(Condition = unique(.exposure$Condition)) %>%
    nest(model = -Condition) %>%
    # Create one update model for condition, assuming that all subjects update the same way
    mutate(
      model = map2(model, Condition, ~ update_NIW_ideal_adaptor_batch(
        prior = .x,
        exposure = .exposure %>% filter(Condition == .y),
        noise_treatment = "no_noise"))) %>%
    left_join(.test, by = "Condition") %>%
    # Nest subject, too, since that makes the sampling more efficient (and currently all subjects
    # in a condition have the same updated model)
    nest(test = c(Subject, all_of(.cues), cue, response)) %>%
    mutate(
      test = map2(test, model, ~
                    .x %>%
                    mutate(
                      Response = get_categorization_from_NIW_ideal_adaptor(
                        x = .x$cue,
                        model = .y,
                        decision_rule = "sampling",
                        simplify = T,
                        noise_treatment = "no_noise",
                        lapse_treatment = "marginalize"))))

  .data <-
    bind_rows(
      .exposure %>%
        crossing(Subject = 1:n_subject) %>%
        select(Phase, Condition, Subject, !!! syms(.cues), cue, category),
      .ia %>%
        unnest(test) %>%
        mutate(Phase = "test") %>%
        select(Phase, Condition, all_of(.cues), cue, Response)) %>%
    mutate(across(c(Phase, Condition, category, Response), factor))

  return(.data)
}

make_data_for_1Dstanfit_with_exposure <- function() {
  .cues <- c("VOT")

  # Make 5 ideal observers to sample EXPOSURE from
  .io <-
    example_MVG_ideal_observer(1) %>%
    mutate(Sigma = map(Sigma, ~ .x))
  .exposure <-
    .io %>%
    mutate(Condition = "baseline") %>%
    bind_rows(
      .io %>%
        mutate(
          mu = map(mu, ~ .x + c(20)),
          Condition = "plus20"),
      .io %>%
        mutate(
          mu = map(mu, ~ .x + c(40)),
          Condition = "plus40"),
      .io %>%
        mutate(
          mu = map(mu, ~ .x + c(-20)),
          Condition = "minus20"),
      .io %>%
        mutate(
          mu = map(mu, ~ .x + c(-40)),
          Condition = "minus40")) %>%
    nest(model = -Condition) %>%
    mutate(data = map(model, ~ sample_data_from_model(.x, Ns = n_exposure_trial, randomize.order = T))) %>%
    unnest(data) %>%
    make_vector_column(cols = .cues, vector_col = "cue") %>%
    mutate(Phase = "exposure")

  # Define a test grid
  .test <-
    crossing(
      !! sym(.cues[1]) := seq(min(.exposure[[.cues[1]]]), max(.exposure[[.cues[1]]]), length.out = n_test_trial),
      response = NA) %>%
    make_vector_column(cols = .cues, vector_col = "cue") %>%
    crossing(
      Condition = unique(.exposure$Condition),
      Subject = 1:n_subject)

  .data <- get_test_responses_after_updating_based_on_exposure(.io, .exposure, .test, .cues)
  return(.data)
}


make_data_for_2Dstanfit_with_exposure <- function() {
  .cues <- c("VOT", "f0_semitones")

  # Make 5 ideal observers to sample EXPOSURE from
  .io <-
    example_MVG_ideal_observer(2) %>%
    mutate(Sigma = map(Sigma, ~ .x))
  .exposure <-
    .io %>%
    mutate(Condition = "baseline") %>%
    bind_rows(
      .io %>%
        mutate(
          mu = map(mu, ~ .x + c(20, 20)),
          Condition = "plus20.20"),
      .io %>%
        mutate(
          mu = map(mu, ~ .x + c(40, 40)),
          Condition = "plus40.40"),
      .io %>%
        mutate(
          mu = map(mu, ~ .x + c(20, 40)),
          Condition = "plus20.40"),
      .io %>%
        mutate(
          mu = map(mu, ~ .x + c(40, 20)),
          Condition = "plus40.20")) %>%
    nest(model = -Condition) %>%
    mutate(data = map(model, ~ sample_data_from_model(.x, Ns = n_exposure_trial, randomize.order = T))) %>%
    unnest(data) %>%
    make_vector_column(cols = .cues, vector_col = "cue") %>%
    mutate(Phase = "exposure")

  # Define a test grid
  .test <-
    crossing(
      !! sym(.cues[1]) := seq(min(.exposure[[.cues[1]]]), max(.exposure[[.cues[1]]]), length.out = n_test_trial),
      !! sym(.cues[2]) := seq(min(.exposure[[.cues[2]]]), max(.exposure[[.cues[2]]]), length.out = n_test_trial),
      response = NA) %>%
    make_vector_column(cols = .cues, vector_col = "cue") %>%
    crossing(
      Condition = unique(.exposure$Condition),
      Subject = 1:n_subject)

  .data <- get_test_responses_after_updating_based_on_exposure(.io, .exposure, .test, .cues)

  p <-
    .data %>%
    filter(Phase == "exposure") %>%
    ggplot(aes(x = !! sym(.cues[1]), y = !! sym(.cues[2]))) +
    stat_ellipse(
      aes(fill = category),
      level = 0.95, geom = "polygon", alpha = 0.3) +
    geom_point(aes(color = category)) +
    facet_wrap(~ Condition, nrow = 1) +
    theme_bw() +
    title("Exposure data")
  plot(p)

  p <-
    .data %>%
    filter(Phase == "test") %>%
    group_by(Condition, !!! syms(.cues)) %>%
    summarise(meanResponse = mean(ifelse(Response == levels(.env$.data$Response)[2], 1, 0))) %>%
    ggplot(aes(x = !! sym(.cues[1]), y = !! sym(.cues[2]))) +
    geom_point(aes(color = meanResponse)) +
    scale_color_gradient(name = as.character(paste("Proportion of", levels(.data$Response)[2], "responses")), aesthetics = c("color", "fill"), low = "pink", high = "cyan") +
    facet_wrap(~ Condition, nrow = 1) +
    theme_bw() + theme(legend.position = "bottom") +
    title("Test data")
  plot(p)

  return(.data)
}

make_data_for_3Dstanfit_with_exposure <- function() {
  .cues <- c("VOT", "f0_semitones", "vowel_duration")

  # Make 5 ideal observers to sample EXPOSURE from
  .io <-
    example_MVG_ideal_observer(3) %>%
    mutate(Sigma = map(Sigma, ~ .x))
  .exposure <-
    .io %>%
    mutate(Condition = "baseline") %>%
    bind_rows(
      .io %>%
        mutate(
          mu = map(mu, ~ .x + c(20, 20, 20)),
          Condition = "plus20.20.20"),
      .io %>%
        mutate(
          mu = map(mu, ~ .x + c(40, 40, 40)),
          Condition = "plus40.40.40"),
      .io %>%
        mutate(
          mu = map(mu, ~ .x + c(20, 40, 60)),
          Condition = "plus20.40.60"),
      .io %>%
        mutate(
          mu = map(mu, ~ .x + c(40, 20, 0)),
          Condition = "plus40.20.0")) %>%
    nest(model = -Condition) %>%
    mutate(data = map(model, ~ sample_data_from_model(.x, Ns = n_exposure_trial, randomize.order = T))) %>%
    unnest(data) %>%
    make_vector_column(cols = .cues, vector_col = "cue") %>%
    mutate(Phase = "exposure")

  # Define a test grid
  .test <-
    crossing(
      !! sym(.cues[1]) := seq(min(.exposure[[.cues[1]]]), max(.exposure[[.cues[1]]]), length.out = n_test_trial),
      !! sym(.cues[2]) := seq(min(.exposure[[.cues[2]]]), max(.exposure[[.cues[2]]]), length.out = n_test_trial),
      !! sym(.cues[3]) := seq(min(.exposure[[.cues[3]]]), max(.exposure[[.cues[3]]]), length.out = n_test_trial),
      response = NA) %>%
    make_vector_column(cols = .cues, vector_col = "cue") %>%
    crossing(
      Condition = unique(.exposure$Condition),
      Subject = 1:n_subject)

  .data <- get_test_responses_after_updating_based_on_exposure(.io, .exposure, .test, .cues)
  return(.data)
}

make_data_for_1Dstanfit_without_exposure <- function() {
  make_data_for_1Dstanfit_with_exposure() %>%
    ungroup() %>%
    filter(Condition != "baseline" | Phase == "test")
}

make_data_for_2Dstanfit_without_exposure <- function() {
  make_data_for_2Dstanfit_with_exposure() %>%
    ungroup() %>%
    filter(Condition != "baseline" | Phase == "test")
}

make_data_for_3Dstanfit_without_exposure <- function() {
  make_data_for_3Dstanfit_with_exposure() %>%
    ungroup() %>%
    filter(Condition != "baseline" | Phase == "test")
}


get_example_staninput <- function(
    example = 1,
    stanmodel = "NIW_ideal_adaptor",
    transform_type = "PCA whiten",
    seed = 42
) {
  .data <- make_data_for_stanfit(example, seed = seed)
  .staninput <-
    make_staninput(
      exposure = .data %>% filter(Phase == "exposure"),
      test = .data %>% filter(Phase == "test"),
      cues =
        if (example %in% 1:3)
        {
          c("VOT", "f0_semitones", "vowel_duration")[1:(example)]
        } else if (example %in% 4:6) {
          c("VOT", "f0_semitones", "vowel_duration")[1:(example - 3)]
        },
      category = "category",
      response = "Response",
      group = "Subject",
      group.unique = "Condition",
      stanmodel = stanmodel,
      transform_type = transform_type)

  if (example %in% c(2, 5)) {
    require(ellipse)
    df <- tibble()

    # Get mean and cov
    L <- .staninput$staninput$transformed$L
    M <- .staninput$staninput$transformed$M
    for (m in 1:M) {
      for (l in 1:L) {
        df %<>%
          bind_rows(
            tibble(
              group = l,
              category = m,
              N = .staninput$staninput$transformed$N_exposure[m, l],
              center = list(.staninput$staninput$transformed$x_mean_exposure[m, l,]),
              uss_matrix = list(.staninput$staninput$transformed$x_ss_exposure[m, l, , ]),
              css_matrix = pmap(list(uss_matrix, N, center), ~ uss2css(..1, n = ..2, mean = ..3)),
              cov_matrix = map2(css_matrix, N, ~ css2cov(.x, .y))) %>%
              mutate(ellipse = map2(cov_matrix, center, ~ ellipse(x = .x, centre = .y))) %>%
              # # This step is necessary since unnest() can't yet unnest lists of matrices
              # # (bug was reported and added as milestone, 11/2019)
              # mutate(ellipse = map(ellipse, ~ as_tibble(.x, .name_repair = "unique"))) %>%
              # select(-c(kappa, nu, m, S, Sigma, lapse_rate)) %>%
              unnest(ellipse)
              )
      }
    }


    eigen_decomp <- eigen(cov_matrix)
    eigenvalues <- eigen_decomp$values
    eigenvectors <- eigen_decomp$vectors

    # 4. Scale ellipse axes by square root of eigenvalues
    scale <- sqrt(eigenvalues)

    # 5. Define rotation angle
    angle <- atan2(eigenvectors[2, 1], eigenvectors[1, 1])

    # 6. Create data for plotting the ellipse
    angles <- seq(0, 2*pi, length.out = 100)
    unit_circle <- cbind(cos(angles), sin(angles))
    rotated_unit_circle <-
      cbind(cos(angle) * unit_circle[, 1] - sin(angle) * unit_circle[, 2],
            sin(angle) * unit_circle[, 1] + cos(angle) * unit_circle[, 2])
    ellipse_df <- data.frame(
      x = center[1] + scale[1] * rotated_unit_circle[, 1],
      y = center[2] + scale[2] * rotated_unit_circle[, 2]
    )

    # 7. Plot the data and the ellipse
    ggplot(df, aes(x = x, y = y)) +
      geom_point() +
      geom_path(data = ellipse_df, aes(x = x, y = y), color = "red")

    p <-
      .staninput$staninput$transformed %>%
      filter(Phase == "exposure") %>%
      ggplot(aes(x = !! sym(.cues[1]), y = !! sym(.cues[2]))) +
      stat_ellipse(
        aes(fill = category),
        level = 0.95, geom = "polygon", alpha = 0.3) +
      geom_point(aes(color = category)) +
      facet_wrap(~ Condition, nrow = 1) +
      theme_bw() +
      title("Exposure statistics in staninput")
    plot(p)

  }

  return(.staninput)
}

get_example_stanfit <- function(
    example = 1,
    silent = 2, refresh = 0, seed = 42,
    file_refit = "on_change",
    stanmodel = "NIW_ideal_adaptor",
    transform_type = "PCA whiten",
    ...
) {
  .staninput <-
    get_example_staninput(
      example = example,
      stanmodel = stanmodel,
      transform_type = transform_type,
      seed = seed)

  filename <-
    paste0(
      "../example-stanfit-",
      paste(
        c(
          if (!is.null(stanmodel)) stanmodel else "",
          example,
          if (!is.null(transform_type)) transform_type else "",
          seed),
        collapse = "-"),
      ".rds")

  fit <-
    fit_ideal_adaptor(
      staninput = .staninput,
      stanmodel = stanmodel,
      file = filename, file_refit = file_refit,
      refresh = refresh,
      silent = silent,
      cores = 4,
      ...)

  return(fit)
}
