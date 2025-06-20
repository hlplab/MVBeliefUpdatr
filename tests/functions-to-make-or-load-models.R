n_subject <- 30
# number of trials in exposure per category per subject
n_exposure_trial <- 50
n_test_trial <- 125

make_data_for_stanfit <- function(example = 1, seed = NULL, verbose = F) {
  require(tidyverse)
  require(magrittr)
  require(MVBeliefUpdatr)

  if (!is.null(seed)) set.seed(seed)

  if (example == 1) {
    return(make_data_for_1Dstanfit_with_exposure(verbose = verbose))
  } else if (example == 2) {
    return(make_data_for_2Dstanfit_with_exposure(verbose = verbose))
  } else if (example == 3) {
    return(make_data_for_3Dstanfit_with_exposure(verbose = verbose))
  } else if (example == 4) {
    return(make_data_for_1Dstanfit_without_exposure(verbose = verbose))
  } else if (example == 5) {
    return(make_data_for_2Dstanfit_without_exposure(verbose = verbose))
  } else if (example == 6) {
    return(make_data_for_3Dstanfit_without_exposure(verbose = verbose))
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
        prior_model = .x,
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
        crossing(Subject = factor(1:n_subject)) %>%
        select(Phase, Condition, Subject, !!! syms(.cues), cue, category),
      .ia %>%
        unnest(test) %>%
        mutate(Phase = "test") %>%
        select(Phase, Condition, Subject, all_of(.cues), cue, Response)) %>%
    mutate(across(c(Phase, Condition, Subject, category, Response), factor))

  return(.data)
}

make_data_for_1Dstanfit_with_exposure <- function(verbose = F) {
  .cues <- c("VOT")

  # Make 5 ideal observers to sample EXPOSURE from
  .io <-
    example_MVG_ideal_observer(1, verbose = verbose) %>%
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
      Subject = factor(1:n_subject))

  .data <- get_test_responses_after_updating_based_on_exposure(.io, .exposure, .test, .cues)
  return(.data)
}


make_data_for_2Dstanfit_with_exposure <- function(verbose = F, plot = F) {
  .cues <- c("VOT", "f0_semitones")

  # Make 5 ideal observers to sample EXPOSURE from
  .io <-
    example_MVG_ideal_observer(2, verbose = verbose) %>%
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
      !! sym(.cues[1]) := seq(min(.exposure[[.cues[1]]]), max(.exposure[[.cues[1]]]), length.out = ceiling(n_test_trial^(1/2))),
      !! sym(.cues[2]) := seq(min(.exposure[[.cues[2]]]), max(.exposure[[.cues[2]]]), length.out = ceiling(n_test_trial^(1/2))),
      response = NA) %>%
    make_vector_column(cols = .cues, vector_col = "cue") %>%
    crossing(
      Condition = unique(.exposure$Condition),
      Subject = factor(1:n_subject))

  .data <- get_test_responses_after_updating_based_on_exposure(.io, .exposure, .test, .cues)

  if (plot) {
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
      ggtitle("Exposure data")
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
      ggtitle("Test data")
    plot(p)
  }

  return(.data)
}

make_data_for_3Dstanfit_with_exposure <- function(verbose = F) {
  .cues <- c("VOT", "f0_semitones", "vowel_duration")

  # Make 5 ideal observers to sample EXPOSURE from
  .io <-
    example_MVG_ideal_observer(3, verbose = verbose) %>%
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
      !! sym(.cues[1]) := seq(min(.exposure[[.cues[1]]]), max(.exposure[[.cues[1]]]), length.out = ceiling(n_test_trial^(1/3))),
      !! sym(.cues[2]) := seq(min(.exposure[[.cues[2]]]), max(.exposure[[.cues[2]]]), length.out = ceiling(n_test_trial^(1/3))),
      !! sym(.cues[3]) := seq(min(.exposure[[.cues[3]]]), max(.exposure[[.cues[3]]]), length.out = ceiling(n_test_trial^(1/3))),
      response = NA) %>%
    make_vector_column(cols = .cues, vector_col = "cue") %>%
    crossing(
      Condition = unique(.exposure$Condition),
      Subject = factor(1:n_subject))

  .data <- get_test_responses_after_updating_based_on_exposure(.io, .exposure, .test, .cues)
  return(.data)
}

make_data_for_1Dstanfit_without_exposure <- function(verbose = F) {
  make_data_for_1Dstanfit_with_exposure() %>%
    ungroup() %>%
    filter(Condition != "baseline" | Phase == "test")
}

make_data_for_2Dstanfit_without_exposure <- function(verbose = F) {
  make_data_for_2Dstanfit_with_exposure() %>%
    ungroup() %>%
    filter(Condition != "baseline" | Phase == "test")
}

make_data_for_3Dstanfit_without_exposure <- function(verbose = F) {
  make_data_for_3Dstanfit_with_exposure() %>%
    ungroup() %>%
    filter(Condition != "baseline" | Phase == "test")
}


get_example_staninput <- function(
    example = 1,
    stanmodel = "NIW_ideal_adaptor",
    lapse_rate = NULL, mu_0 = NULL, Sigma_0 = NULL,
    control = control_staninput(),
    seed = 42, verbose = F
) {
  .data <- make_data_for_stanfit(example, seed = seed, verbose = verbose)
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
      lapse_rate = lapse_rate, mu_0 = mu_0, Sigma_0 = Sigma_0,
      control = control,
      stanmodel = stanmodel,
      verbose = verbose)

  # # For debugging:
  # if (example %in% c(2, 5)) {
  #   require(ellipse)
  #   df.transformed <- df.untransformed <- tibble()
  #
  #   # Get mean and cov
  #   L <- .staninput$staninput$transformed$L
  #   M <- .staninput$staninput$transformed$M
  #   for (m in 1:M) {
  #     for (l in 1:L) {
  #       df.transformed %<>%
  #         bind_rows(
  #           tibble(
  #             Condition = l,
  #             category = m,
  #             N = .staninput$staninput$transformed$N_exposure[m, l],
  #             center = list(.staninput$staninput$transformed$x_mean_exposure[m, l,]),
  #             uss_matrix = list(.staninput$staninput$transformed$x_ss_exposure[m, l, , ]),
  #             css_matrix = pmap(list(uss_matrix, N, center), ~ uss2css(..1, n = ..2, mean = ..3)),
  #             cov_matrix = map2(css_matrix, N, ~ css2cov(.x, .y))) %>%
  #             mutate(ellipse = map2(cov_matrix, center, ~ ellipse(x = .x, centre = .y))) %>%
  #             # # This step is necessary since unnest() can't yet unnest lists of matrices
  #             # # (bug was reported and added as milestone, 11/2019)
  #             mutate(ellipse = map(ellipse, ~ as_tibble(.x, .name_repair = "unique"))) %>%
  #             unnest(ellipse))
  #
  #       df.untransformed %<>%
  #         bind_rows(
  #           tibble(
  #             Condition = l,
  #             category = m,
  #             N = .staninput$staninput$untransformed$N_exposure[m, l],
  #             center = list(.staninput$staninput$untransformed$x_mean_exposure[m, l,]),
  #             uss_matrix = list(.staninput$staninput$untransformed$x_ss_exposure[m, l, , ]),
  #             css_matrix = pmap(list(uss_matrix, N, center), ~ uss2css(..1, n = ..2, mean = ..3)),
  #             cov_matrix = map2(css_matrix, N, ~ css2cov(.x, .y))) %>%
  #             mutate(ellipse = map2(cov_matrix, center, ~ ellipse(x = .x, centre = .y))) %>%
  #             # # This step is necessary since unnest() can't yet unnest lists of matrices
  #             # # (bug was reported and added as milestone, 11/2019)
  #             mutate(ellipse = map(ellipse, ~ as_tibble(.x, .name_repair = "unique"))) %>%
  #             unnest(ellipse))
  #     }
  #   }
  #
  #   # Plot the data and the ellipse
  #   p <-
  #     df.transformed %>%
  #     mutate(across(c(Condition, category), factor)) %>%
  #     ggplot(aes(x = x, y = y)) +
  #     geom_path(aes(x = x, y = y, color = category)) +
  #     facet_wrap(~ Condition) +
  #     theme_bw() +
  #     ggtitle("Exposure statistics in staninput (transformed)")
  #   plot(p)
  #
  #   p <-
  #     df.untransformed %<>%
  #     mutate(across(c(Condition, category), factor)) %>%
  #     ggplot(aes(x = x, y = y)) +
  #     geom_path(aes(x = x, y = y, color = category)) +
  #     facet_wrap(~ Condition) +
  #     theme_bw() +
  #     ggtitle("Exposure statistics in staninput (untransformed)")
  #   plot(p)
  # }

  return(.staninput)
}

get_example_stanfit <- function(
    example = 1,
    silent = 2, refresh = 0, seed = 42, verbose = F,
    file_refit = "on_change",
    stanmodel = "NIW_ideal_adaptor",
    lapse_rate = NULL, mu_0 = NULL, Sigma_0 = NULL,
    control = control_staninput(),
    filename = NULL,
    ...
) {
  transform_type <- control$transform_type

  if (is.null(filename))
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
  if (file.exists(filename) && file_refit == "never") {
    if (verbose) message("File already exists and file_refit is set to 'never'. Loading existing model from file.")
    return(MVBeliefUpdatr:::read_ideal_adaptor_stanfit(filename))
  }

  .staninput <-
    get_example_staninput(
      example = example,
      stanmodel = stanmodel,
      control = control,
      lapse_rate = lapse_rate, mu_0 = mu_0, Sigma_0 = Sigma_0,
      seed = seed, verbose = verbose)

  fit <-
    fit_ideal_adaptor(
      staninput = .staninput,
      stanmodel = stanmodel,
      file = filename, file_refit = file_refit,
      refresh = refresh,
      silent = silent, verbose = verbose,
      ...)

  return(fit)
}
