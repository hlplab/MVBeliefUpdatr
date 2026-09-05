# Helper files are sourced into an environment whose parent is package:MVBeliefUpdatr,
# so packages attached after it are not on this file's lookup path. Base R is used
# throughout to keep this file independent of which packages happen to be attached.

n_subject <- 3
# number of trials in exposure per category per subject
n_exposure_trial <- 50
n_test_trial <- 125

make_data_for_stanfit <- function(example = 1, seed = NULL, verbose = F) {
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

# Shifting sampled cues is equivalent to shifting the category means of the generating model.
make_shifted_exposure <- function(.io, .cues, shifts) {
  parts <- lapply(names(shifts), function(condition) {
    d <- as.data.frame(sample_observations(.io, Ns = n_exposure_trial, randomize.order = TRUE))
    d[.cues] <- sweep(as.matrix(d[.cues]), 2, shifts[[condition]], "+")
    d$Condition <- condition
    d
  })

  .exposure <- do.call(rbind, parts)
  .exposure <- make_vector_column(.exposure, cols = .cues, vector_col = "cue")
  .exposure$Phase <- "exposure"
  rownames(.exposure) <- NULL

  .exposure
}

make_test_grid <- function(.exposure, .cues, n_per_cue) {
  cue_values <- lapply(.cues, function(cue) {
    seq(min(.exposure[[cue]]), max(.exposure[[cue]]), length.out = n_per_cue)
  })
  names(cue_values) <- .cues

  grid <- expand.grid(cue_values, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid$response <- NA
  grid <- make_vector_column(grid, cols = .cues, vector_col = "cue")

  index <- expand.grid(
    row = seq_len(nrow(grid)),
    Condition = unique(.exposure$Condition),
    Subject = seq_len(n_subject),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE)

  .test <- grid[index$row, , drop = FALSE]
  .test$Condition <- index$Condition
  .test$Subject <- factor(paste0(index$Condition, "_", index$Subject))
  rownames(.test) <- NULL

  .test
}

get_test_responses_after_updating_based_on_exposure <- function(.io, .exposure, .test, .cues) {
  .kappa <- 5
  .nu <- 100

  prior_model <- as_niw_ideal_adaptor(.io, kappa = .kappa, nu = .nu)

  # One updated model per condition, assuming all subjects update the same way.
  parts <- lapply(unique(.exposure$Condition), function(condition) {
    updated_model <-
      update_template(
        prior_model,
        .exposure[.exposure$Condition == condition, , drop = FALSE],
        noise_treatment = "no_noise")

    .test_condition <- .test[.test$Condition == condition, , drop = FALSE]
    .test_condition$Response <-
      as.character(
        categorize(
          updated_model,
          do.call(rbind, .test_condition$cue),
          decision_rule = "sampling")$category)
    .test_condition
  })

  .test_with_responses <- do.call(rbind, parts)
  .test_with_responses$Phase <- "test"
  .test_with_responses$category <- NA_character_

  # Exposure is shared across subjects within each condition.
  index <- expand.grid(
    row = seq_len(nrow(.exposure)),
    Subject = seq_len(n_subject),
    KEEP.OUT.ATTRS = FALSE)
  .exposure_by_subject <- .exposure[index$row, , drop = FALSE]
  .exposure_by_subject$Subject <- factor(paste0(.exposure$Condition[index$row], "_", index$Subject))
  .exposure_by_subject$category <- as.character(.exposure_by_subject$category)
  .exposure_by_subject$Response <- NA_character_

  columns <- c("Phase", "Condition", "Subject", .cues, "cue", "category", "Response")
  .data <- rbind(.exposure_by_subject[columns], .test_with_responses[columns])

  for (column in c("Phase", "Condition", "Subject", "category", "Response")) {
    .data[[column]] <- factor(.data[[column]])
  }
  rownames(.data) <- NULL

  .data
}

make_data_for_1Dstanfit_with_exposure <- function(verbose = F) {
  .cues <- c("VOT")

  # Make 5 ideal observers to sample EXPOSURE from
  .io <- example_mvg_ideal_observer(n_cues = 1)
  .exposure <-
    make_shifted_exposure(
      .io, .cues,
      list(
        baseline = 0,
        plus20 = 20,
        plus40 = 40,
        minus20 = -20,
        minus40 = -40))

  # Define a test grid
  .test <- make_test_grid(.exposure, .cues, n_test_trial)

  .data <- get_test_responses_after_updating_based_on_exposure(.io, .exposure, .test, .cues)
  return(.data)
}


make_data_for_2Dstanfit_with_exposure <- function(verbose = F) {
  .cues <- c("VOT", "f0_semitones")

  # Make 5 ideal observers to sample EXPOSURE from
  .io <- example_mvg_ideal_observer(n_cues = 2)
  .exposure <-
    make_shifted_exposure(
      .io, .cues,
      list(
        baseline = c(0, 0),
        plus20.20 = c(20, 20),
        plus40.40 = c(40, 40),
        plus20.40 = c(20, 40),
        plus40.20 = c(40, 20)))

  # Define a test grid
  .test <- make_test_grid(.exposure, .cues, ceiling(n_test_trial^(1/2)))

  .data <- get_test_responses_after_updating_based_on_exposure(.io, .exposure, .test, .cues)

  return(.data)
}

make_data_for_3Dstanfit_with_exposure <- function(verbose = F) {
  .cues <- c("VOT", "f0_semitones", "vowel_duration")

  # Make 5 ideal observers to sample EXPOSURE from
  .io <- example_mvg_ideal_observer(n_cues = 3)
  .exposure <-
    make_shifted_exposure(
      .io, .cues,
      list(
        baseline = c(0, 0, 0),
        plus20.20.20 = c(20, 20, 20),
        plus40.40.40 = c(40, 40, 40),
        plus20.40.60 = c(20, 40, 60),
        plus40.20.0 = c(40, 20, 0)))

  # Define a test grid
  .test <- make_test_grid(.exposure, .cues, ceiling(n_test_trial^(1/3)))

  .data <- get_test_responses_after_updating_based_on_exposure(.io, .exposure, .test, .cues)
  return(.data)
}

make_data_for_1Dstanfit_without_exposure <- function(verbose = F) {
  .data <- make_data_for_1Dstanfit_with_exposure()
  .data[.data$Condition != "baseline" | .data$Phase == "test", , drop = FALSE]
}

make_data_for_2Dstanfit_without_exposure <- function(verbose = F) {
  .data <- make_data_for_2Dstanfit_with_exposure()
  .data[.data$Condition != "baseline" | .data$Phase == "test", , drop = FALSE]
}

make_data_for_3Dstanfit_without_exposure <- function(verbose = F) {
  .data <- make_data_for_3Dstanfit_with_exposure()
  .data[.data$Condition != "baseline" | .data$Phase == "test", , drop = FALSE]
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
    new_ideal_adaptor_stanfit_input(
      exposure = .data[.data$Phase == "exposure", , drop = FALSE],
      test = .data[.data$Phase == "test", , drop = FALSE],
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
      fixed_parameters = list(lapse_rate = lapse_rate, mu_0 = mu_0, Sigma_0 = Sigma_0),
      control = control,
      stanmodel = stanmodel,
      verbose = verbose)

  return(.staninput)
}

example_stanfit_path <- function(example, stanmodel = "NIW_ideal_adaptor", seed = 42, control = control_staninput()) {
  transform_type <- control$transform_type
  file.path(
    .get_ideal_adaptor_fit_model_dir(),
    paste0(
      "example-stanfit-",
      paste(
        c(
          if (!is.null(stanmodel)) stanmodel else "",
          example,
          if (!is.null(transform_type)) transform_type else "",
          seed),
        collapse = "-"),
      ".rds"))
}

get_example_stanfit <- function(
    example = 1,
    silent = 2, refresh = 0, seed = 42, verbose = F,
    file_refit = "on_change",
    stanmodel = "NIW_ideal_adaptor",
    lapse_rate = NULL, mu_0 = NULL, Sigma_0 = NULL,
    control = control_staninput(),
    transform_type = NULL,
    filename = NULL,
    ...
) {
  if (!is.null(transform_type)) {
    control$transform_type <- transform_type
  }
  transform_type <- control$transform_type

  if (is.null(filename))
    filename <- example_stanfit_path(example, stanmodel = stanmodel, seed = seed, control = control)
  if (file.exists(filename) && file_refit == "never") {
    if (verbose) message("File already exists and file_refit is set to 'never'. Loading existing model from file.")
    return(read_stanfit(filename))
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
      stanfit_input = .staninput,
      stanmodel = stanmodel,
      file = filename, file_refit = file_refit,
      refresh = refresh,
      silent = silent, verbose = verbose,
      ...)

  return(fit)
}
