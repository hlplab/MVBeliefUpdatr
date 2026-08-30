make_minimal_staninput_data <- function(
  cues = c("cue1"),
  n_obs_exposure = 8L,
  n_obs_test = 8L,
  n_group = 2L,
  n_category = 2L
) {
  cue_names <- as.character(cues)
  category_levels <- paste0("A", seq_len(n_category))
  group_levels <- paste0("g", seq_len(n_group))

  exposure <- data.frame(
    group = factor(rep(group_levels, length.out = n_obs_exposure), levels = group_levels),
    category = factor(rep(category_levels, length.out = n_obs_exposure), levels = category_levels),
    response = factor(rep(category_levels, length.out = n_obs_exposure), levels = category_levels)
  )
  if (n_obs_exposure == 0L) {
    exposure <- exposure[0, , drop = FALSE]
  }

  if (length(cue_names) > 0) {
    base_values <- c(0.12, 0.82, 0.24, 0.76, 0.36, 0.64, 0.48, 0.88)
    for (i in seq_along(cue_names)) {
      cue_name <- cue_names[[i]]
      if (n_obs_exposure > 0L) {
        cue_values <- pmin(0.98, rep(base_values, length.out = n_obs_exposure) + 0.03 * (i - 1))
      } else {
        cue_values <- numeric(0)
      }
      exposure[[cue_name]] <- cue_values
    }
  }

  test <- data.frame(
    group = factor(rep(group_levels, length.out = n_obs_test), levels = group_levels),
    response = factor(rep(category_levels, length.out = n_obs_test), levels = category_levels)
  )
  if (n_obs_test == 0L) {
    test <- test[0, , drop = FALSE]
  }

  test_values <- c(0.15, 0.85, 0.30, 0.70, 0.45, 0.55, 0.60, 0.90)
  for (i in seq_along(cue_names)) {
    cue_name <- cue_names[[i]]
    if (n_obs_test > 0L) {
      test[[cue_name]] <- pmin(0.98, rep(test_values, length.out = n_obs_test) + 0.01 * (i - 1))
    } else {
      test[[cue_name]] <- numeric(0)
    }
  }

  list(exposure = exposure, test = test)
}

build_shifted_prior_fit_example <- function(
  cues,
  n_exposure_per_condition = 40L,
  n_test_per_condition = 30L,
  seed = 123L
) {
  cue_names <- as.character(cues)
  conditions <- c("low", "high")
  categories <- c("A", "B")
  set.seed(seed)

  make_cue_matrix <- function(category_labels, cue_names) {
    n <- length(category_labels)
    cue_count <- length(cue_names)
    cue_matrix <- matrix(NA_real_, nrow = n, ncol = cue_count)

    if (cue_count == 1L) {
      cue_matrix[, 1L] <- ifelse(
        category_labels == "A",
        stats::rnorm(n, mean = 0.20, sd = 0.08),
        stats::rnorm(n, mean = 0.80, sd = 0.08)
      )
    } else {
      for (i in seq_len(cue_count)) {
        cue_mean <- if (i == 1L) 0.20 else 0.30
        alt_mean <- if (i == 1L) 0.80 else 0.70
        cue_matrix[, i] <- ifelse(
          category_labels == "A",
          stats::rnorm(n, mean = cue_mean, sd = 0.08),
          stats::rnorm(n, mean = alt_mean, sd = 0.08)
        )
      }
    }

    cue_matrix <- matrix(
      pmin(0.98, pmax(0.02, cue_matrix)),
      nrow = n,
      ncol = length(cue_names)
    )
    colnames(cue_matrix) <- cue_names
    cue_matrix
  }

  n_exposure <- length(conditions) * n_exposure_per_condition
  n_test <- length(conditions) * n_test_per_condition

  exposure_condition <- factor(character(n_exposure), levels = conditions)
  exposure_category <- factor(character(n_exposure), levels = categories)
  exposure_group <- factor(character(n_exposure), levels = conditions)
  exposure_cues <- matrix(0, nrow = n_exposure, ncol = length(cue_names))
  colnames(exposure_cues) <- cue_names

  test_condition <- factor(character(n_test), levels = conditions)
  test_response <- factor(character(n_test), levels = categories)
  test_group <- factor(character(n_test), levels = conditions)
  test_cues <- matrix(0, nrow = n_test, ncol = length(cue_names))
  colnames(test_cues) <- cue_names

  exposure_idx <- 1L
  test_idx <- 1L

  for (condition in conditions) {
    p_category_a <- if (condition == "low") 0.25 else 0.75
    category_labels <- sample(categories, size = n_exposure_per_condition, replace = TRUE, prob = c(p_category_a, 1 - p_category_a))
    cues_matrix <- make_cue_matrix(category_labels, cue_names)

    for (i in seq_len(n_exposure_per_condition)) {
      exposure_condition[exposure_idx] <- condition
      exposure_category[exposure_idx] <- category_labels[[i]]
      exposure_group[exposure_idx] <- condition
      exposure_cues[exposure_idx, ] <- cues_matrix[i, , drop = FALSE]
      exposure_idx <- exposure_idx + 1L
    }

    test_category_labels <- sample(categories, size = n_test_per_condition, replace = TRUE, prob = c(p_category_a, 1 - p_category_a))
    test_cues_matrix <- make_cue_matrix(test_category_labels, cue_names)

    for (i in seq_len(n_test_per_condition)) {
      test_condition[test_idx] <- condition
      test_response[test_idx] <- test_category_labels[[i]]
      test_group[test_idx] <- condition
      test_cues[test_idx, ] <- test_cues_matrix[i, , drop = FALSE]
      test_idx <- test_idx + 1L
    }
  }

  exposure <- data.frame(
    Condition = exposure_condition,
    category = exposure_category,
    group = exposure_group,
    exposure_cues,
    stringsAsFactors = FALSE
  )
  exposure$Condition <- factor(exposure$Condition, levels = conditions)
  exposure$category <- factor(exposure$category, levels = categories)
  exposure$group <- factor(exposure$group, levels = conditions)

  test <- data.frame(
    Condition = test_condition,
    response = test_response,
    group = test_group,
    test_cues,
    stringsAsFactors = FALSE
  )
  test$Condition <- factor(test$Condition, levels = conditions)
  test$response <- factor(test$response, levels = categories)
  test$group <- factor(test$group, levels = "g1")

  list(exposure = exposure, test = test)
}

expect_staninput_structure <- function(input, expected_class, required_names, forbidden_names = character()) {
  expect_true(S7::S7_inherits(input, IdealAdaptorStanfitInput))
  expect_true(S7::S7_inherits(input@staninput, expected_class))
  expect_true(S7::S7_inherits(input@staninput, IdealAdaptorStaninput))
  expect_true(is.list(input@staninput@values))
  expect_true(all(required_names %in% names(input@staninput@values)))
  expect_false(any(forbidden_names %in% names(input@staninput@values)))
}

.get_ideal_adaptor_fit_model_dir <- function() {
  root <- if (exists("pkg_root", inherits = TRUE) && length(pkg_root) > 0L && is.character(pkg_root)) {
    pkg_root[1]
  } else {
    "."
  }
  path <- file.path(root, "tests", "testthat", "models")
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

save_ideal_adaptor_fit_model <- function(fit, name) {
  model_dir <- .get_ideal_adaptor_fit_model_dir()
  model_path <- file.path(model_dir, paste0(name, ".rds"))
  saveRDS(fit, model_path, compress = TRUE)
  model_path
}

load_ideal_adaptor_fit_model <- function(name) {
  model_path <- file.path(.get_ideal_adaptor_fit_model_dir(), paste0(name, ".rds"))
  if (!file.exists(model_path)) {
    stop("Model file not found: ", model_path)
  }
  readRDS(model_path)
}

run_fixed_param_stan_program <- function(stan_file, input, model_name) {
  skip_if_not_installed("rstan")

  model_key <- sub("_compat$", "", model_name)
  model <- NULL
  if (exists("stanmodels", envir = asNamespace("MVBeliefUpdatr"), inherits = FALSE)) {
    pkg_models <- get("stanmodels", envir = asNamespace("MVBeliefUpdatr"))
    if (!is.null(pkg_models[[model_name]])) {
      model <- pkg_models[[model_name]]
    } else if (!is.null(pkg_models[[model_key]])) {
      model <- pkg_models[[model_key]]
    } else if (!is.null(pkg_models[[paste0(model_key, "_ideal_adaptor")]])) {
      model <- pkg_models[[paste0(model_key, "_ideal_adaptor")]]
    }
  }

  if (is.null(model)) {
    model <- rstan::stan_model(
      file = stan_file,
      model_name = model_name,
      auto_write = FALSE
    )
  }

  fit <- rstan::sampling(
    object = model,
    data = input@staninput@values,
    chains = 1,
    iter = 1,
    warmup = 0,
    refresh = 0,
    seed = 123,
    algorithm = "Fixed_param"
  )

  expect_s4_class(fit, "stanfit")
}

