#' @include S7-core-classes.R
#' @include S7-transform-information.R

#' @name MVBU-Staninput-class
#' @title S7 base class for processed Stan input
#' @docType class
#' @slot values A named list containing the Stan input values.
#' @export
MVBU_Staninput <- S7::new_class(
  "MVBU_Staninput",
  package = NULL,
  parent = MVBU_Object,
  properties = list(
    values = S7::class_list
  ),
  constructor = function(values = list()) {
    S7::new_object(
      MVBU_Object(),
      values = as.list(values)
    )
  },
  validator = function(self) {
    if (!is.list(self@values)) {
      return("`values` must be a list")
    }
    NULL
  }
)

#' S7 class for ideal-adaptor Stan input.
#'
#' @name IdealAdaptorStaninput-class
#' @docType class
#' @export
IdealAdaptorStaninput <- S7::new_class(
  "IdealAdaptorStaninput",
  package = NULL,
  parent = MVBU_Staninput,
  constructor = function(values = list()) {
    S7::new_object(
      MVBU_Staninput(values = as.list(values))
    )
  }
)

NIX_IdealAdaptorStaninput <- S7::new_class(
  "NIX_IdealAdaptorStaninput",
  package = NULL,
  parent = IdealAdaptorStaninput,
  constructor = function(values = list()) {
    S7::new_object(
      IdealAdaptorStaninput(values = as.list(values))
    )
  }
)

MNIX_IdealAdaptorStaninput <- S7::new_class(
  "MNIX_IdealAdaptorStaninput",
  package = NULL,
  parent = IdealAdaptorStaninput,
  constructor = function(values = list()) {
    S7::new_object(
      IdealAdaptorStaninput(values = as.list(values))
    )
  }
)

NIW_IdealAdaptorStaninput <- S7::new_class(
  "NIW_IdealAdaptorStaninput",
  package = NULL,
  parent = IdealAdaptorStaninput,
  constructor = function(values = list()) {
    S7::new_object(
      IdealAdaptorStaninput(values = as.list(values))
    )
  }
)

#' Ensure scalar control values are passed to Stan as the correct R object type.
#'
#' Internal helper that preserves vector semantics for single-cue models so
#' Stan receives arrays rather than bare scalars for vector-valued data.
#'
#' @keywords internal
#' @noRd
.make_stan_tauscale <- function(tau_scale, n_cues) {
  tau_scale <- as.numeric(tau_scale)
  if (length(tau_scale) == 1) {
    tau_scale <- rep(tau_scale, n_cues)
  }
  array(tau_scale, dim = length(tau_scale))
}

.make_stan_shift <- function(shift, n_cues) {
  shift <- as.numeric(shift)
  if (length(shift) == 1) {
    shift <- rep(shift, n_cues)
  }
  array(shift, dim = length(shift))
}

.make_stan_inv_scale <- function(inv_scale, n_cues) {
  if (length(inv_scale) == 1 && n_cues == 1) {
    return(matrix(as.numeric(inv_scale), nrow = 1, ncol = 1))
  }
  if (is.matrix(inv_scale)) {
    return(inv_scale)
  }
  matrix(as.numeric(inv_scale), nrow = n_cues, ncol = n_cues)
}


#' Build the typed Stan input object for the NIX ideal-adaptor model.
#'
#' Internal helper that summarizes exposure and test data and assembles the
#' model-specific list expected by the NIX Stan program, then wraps it in an
#' S7 staninput object.
#'
#' @keywords internal
#' @noRd
new_nix_staninput <- function(
  exposure,
  test,
  cues,
  category,
  response,
  group,
  category_levels,
  group_levels,
  tau_scale,
  L_omega_eta,
  split_loglik_per_observation,
  lapse_rate,
  mu_0,
  Sigma_0,
  transform,
  n_cues,
  n_categories,
  n_groups
) {
  exposure_summary <- .summarize_exposure(
    exposure,
    cues,
    category,
    group,
    category_levels,
    group_levels,
    model = "NIX"
  )
  test_summary <- .summarize_test(
    test,
    cues,
    response,
    group,
    category_levels,
    group_levels,
    model = "NIX"
  )

  staninput <- list(
    K = n_cues,
    M = n_categories,
    L = n_groups,
    tau_scale = if (n_cues == 1) {
      as.numeric(tau_scale[1])
    } else {
      .make_stan_tauscale(tau_scale, n_cues)
    },
    L_omega_eta = as.numeric(L_omega_eta),
    split_loglik_per_observation = as.numeric(split_loglik_per_observation),
    lapse_rate_known = if (is.null(lapse_rate)) 0 else 1,
    N_exposure = exposure_summary$N_exposure,
    x_mean_exposure = exposure_summary$x_mean_exposure,
    x_sd_exposure = exposure_summary$x_sd_exposure,
    x_test = test_summary$x_test,
    y_test = test_summary$y_test,
    z_test_counts = test_summary$z_test_counts,
    N_test = test_summary$N_test,
    mu_0_known = if (is.null(mu_0)) 0 else 1,
    Sigma_0_known = if (is.null(Sigma_0)) 0 else 1,
    shift = if (n_cues == 1) {
      as.numeric(transform$transform.parameters[["shift"]][1])
    } else {
      .make_stan_shift(transform$transform.parameters[["shift"]], n_cues)
    },
    INV_SCALE = if (n_cues == 1) {
      as.numeric(transform$transform.parameters[["INV_SCALE"]])
    } else {
      .make_stan_inv_scale(transform$transform.parameters[["INV_SCALE"]], n_cues)
    }
  )

  if (!is.null(lapse_rate)) {
    staninput$lapse_rate_data <- as.numeric(lapse_rate)
  } else {
    staninput$lapse_rate_data <- numeric(0)
  }
  if (!is.null(mu_0)) {
    staninput$mu_0_data <- to_array(
      vapply(mu_0, function(x) as.numeric(x[1]), numeric(1)),
      inner_dims = length(mu_0),
      outer_dims = NULL,
      simplify = FALSE
    )
  } else {
    staninput$mu_0_data <- numeric(0)
  }
  if (!is.null(Sigma_0)) {
    staninput$Sigma_0_data <- to_array(
      vapply(Sigma_0, function(x) as.numeric(x[1, 1]), numeric(1)),
      inner_dims = length(Sigma_0),
      outer_dims = NULL,
      simplify = FALSE
    )
  } else {
    staninput$Sigma_0_data <- numeric(0)
  }

  NIX_IdealAdaptorStaninput(values = staninput)
}

#' Build the typed Stan input object for the NIW ideal-adaptor model.
#'
#' Internal helper that summarizes exposure and test data and assembles the
#' model-specific list expected by the NIW Stan program, then wraps it in an
#' S7 staninput object.
#'
#' @keywords internal
#' @noRd
new_niw_staninput <- function(
  exposure,
  test,
  cues,
  category,
  response,
  group,
  category_levels,
  group_levels,
  tau_scale,
  L_omega_eta,
  split_loglik_per_observation,
  lapse_rate,
  mu_0,
  Sigma_0,
  transform,
  n_cues,
  n_categories,
  n_groups
) {
  exposure_summary <- .summarize_exposure(
    exposure,
    cues,
    category,
    group,
    category_levels,
    group_levels,
    model = "NIW"
  )
  test_summary <- .summarize_test(
    test,
    cues,
    response,
    group,
    category_levels,
    group_levels,
    model = "NIW"
  )

  staninput <- list(
    K = n_cues,
    M = n_categories,
    L = n_groups,
    tau_scale = .make_stan_tauscale(tau_scale, n_cues),
    L_omega_eta = as.numeric(L_omega_eta),
    split_loglik_per_observation = as.numeric(split_loglik_per_observation),
    lapse_rate_known = if (is.null(lapse_rate)) 0 else 1,
    N_exposure = exposure_summary$N_exposure,
    x_mean_exposure = exposure_summary$x_mean_exposure,
    x_ss_exposure = exposure_summary$x_ss_exposure,
    x_test = test_summary$x_test,
    y_test = test_summary$y_test,
    z_test_counts = test_summary$z_test_counts,
    N_test = test_summary$N_test,
    mu_0_known = if (is.null(mu_0)) 0 else 1,
    Sigma_0_known = if (is.null(Sigma_0)) 0 else 1,
    shift = .make_stan_shift(transform$transform.parameters[["shift"]], n_cues),
    INV_SCALE = .make_stan_inv_scale(
      transform$transform.parameters[["INV_SCALE"]],
      n_cues
    )
  )

  if (!is.null(lapse_rate)) {
    staninput$lapse_rate_data <- as.numeric(lapse_rate)
  } else {
    staninput$lapse_rate_data <- numeric(0)
  }
  if (!is.null(mu_0)) {
    mu_0_data <- to_array(
      lapply(mu_0, function(x) as.numeric(x)),
      inner_dims = c(n_cues),
      outer_dims = c(length(mu_0)),
      simplify = FALSE
    )
    staninput$mu_0_data <- mu_0_data
  } else {
    staninput$mu_0_data <- array(0, dim = c(0, 0))
  }
  if (!is.null(Sigma_0)) {
    Sigma_0_data <- to_array(
      lapply(Sigma_0, function(x) as.matrix(x)),
      inner_dims = c(n_cues, n_cues),
      outer_dims = c(length(Sigma_0)),
      simplify = FALSE
    )
    staninput$Sigma_0_data <- Sigma_0_data
  } else {
    staninput$Sigma_0_data <- array(0, dim = c(0, 0, 0))
  }

  NIW_IdealAdaptorStaninput(values = staninput)
}

#' Build the typed Stan input object for the MNIX ideal-adaptor model.
#'
#' Internal helper that summarizes exposure and test data and assembles the
#' model-specific list expected by the MNIX Stan program, then wraps it in an
#' S7 staninput object.
#'
#' @keywords internal
#' @noRd
new_mnix_staninput <- function(
  exposure,
  test,
  cues,
  category,
  response,
  group,
  category_levels,
  group_levels,
  tau_scale,
  L_omega_eta,
  split_loglik_per_observation,
  lapse_rate,
  mu_0,
  Sigma_0,
  transform,
  n_cues,
  n_categories,
  n_groups
) {
  exposure_summary <- .summarize_exposure(
    exposure,
    cues,
    category,
    group,
    category_levels,
    group_levels,
    model = "MNIX"
  )
  test_summary <- .summarize_test(
    test,
    cues,
    response,
    group,
    category_levels,
    group_levels,
    model = "MNIX"
  )

  staninput <- list(
    K = n_cues,
    M = n_categories,
    L = n_groups,
    tau_scale = .make_stan_tauscale(tau_scale, n_cues),
    L_omega_eta = as.numeric(L_omega_eta),
    split_loglik_per_observation = as.numeric(split_loglik_per_observation),
    lapse_rate_known = if (is.null(lapse_rate)) 0 else 1,
    p_cat = array(rep(1 / n_categories, n_categories), dim = c(n_categories)),
    N_exposure = exposure_summary$N_exposure,
    x_mean_exposure = exposure_summary$x_mean_exposure,
    x_ss_exposure = exposure_summary$x_ss_exposure,
    x_test = test_summary$x_test,
    y_test = test_summary$y_test,
    z_test_counts = test_summary$z_test_counts,
    N_test = test_summary$N_test,
    mu_0_known = if (is.null(mu_0)) 0 else 1,
    Sigma_0_known = if (is.null(Sigma_0)) 0 else 1,
    shift = .make_stan_shift(transform$transform.parameters[["shift"]], n_cues),
    INV_SCALE = .make_stan_inv_scale(
      transform$transform.parameters[["INV_SCALE"]],
      n_cues
    )
  )

  if (!is.null(lapse_rate)) {
    staninput$lapse_rate_data <- as.numeric(lapse_rate)
  } else {
    staninput$lapse_rate_data <- numeric(0)
  }
  if (!is.null(mu_0)) {
    mu_0_data <- to_array(
      lapply(mu_0, function(x) as.numeric(x)),
      inner_dims = c(n_cues),
      outer_dims = c(length(mu_0)),
      simplify = FALSE
    )
    staninput$mu_0_data <- mu_0_data
  } else {
    staninput$mu_0_data <- array(0, dim = c(0, 0))
  }
  if (!is.null(Sigma_0)) {
    Sigma_0_data <- to_array(
      lapply(Sigma_0, function(x) as.matrix(x)),
      inner_dims = c(n_cues),
      outer_dims = c(length(Sigma_0)),
      simplify = FALSE
    )
    staninput$Sigma_0_data <- Sigma_0_data
  } else {
    staninput$Sigma_0_data <- array(0, dim = c(0, 0))
  }

  MNIX_IdealAdaptorStaninput(values = staninput)
}

#' Summarize exposure data for Stanfit input construction.
#'
#' Internal helper that computes exposure counts and category- or group-specific
#' cue summaries needed by the model-specific Stan input builders.
#'
#' @keywords internal
#' @noRd
.summarize_exposure <- function(exposure, cues, category, group, category_levels, group_levels, model) {
  n_categories <- length(category_levels)
  n_groups <- length(group_levels)
  n_cues <- length(cues)

  N_exposure <- matrix(0L, nrow = n_categories, ncol = n_groups)
  if (model == "NIX") {
    x_mean_exposure <- matrix(0, nrow = n_categories, ncol = n_groups)
    x_sd_exposure <- matrix(0, nrow = n_categories, ncol = n_groups)
  } else if (model == "MNIX") {
    x_mean_exposure <- array(0, dim = c(n_categories, n_groups, n_cues))
    x_ss_exposure <- array(0, dim = c(n_categories, n_groups, n_cues))
  } else if (model == "NIW") {
    x_mean_exposure <- array(0, dim = c(n_categories, n_groups, n_cues))
    x_ss_exposure <- array(0, dim = c(n_categories, n_groups, n_cues, n_cues))
  } else {
    stop("Unknown model type: ", model)
  }

  for (i in seq_along(category_levels)) {
    cat_level <- category_levels[i]
    for (j in seq_along(group_levels)) {
      group_level <- group_levels[j]
      subset <- exposure[exposure[[category]] == cat_level & exposure[[group]] == group_level, , drop = FALSE]
      n_obs <- nrow(subset)
      N_exposure[i, j] <- n_obs
      if (n_obs > 0) {
        cue_values <- as.matrix(subset[, cues, drop = FALSE])
        if (model == "NIX") {
          x_mean_exposure[i, j] <- mean(cue_values[, 1])
          x_sd_exposure[i, j] <- if (n_obs > 1) stats::sd(cue_values[, 1]) else 0
        } else if (model == "MNIX")  {
          x_mean_exposure[i, j, ] <- .colMeans(cue_values)
          if (n_obs > 1) {
            centered <- sweep(cue_values, 2, .colMeans(cue_values), "-")
            x_ss_exposure[i, j, ] <- colSums(centered^2)
          } else {
            x_ss_exposure[i, j, ] <- rep(0, n_cues)
          }
        } else if (model == "NIW") {
          x_mean_exposure[i, j, ] <- .colMeans(cue_values)
          if (n_obs > 1) {
            x_ss_exposure[i, j, , ] <- crossprod(cue_values)
          } else {
            # The Stan update of S handles n_obs == 1 as (different) special cases, so any positive-definite placeholder is harmless.
            x_ss_exposure[i, j, , ] <- diag(n_cues)
          }
        } 
      } else {
        if (model == "NIX") {
          x_mean_exposure[i, j] <- 0
          x_sd_exposure[i, j] <- 1
        } else if (model == "MNIX")  {
          x_mean_exposure[i, j, ] <- rep(0, n_cues)
          # For n_obs == 0, no update occurs in Stan, so a zero placeholder is harmless.
          x_ss_exposure[i, j, ] <- rep(0, n_cues)
        } else if (model == "NIW") {
          x_mean_exposure[i, j, ] <- rep(0, n_cues)
          # For n_obs == 0, no update occurs in Stan, so a positive-definite placeholder is harmless.
          x_ss_exposure[i, j, , ] <- diag(n_cues)
        } 
      }
    }
  }

  if (model == "NIX") {
    list(N_exposure = N_exposure, x_mean_exposure = x_mean_exposure, x_sd_exposure = x_sd_exposure)
  } else if (model == "MNIX") {
    list(N_exposure = N_exposure, x_mean_exposure = x_mean_exposure, x_ss_exposure = x_ss_exposure)
  } else if (model == "NIW") {
    list(N_exposure = N_exposure, x_mean_exposure = x_mean_exposure, x_ss_exposure = x_ss_exposure)
  } else {
    stop("Unknown model type: ", model)
  }
}

#' Summarize test data for Stanfit input construction.
#'
#' Internal helper that aggregates the test responses into the counts and cue
#' vectors expected by the Stan input builders.
#'
#' @keywords internal
#' @noRd
.summarize_test <- function(
  test,
  cues,
  response,
  group,
  category_levels,
  group_levels,
  model
) {
  n_cues <- length(cues)
  n_categories <- length(category_levels)
  n_groups <- length(group_levels)

  unique_cols <- unique(c(group, cues))
  unique_rows <- unique(test[, unique_cols, drop = FALSE])
  n_test <- nrow(unique_rows)
  if (n_test == 0) {
    x_test <- if (model == "NIX") {
      array(0, dim = c(0))
    } else {
      matrix(0, nrow = 0, ncol = n_cues)
    }
    y_test <- array(integer(0), dim = c(0))
    z_test_counts <- matrix(0L, nrow = 0, ncol = n_categories)
    return(
      list(
        x_test = x_test,
        y_test = y_test,
        z_test_counts = z_test_counts,
        N_test = 0
      )
    )
  }

  x_test <- if (model == "NIX") {
    array(0, dim = c(n_test))
  } else {
    matrix(0, nrow = n_test, ncol = n_cues)
  }
  y_test <- array(0L, dim = c(n_test))
  z_test_counts <- matrix(0L, nrow = n_test, ncol = n_categories)

  for (i in seq_len(n_test)) {
    row <- unique_rows[i, , drop = FALSE]
    matching_rows <- which(
      test[[group]] == row[[group]] &
        apply(test[, cues, drop = FALSE], 1, function(x) all(x == row[1, cues]))
    )
    if (length(matching_rows) == 0) {
      next
    }
    response_values <- factor(
      test[[response]][matching_rows],
      levels = category_levels
    )
    counts <- tabulate(as.integer(response_values), nbins = n_categories)
    z_test_counts[i, ] <- counts
    y_test[i] <- match(as.character(row[[group]]), group_levels)
    if (model == "NIX") {
      x_test[i] <- as.numeric(row[[cues[1]]])
    } else {
      x_test[i, ] <- as.numeric(as.matrix(row[, cues, drop = FALSE]))
    }
  }

  list(
    x_test = x_test,
    y_test = y_test,
    z_test_counts = z_test_counts,
    N_test = n_test
  )
}

