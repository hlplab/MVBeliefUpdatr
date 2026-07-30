#' Specify control parameters for new_ideal_adaptor_staninput()
#'
#' This function is used to specify control parameters for the `new_ideal_adaptor_staninput()` function, and to provide
#' reasonable defaults for any of the unspecified parameters.
#'
#' @param tau_scale A vector of scales for the Cauchy priors for each cue's standard deviations. Used in
#'   both the prior for m_0 and the prior for S_0. (default: vector of `5`s, assuming that the data are standardized).
#' @param L_omega_eta A vector of etas of the LKJ prior for the correlations of the covariance matrix of \code{mu_0}. Only used for
#'   models with multivariate categories (e.g., NIW_ideal_adaptor). (default: `1`,
#'   which corresponds to a uniform prior of correlation matrices)
#' @param split_loglik_per_observation Optionally, split the log likelihood per observation. This can be helpful of leave-one-out
#'   estimation in order to avoid high Pareto k, but it also makes the stored stanfit object much larger. (default: `0`)
#' @param transform_type An affine transformation that can be applied to the data. See `type` in \code{\link{get_affine_transform}}
#'    for details. (default: "standardize", which standardizes each cue separately)
#'
#' @return A list of control parameters that can be passed to \code{\link{new_ideal_adaptor_staninput}}.
#'
#' @export
control_staninput <- function(
    tau_scale = 5,
    L_omega_eta = 1,
    split_loglik_per_observation = 0,
    # THIS IS KEPT HERE JUST FOR NOW UNTIL I HAVE DETERMINED WHICH TRANSFORM IS BEST SUITED FOR FITTING.
    # (also remove the documentation for this once it's no longer needed AND remove transform_information
    # from the returned information AND change the documentation for the returned object above.)
    transform_type = c("identity", "center", "standardize", "PCA whiten", "ZCA whiten")[3]
) {
  list(
    tau_scale = tau_scale,
    L_omega_eta = L_omega_eta,
    split_loglik_per_observation = split_loglik_per_observation,
    transform_type = transform_type)
}


#' Construct a streamlined ideal-adaptor Stan input object.
#'
#' This function composes the exposure and test data in a compact form that can
#' be used as input to the ideal-adaptor Stan models. It accepts empty exposure
#' data as long as test data are present, and returns a compact list with the
#' transformed Stan input, the untransformed data object, and transformation
#' metadata.
#'
#' Exposure and test data are checked for the required cue, group, category,
#' and response columns. The exposure and test data are transformed with an affine
#' transformation to facilitate model-fitting (see \code{\link{control_staninput}}). 
#' This transformation is undone in in the stanfit object returned by the Stan model, 
#' so that the model parameters are returned in the original cue space.
#'
#' It is important to use `group` to identify individuals that had a specific
#' exposure (or no exposure at all) and specific test trials. You should not
#' use `group` to identify exposure conditions. Setting `group` to an exposure
#' condition results in exposure that concatenates the exposure observations
#' from all subjects in that condition. Typically, this is not what users
#' intend, as it models exposure to the combination of exposure tokens across
#' all subjects, rather than exposure to one set of those exposure tokens. To
#' achieve this intended outcome, use `group.unique` to identify groups with
#' identical exposure. This will correctly use only one unique instance of the
#' observations that any level of `group` receives during exposure.
#'
#' @param exposure A data frame or tibble with exposure data. Each row is assumed to contain one observation.
#' @param test A data frame or tibble with test data. Each row is assumed to contain one observation.
#' @param cues Names of columns with cue values. These columns must be present in both `exposure` and `test`.
#' @param category Name of the column in `exposure` that stores the category information. (default: "category")
#' @param response Name of the column in `test` that stores the response information. (default: "response")
#' @param group Name of the column in `exposure` and `test` that stores the grouping information. (default: "group")
#' @param group.unique Optional grouping column used to collapse repeated exposure conditions before constructing the Stan input.
#' @param fix_parameters Optional list of fixed model parameters. Currently recognized are the following parameters:  
#' @param fix_parameters Optional list of fixed model parameters. Currently recognized are the following parameters:  
#' `lapse_rate` Optional scalar lapse-rate value. If provided, must be between 0 and 1.
#' `mu_0` Optional prior mean(s) for the category means. The shape should match the model requirements.
#' `Sigma_0` Optional prior covariance matrix or matrices for the category covariance. The shape should match the model requirements.
#' @param control A list of control parameters for the constructor, typically produced by `control_staninput()`.
#' @param stanmodel Character string naming the Stan model. Must be one of `NIX_ideal_adaptor`, `NIW_ideal_adaptor`, or `MNIX_ideal_adaptor`.
#' @param verbose Should verbose output be provided? (default: `FALSE`)
#' @return A list with components `staninput`, `data`, and `transform_information`.
#' @export
new_ideal_adaptor_staninput <- function(
    exposure,
    test,
    cues,
    category = "category",
    response = "response",
    group = "group",
    group.unique = NULL,
    fix_parameters = NULL,
    control = control_staninput(),
    stanmodel = "NIW_ideal_adaptor",
    verbose = FALSE
) {
  if (!is.list(control)) {
    stop("control must be a list.")
  }

  expected_control <- c("tau_scale", "L_omega_eta", "split_loglik_per_observation", "transform_type")
  if (!all(expected_control %in% names(control))) {
    stop("control must contain tau_scale, L_omega_eta, split_loglik_per_observation, and transform_type.")
  }

  if (!is.character(cues) || length(cues) < 1) {
    stop("cues must be a non-empty character vector.")
  }
  cues <- unique(cues)

  if (!is.data.frame(exposure)) exposure <- as.data.frame(exposure)
  if (!is.data.frame(test)) test <- as.data.frame(test)

  if (!all(cues %in% names(exposure))) {
    stop("All cue columns must be present in exposure.")
  }
  if (!all(cues %in% names(test))) {
    stop("All cue columns must be present in test.")
  }

  if (!group %in% names(exposure) || !group %in% names(test)) {
    stop("group column must be present in both exposure and test.")
  }
  if (!is.null(category) && !category %in% names(exposure)) {
    stop("category column must be present in exposure.")
  }
  if (!is.null(response) && !response %in% names(test)) {
    stop("response column must be present in test.")
  }

  if (nrow(test) < 1) {
    stop("new_ideal_adaptor_staninput requires non-empty test data.")
  }

  if (is.null(fix_parameters)) {
    fix_parameters <- list()
  } else if (!is.list(fix_parameters)) {
    stop("fix_parameters must be a list.")
  }

  if (!is.null(fix_parameters$lapse_rate)) {
    lapse_rate <- fix_parameters$lapse_rate
    if (!is.numeric(lapse_rate) || length(lapse_rate) != 1 || is.na(lapse_rate) || lapse_rate < 0 || lapse_rate > 1) {
      stop("lapse_rate must be a numeric value between 0 and 1.")
    }
  } else {
    lapse_rate <- NULL
  }

  mu_0 <- fix_parameters$mu_0
  Sigma_0 <- fix_parameters$Sigma_0

  stanmodel <- match.arg(stanmodel, c("NIX_ideal_adaptor", "NIW_ideal_adaptor", "MNIX_ideal_adaptor"))

  if (stanmodel == "NIX_ideal_adaptor" && length(cues) != 1) {
    stop("NIX_ideal_adaptor requires exactly one cue.")
  }
  if (stanmodel == "MNIX_ideal_adaptor" && length(cues) < 2) {
    stop("MNIX_ideal_adaptor requires at least two cues.")
  }

  tau_scale <- control$tau_scale
  if (length(tau_scale) == 1) tau_scale <- rep(tau_scale, length(cues))
  tau_scale <- as.numeric(tau_scale)

  if (length(tau_scale) != length(cues)) {
    stop("tau_scale must have length 1 or length(cues).")
  }

  transform_type <- control$transform_type
  if (!is.character(transform_type) || length(transform_type) != 1) {
    stop("transform_type must be a single character value.")
  }
  if (!transform_type %in% c("identity", "center", "standardize", "PCA whiten", "ZCA whiten")) {
    stop("transform_type must be one of identity, center, standardize, PCA whiten, or ZCA whiten.")
  }

  exposure <- .prepare_staninput_frame(
    exposure,
    cues = cues,
    category = category,
    response = NULL,
    group = group,
    group.unique = group.unique,
    verbose = verbose
  )
  test <- .prepare_staninput_frame(
    test,
    cues = cues,
    category = NULL,
    response = response,
    group = group,
    group.unique = group.unique,
    verbose = verbose
  )

  if (!is.null(group.unique)) {
    if (!group.unique %in% names(exposure)) {
      stop("group.unique column must be present in exposure.")
    }
    exposure[[group.unique]] <- factor(exposure[[group.unique]])
    key <- interaction(exposure[[group.unique]], exposure[[category]], do.call("interaction", exposure[, cues, drop = FALSE]), drop = TRUE)
    keep_idx <- !duplicated(key)
    exposure <- exposure[keep_idx, , drop = FALSE]
  }

  if (!is.null(category) && !is.null(response)) {
    shared_levels <- union(levels(exposure[[category]]), levels(test[[response]]))
    exposure[[category]] <- factor(exposure[[category]], levels = shared_levels)
    test[[response]] <- factor(test[[response]], levels = shared_levels)
  }
  shared_group_levels <- union(levels(exposure[[group]]), levels(test[[group]]))
  exposure[[group]] <- factor(exposure[[group]], levels = shared_group_levels)
  test[[group]] <- factor(test[[group]], levels = shared_group_levels)

  # Use transformed and untransformed data copies.
  exposure_untransformed <- exposure
  test_untransformed <- test

  transform <- get_affine_transform(exposure, cues, transform_type)
  exposure <- transform$transform.function(exposure, return_type = "replace")
  test <- transform$transform.function(test, return_type = "replace")

  if (!is.null(mu_0)) {
    mu_0 <- .validate_and_transform_prior(mu_0, exposure[[category]], cues, transform, which = "mu_0")
  }
  if (!is.null(Sigma_0)) {
    Sigma_0 <- .validate_and_transform_prior(Sigma_0, exposure[[category]], cues, transform, which = "Sigma_0")
  }

  n_cues <- length(cues)
  n_categories <- nlevels(exposure[[category]])
  n_groups <- nlevels(exposure[[group]])
  category_levels <- levels(exposure[[category]])
  group_levels <- levels(exposure[[group]])

  if (is.null(category_levels) || length(category_levels) == 0) {
    category_levels <- levels(test[[response]])
  }
  if (is.null(group_levels) || length(group_levels) == 0) {
    group_levels <- levels(test[[group]])
  }

  if (stanmodel == "NIX_ideal_adaptor") {
    transformed_staninput <- .make_nix_staninput(
      exposure = exposure,
      test = test,
      cues = cues,
      category = category,
      response = response,
      group = group,
      category_levels = category_levels,
      group_levels = group_levels,
      tau_scale = tau_scale,
      L_omega_eta = control$L_omega_eta,
      split_loglik_per_observation = control$split_loglik_per_observation,
      lapse_rate = lapse_rate,
      mu_0 = mu_0,
      Sigma_0 = Sigma_0,
      transform = transform,
      n_cues = n_cues,
      n_categories = n_categories,
      n_groups = n_groups
    )
  } else if (stanmodel == "NIW_ideal_adaptor") {
    transformed_staninput <- .make_niw_staninput(
      exposure = exposure,
      test = test,
      cues = cues,
      category = category,
      response = response,
      group = group,
      category_levels = category_levels,
      group_levels = group_levels,
      tau_scale = tau_scale,
      L_omega_eta = control$L_omega_eta,
      split_loglik_per_observation = control$split_loglik_per_observation,
      lapse_rate = lapse_rate,
      mu_0 = mu_0,
      Sigma_0 = Sigma_0,
      transform = transform,
      n_cues = n_cues,
      n_categories = n_categories,
      n_groups = n_groups
    )
  } else {
    transformed_staninput <- .make_mnix_staninput(
      exposure = exposure,
      test = test,
      cues = cues,
      category = category,
      response = response,
      group = group,
      category_levels = category_levels,
      group_levels = group_levels,
      tau_scale = tau_scale,
      L_omega_eta = control$L_omega_eta,
      split_loglik_per_observation = control$split_loglik_per_observation,
      lapse_rate = lapse_rate,
      mu_0 = mu_0,
      Sigma_0 = Sigma_0,
      transform = transform,
      n_cues = n_cues,
      n_categories = n_categories,
      n_groups = n_groups
    )
  }

  transform_information <- .make_transform_information(transform)

  data <- .build_staninput_data(
    exposure = exposure_untransformed,
    test = test_untransformed,
    category = category,
    response = response,
    group = group,
    group.unique = group.unique,
    cues = cues
  )

  list(
    staninput = list(
      transformed = transformed_staninput,
      untransformed = transformed_staninput
    ),
    data = data,
    transform_information = transform_information
  )
}

#' Prepare a data frame for Stan input construction.
#'
#' Internal helper that retains the required columns, coerces cue values to
#' numeric data, removes rows with missing values, and converts the required
#' columns to factors before downstream processing.
#'
#' @param data Input data frame.
#' @param cues Character vector of cue columns.
#' @param category Optional category column.
#' @param response Optional response column.
#' @param group Grouping column.
#' @param group.unique Optional grouping column used to collapse repeated exposure conditions.
#' @param verbose Logical flag for verbose output.
#' @keywords internal
#' @noRd
.prepare_staninput_frame <- function(data, cues, category, response, group, group.unique = NULL, verbose = FALSE) {
  data <- as.data.frame(data)
  required_cols <- c(group)
  if (!is.null(category)) required_cols <- c(required_cols, category)
  if (!is.null(response)) required_cols <- c(required_cols, response)
  if (!is.null(group.unique)) required_cols <- c(required_cols, group.unique)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing columns in data: %s", paste(missing_cols, collapse = ", ")))
  }
  missing_cues <- setdiff(cues, names(data))
  if (length(missing_cues) > 0) {
    stop(sprintf("Missing cue columns in data: %s", paste(missing_cues, collapse = ", ")))
  }

  keep_cols <- unique(c(group, cues, required_cols))
  data <- data[, keep_cols, drop = FALSE]

  factor_cols <- required_cols
  preserved_levels <- lapply(data[, factor_cols, drop = FALSE], function(x) {
    if (is.factor(x)) levels(x) else NULL
  })
  data[factor_cols] <- lapply(data[factor_cols], function(x) {
    values <- as.character(x)
    values[is.na(values)] <- NA
    values
  })
  data[, cues] <- lapply(data[, cues, drop = FALSE], function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  })
  complete_rows <- stats::complete.cases(data[, c(cues, required_cols), drop = FALSE])
  data <- data[complete_rows, , drop = FALSE]

  for (col_name in factor_cols) {
    current_levels <- preserved_levels[[col_name]]
    if (!is.null(current_levels) && length(current_levels) > 0) {
      data[[col_name]] <- factor(data[[col_name]], levels = current_levels)
    } else {
      data[[col_name]] <- factor(data[[col_name]])
    }
  }

  data
}

#' Align factor levels across exposure and test data.
#'
#' Internal helper that ensures the grouping and category factors are aligned
#' with the levels present in the reference data.
#'
#' @param data Input data frame to align.
#' @param reference Reference data frame providing the target factor levels.
#' @param category Optional category column.
#' @param response Optional response column.
#' @param group Grouping column.
#' @keywords internal
#' @noRd
.align_staninput_factors <- function(data, reference, category, response, group) {
  if (!is.null(category) && !is.null(response)) {
    category_levels <- levels(reference[[category]])
    if (length(category_levels) == 0) {
      category_levels <- levels(reference[[response]])
    }
    if (length(category_levels) > 0) {
      data[[category]] <- factor(data[[category]], levels = category_levels)
      data[[response]] <- factor(data[[response]], levels = category_levels)
    }
  }
  group_levels <- levels(reference[[group]])
  if (length(group_levels) == 0) {
    group_levels <- levels(data[[group]])
  }
  if (length(group_levels) > 0) {
    data[[group]] <- factor(data[[group]], levels = group_levels)
  }
  data
}

#' Validate and transform prior information.
#'
#' Internal helper that checks the dimensionality of supplied priors and
#' applies the affine transform to the category means and covariance matrices.
#'
#' @param prior Prior values supplied by the user.
#' @param category_factor Factor of category labels.
#' @param cues Character vector of cue columns.
#' @param transform Affine transform object.
#' @param which Whether the prior is for `mu_0` or `Sigma_0`.
#' @keywords internal
#' @noRd
.validate_and_transform_prior <- function(prior, category_factor, cues, transform, which = c("mu_0", "Sigma_0")) {
  which <- match.arg(which)
  category_levels <- levels(category_factor)
  if (length(category_levels) == 1) {
    if (which == "mu_0") {
      if (!is.vector(prior)) {
        stop("mu_0 must be a vector when there is only one category.")
      }
      prior <- list(as.numeric(prior))
    } else if (which == "Sigma_0") {
      if (!is.array(prior) && !is.matrix(prior)) {
        stop("Sigma_0 must be a matrix when there is only one category.")
      }
      prior <- list(as.matrix(prior))
    }
  } else {
    if (!is.list(prior) || length(prior) != length(category_levels)) {
      stop(sprintf("%s must be a list with one entry per category.", which))
    }
  }

  if (which == "mu_0") {
    prior <- lapply(prior, function(x) transform_category_mean(as.numeric(x), transform))
  } else {
    prior <- lapply(prior, function(x) transform_category_cov(as.matrix(x), transform))
  }
  prior
}

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

#' Build the Stan input list for the NIX ideal-adaptor model.
#'
#' Internal helper that summarizes exposure and test data and assembles the
#' model-specific list expected by the NIX Stan program.
#'
#' @keywords internal
#' @noRd
.make_nix_staninput <- function(exposure, test, cues, category, response, group, category_levels, group_levels, tau_scale, L_omega_eta, split_loglik_per_observation, lapse_rate, mu_0, Sigma_0, transform, n_cues, n_categories, n_groups) {
  exposure_summary <- .summarize_exposure(exposure, cues, category, group, category_levels, group_levels, model = "NIX")
  test_summary <- .summarize_test(test, cues, response, group, category_levels, group_levels, model = "NIX")

  staninput <- list(
    K = n_cues,
    M = n_categories,
    L = n_groups,
    tau_scale = if (n_cues == 1) as.numeric(tau_scale[1]) else .make_stan_tauscale(tau_scale, n_cues),
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
    shift = if (n_cues == 1) as.numeric(transform$transform.parameters[["shift"]][1]) else .make_stan_shift(transform$transform.parameters[["shift"]], n_cues),
    INV_SCALE = if (n_cues == 1) as.numeric(transform$transform.parameters[["INV_SCALE"]]) else .make_stan_inv_scale(transform$transform.parameters[["INV_SCALE"]], n_cues)
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

  staninput
}

#' Build the Stan input list for the NIW ideal-adaptor model.
#'
#' Internal helper that summarizes exposure and test data and assembles the
#' model-specific list expected by the NIW Stan program.
#'
#' @keywords internal
#' @noRd
.make_niw_staninput <- function(exposure, test, cues, category, response, group, category_levels, group_levels, tau_scale, L_omega_eta, split_loglik_per_observation, lapse_rate, mu_0, Sigma_0, transform, n_cues, n_categories, n_groups) {
  exposure_summary <- .summarize_exposure(exposure, cues, category, group, category_levels, group_levels, model = "NIW")
  test_summary <- .summarize_test(test, cues, response, group, category_levels, group_levels, model = "NIW")

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
    INV_SCALE = .make_stan_inv_scale(transform$transform.parameters[["INV_SCALE"]], n_cues)
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

  staninput
}

#' Build the Stan input list for the MNIX ideal-adaptor model.
#'
#' Internal helper that summarizes exposure and test data and assembles the
#' model-specific list expected by the MNIX Stan program.
#'
#' @keywords internal
#' @noRd
.make_mnix_staninput <- function(exposure, test, cues, category, response, group, category_levels, group_levels, tau_scale, L_omega_eta, split_loglik_per_observation, lapse_rate, mu_0, Sigma_0, transform, n_cues, n_categories, n_groups) {
  exposure_summary <- .summarize_exposure(exposure, cues, category, group, category_levels, group_levels, model = "MNIX")
  test_summary <- .summarize_test(test, cues, response, group, category_levels, group_levels, model = "MNIX")

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
    INV_SCALE = .make_stan_inv_scale(transform$transform.parameters[["INV_SCALE"]], n_cues)
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

  staninput
}

#' Summarize exposure data for Stan input construction.
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
          x_mean_exposure[i, j, ] <- colMeans(cue_values)
          if (n_obs > 1) {
            centered <- sweep(cue_values, 2, colMeans(cue_values), "-")
            x_ss_exposure[i, j, ] <- colSums(centered^2)
          } else {
            x_ss_exposure[i, j, ] <- rep(0, n_cues)
          }
        } else if (model == "NIW") {
          x_mean_exposure[i, j, ] <- colMeans(cue_values)
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
  } else if (model == "NIW") {
    list(N_exposure = N_exposure, x_mean_exposure = x_mean_exposure, x_ss_exposure = x_ss_exposure)
  } else {
    list(N_exposure = N_exposure, x_mean_exposure = x_mean_exposure, x_ss_exposure = x_ss_exposure)
  }
}

#' Summarize test data for Stan input construction.
#'
#' Internal helper that aggregates the test responses into the counts and cue
#' vectors expected by the Stan input builders.
#'
#' @keywords internal
#' @noRd
.summarize_test <- function(test, cues, response, group, category_levels, group_levels, model) {
  n_cues <- length(cues)
  n_categories <- length(category_levels)
  n_groups <- length(group_levels)
  group_codes <- match(levels(test[[group]]), group_levels)
  if (length(group_codes) == 0) {
    group_codes <- integer(0)
  }

  unique_rows <- unique(test[, c(group, cues), drop = FALSE])
  n_test <- nrow(unique_rows)
  if (n_test == 0) {
    x_test <- if (model == "NIX") {
      array(0, dim = c(0))
    } else {
      matrix(0, nrow = 0, ncol = n_cues)
    }
    y_test <- array(integer(0), dim = c(0))
    z_test_counts <- matrix(0L, nrow = 0, ncol = n_categories)
    return(list(x_test = x_test, y_test = y_test, z_test_counts = z_test_counts, N_test = 0))
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
    matching_rows <- which(test[[group]] == row[[group]] & apply(test[, cues, drop = FALSE], 1, function(x) all(x == row[1, cues])))
    if (length(matching_rows) == 0) {
      next
    }
    response_values <- factor(test[[response]][matching_rows], levels = category_levels)
    counts <- tabulate(as.integer(response_values), nbins = n_categories)
    z_test_counts[i, ] <- counts
    y_test[i] <- match(as.character(row[[group]]), group_levels)
    if (model == "NIX") {
      x_test[i] <- as.numeric(row[[cues[1]]])
    } else {
      x_test[i, ] <- as.numeric(as.matrix(row[, cues, drop = FALSE]))
    }
  }

  list(x_test = x_test, y_test = y_test, z_test_counts = z_test_counts, N_test = n_test)
}

#' Package the transformation metadata.
#'
#' Internal helper that exposes the affine transform parameters and functions
#' in a simple list suitable for the returned Stan-input object.
#'
#' @keywords internal
#' @noRd
.make_transform_information <- function(transform) {
  list(
    transform.parameters = transform$transform.parameters,
    transform.function = transform$transform.function,
    untransform.function = transform$untransform.function
  )
}

#' Build the combined exposure and test data frame.
#'
#' Internal helper that stacks the processed exposure and test rows into a
#' single data frame annotated with a phase column.
#'
#' @keywords internal
#' @noRd
.build_staninput_data <- function(exposure, test, category, response, group, group.unique, cues) {
  keep_exposure <- c(group, category, cues)
  keep_test <- c(group, response, cues)
  if (!is.null(group.unique)) {
    keep_exposure <- c(group.unique, keep_exposure)
    keep_test <- c(group.unique, keep_test)
  }
  keep_exposure <- keep_exposure[keep_exposure %in% names(exposure)]
  keep_test <- keep_test[keep_test %in% names(test)]

  exposure_data <- exposure[, keep_exposure, drop = FALSE]
  test_data <- test[, keep_test, drop = FALSE]

  all_cols <- unique(c(names(exposure_data), names(test_data)))
  for (col_name in setdiff(all_cols, names(exposure_data))) {
    exposure_data[[col_name]] <- rep(NA, nrow(exposure_data))
  }
  for (col_name in setdiff(all_cols, names(test_data))) {
    test_data[[col_name]] <- rep(NA, nrow(test_data))
  }

  exposure_data <- exposure_data[, all_cols, drop = FALSE]
  test_data <- test_data[, all_cols, drop = FALSE]

  data <- rbind(exposure_data, test_data)
  data$Phase <- c(rep("exposure", nrow(exposure_data)), rep("test", nrow(test_data)))
  data$Phase <- factor(data$Phase, levels = c("exposure", "test"))
  data
}



#' Prepare long data from incremental exposure-test design for input to Stan
#'
#' Takes \code{data.frame} or \code{tibble} that contains the exposure and test data from an incremental
#' exposure-test design in long format, and prepares it for input to the \code{\link{ideal_adaptor_stanfit}}
#' Stan programs. This is done by pretending that each incremental test block (and its preceding exposure)
#' constitute a separate between-participant condition. Note that this does not capture the dependency
#' between test responses of participants in the same between-participant conditions, but such dependencies
#' are not modeled by current `MVBeliefUpdatr` Stan programs anyway (which do not include random effects by
#' participants).
#'
#' @param data Data frame or tibble to be sliced. Each row should be a single exposure or test observation.
#' @param group Character string indicating the name of the column that contains the information about
#'   the between-participant condition. (default: "Group")
#' @param phase Character string indicating the name of the column that contains the information about
#'   whether an observation is part of "exposure" or "test". This column must contain the values "exposure"
#'   and "test". Observation with other values will be ignored. (default: "Phase")
#' @param block Character string indicating the name of the column that contains the information about the
#'   incremental exposure and test blocks. Must be a factor with the levels indicating the order of the blocks.
#'   (default: "Block")
#' @param join_adjacent_test_blocks Logical indicating whether adjacent test blocks without intervening
#'   exposure blocks should be joined into a single test block. This will speed up \code{\link{fit_ideal_adaptor}}
#'   since there will be fewer conditions to iterate over but also means that the default plotting functions
#'   won't be able to plot the results of the different test blocks separately. (default: `FALSE`)
#' @param verbose Should verbose output be provided? (default: `FALSE`)
#'
#' @return A data frame or tibble in long format with a new column "ExposureGroup" that contains a unique
#'   label for each unique combination of `group` and `block`.
#'
#' @export
reshape_incremental_design_into_unique_exposure_test_combinations <- function(
    data,
    group = "Group",
    phase = "Phase",
    block = "Block",
    join_adjacent_test_blocks = FALSE,
    verbose = FALSE
) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  stopifnot(all(c(phase, group, block) %in% names(data)))
  stopifnot(all(c("exposure", "test") %in% unique(as.character(data[[phase]]))))
  if (!is.factor(data[[block]])) {
    data[[block]] <- factor(data[[block]])
  }

  exposure_blocks <- unique(as.character(data[[block]][data[[phase]] == "exposure"]))
  test_blocks <- unique(as.character(data[[block]][data[[phase]] == "test"]))
  if (any(exposure_blocks %in% test_blocks)) {
    stop2("The levels of the block variable in the exposure phase must not overlap with those in the test phase. Please check your data.")
  }

  if (verbose && length(setdiff(unique(as.character(data[[phase]])), c("exposure", "test"))) > 0) {
    message(
      paste("The following values in the", phase, "column are not recognized as exposure or test and thus removed:",
            paste(setdiff(unique(as.character(data[[phase]])), c("exposure", "test")), collapse = ", ")))
  }

  keep_rows <- data[[phase]] %in% c("exposure", "test")
  data <- data[keep_rows, , drop = FALSE]
  data[["..block_order"]] <- as.numeric(data[[block]])

  block_levels <- levels(data[[block]])
  phase_table <- unique(data[, c(phase, block, "..block_order"), drop = FALSE])
  phase_table <- phase_table[order(phase_table[["..block_order"]]), , drop = FALSE]
  phase_levels <- as.character(phase_table[[phase]])

  if (join_adjacent_test_blocks) {
    for (b in seq_len(length(block_levels) - 1)) {
      if (all(phase_levels[b:(b + 1)] == "test")) {
        if (verbose) {
          message("Joining adjacent test blocks ", block_levels[b], " and ", block_levels[b + 1], " into a single test block.")
        }

        block_levels[b] <- paste(block_levels[b], block_levels[b + 1], sep = "_")
        block_levels <- block_levels[-(b + 1)]

        data[["..block_order"]] <- ifelse(data[["..block_order"]] == b + 1, b, data[["..block_order"]])
        data[[block]] <- factor(
          ifelse(data[["..block_order"]] %in% c(b, b + 1), block_levels[b], as.character(data[[block]])),
          levels = block_levels
        )
      }
    }
  }

  if (verbose) {
    message("Inferred block order: ", paste(block_levels, collapse = ", "))
  }

  testblock_order <- sort(unique(as.numeric(data[["..block_order"]][data[[phase]] == "test"])))
  if (verbose) {
    message("Inferred test block order: ", paste(testblock_order, collapse = ", "))
  }

  rows <- list()
  group_levels <- unique(as.character(data[[group]]))
  idx <- 1L
  for (g in group_levels) {
    for (b in testblock_order) {
      subset <- data[
        data[[group]] == g & data[["..block_order"]] <= b & (data[["..block_order"]] == b | data[[phase]] != "test"),
        , drop = FALSE
      ]
      subset[["ExposureGroup"]] <- if (b == 1) "no exposure" else paste0("Group ", g, "_up to block ", block_levels[b])
      rows[[idx]] <- subset
      idx <- idx + 1L
    }
  }

  if (length(rows) == 0L) {
    df.new <- data.frame(ExposureGroup = character(0), stringsAsFactors = FALSE)
  } else {
    df.new <- do.call(rbind, rows)
  }

  keep_cols <- c("ExposureGroup", group, phase, block, setdiff(names(df.new), c("ExposureGroup", group, phase, block)))
  df.new[, keep_cols, drop = FALSE]
}
