#' @include asserts.R
#' @include S7-core-classes.R
#' @include S7-transform-information.R
#' @include S7-staninput.R
NULL

#' Specify control parameters for new_ideal_adaptor_stanfit_input()
#'
#' This function is used to specify control parameters for the `new_ideal_adaptor_stanfit_input()` function, and to provide
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
#' @return A list of control parameters that can be passed to \code{\link{new_ideal_adaptor_stanfit_input}}.
#'
#' @export
control_staninput <- function(
    tau_scale = 5,
    L_omega_eta = 1,
    split_loglik_per_observation = 0,
    transform_type = c("identity", "center", "standardize", "PCA whiten", "ZCA whiten")[3]
) {
  list(
    tau_scale = tau_scale,
    L_omega_eta = L_omega_eta,
    split_loglik_per_observation = split_loglik_per_observation,
    transform_type = transform_type
  )
}

#' An S7 class for ideal-adaptor fit-input objects.
#'
#' @name IdealAdaptorStanfitInput-class
#' @docType class
#' @export
IdealAdaptorStanfitInput <- S7::new_class(
  "IdealAdaptorStanfitInput",
  package = NULL,
  parent = MVBU_Object,
  properties = list(
    data = S7::new_S3_class("data.frame"),
    staninput = IdealAdaptorStaninput,
    transform_information = MVBU_TransformInformation
  ),
  constructor = function(
    data = data.frame(),
    staninput = NULL,
    transform_information = NULL
  ) {
    if (is.null(transform_information)) {
      transform_information <- MVBU_TransformInformation()
    }

    if (is.null(staninput)) {
      staninput <- IdealAdaptorStaninput(values = list())
    }

    S7::new_object(
      MVBU_Object(),
      data = as.data.frame(data),
      staninput = staninput,
      transform_information = transform_information
    )
  },
  validator = function(self) {
    if (!is.data.frame(self@data)) {
      return("`data` must be a data.frame")
    }
    if (!is.null(self@staninput) && !S7::S7_inherits(self@staninput, IdealAdaptorStaninput)) {
      return("`staninput` must be NULL or an S7 object inheriting from IdealAdaptorStaninput")
    }
    if (!is.null(self@transform_information) && !S7::S7_inherits(self@transform_information, MVBU_TransformInformation)) {
      return("`transform_information` must inherit from MVBU_TransformInformation")
    }
    NULL
  }
)

#' Construct a streamlined ideal-adaptor Stanfit input object.
#'
#' This function composes the exposure and test data in a compact form that can
#' be used as input to the ideal-adaptor Stan models. It accepts empty exposure
#' data as long as test data are present, and returns an S7 fit-input object with
#' the prepared data, typed Stan input, and transformation metadata.
#'
#' Exposure and test data are checked for the required cue, group, category,
#' and response columns. The exposure and test data are transformed to
#' facilitate model-fitting (see \code{\link{control_staninput}}). This
#' transformation is undone in in the stanfit object returned by the Stan model,
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
#' @return An object of class `IdealAdaptorStanfitInput` containing the prepared data, typed Stan input, and transform metadata.
#' @export
new_ideal_adaptor_stanfit_input <- function(
    exposure,
    test,
    cues,
    category = "category",
    response = "response",
    group = "group",
    group.unique = NULL,
    fixed_parameters = NULL,
    control = control_staninput(),
    stanmodel = "NIW_ideal_adaptor",
    verbose = FALSE
) {
  .assert_list(control)

  expected_control <- c("tau_scale", "L_omega_eta", "split_loglik_per_observation", "transform_type")
  .assert_all(
    expected_control %in% names(control),
    msg = "control must contain tau_scale, L_omega_eta, split_loglik_per_observation, and transform_type."
  )

  .assert_non_NA_character(cues)
  cues <- unique(cues)

  if (!is.data.frame(exposure)) exposure <- as.data.frame(exposure)
  if (!is.data.frame(test)) test <- as.data.frame(test)

  .assert_data_contains_cols(exposure, c(cues, group), msg = "exposure data must contain all cues and the group column.")
  .assert_data_contains_cols(test, c(cues, group), msg = "test data must contain all cues and the group column.")

  .assert_true(
    is.null(category) || category %in% names(exposure),
    msg = "category column must be present in exposure."
  )
  .assert_true(
    is.null(response) || response %in% names(test),
    msg = "response column must be present in test."
  )

  .assert_true(nrow(test) >= 1, msg = "new_ideal_adaptor_stanfit_input requires non-empty test data.")

  if (is.null(fixed_parameters)) {
    fixed_parameters <- list()
  } else {
    .assert_list(fixed_parameters)
  }

  if (!is.null(fixed_parameters$lapse_rate)) {
    lapse_rate <- fixed_parameters$lapse_rate
    .assert_all(
      is.numeric(lapse_rate),
      length(lapse_rate) == 1L,
      !is.na(lapse_rate),
      lapse_rate >= 0,
      lapse_rate <= 1,
      msg = "lapse_rate must be a numeric value between 0 and 1."
    )
  } else {
    lapse_rate <- NULL
  }

  mu_0 <- fixed_parameters$mu_0
  Sigma_0 <- fixed_parameters$Sigma_0

  stanmodel <- match.arg(stanmodel, c("NIX_ideal_adaptor", "NIW_ideal_adaptor", "MNIX_ideal_adaptor"))

  n_cues <- length(cues)
  .assert_any(
    stanmodel != "NIX_ideal_adaptor",
    n_cues == 1,
    msg = "NIX_ideal_adaptor requires exactly one cue."
  )
  .assert_any(
    stanmodel != "MNIX_ideal_adaptor",
    n_cues >= 2,
    msg = "MNIX_ideal_adaptor requires at least two cues."
  )

  tau_scale <- control$tau_scale
  if (length(tau_scale) == 1) tau_scale <- rep(tau_scale, length(cues))
  tau_scale <- as.numeric(tau_scale)

  .assert_true(length(tau_scale) == length(cues), msg = "tau_scale must have length 1 or length(cues).")

  transform_type <- control$transform_type
  .assert_all(
    is.character(transform_type),
    length(transform_type) == 1L,
    msg = "transform_type must be a single character value."
  )
  .assert_true(transform_type %in% c("identity", "center", "standardize", "PCA whiten", "ZCA whiten"), msg = "transform_type must be one of identity, center, standardize, PCA whiten, or ZCA whiten.")

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
    .assert_data_contains_cols(exposure, group.unique)
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

  exposure_untransformed <- exposure
  test_untransformed <- test

  transform <- get_affine_transform(exposure, cues, transform_type)
  exposure <- transform$transform.function(exposure, return_type = "replace")
  test <- transform$transform.function(test, return_type = "replace")

  if (!is.null(mu_0)) {
    mu_0 <- .validate_and_transform_prior_likelihood(mu_0, exposure[[category]], n_cues = n_cues, transform, which = "mu_0")
  }
  if (!is.null(Sigma_0)) {
    Sigma_0 <- .validate_and_transform_prior_likelihood(Sigma_0, exposure[[category]], n_cues = n_cues, transform, which = "Sigma_0")
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
    staninput <- new_nix_staninput(
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
  } else if (stanmodel == "MNIX_ideal_adaptor") {
    staninput <- new_mnix_staninput(
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
    staninput <- new_niw_staninput(
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

  transform_information <- new_transform_information(transform)

  data <- .build_stanfit_input_data(
    exposure = exposure_untransformed,
    test = test_untransformed,
    category = category,
    response = response,
    group = group,
    group.unique = group.unique,
    cues = cues
  )

  attr(data, "category") <- category
  attr(data, "group") <- group
  attr(data, "response") <- response
  attr(data, "group.unique") <- group.unique
  attr(data, "cues") <- cues

  IdealAdaptorStanfitInput(
    data = data,
    staninput = staninput,
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
  .assert_non_NA_character(cues)
  required_cols <- c(group)
  if (!is.null(category)) required_cols <- c(required_cols, category)
  if (!is.null(response)) required_cols <- c(required_cols, response)
  if (!is.null(group.unique)) required_cols <- c(required_cols, group.unique)
  required_cols <- unique(required_cols)
  missing_cols <- setdiff(required_cols, names(data))
  .assert_true(length(missing_cols) == 0, msg = sprintf("Missing columns in data: %s", paste(missing_cols, collapse = ", ")))
  missing_cues <- setdiff(cues, names(data))
  .assert_true(length(missing_cues) == 0, msg = sprintf("Missing cue columns in data: %s", paste(missing_cues, collapse = ", ")))

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

#' Validate and transform prior likelihood information.
#'
#' Internal helper that checks the dimensionality of supplied priors and
#' applies the affine transform to the category means and covariance matrices.
#'
#' @param prior Prior likelihood values supplied by the user.
#' @param category_var Factor of category labels.
#' @param n_cues Number of cue dimensions expected for the supplied prior values.
#' @param transform Affine transform object.
#' @param which Whether the prior is for `mu_0` or `Sigma_0`.
#' @keywords internal
#' @noRd
.validate_and_transform_prior_likelihood <- function(prior, category_var, n_cues, transform, which = c("mu_0", "Sigma_0")) {
  .assert_numeric_scalar(n_cues)
  n_cues <- as.integer(n_cues)

  which <- match.arg(which)
  category_levels <- levels(category_var)

  if (length(category_levels) == 1) {
    if (which == "mu_0") {
      .assert_true(is.vector(prior), msg = "mu_0 must be a vector when there is only one category.")
      prior <- list(as.numeric(prior))
    } else if (which == "Sigma_0") {
      .assert_true(is.array(prior) || is.matrix(prior), msg = "Sigma_0 must be a matrix when there is only one category.")
      prior <- list(as.matrix(prior))
    }
  } else {
    .assert_list(prior)

    .assert_true(length(prior) == length(category_levels), msg = sprintf("%s must be a named list with names matching the categories present in the exposure and test data.", which))

    prior_names <- names(prior)
    .assert_true(!is.null(prior_names) && setequal(prior_names, category_levels), msg = sprintf("%s must be a named list with names matching the categories present in the exposure and test data.", which))
  }

  if (which == "mu_0") {
    .assert_true(
      all(vapply(prior, function(x) is.numeric(x) && length(as.numeric(x)) == n_cues, logical(1))),
      msg = sprintf("mu_0 must be a named list with each entry a numeric vector of length %s", n_cues)
    )
    prior <- lapply(prior, function(x) transform_category_mean(as.numeric(x), transform))
  } else {
    .assert_true(
      all(vapply(prior, function(x) {
        x <- as.matrix(x)
        is.numeric(x) && nrow(x) == n_cues && ncol(x) == n_cues
      }, logical(1))),
      msg = sprintf("Sigma_0 must be a named list with each entry a square matrix of dimension %s x %s", n_cues, n_cues)
    )
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

#' Build the combined exposure and test data frame.
#'
#' Internal helper that stacks the processed exposure and test rows into a
#' single data frame annotated with a phase column.
#'
#' @keywords internal
#' @noRd
.build_stanfit_input_data <- function(exposure, test, category, response, group, group.unique, cues) {
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
