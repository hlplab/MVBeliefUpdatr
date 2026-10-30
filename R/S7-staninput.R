#' @include S7-class.R
#' @include S7-transform-information.R
#' @include S7-family.R
#' @include S7-stanfit-input.R
#' @include S7-prior.R

#' @title MVBeliefUpdatr Stan Input Classes and Constructors
#' @param values A named list containing the Stan input values.
#'   Default: \code{list()}.
#' @param family Character string indicating the model family name.
#'   Default: \code{character(0)}.
#'
#' @docType class
#' @name MVBU_Staninput
#' @rdname MVBU_Staninput
#' @export
MVBU_Staninput <- S7::new_class(
  "MVBU_Staninput",
  package = NULL,
  parent = MVBU_Object,
  properties = list(
    values = S7::class_list,
    family = S7::class_character
  ),
  constructor = function(values = list(), family = character(0)) {
    if (is.null(family)) {
      family <- character(0)
    }
    S7::new_object(
      MVBU_Object(),
      values = as.list(values),
      family = as.character(family)
    )
  },
  validator = function(self) {
    if (!is.list(self@values)) {
      return("`values` must be a list")
    }
    if (length(self@family) > 1) {
      return("`family` must have length at most 1")
    }
    NULL
  }
)

#' @docType class
#' @rdname MVBU_Staninput
#' @export
IdealAdaptorStaninput <- S7::new_class(
  "IdealAdaptorStaninput",
  package = NULL,
  parent = MVBU_Staninput,
  constructor = function(values = list(), family = character(0)) {
    if (is.null(family)) {
      family <- character(0)
    }
    S7::new_object(
      MVBU_Staninput(values = as.list(values), family = as.character(family))
    )
  },
  validator = function(self) {
    if (length(self@family) != 1L || !nzchar(self@family)) {
      return(paste0(
        "`family` must be a non-empty character string specifying a ",
        "registered Stan model family."
      ))
    }
    if (!.is_registered_stan_family(self@family)) {
      reg_stan_fams <- sort(names(.mvbu_family_registry$family_to_stanmodel))
      return(sprintf(
        paste0(
          "`family` '%s' is not a registered Stan model family. ",
          "Available Stan model families: %s."
        ),
        self@family,
        paste(reg_stan_fams, collapse = ", ")
      ))
    }
    NULL
  }
)

#' Ensure scalar control values are passed to Stan as the correct R object type
#'
#' @keywords internal
#' @noRd

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


#' Build the typed Stan input object for ideal-adaptor models
#'
#' Internal helper that summarizes exposure and test data from an
#' \code{\link{IdealAdaptorStanfitInput}} object, assembles the model-specific
#' list expected by Stan programs, and wraps it in an
#' \code{\link{IdealAdaptorStaninput}} object.
#'
#' @param stanfit_input An \code{\link{IdealAdaptorStanfitInput}} object.
#' @param stanmodel Character string naming the Stan model.
#' @param fixed_parameters_transformed Named list of transformed fixed
#'   parameters.
#' @param priors_transformed List of transformed priors.
#' @param control List of Stan control parameters (e.g.,
#'    \code{split_loglik_per_observation}).
#'
#' @return An \code{\link{IdealAdaptorStaninput}} object.
#' @keywords internal
#' @noRd
.new_staninput <- function(
  stanfit_input,
  stanmodel = "NIW_ideal_adaptor",
  fixed_parameters_transformed = list(),
  priors_transformed = list(),
  control = list()
) {
  .assert_true(
    S7::S7_inherits(stanfit_input, IdealAdaptorStanfitInput),
    msg = "stanfit_input must be an IdealAdaptorStanfitInput object."
  )

  data_transformed <- get_data(stanfit_input, space = "transformed")
  cues <- get_cue_labels(stanfit_input)
  categories <- get_category_labels(stanfit_input)
  groups <- get_group_labels(stanfit_input, include_prior = FALSE)
  transform <- stanfit_input@transform_information

  n_cues <- length(cues)
  n_categories <- length(categories)
  n_groups <- length(groups)

  .assert_any(
    stanmodel != "NIX_ideal_adaptor",
    n_cues == 1L,
    msg = "NIX_ideal_adaptor requires exactly one cue."
  )
  .assert_any(
    !startsWith(stanmodel, "MNIX_ideal_adaptor"),
    n_cues >= 2L,
    msg = sprintf("%s requires at least two cues.", stanmodel)
  )
  .assert_any(
    !grepl("MNIX_ideal_adaptor$", stanmodel),
    n_cues >= 2L,
    msg = sprintf("%s requires at least two cues.", stanmodel)
  )

  exposure <- data_transformed[
    data_transformed$Phase == "exposure", ,
    drop = FALSE
  ]
  test <- data_transformed[
    data_transformed$Phase == "test", ,
    drop = FALSE
  ]

  is_1d_nix <- (stanmodel == "NIX_ideal_adaptor")
  is_late_mnix <- startsWith(stanmodel, "MNIX_ideal_adaptor")

  exp_model <- if (is_1d_nix) "NIX" else if (is_late_mnix) "MNIX" else "NIW"
  test_model <- if (is_1d_nix) "NIX" else "NIW"

  orig_names <- get_original_variable_names(stanfit_input)
  has_unique_groups <- !is.null(orig_names$group_unique) &&
    !identical(orig_names$group_unique, orig_names$group) &&
    "group_unique" %in% names(data_transformed)

  if (has_unique_groups) {
    first_groups <- tapply(
      as.character(exposure$group),
      exposure$group_unique,
      `[`,
      1
    )
    exposure_subset <- exposure[
      exposure$group %in% first_groups, ,
      drop = FALSE
    ]
    exposure_summary <- .summarize_exposure(
      exposure = exposure_subset,
      cues = cues,
      categories = categories,
      groups = groups,
      category = "category",
      group = "group_unique",
      model = exp_model
    )
    test_summary <- .summarize_test(
      test = test,
      cues = cues,
      categories = categories,
      groups = groups,
      response = "response_category",
      group = "group_unique",
      model = test_model
    )
  } else {
    exposure_summary <- .summarize_exposure(
      exposure = exposure,
      cues = cues,
      categories = categories,
      groups = groups,
      category = "category",
      group = "group",
      model = exp_model
    )
    test_summary <- .summarize_test(
      test = test,
      cues = cues,
      categories = categories,
      groups = groups,
      response = "response_category",
      group = "group",
      model = test_model
    )
  }

  shift_param <- transform@shift
  inv_scale_param <- if (length(shift_param) == 1L) {
    1 / as.numeric(transform@SCALE)
  } else {
    solve(transform@SCALE)
  }

  split_loglik_per_observation <- if (!is.null(
    control$split_loglik_per_observation
  )) {
    as.numeric(control$split_loglik_per_observation)
  } else {
    0
  }

  struct <- fixed_parameters_transformed$structured
  if (is.null(struct)) {
    struct <- .resolve_fixed_parameters(
      fixed = MVBU_FixedParameter(),
      categories = categories,
      cues = cues,
      transform = transform
    )
  }

  stan_priors <- priors_transformed$stan
  if (is.null(stan_priors)) {
    stan_priors <- .resolve_priors(
      priors = MVBU_Prior(),
      stanmodel = stanmodel,
      cues = cues,
      category_labels = categories,
      transform = transform
    )$stan
  }

  # Build common data elements shared across NIW, NIX, and MNIX models
  base_values <- list(
    L = n_groups,
    M = n_categories,
    n_exposure = exposure_summary$n_exposure,
    mean_exposure = exposure_summary$mean_exposure,
    css_exposure = exposure_summary$css_exposure,
    n_test = test_summary$n_test,
    y_test = test_summary$y_test,
    z_test_counts = test_summary$z_test_counts,
    category_prior_known = as.integer(struct$category_prior$known),
    category_prior_data = matrix(
      as.numeric(struct$category_prior$data),
      nrow = struct$category_prior$known,
      ncol = n_categories
    ),
    lapse_rate_known = as.integer(struct$lapse_rate$known),
    lapse_rate_data = array(
      as.numeric(struct$lapse_rate$data),
      dim = struct$lapse_rate$known
    ),
    lapse_bias_known = as.integer(struct$lapse_bias$known),
    lapse_bias_data = matrix(
      as.numeric(struct$lapse_bias$data),
      nrow = struct$lapse_bias$known,
      ncol = n_categories
    ),
    M_m_0_known = as.integer(struct$m_0$M_known),
    idx_m_0_known = as.array(as.integer(struct$m_0$idx_known)),
    M_S_0_known = as.integer(struct$S_0$M_known),
    idx_S_0_known = as.array(as.integer(struct$S_0$idx_known)),
    mode_S_0_known = as.array(as.integer(struct$S_0$mode_known)),
    kappa_0_category_specific = as.integer(
      struct$kappa_0$is_category_specific
    ),
    M_kappa_0_known = as.integer(struct$kappa_0$M_known),
    idx_kappa_0_known = as.array(as.integer(struct$kappa_0$idx_known)),
    kappa_0_known_data = as.array(as.numeric(struct$kappa_0$data)),
    nu_0_category_specific = as.integer(
      struct$nu_0$is_category_specific
    ),
    M_nu_0_known = as.integer(struct$nu_0$M_known),
    idx_nu_0_known = as.array(as.integer(struct$nu_0$idx_known)),
    nu_0_known_data = as.array(as.numeric(struct$nu_0$data)),
    prior_alphas_of_category_prior = as.array(
      as.numeric(stan_priors$prior_alphas_of_category_prior)
    ),
    prior_alpha_of_lapse_rate = as.numeric(
      stan_priors$prior_alpha_of_lapse_rate
    ),
    prior_beta_of_lapse_rate = as.numeric(
      stan_priors$prior_beta_of_lapse_rate
    ),
    prior_alphas_of_lapse_bias = as.array(
      as.numeric(stan_priors$prior_alphas_of_lapse_bias)
    ),
    prior_loc_of_kappa = as.array(as.numeric(
      if (struct$kappa_0$is_category_specific == 1L) {
        stan_priors$prior_loc_of_kappa
      } else {
        stan_priors$prior_loc_of_kappa[1]
      }
    )),
    prior_scale_of_kappa = as.array(as.numeric(
      if (struct$kappa_0$is_category_specific == 1L) {
        stan_priors$prior_scale_of_kappa
      } else {
        stan_priors$prior_scale_of_kappa[1]
      }
    )),
    prior_loc_of_nu = as.array(as.numeric(
      if (struct$nu_0$is_category_specific == 1L) {
        stan_priors$prior_loc_of_nu
      } else {
        stan_priors$prior_loc_of_nu[1]
      }
    )),
    prior_scale_of_nu = as.array(as.numeric(
      if (struct$nu_0$is_category_specific == 1L) {
        stan_priors$prior_scale_of_nu
      } else {
        stan_priors$prior_scale_of_nu[1]
      }
    )),
    is_mu_0_known = as.array(as.integer(struct$mu_0$is_known)),
    is_Sigma_0_known = as.array(as.integer(struct$Sigma_0$is_known)),
    split_loglik_per_observation = as.integer(
      split_loglik_per_observation
    )
  )

  if (stanmodel == "NIX_ideal_adaptor") {
    m_0_data_1d <- if (struct$m_0$M_known > 0) {
      as.numeric(struct$m_0$data[, 1])
    } else {
      numeric(0)
    }
    S_0_data_1d <- if (struct$S_0$M_known > 0) {
      as.numeric(struct$S_0$data[, 1, 1])
    } else {
      numeric(0)
    }

    staninput_nix <- c(
      base_values,
      list(
        x_test = as.numeric(test_summary$x_test),
        m_0_known_data = as.array(m_0_data_1d),
        S_0_known_data = as.array(S_0_data_1d),
        mu_0_known_data = as.array(as.numeric(struct$mu_0$data[, 1])),
        Sigma_0_known_data = as.array(
          as.numeric(struct$Sigma_0$data[, 1, 1])
        ),
        prior_mu_of_m = as.array(as.numeric(stan_priors$prior_mu_of_m[, 1])),
        prior_loc_of_tau = as.numeric(stan_priors$prior_loc_of_tau[1]),
        prior_scale_of_tau = as.numeric(stan_priors$prior_scale_of_tau[1]),
        INV_SCALE = as.numeric(inv_scale_param),
        shift = as.numeric(shift_param[1])
      )
    )
    return(IdealAdaptorStaninput(values = staninput_nix, family = "NIX"))
  }

  if (stanmodel == "NIW_ideal_adaptor") {
    staninput_niw <- c(
      list(K = n_cues),
      base_values,
      list(
        x_test = test_summary$x_test,
        m_0_known_data = matrix(
          as.numeric(struct$m_0$data),
          nrow = struct$m_0$M_known,
          ncol = n_cues
        ),
        S_0_known_data = array(
          as.numeric(struct$S_0$data),
          dim = c(struct$S_0$M_known, n_cues, n_cues)
        ),
        mu_0_known_data = matrix(
          as.numeric(struct$mu_0$data),
          nrow = n_categories,
          ncol = n_cues
        ),
        Sigma_0_known_data = array(
          as.numeric(struct$Sigma_0$data),
          dim = c(n_categories, n_cues, n_cues)
        ),
        prior_mu_of_m = matrix(
          as.numeric(stan_priors$prior_mu_of_m),
          nrow = n_categories,
          ncol = n_cues
        ),
        prior_loc_of_tau = as.array(as.numeric(stan_priors$prior_loc_of_tau)),
        prior_scale_of_tau = as.array(
          as.numeric(stan_priors$prior_scale_of_tau)
        ),
        prior_eta_of_cor = as.numeric(stan_priors$prior_eta_of_cor),
        INV_SCALE = .make_stan_inv_scale(inv_scale_param, n_cues),
        shift = .make_stan_shift(shift_param, n_cues)
      )
    )
    return(IdealAdaptorStaninput(values = staninput_niw, family = "NIW"))
  }

  # MNIX models
  staninput_mnix <- c(
    list(K = n_cues),
    base_values,
    list(
      x_test = test_summary$x_test,
      m_0_known_data = matrix(
        as.numeric(struct$m_0$data),
        nrow = struct$m_0$M_known,
        ncol = n_cues
      ),
      S_0_known_data = array(
        as.numeric(struct$S_0$data),
        dim = c(struct$S_0$M_known, n_cues, n_cues)
      ),
      mu_0_known_data = matrix(
        as.numeric(struct$mu_0$data),
        nrow = n_categories,
        ncol = n_cues
      ),
      Sigma_0_known_data = array(
        as.numeric(struct$Sigma_0$data),
        dim = c(n_categories, n_cues, n_cues)
      ),
      prior_mu_of_m = matrix(
        as.numeric(stan_priors$prior_mu_of_m),
        nrow = n_categories,
        ncol = n_cues
      ),
      prior_loc_of_tau = as.array(as.numeric(stan_priors$prior_loc_of_tau)),
      prior_scale_of_tau = as.array(
        as.numeric(stan_priors$prior_scale_of_tau)
      ),
      INV_SCALE = .make_stan_inv_scale(inv_scale_param, n_cues),
      shift = .make_stan_shift(shift_param, n_cues)
    )
  )

  is_global_weighted <- stanmodel %in% c(
    "MNIX_ideal_adaptor_weighted_integrator",
    "weighted_integrator_MNIX_ideal_adaptor"
  )
  is_cat_weighted <- stanmodel %in% c(
    "MNIX_ideal_adaptor_category_weighted_integrator",
    "category_weighted_integrator_MNIX_ideal_adaptor"
  )
  is_grp_weighted <- stanmodel %in% c(
    "MNIX_ideal_adaptor_group_weighted_integrator",
    "group_weighted_integrator_MNIX_ideal_adaptor"
  )
  is_grpc_weighted <- stanmodel %in% c(
    "MNIX_ideal_adaptor_group_category_weighted_integrator",
    "group_category_weighted_integrator_MNIX_ideal_adaptor"
  )

  user_weights <- if (
    !is.null(struct$weights) && struct$weights$known == 1L
  ) {
    struct$weights$data
  } else {
    fixed_parameters_transformed$weights
  }

  if (is_global_weighted) {
    staninput_mnix$prior_alphas_of_weights <- as.array(
      as.numeric(rep(1, n_cues))
    )
    if (!is.null(user_weights)) {
      staninput_mnix$weights_known <- 1L
      staninput_mnix$weights_data <- matrix(
        as.numeric(user_weights),
        nrow = 1, ncol = n_cues
      )
      staninput_mnix$weights <- as.numeric(user_weights)
    } else {
      staninput_mnix$weights_known <- 0L
      staninput_mnix$weights_data <- matrix(0, nrow = 0, ncol = n_cues)
    }
  } else if (is_cat_weighted) {
    staninput_mnix$prior_alphas_of_weights <- as.array(
      as.numeric(rep(1, n_cues))
    )
    if (!is.null(user_weights)) {
      staninput_mnix$weights_known <- 1L
      w_mat <- if (is.matrix(user_weights)) {
        user_weights
      } else {
        do.call(rbind, lapply(user_weights, as.numeric))
      }
      staninput_mnix$weights_data <- matrix(
        as.numeric(w_mat),
        nrow = n_categories, ncol = n_cues
      )
      staninput_mnix$weights <- w_mat
    } else {
      staninput_mnix$weights_known <- 0L
      staninput_mnix$weights_data <- matrix(0, nrow = 0, ncol = n_cues)
    }
  } else if (is_grp_weighted) {
    staninput_mnix$prior_alphas_of_weights <- as.array(
      as.numeric(rep(1, n_cues))
    )
    if (!is.null(user_weights)) {
      staninput_mnix$weights_known <- 1L
      w_mat <- if (is.matrix(user_weights)) {
        user_weights
      } else {
        do.call(rbind, lapply(user_weights, as.numeric))
      }
      staninput_mnix$weights_data <- matrix(
        as.numeric(w_mat),
        nrow = n_groups, ncol = n_cues
      )
      staninput_mnix$weights <- w_mat
    } else {
      staninput_mnix$weights_known <- 0L
      staninput_mnix$weights_data <- matrix(0, nrow = 0, ncol = n_cues)
    }
  } else if (is_grpc_weighted) {
    staninput_mnix$prior_alphas_of_weights <- as.array(
      as.numeric(rep(1, n_cues))
    )
    if (!is.null(user_weights)) {
      staninput_mnix$weights_known <- 1L
      staninput_mnix$weights_data <- array(
        as.numeric(user_weights),
        dim = c(n_groups, n_categories, n_cues)
      )
      staninput_mnix$weights <- user_weights
    } else {
      staninput_mnix$weights_known <- 0L
      staninput_mnix$weights_data <- array(0, dim = c(0, 0, n_cues))
    }
  }

  fam <- .to_family(stanmodel)
  return(IdealAdaptorStaninput(values = staninput_mnix, family = fam))
}

#' Summarize exposure data for Stanfit input construction
#'
#' Internal helper that computes exposure counts and category- or group-specific
#' cue summaries needed by the model-specific Stan input builders.
#'
#' @keywords internal
#' @noRd
.summarize_exposure <- function(
  exposure, cues, categories, groups, category, group, model
) {
  n_categories <- length(categories)
  n_groups <- length(groups)
  n_cues <- length(cues)

  n_exposure <- matrix(0L, nrow = n_categories, ncol = n_groups)
  if (model == "NIX") {
    mean_exposure <- matrix(0, nrow = n_categories, ncol = n_groups)
    css_exposure <- matrix(0, nrow = n_categories, ncol = n_groups)
  } else if (model == "MNIX") {
    mean_exposure <- array(0, dim = c(n_categories, n_groups, n_cues))
    css_exposure <- array(0, dim = c(n_categories, n_groups, n_cues))
  } else if (model == "NIW") {
    mean_exposure <- array(0, dim = c(n_categories, n_groups, n_cues))
    css_exposure <- array(0, dim = c(n_categories, n_groups, n_cues, n_cues))
  } else {
    stop("Unknown model type: ", model)
  }

  for (i in seq_along(categories)) {
    cat_label <- categories[i]
    for (j in seq_along(groups)) {
      group_label <- groups[j]
      subset <- exposure[
        exposure[[category]] == cat_label & exposure[[group]] == group_label, ,
        drop = FALSE
      ]
      n_obs <- nrow(subset)
      n_exposure[i, j] <- n_obs
      if (n_obs > 0) {
        cue_values <- as.matrix(subset[, cues, drop = FALSE])
        if (model == "NIX") {
          mean_exposure[i, j] <- mean(cue_values[, 1])
          if (n_obs > 1) {
            css_exposure[i, j] <- sum(
              (cue_values[, 1] - mean_exposure[i, j])^2
            )
          }
        } else {
          col_means <- .colMeans(cue_values)
          mean_exposure[i, j, ] <- col_means
          if (n_obs > 1) {
            centered <- sweep(cue_values, 2, col_means, "-")
            if (model == "MNIX") {
              css_exposure[i, j, ] <- colSums(centered^2)
            } else {
              css_exposure[i, j, , ] <- crossprod(centered)
            }
          }
        }
      }
    }
  }

  list(
    n_exposure = n_exposure,
    mean_exposure = mean_exposure,
    css_exposure = css_exposure
  )
}

#' Summarize test data for Stanfit input construction
#'
#' Internal helper that aggregates the test responses into the counts and cue
#' vectors expected by the Stan input builders.
#'
#' @keywords internal
#' @noRd
.summarize_test <- function(
  test,
  cues,
  categories,
  groups,
  response,
  group,
  model
) {
  n_cues <- length(cues)
  n_categories <- length(categories)
  n_groups <- length(groups)

  unique_cols <- unique(c(group, cues))
  unique_rows <- unique(test[, unique_cols, drop = FALSE])
  n_test <- nrow(unique_rows)

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
      levels = categories
    )
    counts <- tabulate(as.integer(response_values), nbins = n_categories)
    z_test_counts[i, ] <- counts
    y_test[i] <- match(as.character(row[[group]]), groups)
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
    n_test = n_test
  )
}
