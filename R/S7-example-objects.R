#' @include S7-make-objects.R
#' @include S7-stanfit-input.R
#' @include S7-stanfit-fitting.R
NULL

#' Example S7 objects built from the bundled Chodroff-Wilson data
#'
#' Helper functions to construct example S7 representations, templates,
#' and cognitive models from the bundled \code{mixer6} dataset.
#'
#' The parameter \code{n_cues} specifies the number of acoustic cues to
#' include (up to a maximum of 3 cues):
#' \itemize{
#'   \item \code{n_cues = 1}: Voice onset time (\code{"VOT"}). Required for
#'     univariate models (UVG and NIX).
#'   \item \code{n_cues = 2}: Voice onset time and fundamental frequency
#'     (\code{"VOT"} and \code{"f0_semitones"}).
#'   \item \code{n_cues = 3} (maximum): Voice onset time, fundamental
#'     frequency, and vowel duration (\code{"VOT"}, \code{"f0_semitones"},
#'     and \code{"vowel_duration"}).
#' }
#'
#' @param n_cues Integer number of cues: 1 (\code{"VOT"}), 2 (\code{"VOT"},
#'   \code{"f0_semitones"}), or 3 (\code{"VOT"}, \code{"f0_semitones"},
#'   \code{"vowel_duration"}; maximum is 3). Default: 1.
#' @param category Character string giving the category label.
#' @param categories Character vector of category labels.
#' @param kappa Strength of belief (pseudocount) about the category mean
#'   \eqn{\mu}. Default: 1.
#' @param nu Degrees of freedom (pseudocount) for covariance/variance beliefs
#'   \eqn{\Sigma}. Default: \code{n_cues + 2}.
#' @param type Character string indicating the model types. Valid choices are:
#'   \itemize{
#'     \item \code{"UVG"}: Univariate Gaussian.
#'     \item \code{"NIX"}: Normal-Inverse-chi^2.
#'     \item \code{"MUVG"}: Multivariate Univariate Gaussian.
#'     \item \code{"MNIX"}: Multivariate Normal-Inverse-chi^2.
#'     \item \code{"MVG"}: Multivariate Gaussian.
#'     \item \code{"NIW"}: Normal-Inverse-Wishart.
#'     \item \code{"EXEMPLAR"}: Exemplar.
#'   }
#' @param ... Additional arguments passed to underlying constructors.
#' @name example-s7-objects
NULL

.example_data <- function(n_cues, categories = NULL) {
  .assert_true(
    .is_scalar_count(n_cues),
    msg = "n_cues must be a positive whole number."
  )
  .assert_true(n_cues %in% 1:3, msg = "n_cues must be one of 1, 2, or 3.")
  d_env <- new.env(parent = emptyenv())
  utils::data("mixer6", package = "MVBeliefUpdatr", envir = d_env)
  mixer6 <- d_env$mixer6
  if (is.null(categories)) categories <- levels(mixer6$stop)
  cue_cols <- c("vot", "f0_semitones", "vowel_duration")[seq_len(n_cues)]
  out <- mixer6[mixer6$stop %in% categories, c("stop", cue_cols)]
  out <- out[stats::complete.cases(out), ]
  colnames(out)[colnames(out) == "stop"] <- "category"
  colnames(out)[colnames(out) == "vot"] <- "VOT"
  as.data.frame(out)
}

#' @rdname example-s7-objects
#' @export
example_uvg_category_representation <- function(
  n_cues = 1,
  category = "/b/"
) {
  .assert_true(n_cues == 1, msg = "UVG examples require n_cues = 1.")
  new_uvg_category_representation_from_data(
    .example_data(n_cues, category),
    cues = "VOT"
  )
}

#' @rdname example-s7-objects
#' @export
example_nix_category_representation <- function(
  n_cues = 1,
  category = "/b/",
  kappa = 1,
  nu = n_cues + 2
) {
  .assert_true(n_cues == 1, msg = "NIX examples require n_cues = 1.")
  new_nix_category_representation_from_data(
    .example_data(n_cues, category),
    cues = "VOT",
    kappa = kappa,
    nu = nu
  )
}

#' @rdname example-s7-objects
#' @export
example_muvg_category_representation <- function(
  n_cues = 2,
  category = "/b/"
) {
  .assert_true(n_cues >= 2, msg = "MUVG examples require n_cues >= 2.")
  new_muvg_category_representation_from_data(
    .example_data(n_cues, category),
    cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)]
  )
}

#' @rdname example-s7-objects
#' @export
example_mnix_category_representation <- function(
  n_cues = 2,
  category = "/b/",
  kappa = 1,
  nu = n_cues + 2
) {
  .assert_true(n_cues >= 2, msg = "MNIX examples require n_cues >= 2.")
  new_mnix_category_representation_from_data(
    .example_data(n_cues, category),
    cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)],
    kappa = kappa,
    nu = nu
  )
}

#' @rdname example-s7-objects
#' @export
example_mvg_category_representation <- function(
  n_cues = 1,
  category = "/b/"
) {
  new_mvg_category_representation_from_data(
    .example_data(n_cues, category),
    cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)]
  )
}

#' @rdname example-s7-objects
#' @export
example_niw_category_representation <- function(
  n_cues = 1,
  category = "/b/",
  kappa = 1,
  nu = n_cues + 2
) {
  new_niw_category_representation_from_data(
    .example_data(n_cues, category),
    cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)],
    kappa = kappa,
    nu = nu
  )
}

#' @rdname example-s7-objects
#' @export
example_exemplar_category_representation <- function(
  n_cues = 1,
  category = "/b/"
) {
  new_exemplar_category_representation_from_data(
    .example_data(n_cues, category),
    cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)]
  )
}

#' @rdname example-s7-objects
#' @export
example_category_representation <- function(
  type,
  n_cues = NULL,
  category = "/b/",
  ...
) {
  type <- toupper(type)
  if (is.null(n_cues)) {
    n_cues <- if (type %in% c("MUVG", "MNIX")) 2 else 1
  }
  constructor <- switch(type,
    UVG = example_uvg_category_representation,
    NIX = example_nix_category_representation,
    MUVG = example_muvg_category_representation,
    MNIX = example_mnix_category_representation,
    MVG = example_mvg_category_representation,
    NIW = example_niw_category_representation,
    EXEMPLAR = example_exemplar_category_representation,
    .stop("type must be one of UVG, NIX, MUVG, MNIX, MVG, NIW, or EXEMPLAR.")
  )
  constructor(n_cues = n_cues, category = category, ...)
}

#' @rdname example-s7-objects
#' @export
example_uvg_category_representation_template <- function(
  n_cues = 1,
  categories = NULL
) {
  .assert_true(n_cues == 1, msg = "UVG examples require n_cues = 1.")
  new_uvg_category_representation_template_from_data(
    .example_data(n_cues, categories),
    cues = "VOT"
  )
}

#' @rdname example-s7-objects
#' @export
example_nix_category_representation_template <- function(
  n_cues = 1,
  categories = NULL,
  kappa = 1,
  nu = n_cues + 2
) {
  .assert_true(n_cues == 1, msg = "NIX examples require n_cues = 1.")
  new_nix_category_representation_template_from_data(
    .example_data(n_cues, categories),
    cues = "VOT",
    kappa = kappa,
    nu = nu
  )
}

#' @rdname example-s7-objects
#' @export
example_muvg_category_representation_template <- function(
  n_cues = 2,
  categories = NULL
) {
  .assert_true(n_cues >= 2, msg = "MUVG examples require n_cues >= 2.")
  new_muvg_category_representation_template_from_data(
    .example_data(n_cues, categories),
    cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)]
  )
}

#' @rdname example-s7-objects
#' @export
example_mnix_category_representation_template <- function(
  n_cues = 2,
  categories = NULL,
  kappa = 1,
  nu = n_cues + 2
) {
  .assert_true(n_cues >= 2, msg = "MNIX examples require n_cues >= 2.")
  new_mnix_category_representation_template_from_data(
    .example_data(n_cues, categories),
    cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)],
    kappa = kappa,
    nu = nu
  )
}

#' @rdname example-s7-objects
#' @export
example_mvg_category_representation_template <- function(
  n_cues = 1,
  categories = NULL
) {
  new_mvg_category_representation_template_from_data(
    .example_data(n_cues, categories),
    cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)]
  )
}

#' @rdname example-s7-objects
#' @export
example_niw_category_representation_template <- function(
  n_cues = 1,
  categories = NULL,
  kappa = 1,
  nu = n_cues + 2
) {
  new_niw_category_representation_template_from_data(
    .example_data(n_cues, categories),
    cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)],
    kappa = kappa,
    nu = nu
  )
}

#' @rdname example-s7-objects
#' @export
example_exemplar_category_representation_template <- function(
  n_cues = 1,
  categories = NULL
) {
  new_exemplar_category_representation_template_from_data(
    .example_data(n_cues, categories),
    cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)]
  )
}

#' @rdname example-s7-objects
#' @export
example_category_representation_template <- function(
  type,
  n_cues = NULL,
  categories = NULL,
  ...
) {
  type <- toupper(type)
  if (is.null(n_cues)) {
    n_cues <- if (type %in% c("MUVG", "MNIX")) 2 else 1
  }
  constructor <- switch(type,
    UVG = example_uvg_category_representation_template,
    NIX = example_nix_category_representation_template,
    MUVG = example_muvg_category_representation_template,
    MNIX = example_mnix_category_representation_template,
    MVG = example_mvg_category_representation_template,
    NIW = example_niw_category_representation_template,
    EXEMPLAR = example_exemplar_category_representation_template,
    .stop("type must be one of UVG, NIX, MUVG, MNIX, MVG, NIW, or EXEMPLAR.")
  )
  constructor(n_cues = n_cues, categories = categories, ...)
}

#' @rdname example-s7-objects
#' @export
example_uvg_ideal_observer <- function(
  n_cues = 1,
  categories = NULL,
  ...
) {
  .assert_true(n_cues == 1, msg = "UVG examples require n_cues = 1.")
  new_uvg_ideal_observer_from_data(
    .example_data(n_cues, categories),
    cues = "VOT",
    ...
  )
}

#' @rdname example-s7-objects
#' @export
example_nix_ideal_adaptor <- function(
  n_cues = 1,
  categories = NULL,
  kappa = 1,
  nu = n_cues + 2,
  ...
) {
  .assert_true(n_cues == 1, msg = "NIX examples require n_cues = 1.")
  new_nix_ideal_adaptor_from_data(
    .example_data(n_cues, categories),
    cues = "VOT",
    kappa = kappa,
    nu = nu,
    ...
  )
}

#' @rdname example-s7-objects
#' @export
example_muvg_ideal_observer <- function(
  n_cues = 2,
  categories = NULL,
  ...
) {
  .assert_true(n_cues >= 2, msg = "MUVG examples require n_cues >= 2.")
  new_muvg_ideal_observer_from_data(
    .example_data(n_cues, categories),
    cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)],
    ...
  )
}

#' @rdname example-s7-objects
#' @export
example_mnix_ideal_adaptor <- function(
  n_cues = 2,
  categories = NULL,
  kappa = 1,
  nu = n_cues + 2,
  ...
) {
  .assert_true(n_cues >= 2, msg = "MNIX examples require n_cues >= 2.")
  new_mnix_ideal_adaptor_from_data(
    .example_data(n_cues, categories),
    cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)],
    kappa = kappa,
    nu = nu,
    ...
  )
}

#' @rdname example-s7-objects
#' @export
example_mvg_ideal_observer <- function(
  n_cues = 1,
  categories = NULL,
  ...
) {
  new_mvg_ideal_observer_from_data(
    .example_data(n_cues, categories),
    cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)],
    ...
  )
}

#' @rdname example-s7-objects
#' @export
example_niw_ideal_adaptor <- function(
  n_cues = 1,
  categories = NULL,
  kappa = 1,
  nu = n_cues + 2,
  ...
) {
  new_niw_ideal_adaptor_from_data(
    .example_data(n_cues, categories),
    cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)],
    kappa = kappa,
    nu = nu,
    ...
  )
}

#' @rdname example-s7-objects
#' @export
example_exemplar_model <- function(
  n_cues = 1,
  categories = NULL,
  ...
) {
  new_exemplar_model_from_data(
    .example_data(n_cues, categories),
    cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)],
    ...
  )
}

#' @rdname example-s7-objects
#' @export
example_model <- function(type, n_cues = NULL, categories = NULL, ...) {
  type <- toupper(type)
  if (is.null(n_cues)) {
    n_cues <- if (type %in% c("MUVG", "MNIX")) 2 else 1
  }
  constructor <- switch(type,
    UVG = example_uvg_ideal_observer,
    NIX = example_nix_ideal_adaptor,
    MUVG = example_muvg_ideal_observer,
    MNIX = example_mnix_ideal_adaptor,
    MVG = example_mvg_ideal_observer,
    NIW = example_niw_ideal_adaptor,
    EXEMPLAR = example_exemplar_model,
    .stop("type must be one of UVG, NIX, MUVG, MNIX, MVG, NIW, or EXEMPLAR.")
  )
  constructor(n_cues = n_cues, categories = categories, ...)
}

#' Example exposure and test experimental data
#'
#' Generates simulated exposure-test experiment datasets for phonetic
#' categorization and belief updating models. Shifts between exposure
#' conditions hold the mean distance between categories constant, while
#' variances remain unchanged. Shift magnitudes are approximately 50\% of the
#' distance between prior category means. A \code{"no_exposure"} baseline
#' condition tests the unshifted prior. Test locations form a regular grid
#' spanning the category continuum and the adaptation shift (forming a
#' parallelogram in 2D space), with simulated categorization responses.
#'
#' @param model_family Character string indicating the model family to
#'   simulate: \code{"NIW"} (default), \code{"NIX"}, or \code{"MNIX"}.
#' @param n_cues Integer number of acoustic cues (1 for NIX, 2-3 for MNIX,
#'   1 to 3 for NIW). Default: 1 (or 2 for MNIX).
#' @param categories Character vector of length 2 naming the categories.
#'   Default: \code{c("/d/", "/t/")}.
#' @param n_exposure_per_category Integer number of exposure observations
#'   per category for each shifted condition. Default: 20L.
#' @param n_test_steps_category Integer number of test continuum steps between
#'   categories. Default: 5L.
#' @param n_test_steps_adaptation Integer number of test continuum steps
#'   between prior and posterior means. Default: 3L.
#' @param n_test_tokens_per_location Integer number of test trials evaluated
#'   at each test location. Default: 40L.
#' @param seed Optional integer random seed for reproducibility. Default: 42L.
#'
#' @return A list with two tibbles:
#'   \describe{
#'     \item{exposure}{Exposure trials containing \code{Phase},
#'       \code{Condition}, \code{group}, \code{Subject}, cue columns, and
#'       \code{category}.}
#'     \item{test}{Test trials containing \code{Phase}, \code{Condition},
#'       \code{group}, \code{Subject}, cue columns, and
#'       \code{response_category}.}
#'   }
#'
#' @seealso \code{\link{example_ideal_adaptor_stanfit_input}},
#'   \code{\link{example_ideal_adaptor_staninput}},
#'   \code{\link{example_model}}
#' @export
example_exposure_test_data <- function(
  model_family = c("NIW", "NIX", "MNIX"),
  n_cues = NULL,
  categories = c("/d/", "/t/"),
  n_exposure_per_category = 20L,
  n_test_steps_category = 5L,
  n_test_steps_adaptation = 3L,
  n_test_tokens_per_location = 40L,
  seed = 42L
) {
  model_family <- match.arg(model_family)
  if (is.null(n_cues)) {
    n_cues <- if (model_family == "MNIX") 2L else 1L
  }
  .assert_true(
    model_family %in% c("NIW", "NIX", "MNIX"),
    msg = "model_family must be NIW, NIX, or MNIX."
  )
  .assert_true(
    length(categories) == 2L,
    msg = "categories must contain exactly 2 category labels."
  )
  .assert_between(n_cues, 1, 3)
  if (model_family == "NIX") {
    .assert_true(n_cues == 1L, msg = "NIX requires n_cues = 1.")
  } else if (model_family == "MNIX") {
    .assert_true(n_cues %in% 2:3, msg = "MNIX requires n_cues >= 2.")
  }
  .assert_true(
    n_exposure_per_category >= 0L,
    msg = "n_exposure_per_category must be at least 0."
  )
  .assert_true(
    n_test_steps_category >= 2L,
    msg = "n_test_steps_category must be at least 2."
  )
  .assert_true(
    n_test_steps_adaptation >= 2L,
    msg = "n_test_steps_adaptation must be at least 2."
  )
  .assert_true(
    n_test_tokens_per_location >= 1L,
    msg = "n_test_tokens_per_location must be at least 1."
  )

  if (!is.null(seed)) {
    set.seed(seed)
  }

  prior_model <- example_model(
    type = model_family,
    n_cues = n_cues,
    categories = categories,
    kappa = 5L,
    nu = 5L
  )

  prior_mus <- get_expected_mu(prior_model)
  prior_sigmas <- get_expected_sigma(prior_model)
  mu1 <- as.numeric(prior_mus[[categories[1]]])
  mu2 <- as.numeric(prior_mus[[categories[2]]])
  sigma1 <- prior_sigmas[[categories[1]]]
  sigma2 <- prior_sigmas[[categories[2]]]
  if (n_cues == 1L) {
    sigma1 <- as.matrix(sigma1)
    sigma2 <- as.matrix(sigma2)
  }
  cue_names <- get_cue_labels(prior_model)
  delta_mu <- mu2 - mu1

  if (n_cues == 1L) {
    conditions <- list(
      no_exposure = 0,
      left_shifted = -0.5 * delta_mu,
      right_shifted = 0.5 * delta_mu
    )
  } else if (n_cues == 2L) {
    conditions <- list(
      no_exposure = c(0, 0),
      top_right = c(0.5 * delta_mu[1], 0.5 * delta_mu[2]),
      bottom_left = c(-0.5 * delta_mu[1], -0.5 * delta_mu[2]),
      top_left = c(-0.5 * delta_mu[1], 0.5 * delta_mu[2]),
      bottom_right = c(0.5 * delta_mu[1], -0.5 * delta_mu[2])
    )
  } else {
    conditions <- list(
      no_exposure = rep(0, n_cues),
      top_right_front = c(
        0.5 * delta_mu[1], 0.5 * delta_mu[2], 0.5 * delta_mu[3]
      ),
      top_right_back = c(
        0.5 * delta_mu[1], 0.5 * delta_mu[2], -0.5 * delta_mu[3]
      ),
      top_left_front = c(
        -0.5 * delta_mu[1], 0.5 * delta_mu[2], 0.5 * delta_mu[3]
      ),
      top_left_back = c(
        -0.5 * delta_mu[1], 0.5 * delta_mu[2], -0.5 * delta_mu[3]
      ),
      bottom_right_front = c(
        0.5 * delta_mu[1], -0.5 * delta_mu[2], 0.5 * delta_mu[3]
      ),
      bottom_right_back = c(
        0.5 * delta_mu[1], -0.5 * delta_mu[2], -0.5 * delta_mu[3]
      ),
      bottom_left_front = c(
        -0.5 * delta_mu[1], -0.5 * delta_mu[2], 0.5 * delta_mu[3]
      ),
      bottom_left_back = c(
        -0.5 * delta_mu[1], -0.5 * delta_mu[2], -0.5 * delta_mu[3]
      )
    )
  }

  exposure_list <- list()
  test_list <- list()

  for (cond_name in names(conditions)) {
    s <- conditions[[cond_name]]
    is_no_exp <- identical(cond_name, "no_exposure")

    if (!is_no_exp && n_exposure_per_category > 0L) {
      mu1_shift <- mu1 + s
      mu2_shift <- mu2 + s
      d1 <- mvtnorm::rmvnorm(
        n_exposure_per_category,
        mean = mu1_shift,
        sigma = sigma1
      )
      d2 <- mvtnorm::rmvnorm(
        n_exposure_per_category,
        mean = mu2_shift,
        sigma = sigma2
      )
      exp_cues <- rbind(d1, d2)
      colnames(exp_cues) <- cue_names
      exp_df <- data.frame(
        Phase = "exposure",
        Condition = cond_name,
        group = cond_name,
        Subject = paste0(cond_name, "_sub1"),
        stringsAsFactors = FALSE
      )
      exp_df <- cbind(exp_df, as.data.frame(exp_cues))
      exp_df$category <- rep(categories, each = n_exposure_per_category)
      exposure_list[[cond_name]] <- exp_df

      updated_model <- update_template(
        prior_model,
        exp_df,
        noise_treatment = "no_noise"
      )
    } else {
      updated_model <- prior_model
    }

    post_mus <- get_expected_mu(updated_model)
    m1 <- as.numeric(post_mus[[categories[1]]])
    m2 <- as.numeric(post_mus[[categories[2]]])

    u_vals <- seq(0, 1, length.out = n_test_steps_category)
    if (is_no_exp || n_exposure_per_category == 0L) {
      test_locations <- t(vapply(u_vals, function(u) {
        (1 - u) * mu1 + u * mu2
      }, numeric(n_cues)))
    } else {
      v_vals <- seq(0, 1, length.out = n_test_steps_adaptation)
      grid_uv <- expand.grid(u = u_vals, v = v_vals)
      test_locations <- t(apply(grid_uv, 1, function(row) {
        u <- row[["u"]]
        v <- row[["v"]]
        (1 - v) * ((1 - u) * mu1 + u * mu2) + v * ((1 - u) * m1 + u * m2)
      }))
    }
    if (n_cues == 1L) {
      test_locations <- matrix(test_locations, ncol = 1L)
    }
    colnames(test_locations) <- cue_names

    n_locs <- nrow(test_locations)
    test_rep_idx <- rep(seq_len(n_locs), each = n_test_tokens_per_location)
    test_cues_df <- as.data.frame(
      test_locations[test_rep_idx, , drop = FALSE]
    )
    rownames(test_cues_df) <- NULL

    resp <- categorize(
      updated_model,
      test_cues_df,
      decision_rule = "sample",
      simplify = TRUE
    )

    test_df <- data.frame(
      Phase = "test",
      Condition = cond_name,
      group = cond_name,
      Subject = paste0(cond_name, "_sub1"),
      stringsAsFactors = FALSE
    )
    test_df <- cbind(test_df, test_cues_df)
    test_df$response_category <- factor(resp, levels = categories)
    test_list[[cond_name]] <- test_df
  }

  exposure_all <- if (length(exposure_list) > 0L) {
    do.call(rbind, exposure_list)
  } else {
    base_cols <- list(
      Phase = character(0),
      Condition = character(0),
      group = character(0),
      Subject = character(0)
    )
    for (cn in cue_names) {
      base_cols[[cn]] <- numeric(0)
    }
    base_cols[["category"]] <- character(0)
    as.data.frame(base_cols, stringsAsFactors = FALSE)
  }
  test_all <- do.call(rbind, test_list)

  rownames(exposure_all) <- NULL
  rownames(test_all) <- NULL

  exposure_all$Phase <- factor(exposure_all$Phase)
  exposure_all$Condition <- factor(
    exposure_all$Condition,
    levels = names(conditions)
  )
  exposure_all$group <- factor(
    exposure_all$group,
    levels = names(conditions)
  )
  exposure_all$category <- factor(
    exposure_all$category,
    levels = categories
  )
  test_all$Phase <- factor(test_all$Phase)
  test_all$Condition <- factor(
    test_all$Condition,
    levels = names(conditions)
  )
  test_all$group <- factor(
    test_all$group,
    levels = names(conditions)
  )
  test_all$response_category <- factor(
    test_all$response_category,
    levels = categories
  )

  list(
    exposure = tibble::as_tibble(exposure_all),
    test = tibble::as_tibble(test_all)
  )
}

#' Example IdealAdaptorStanfitInput and IdealAdaptorStaninput objects
#'
#' Constructs example \code{\link{IdealAdaptorStanfitInput}} and
#' \code{\link{IdealAdaptorStaninput}} objects from data generated by
#' \code{\link{example_exposure_test_data}}.
#'
#' \code{example_ideal_adaptor_stanfit_input} constructs the complete high-level
#' \code{\link{IdealAdaptorStanfitInput}} object (including data, transform
#' information, and metadata). \code{example_ideal_adaptor_staninput} extracts
#' the underlying typed Stan-ready input (\code{\link{IdealAdaptorStaninput}})
#' from its \code{@staninput} property.
#'
#' @param model_family Character string indicating the model family:
#'   \code{"NIW"} (default), \code{"NIX"}, or \code{"MNIX"}.
#' @param n_cues Integer number of acoustic cues (1 or 2). Default: 1 (or 2
#'   for \code{"MNIX"}).
#' @param categories Character vector of length 2 naming the categories.
#'   Default: \code{c("/d/", "/t/")}.
#' @param control Named list of Stan input control parameters created by
#'   \code{\link{control_staninput}}. Default: \code{control_staninput()}.
#' @param data Optional list of \code{exposure} and \code{test} tibbles as
#'   returned by \code{\link{example_exposure_test_data}}. If \code{NULL}
#'   (default), data are generated via \code{\link{example_exposure_test_data}}.
#' @param ... Additional arguments passed to
#'   \code{\link{example_exposure_test_data}} when \code{data} is \code{NULL}.
#'
#' @return For \code{example_ideal_adaptor_stanfit_input}, an object of class
#'   \code{\link{IdealAdaptorStanfitInput}}. For
#'   \code{example_ideal_adaptor_staninput}, an object of class
#'   \code{\link{IdealAdaptorStaninput}}.
#'
#' @seealso \code{\link{example_exposure_test_data}},
#'   \code{\link{new_ideal_adaptor_stanfit_input}}
#' @rdname example_ideal_adaptor_stanfit_input
#' @export
example_ideal_adaptor_stanfit_input <- function(
  model_family = c("NIW", "NIX", "MNIX"),
  n_cues = NULL,
  categories = c("/d/", "/t/"),
  control = control_staninput(),
  data = NULL,
  ...
) {
  model_family <- match.arg(model_family)
  if (is.null(n_cues)) {
    n_cues <- if (model_family == "MNIX") 2L else 1L
  }
  if (is.null(data)) {
    data <- example_exposure_test_data(
      model_family = model_family,
      n_cues = n_cues,
      categories = categories,
      ...
    )
  }
  cues <- switch(as.character(n_cues),
    "1" = "VOT",
    "2" = c("VOT", "f0_semitones"),
    "3" = c("VOT", "f0_semitones", "vowel_duration")
  )
  stanmodel <- paste0(model_family, "_ideal_adaptor")
  new_ideal_adaptor_stanfit_input(
    exposure = data$exposure,
    test = data$test,
    cues = cues,
    category = "category",
    response = "response_category",
    group = "group",
    control = control,
    stanmodel = stanmodel
  )
}

#' @rdname example_ideal_adaptor_stanfit_input
#' @export
example_ideal_adaptor_staninput <- function(
  model_family = c("NIW", "NIX", "MNIX"),
  n_cues = NULL,
  categories = c("/d/", "/t/"),
  ...
) {
  stanfit_input <- example_ideal_adaptor_stanfit_input(
    model_family = model_family,
    n_cues = n_cues,
    categories = categories,
    ...
  )
  stanfit_input@staninput
}

#' Example IdealAdaptorStanfit object
#'
#' Fits or loads a pre-fitted example \code{\link{IdealAdaptorStanfit}} model.
#' If \code{file} is specified, the model is loaded from or saved to disk
#' using \code{\link{fit_ideal_adaptor}} according to \code{file_refit}.
#'
#' @param model_family Character string indicating the model family:
#'   \code{"NIW"} (default), \code{"NIX"}, or \code{"MNIX"}.
#' @param n_cues Integer number of acoustic cues: 1 (\code{"VOT"}),
#'   2 (\code{"VOT"} and \code{"f0_semitones"}), or 3 (\code{"VOT"},
#'   \code{"f0_semitones"}, and \code{"vowel_duration"}).
#'   Default: 1 (or 2 for \code{"MNIX"}).
#' @param categories Character vector of length 2 naming the categories.
#'   Default: \code{c("/d/", "/t/")}.
#' @param seed Random seed for reproducible data generation and sampling.
#'   Default: 42L.
#' @param file Either \code{NULL} or a character string naming the file
#'   (without or with \code{.rds} extension) to load or save the model.
#' @param file_refit Character string controlling cache behavior when
#'   \code{file} is specified: \code{"on_change"} (default),
#'   \code{"never"}, or \code{"always"}. See \code{\link{fit_ideal_adaptor}}.
#' @param staninput_control Named list of Stan input control parameters
#'   created by \code{\link{control_staninput}}. Default:
#'   \code{control_staninput()}.
#' @param ... Additional arguments passed to \code{\link{fit_ideal_adaptor}}
#'   (e.g., sampler arguments such as \code{chains}, \code{iter},
#'   \code{control = list(adapt_delta = 0.99)}).
#'
#' @return An object of class \code{\link{IdealAdaptorStanfit}}.
#'
#' @seealso \code{\link{example_ideal_adaptor_stanfit_input}},
#'   \code{\link{fit_ideal_adaptor}}
#' @rdname example_ideal_adaptor_stanfit
#' @export
example_ideal_adaptor_stanfit <- function(
  model_family = c("NIW", "NIX", "MNIX"),
  n_cues = NULL,
  categories = c("/d/", "/t/"),
  seed = 42L,
  staninput_control = MVBeliefUpdatr::control_staninput(),
  file = NULL,
  file_refit = c("on_change", "never", "always"),
  ...
) {
  model_family <- match.arg(model_family)
  file_refit <- match.arg(file_refit)
  if (is.null(n_cues)) {
    n_cues <- if (model_family == "MNIX") 2L else 1L
  }
  stanfit_input <- example_ideal_adaptor_stanfit_input(
    model_family = model_family,
    n_cues = n_cues,
    categories = categories,
    seed = seed,
    control = staninput_control
  )
  stanmodel <- paste0(model_family, "_ideal_adaptor")
  fit_ideal_adaptor(
    stanfit_input = stanfit_input,
    file = file,
    file_refit = file_refit,
    seed = seed,
    stanmodel = stanmodel,
    ...
  )
}
