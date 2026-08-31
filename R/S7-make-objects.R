#' @include S7-family-coercion.R
NULL

# Construct S7 category-representation templates and full cognitive models
# directly from data. These are the modern replacements for the deprecated
# make_*_from_data()/lift_*() functions in make-objects.R. Grouped data is
# not (yet) supported: construct one template/model per group by calling
# these functions once per group's subset of the data.

#' @rdname family-uvg
#' @export
new_uvg_category_representation_from_data <- function(
  data,
  category = "category",
  cues
) {
  .assert_true(
    length(cues) == 1L,
    msg = "UVG representations require exactly one cue."
  )
  .assert_data_frame_like(data)
  .assert_non_NA_scalar_character(
    category,
    msg = "category must be a non-empty scalar character value."
  )
  .assert_data_contains_cols(data, category)
  .assert_data_contains_cols(data, cues)
  category_labels <- unique(as.character(data[[category]]))
  .assert_true(
    length(category_labels) == 1L,
    msg = "data must contain exactly one category."
  )

  cue_values <- as.numeric(data[[cues]])
  new_uvg_category_representation(
    category_labels = category_labels,
    cue_labels = cues,
    mu = mean(cue_values),
    sigma2 = stats::var(cue_values)
  )
}

#' @rdname family-nix
#' @export
new_nix_category_representation_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = 3
) {
  .assert_non_NA_scalar_numeric(
    kappa,
    msg = "kappa must be a non-NA scalar numeric value."
  )
  .assert_non_NA_scalar_numeric(
    nu,
    msg = "nu must be a non-NA scalar numeric value."
  )
  .assert_true(
    nu > 2,
    msg = "nu must be greater than 2 for a univariate NIX representation."
  )

  uvg <- new_uvg_category_representation_from_data(
    data,
    category = category,
    cues = cues
  )
  as_nix_category_representation(uvg, kappa = kappa, nu = nu)
}

#' @rdname family-muvg
#' @export
new_muvg_category_representation_from_data <- function(
  data,
  category = "category",
  cues
) {
  .assert_data_frame_like(data)
  .assert_non_NA_scalar_character(
    category,
    msg = "category must be a non-empty scalar character value."
  )
  .assert_true(
    is.character(cues) && length(cues) > 0,
    msg = "cues must be a non-empty character vector."
  )
  .assert_data_contains_cols(data, category)
  .assert_data_contains_cols(data, cues)
  category_labels <- unique(as.character(data[[category]]))
  .assert_true(
    length(category_labels) == 1L,
    msg = "data must contain exactly one category."
  )

  cue_values <- as.matrix(data[, cues, drop = FALSE])
  new_muvg_category_representation(
    category_labels = category_labels,
    cue_labels = cues,
    component_mu = .colMeans(cue_values),
    component_sigma2 = diag(.cov(cue_values))
  )
}

#' @rdname family-mnix
#' @export
new_mnix_category_representation_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = 3
) {
  .assert_non_NA_scalar_numeric(
    kappa,
    msg = "kappa must be a non-NA scalar numeric value."
  )
  .assert_non_NA_scalar_numeric(
    nu,
    msg = "nu must be a non-NA scalar numeric value."
  )
  .assert_true(
    nu > 2,
    msg = "nu must be greater than 2 for each MNIX component."
  )

  muvg <- new_muvg_category_representation_from_data(
    data,
    category = category,
    cues = cues
  )
  as_mnix_category_representation(muvg, kappa = kappa, nu = nu)
}

#' @rdname family-mvg
#' @export
new_mvg_category_representation_from_data <- function(
  data,
  category = "category",
  cues
) {
  .assert_data_frame_like(data)
  .assert_non_NA_scalar_character(
    category,
    msg = "category must be a non-empty scalar character value."
  )
  .assert_true(
    is.character(cues) && length(cues) > 0,
    msg = "cues must be a non-empty character vector."
  )
  .assert_data_contains_cols(data, category)
  .assert_data_contains_cols(data, cues)
  category_labels <- unique(as.character(data[[category]]))
  .assert_true(
    length(category_labels) == 1L,
    msg = "data must contain exactly one category."
  )

  cue_values <- as.matrix(data[, cues, drop = FALSE])
  new_mvg_category_representation(
    category_labels = category_labels,
    cue_labels = cues,
    mu = .colMeans(cue_values),
    Sigma = .cov(cue_values)
  )
}

#' @rdname family-niw
#' @export
new_niw_category_representation_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = length(cues) + 2
) {
  .assert_non_NA_scalar_numeric(
    kappa,
    msg = "kappa must be a non-NA scalar numeric value."
  )
  .assert_non_NA_scalar_numeric(
    nu,
    msg = "nu must be a non-NA scalar numeric value."
  )
  .assert_true(
    nu > length(cues) + 1,
    msg = paste0(
      "nu must be larger than dimensionality of cues + 1 (>",
      length(cues) + 1, ")."
    )
  )

  mvg <- new_mvg_category_representation_from_data(
    data,
    category = category,
    cues = cues
  )
  as_niw_category_representation(mvg, kappa = kappa, nu = nu)
}

#' @rdname family-exemplar
#' @export
new_exemplar_category_representation_from_data <- function(
  data,
  category = "category",
  cues,
  c = NULL
) {
  .assert_data_frame_like(data)
  .assert_non_NA_scalar_character(
    category,
    msg = "category must be a non-empty scalar character value."
  )
  .assert_true(
    is.character(cues) && length(cues) > 0,
    msg = "cues must be a non-empty character vector."
  )
  .assert_data_contains_cols(data, category)
  .assert_data_contains_cols(data, cues)
  category_labels <- unique(as.character(data[[category]]))
  .assert_true(
    length(category_labels) == 1L,
    msg = "data must contain exactly one category."
  )

  new_exemplar_category_representation(
    category_labels = category_labels,
    cue_labels = cues,
    exemplars = as.matrix(data[, cues, drop = FALSE]),
    c = c
  )
}

#' Construct S7 Objects from Data by Family Type
#'
#' Generic dynamic-dispatch constructors that create category representations,
#' representation templates, or full cognitive models from empirical data by
#' specifying the model family type string.
#'
#' @name new-objects-from-data
#' @aliases new-objects-from-data
#' @param data A data frame or tibble containing category and cue
#'   observations. For category representations, data must contain
#'   observations from exactly one category.
#' @param type Family type string: `"UVG"`, `"NIX"`, `"MUVG"`, `"MNIX"`,
#'   `"MVG"`, `"NIW"`, or `"EXEMPLAR"`.
#' @param category Name of the column in `data` containing category labels.
#'   Defaults to `"category"`.
#' @param cues Column name(s) in `data` containing cue observations.
#' @param ... Additional arguments forwarded to the specific family
#'   constructor (e.g., `kappa`, `nu`, `c`, `decision_rule`, `Sigma_noise`).
#' @return An S7 object of the requested type (category representation,
#'   template, or cognitive model).
#' @seealso [UVG_CategoryRepresentation], [NIX_CategoryRepresentation],
#'   [MUVG_CategoryRepresentation], [MNIX_CategoryRepresentation],
#'   [MVG_CategoryRepresentation], [NIW_CategoryRepresentation],
#'   [Exemplar_CategoryRepresentation]
#' @rdname new-objects-from-data
#' @export
new_category_representation_from_data <- function(
  data,
  type,
  category = "category",
  cues,
  ...
) {
  type <- toupper(type)
  constructor <- switch(
    type,
    UVG = new_uvg_category_representation_from_data,
    NIX = new_nix_category_representation_from_data,
    MUVG = new_muvg_category_representation_from_data,
    MNIX = new_mnix_category_representation_from_data,
    MVG = new_mvg_category_representation_from_data,
    NIW = new_niw_category_representation_from_data,
    EXEMPLAR = new_exemplar_category_representation_from_data,
    .stop("type must be one of UVG, NIX, MUVG, MNIX, MVG, NIW, or EXEMPLAR.")
  )
  constructor(data = data, category = category, cues = cues, ...)
}

#' @rdname family-uvg
#' @export
new_uvg_category_representation_template_from_data <- function(
  data,
  category = "category",
  cues,
  verbose = FALSE
) {
  .assert_true(
    length(cues) == 1L,
    msg = "UVG templates require exactly one cue."
  )
  category_labels <- sort(unique(as.character(data[[category]])))
  representations <- lapply(category_labels, function(label) {
    new_uvg_category_representation_from_data(
      data[data[[category]] == label, , drop = FALSE],
      category = category,
      cues = cues
    )
  })
  names(representations) <- category_labels
  if (verbose) {
    message(
      "Constructed a UVG category-representation template with ",
      length(category_labels), " categories."
    )
  }
  new_category_representation_template(representations)
}

#' @rdname family-nix
#' @export
new_nix_category_representation_template_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = 3,
  verbose = FALSE
) {
  .assert_true(
    length(cues) == 1L,
    msg = "NIX templates require exactly one cue."
  )
  category_labels <- sort(unique(as.character(data[[category]])))
  representations <- lapply(category_labels, function(label) {
    new_nix_category_representation_from_data(
      data[data[[category]] == label, , drop = FALSE],
      category = category,
      cues = cues,
      kappa = kappa,
      nu = nu
    )
  })
  names(representations) <- category_labels
  if (verbose) {
    message(
      "Constructed a NIX category-representation template with ",
      length(category_labels), " categories."
    )
  }
  new_category_representation_template(representations)
}

#' @rdname family-muvg
#' @export
new_muvg_category_representation_template_from_data <- function(
  data,
  category = "category",
  cues,
  verbose = FALSE
) {
  category_labels <- sort(unique(as.character(data[[category]])))
  representations <- lapply(category_labels, function(label) {
    new_muvg_category_representation_from_data(
      data[data[[category]] == label, , drop = FALSE],
      category = category,
      cues = cues
    )
  })
  names(representations) <- category_labels
  if (verbose) {
    message(
      "Constructed a MUVG category-representation template with ",
      length(category_labels), " categories and ",
      length(cues), " cue(s)."
    )
  }
  new_category_representation_template(representations)
}

#' @rdname family-mnix
#' @export
new_mnix_category_representation_template_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = 3,
  verbose = FALSE
) {
  category_labels <- sort(unique(as.character(data[[category]])))
  representations <- lapply(category_labels, function(label) {
    new_mnix_category_representation_from_data(
      data[data[[category]] == label, , drop = FALSE],
      category = category,
      cues = cues,
      kappa = kappa,
      nu = nu
    )
  })
  names(representations) <- category_labels
  if (verbose) {
    message(
      "Constructed a MNIX category-representation template with ",
      length(category_labels), " categories and ",
      length(cues), " cue(s)."
    )
  }
  new_category_representation_template(representations)
}

#' @rdname family-mvg
#' @export
new_mvg_category_representation_template_from_data <- function(
  data,
  category = "category",
  cues,
  verbose = FALSE
) {
  .assert_data_frame_like(data)
  .assert_non_NA_scalar_character(
    category,
    msg = "category must be a non-empty scalar character value."
  )
  .assert_true(
    is.character(cues) && length(cues) > 0,
    msg = "cues must be a non-empty character vector."
  )
  .assert_data_contains_cols(data, category)
  .assert_data_contains_cols(data, cues)

  category_labels <- sort(unique(as.character(data[[category]])))
  representations <- vector("list", length(category_labels))
  names(representations) <- category_labels

  for (label in category_labels) {
    representations[[label]] <- new_mvg_category_representation_from_data(
      data[data[[category]] == label, , drop = FALSE],
      category = category,
      cues = cues
    )
  }

  if (verbose) {
    message(
      "Constructed an MVG category-representation template with ",
      length(category_labels), " categories and ",
      length(cues), " cue(s)."
    )
  }

  new_category_representation_template(representations)
}

#' @rdname family-niw
#' @export
new_niw_category_representation_template_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = length(cues) + 2,
  verbose = FALSE
) {
  .assert_non_NA_scalar_numeric(
    kappa,
    msg = "kappa must be a non-NA scalar numeric value."
  )
  .assert_non_NA_scalar_numeric(
    nu,
    msg = "nu must be a non-NA scalar numeric value."
  )
  .assert_true(
    nu > length(cues) + 1,
    msg = paste0(
      "nu must be larger than dimensionality of cues + 1 (>",
      length(cues) + 1, ")."
    )
  )

  mvg_template <- new_mvg_category_representation_template_from_data(
    data,
    category = category,
    cues = cues,
    verbose = FALSE
  )

  representations <- lapply(
    names(mvg_template@representations),
    function(label) {
      new_niw_category_representation_from_data(
        data[data[[category]] == label, , drop = FALSE],
        category = category,
        cues = cues,
        kappa = kappa,
        nu = nu
      )
    }
  )
  names(representations) <- names(mvg_template@representations)

  if (verbose) {
    message(
      "S is set so that the expected category covariance matrix Sigma ",
      "matches the category covariance in the sample (given nu). ",
      "It might be safer to fit an Inverse-Wishart distribution to the ",
      "entire set of covariance matrices."
    )
  }

  new_category_representation_template(representations)
}

#' @rdname family-exemplar
#' @export
new_exemplar_category_representation_template_from_data <- function(
  data,
  category = "category",
  cues,
  c = NULL,
  verbose = FALSE
) {
  .assert_data_frame_like(data)
  .assert_non_NA_scalar_character(
    category,
    msg = "category must be a non-empty scalar character value."
  )
  .assert_true(
    is.character(cues) && length(cues) > 0,
    msg = "cues must be a non-empty character vector."
  )
  .assert_data_contains_cols(data, category)
  .assert_data_contains_cols(data, cues)

  category_labels <- sort(unique(as.character(data[[category]])))
  representations <- vector("list", length(category_labels))
  names(representations) <- category_labels

  for (label in category_labels) {
    representations[[label]] <- new_exemplar_category_representation_from_data(
      data[data[[category]] == label, , drop = FALSE],
      category = category,
      cues = cues,
      c = c
    )
  }

  if (verbose) {
    message(
      "Constructed an exemplar category-representation template with ",
      length(category_labels), " categories and ", length(cues), " cue(s)."
    )
  }

  new_category_representation_template(representations)
}

#' @rdname new-objects-from-data
#' @export
new_category_representation_template_from_data <- function(
  data,
  type,
  category = "category",
  cues,
  ...
) {
  type <- toupper(type)
  constructor <- switch(
    type,
    UVG = new_uvg_category_representation_template_from_data,
    NIX = new_nix_category_representation_template_from_data,
    MUVG = new_muvg_category_representation_template_from_data,
    MNIX = new_mnix_category_representation_template_from_data,
    MVG = new_mvg_category_representation_template_from_data,
    NIW = new_niw_category_representation_template_from_data,
    EXEMPLAR = new_exemplar_category_representation_template_from_data,
    .stop("type must be one of UVG, NIX, MUVG, MNIX, MVG, NIW, or EXEMPLAR.")
  )
  constructor(data = data, category = category, cues = cues, ...)
}

#' @rdname family-uvg
#' @export
new_uvg_ideal_observer_from_data <- function(
  data,
  category = "category",
  cues,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = "no_noise",
  lapse_treatment = "no_lapses",
  verbose = FALSE
) {
  template <- new_uvg_category_representation_template_from_data(
    data,
    category = category,
    cues = cues,
    verbose = verbose
  )
  new_uvg_ideal_observer(
    category_template = template,
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment
  )
}

#' @rdname family-nix
#' @export
new_nix_ideal_adaptor_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = 3,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = "no_noise",
  lapse_treatment = "no_lapses",
  verbose = FALSE
) {
  template <- new_nix_category_representation_template_from_data(
    data,
    category = category,
    cues = cues,
    kappa = kappa,
    nu = nu,
    verbose = verbose
  )
  new_nix_ideal_adaptor(
    category_template = template,
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment
  )
}

#' @rdname family-muvg
#' @export
new_muvg_ideal_observer_from_data <- function(
  data,
  category = "category",
  cues,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = "no_noise",
  lapse_treatment = "no_lapses",
  verbose = FALSE
) {
  template <- new_muvg_category_representation_template_from_data(
    data,
    category = category,
    cues = cues,
    verbose = verbose
  )
  new_muvg_ideal_observer(
    category_template = template,
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment
  )
}

#' @rdname family-mnix
#' @export
new_mnix_ideal_adaptor_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = 3,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = "no_noise",
  lapse_treatment = "no_lapses",
  verbose = FALSE
) {
  template <- new_mnix_category_representation_template_from_data(
    data,
    category = category,
    cues = cues,
    kappa = kappa,
    nu = nu,
    verbose = verbose
  )
  new_mnix_ideal_adaptor(
    category_template = template,
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment
  )
}

#' @rdname family-mvg
#' @export
new_mvg_ideal_observer_from_data <- function(
  data,
  category = "category",
  cues,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = "no_noise",
  lapse_treatment = "no_lapses",
  verbose = FALSE
) {
  template <- new_mvg_category_representation_template_from_data(
    data,
    category = category,
    cues = cues,
    verbose = verbose
  )
  new_mvg_ideal_observer(
    category_template = template,
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment
  )
}

#' @rdname family-niw
#' @export
new_niw_ideal_adaptor_from_data <- function(
  data,
  category = "category",
  cues,
  kappa = nu,
  nu = length(cues) + 2,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = "no_noise",
  lapse_treatment = "no_lapses",
  verbose = FALSE
) {
  template <- new_niw_category_representation_template_from_data(
    data,
    category = category,
    cues = cues,
    kappa = kappa,
    nu = nu,
    verbose = verbose
  )
  new_niw_ideal_adaptor(
    category_template = template,
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment
  )
}

#' @rdname family-exemplar
#' @export
new_exemplar_model_from_data <- function(
  data,
  category = "category",
  cues,
  c = NULL,
  decision_rule = "sampling",
  category_prior = NULL,
  lapse_rate = 0,
  lapse_bias = NULL,
  Sigma_noise = NULL,
  noise_treatment = "no_noise",
  lapse_treatment = "no_lapses",
  verbose = FALSE
) {
  template <- new_exemplar_category_representation_template_from_data(
    data,
    category = category,
    cues = cues,
    c = c,
    verbose = verbose
  )
  new_exemplar_model(
    category_template = template,
    decision_rule = decision_rule,
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    Sigma_noise = Sigma_noise,
    noise_treatment = noise_treatment,
    lapse_treatment = lapse_treatment
  )
}

#' @rdname new-objects-from-data
#' @export
new_model_from_data <- function(
  data,
  type,
  category = "category",
  cues,
  ...
) {
  type <- toupper(type)
  constructor <- switch(
    type,
    UVG = new_uvg_ideal_observer_from_data,
    NIX = new_nix_ideal_adaptor_from_data,
    MUVG = new_muvg_ideal_observer_from_data,
    MNIX = new_mnix_ideal_adaptor_from_data,
    MVG = new_mvg_ideal_observer_from_data,
    NIW = new_niw_ideal_adaptor_from_data,
    EXEMPLAR = new_exemplar_model_from_data,
    .stop("type must be one of UVG, NIX, MUVG, MNIX, MVG, NIW, or EXEMPLAR.")
  )
  constructor(data = data, category = category, cues = cues, ...)
}

#' Aggregate a list of cognitive models into a single model
#'
#' Takes a list of S7 cognitive models of the same class and averages their
#' parameters and category representations into a single aggregated model.
#'
#' @param models A non-empty list of S7 cognitive model objects (all of the
#'   same class).
#' @param weights An optional numeric vector of positive weights corresponding
#'   to each model in `models`. If `NULL`, uniform weights are used. Weights
#'   are automatically normalized to sum to 1.
#' @param ... Additional arguments (currently unused).
#'
#' @return An S7 cognitive model of the same class as the elements in `models`,
#'   with averaged parameters.
#' @seealso \code{\link{new_mvg_ideal_observer}},
#'   \code{\link{new_niw_ideal_adaptor}}
#' @export
aggregate_models <- function(models, weights = NULL, ...) {
  .assert_true(
    is.list(models) && length(models) > 0,
    msg = "models must be a non-empty list of cognitive model objects."
  )

  first_model <- models[[1]]
  .assert_true(
    S7::S7_inherits(first_model, MVBU_CognitiveModel),
    msg = "All elements in models must be S7 MVBU_CognitiveModel objects."
  )

  target_class <- S7::S7_class(first_model)
  for (i in seq_along(models)) {
    .assert_true(
      S7::S7_inherits(models[[i]], target_class),
      msg = paste0(
        "All elements in models must have the same S7 class (expected ",
        target_class@name, ")."
      )
    )
  }

  n_models <- length(models)
  if (is.null(weights)) {
    weights <- rep(1 / n_models, n_models)
  } else {
    .assert_true(
      is.numeric(weights) && length(weights) == n_models &&
        all(weights >= 0) && sum(weights) > 0,
      msg = paste0(
        "weights must be a numeric vector of positive weights ",
        "with the same length as models."
      )
    )
    weights <- weights / sum(weights)
  }

  if (n_models == 1) {
    return(first_model)
  }

  # 1. Aggregate model-level parameters
  # category_prior
  priors_list <- lapply(models, function(m) get_category_prior(m))
  first_prior <- priors_list[[1]]
  agg_category_prior <- if (!is.null(first_prior)) {
    prior_names <- names(first_prior)
    weighted_prior <- Reduce(
      `+`,
      mapply(function(p, w) p * w, priors_list, weights, SIMPLIFY = FALSE)
    )
    if (!is.null(prior_names)) names(weighted_prior) <- prior_names
    weighted_prior
  } else {
    NULL
  }

  # lapse_rate
  lapse_rates <- vapply(models, function(m) get_lapse_rate(m), numeric(1))
  agg_lapse_rate <- sum(lapse_rates * weights)

  # lapse_bias
  lapse_biases_list <- lapply(models, function(m) get_lapse_bias(m))
  agg_lapse_bias <- if (!is.null(lapse_biases_list[[1]])) {
    lb_names <- names(lapse_biases_list[[1]])
    weighted_lb <- Reduce(
      `+`,
      mapply(
        function(lb, w) lb * w,
        lapse_biases_list,
        weights,
        SIMPLIFY = FALSE
      )
    )
    if (!is.null(lb_names)) names(weighted_lb) <- lb_names
    weighted_lb
  } else {
    NULL
  }

  # Sigma_noise
  sigma_noise_list <- lapply(models, function(m) m@noise_behavior$Sigma_noise)
  agg_sigma_noise <- if (!is.null(sigma_noise_list[[1]])) {
    Reduce(
      `+`,
      mapply(
        function(sn, w) sn * w,
        sigma_noise_list,
        weights,
        SIMPLIFY = FALSE
      )
    )
  } else {
    NULL
  }

  # 2. Aggregate representations / category templates
  first_reps <- get_category_representations(first_model)
  categories <- names(first_reps)

  agg_reps <- list()
  for (cat in categories) {
    cat_reps <- lapply(models, function(m) {
      get_category_representations(m)[[cat]]
    })
    first_cat_rep <- cat_reps[[1]]
    cat_rep_class <- S7::S7_class(first_cat_rep)

    cat_labels <- get_category_labels(first_cat_rep)
    cue_labels <- get_cue_labels(first_cat_rep)

    if (S7::S7_inherits(first_cat_rep, MVG_CategoryRepresentation)) {
      mu_weighted <- Reduce(
        `+`,
        mapply(function(r, w) r@mu * w, cat_reps, weights, SIMPLIFY = FALSE)
      )
      Sigma_weighted <- Reduce(
        `+`,
        mapply(function(r, w) r@Sigma * w, cat_reps, weights, SIMPLIFY = FALSE)
      )
      agg_reps[[cat]] <- new_mvg_category_representation(
        category_labels = cat_labels,
        cue_labels = cue_labels,
        mu = mu_weighted,
        Sigma = Sigma_weighted
      )
    } else if (S7::S7_inherits(first_cat_rep, NIW_CategoryRepresentation)) {
      m_weighted <- Reduce(
        `+`,
        mapply(function(r, w) r@m * w, cat_reps, weights, SIMPLIFY = FALSE)
      )
      S_weighted <- Reduce(
        `+`,
        mapply(function(r, w) r@S * w, cat_reps, weights, SIMPLIFY = FALSE)
      )
      kappas <- vapply(cat_reps, function(r) r@kappa, numeric(1))
      nus <- vapply(cat_reps, function(r) r@nu, numeric(1))
      agg_reps[[cat]] <- new_niw_category_representation(
        category_labels = cat_labels,
        cue_labels = cue_labels,
        m = m_weighted,
        S = S_weighted,
        kappa = sum(kappas * weights),
        nu = sum(nus * weights)
      )
    } else if (S7::S7_inherits(first_cat_rep, UVG_CategoryRepresentation)) {
      mu_weighted <- sum(
        vapply(cat_reps, function(r) r@mu, numeric(1)) * weights
      )
      sigma2_weighted <- sum(
        vapply(cat_reps, function(r) r@sigma2, numeric(1)) * weights
      )
      agg_reps[[cat]] <- new_uvg_category_representation(
        category_labels = cat_labels,
        cue_labels = cue_labels,
        mu = mu_weighted,
        sigma2 = sigma2_weighted
      )
    } else if (S7::S7_inherits(first_cat_rep, NIX_CategoryRepresentation)) {
      m_weighted <- sum(
        vapply(cat_reps, function(r) r@m, numeric(1)) * weights
      )
      S_weighted <- sum(
        vapply(cat_reps, function(r) r@S, numeric(1)) * weights
      )
      kappas <- vapply(cat_reps, function(r) r@kappa, numeric(1))
      nus <- vapply(cat_reps, function(r) r@nu, numeric(1))
      agg_reps[[cat]] <- new_nix_category_representation(
        category_labels = cat_labels,
        cue_labels = cue_labels,
        m = m_weighted,
        S = S_weighted,
        kappa = sum(kappas * weights),
        nu = sum(nus * weights)
      )
    } else if (S7::S7_inherits(first_cat_rep, MNIX_CategoryRepresentation)) {
      m_weighted <- Reduce(
        `+`,
        mapply(function(r, w) r@m * w, cat_reps, weights, SIMPLIFY = FALSE)
      )
      S_weighted <- Reduce(
        `+`,
        mapply(function(r, w) r@S * w, cat_reps, weights, SIMPLIFY = FALSE)
      )
      kappas <- Reduce(
        `+`,
        mapply(function(r, w) r@kappa * w, cat_reps, weights, SIMPLIFY = FALSE)
      )
      nus <- Reduce(
        `+`,
        mapply(function(r, w) r@nu * w, cat_reps, weights, SIMPLIFY = FALSE)
      )
      agg_reps[[cat]] <- new_mnix_category_representation(
        category_labels = cat_labels,
        cue_labels = cue_labels,
        m = m_weighted,
        S = S_weighted,
        kappa = kappas,
        nu = nus
      )
    } else {
      .stop(
        "Aggregation is not implemented for category representation class: ",
        cat_rep_class@name
      )
    }

  }

  agg_template <- new_category_representation_template(
    representations = agg_reps
  )

  # 3. Construct aggregated model
  args <- list(
    category_template = agg_template,
    category_prior = agg_category_prior,
    lapse_rate = agg_lapse_rate,
    lapse_bias = agg_lapse_bias,
    Sigma_noise = agg_sigma_noise,
    decision_rule = first_model@decision_rule,
    noise_treatment = first_model@noise_behavior$noise_treatment,
    lapse_treatment = first_model@lapse_behavior$lapse_treatment,
    metadata = first_model@metadata
  )

  constructor <- if (S7::S7_inherits(first_model, MVG_IdealObserver)) {
    new_mvg_ideal_observer
  } else if (S7::S7_inherits(first_model, NIW_IdealAdaptor)) {
    new_niw_ideal_adaptor
  } else if (S7::S7_inherits(first_model, UVG_IdealObserver)) {
    new_uvg_ideal_observer
  } else if (S7::S7_inherits(first_model, NIX_IdealAdaptor)) {
    new_nix_ideal_adaptor
  } else if (S7::S7_inherits(first_model, MUVG_IdealObserver)) {
    new_muvg_ideal_observer
  } else if (S7::S7_inherits(first_model, MNIX_IdealAdaptor)) {
    new_mnix_ideal_adaptor
  } else if (S7::S7_inherits(first_model, Exemplar_Model)) {
    new_exemplar_model
  } else {
    .stop("Aggregation is not supported for model class: ", target_class@name)
  }

  do.call(constructor, args)
}


