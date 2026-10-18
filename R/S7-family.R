#' @include internal-asserts.R
#' @include internal-globals.R
#' @importFrom S7 new_class new_property new_union class_list class_function
#' @importFrom S7 class_any class_character S7_inherits
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom stats cov cov.wt dnorm dt rnorm rt
NULL

# ==============================================================================
# Model Family Framework
# ==============================================================================
#
# Encapsulates the first-class S7 model family representation (MVBU_Family),
# the built-in family definitions, global registry environment, registration
# APIs, bidirectional family <-> Stan model resolution, and computational
# density/moment primitives.

# ------------------------------------------------------------------------------
# 1. Model Family Class & Constructor
# ------------------------------------------------------------------------------

#' Model Family Class
#'
#' Encapsulates the mathematical primitives, moment extraction, and naming
#' conventions for a distributional model family in `MVBeliefUpdatr`.
#'
#' @param family_names A list containing `name_short`, `name_long`, and
#'   `name_details`. Default: \code{list()}.
#' @param density_function Function `f(x, rep, log = FALSE)` computing density.
#' @param sampling_function Function `r(rep, n = 1L, with_replacement = TRUE)`
#'   drawing sample observations from the category representation.
#' @param expected_moments Function `f(rep)` returning a list with `mu`,
#'   `Sigma`, and optional component `weights`.
#' @param marginal_moments Optional function `f(rep)` returning marginal
#'   predictive moments for adaptors. Default: \code{NULL}.
#' @param sufficient_statistics_function Optional function extracting
#'   sufficient data statistics. Default: \code{NULL}.
#' @param stanmodel Optional list bundling `stanmodel_name` and `stanfit_class`
#'   for families associated with a Stan model, or \code{NULL}.
#'   Default: \code{NULL}.
#'
#' @name MVBU_Family
#' @seealso [new_mvbu_family()], [register_model_family()], [get_model_family()]
#' @export
MVBU_Family <- S7::new_class(
  "MVBU_Family",
  properties = list(
    family_names = S7::class_list,
    density_function = S7::class_function,
    sampling_function = S7::class_function,
    expected_moments = S7::class_function,
    marginal_moments = S7::class_any,
    sufficient_statistics_function = S7::class_any,
    stanmodel = S7::new_property(
      class = S7::new_union(NULL, S7::class_list),
      default = quote(NULL)
    )
  ),
  validator = function(self) {
    if (is.null(self@family_names$name_short) ||
      !nzchar(self@family_names$name_short)) {
      return(
        "family_names must contain a non-empty name_short character string."
      )
    }
    if (!is.null(self@stanmodel)) {
      if (!is.list(self@stanmodel)) {
        return("stanmodel must be NULL or a list.")
      }
      sm_name <- self@stanmodel$stanmodel_name
      if (is.null(sm_name) || !is.character(sm_name) ||
        length(sm_name) != 1L || !nzchar(sm_name)) {
        return(paste0(
          "stanmodel must be a list containing a non-empty scalar ",
          "character string 'stanmodel_name'."
        ))
      }
      sfc <- self@stanmodel$stanfit_class
      if (!is.null(sfc)) {
        if (!is.character(sfc) || length(sfc) != 1L || !nzchar(sfc)) {
          return(paste0(
            "stanmodel$stanfit_class must be a non-empty scalar ",
            "character string or NULL."
          ))
        }
      }
    }
    NULL
  }
)

#' Construct a new MVBU_Family object
#'
#' @param name_short Short uppercase identifier for the family (e.g., "UVG").
#' @param name_long Long descriptive name. Default: \code{name_short}.
#' @param name_details Optional list specifying tier class/type names:
#'   `CategoryRepresentation_name`, `CategoryTemplate_name`,
#'   `CognitiveModel_name`, `CognitiveModelList_name`. Default: \code{NULL}.
#' @param density_function Density function `function(x, rep, log = FALSE)`.
#' @param sampling_function Sampling function
#'   `function(rep, n = 1L, with_replacement = TRUE, ...)`.
#' @param expected_moments Function `function(rep)` returning
#'   `list(mu, Sigma, weights)`.
#' @param marginal_moments Optional function returning marginal predictive
#'   moments. Default: \code{NULL}.
#' @param sufficient_statistics_function Optional function returning
#'   sufficient statistics. Default: \code{NULL}.
#' @param stanmodel Optional list bundling `stanmodel_name` and `stanfit_class`
#'   for Stan model families, or \code{NULL}. Default: \code{NULL}.
#' @return An [MVBU_Family] object.
#' @seealso [MVBU_Family], [register_model_family()], [get_model_family()]
#' @export
new_mvbu_family <- function(
  name_short,
  name_long = NULL,
  name_details = NULL,
  density_function,
  sampling_function,
  expected_moments,
  marginal_moments = NULL,
  sufficient_statistics_function = NULL,
  stanmodel = NULL
) {
  .assert_character_scalar(name_short)
  name_short <- trimws(name_short)
  if (!grepl("_", name_short)) {
    name_short <- toupper(name_short)
  }
  name_long <- if (!is.null(name_long)) {
    .assert_character_scalar(name_long)
    trimws(name_long)
  } else {
    name_short
  }

  if (!is.null(stanmodel)) {
    .assert_true(is.list(stanmodel), msg = "stanmodel must be a list or NULL.")
  }

  defaults_details <- list(
    CategoryRepresentation_name =
      paste0(name_short, "_CategoryRepresentation"),
    CategoryTemplate_name =
      paste0(name_short, "_CategoryRepresentationTemplate"),
    CognitiveModel_name = paste0(name_short, "_CognitiveModel"),
    CognitiveModelList_name = paste0(name_short, "_CognitiveModelList")
  )
  if (is.null(name_details)) {
    name_details <- defaults_details
  } else {
    for (k in names(defaults_details)) {
      if (is.null(name_details[[k]])) {
        name_details[[k]] <- defaults_details[[k]]
      }
    }
  }

  MVBU_Family(
    family_names = list(
      name_short = name_short,
      name_long = name_long,
      name_details = name_details
    ),
    density_function = density_function,
    sampling_function = sampling_function,
    expected_moments = expected_moments,
    marginal_moments = marginal_moments,
    sufficient_statistics_function = sufficient_statistics_function,
    stanmodel = stanmodel
  )
}

# ------------------------------------------------------------------------------
# 2. Registry Environment & Built-in Families
# ------------------------------------------------------------------------------

.mvbu_family_registry <- new.env(parent = emptyenv())

#' Build the list of built-in model families
#'
#' @return A named list of [MVBU_Family] objects.
#' @keywords internal
#' @noRd
.build_builtin_families <- function() {
  list(
    UVG = new_mvbu_family(
      name_short = "UVG",
      name_long = "Univariate Gaussian Ideal Observer",
      stanmodel = NULL,
      name_details = list(
        CategoryRepresentation_name = "UVG_CategoryRepresentation",
        CategoryTemplate_name = "UVG_CategoryRepresentationTemplate",
        CognitiveModel_name = "UVG_IdealObserver",
        CognitiveModelList_name = "UVG_CognitiveModelList"
      ),
      density_function = function(x, rep, log = FALSE,
                                  noise_treatment = "no_noise",
                                  Sigma_noise = NULL) {
        x <- .as_observation_matrix(x, d = 1, arg_name = "x")
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(
            n = nrow(x),
            mean = rep(0, ncol(x)),
            sigma = Sigma_noise
          )
        }
        noise_variance <- if (!is.null(Sigma_noise) &&
          (identical(noise_treatment, "sample") ||
            identical(noise_treatment, "marginalize"))) {
          as.numeric(Sigma_noise[1, 1])
        } else {
          0
        }
        stats::dnorm(
          x[, 1],
          mean = rep@mu,
          sd = sqrt(rep@sigma2 + noise_variance),
          log = log
        )
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        matrix(stats::rnorm(n, mean = rep@mu, sd = sqrt(rep@sigma2)), ncol = 1)
      },
      expected_moments = function(rep) {
        list(mu = rep@mu, Sigma = matrix(rep@sigma2, 1, 1), weights = 1)
      }
    ),
    NIX = new_mvbu_family(
      name_short = "NIX",
      name_long = "Normal-Inverse-Chi-Squared Ideal Adaptor",
      stanmodel = list(
        stanmodel_name = "NIX_ideal_adaptor",
        stanfit_class = "NIX_IdealAdaptorStanfit"
      ),
      name_details = list(
        CategoryRepresentation_name = "NIX_CategoryRepresentation",
        CategoryTemplate_name = "NIX_CategoryRepresentationTemplate",
        CognitiveModel_name = "NIX_IdealAdaptor",
        CognitiveModelList_name = "NIX_CognitiveModelList"
      ),
      density_function = function(x, rep, log = FALSE,
                                  noise_treatment = "no_noise",
                                  Sigma_noise = NULL) {
        x <- .as_observation_matrix(x, d = 1, arg_name = "x")
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(
            n = nrow(x),
            mean = rep(0, ncol(x)),
            sigma = Sigma_noise
          )
        }
        noise_variance <- if (!is.null(Sigma_noise) &&
          (identical(noise_treatment, "sample") ||
            identical(noise_treatment, "marginalize"))) {
          as.numeric(Sigma_noise[1, 1])
        } else {
          0
        }
        scale <- sqrt((rep@sigma2 + noise_variance) *
          (rep@kappa + 1) / rep@kappa)
        d <- stats::dt((x[, 1] - rep@m) / scale, df = rep@nu, log = log)
        if (log) d - log(scale) else d / scale
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        scale <- sqrt(rep@sigma2 * (rep@kappa + 1) / rep@kappa)
        matrix(rep@m + scale * stats::rt(n, df = rep@nu), ncol = 1)
      },
      expected_moments = function(rep) {
        list(
          mu = rep@m, Sigma = matrix(get_expected_sigma(rep), 1, 1),
          weights = 1
        )
      },
      marginal_moments = function(rep) {
        list(
          mu = rep@m, Sigma = matrix(get_marginal_sigma(rep), 1, 1),
          weights = 1
        )
      }
    ),
    MUVG = new_mvbu_family(
      name_short = "MUVG",
      name_long = "Multivariate Independent Gaussian Ideal Observer",
      stanmodel = NULL,
      name_details = list(
        CategoryRepresentation_name = "MUVG_CategoryRepresentation",
        CategoryTemplate_name = "MUVG_CategoryRepresentationTemplate",
        CognitiveModel_name = "MUVG_IdealObserver",
        CognitiveModelList_name = "MUVG_CognitiveModelList"
      ),
      density_function = function(x, rep, log = FALSE,
                                  noise_treatment = "no_noise",
                                  Sigma_noise = NULL) {
        k <- length(rep@mu)
        x <- .as_observation_matrix(x, d = k, arg_name = "x")
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(
            n = nrow(x),
            mean = rep(0, ncol(x)),
            sigma = Sigma_noise
          )
        }
        noise_diag <- if (!is.null(Sigma_noise) &&
          (identical(noise_treatment, "sample") ||
            identical(noise_treatment, "marginalize"))) {
          diag(Sigma_noise)
        } else {
          rep(0, k)
        }
        d_mat <- vapply(seq_len(k), function(i) {
          stats::dnorm(
            x[, i],
            mean = rep@mu[i],
            sd = sqrt(rep@sigma2[i] + noise_diag[i]),
            log = TRUE
          )
        }, numeric(nrow(x)))
        tot <- rowSums(matrix(d_mat, nrow = nrow(x), ncol = k))
        if (log) tot else exp(tot)
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        k <- length(rep@mu)
        out <- vapply(seq_len(k), function(i) {
          stats::rnorm(n, mean = rep@mu[i], sd = sqrt(rep@sigma2[i]))
        }, numeric(n))
        matrix(out, nrow = n, ncol = k)
      },
      expected_moments = function(rep) {
        list(
          mu = rep@mu, Sigma = diag(rep@sigma2, nrow = length(rep@sigma2)),
          weights = rep@weights
        )
      }
    ),
    MNIX = new_mvbu_family(
      name_short = "MNIX",
      name_long = paste(
        "Multivariate Independent Normal-Inverse-Chi-Squared",
        "Ideal Adaptor"
      ),
      stanmodel = list(
        stanmodel_name = "MNIX_ideal_adaptor_ideal_integrator",
        stanfit_class = "MNIX_IdealAdaptor_IdealIntegratorStanfit"
      ),
      name_details = list(
        CategoryRepresentation_name = "MNIX_CategoryRepresentation",
        CategoryTemplate_name = "MNIX_CategoryRepresentationTemplate",
        CognitiveModel_name = "MNIX_IdealAdaptor",
        CognitiveModelList_name = "MNIX_CognitiveModelList"
      ),
      density_function = function(x, rep, log = FALSE,
                                  noise_treatment = "no_noise",
                                  Sigma_noise = NULL) {
        if (length(rep@kappa) == 1L) {
          .mnix_early_density_function(
            x, rep,
            log = log, noise_treatment = noise_treatment,
            Sigma_noise = Sigma_noise
          )
        } else {
          .mnix_late_density_function(
            x, rep,
            weights = NULL, log = log,
            noise_treatment = noise_treatment, Sigma_noise = Sigma_noise
          )
        }
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        k <- length(rep@m)
        kappa_vec <- if (length(rep@kappa) == 1L) {
          rep(rep@kappa, k)
        } else {
          rep@kappa
        }
        nu_vec <- if (length(rep@nu) == 1L) rep(rep@nu, k) else rep@nu
        out <- vapply(seq_len(k), function(i) {
          scale <- sqrt(rep@sigma2[i] * (kappa_vec[i] + 1) / kappa_vec[i])
          rep@m[i] + scale * stats::rt(n, df = nu_vec[i])
        }, numeric(n))
        matrix(out, nrow = n, ncol = k)
      },
      expected_moments = function(rep) {
        list(mu = rep@m, Sigma = get_expected_sigma(rep), weights = rep@weights)
      },
      marginal_moments = function(rep) {
        list(mu = rep@m, Sigma = get_marginal_sigma(rep), weights = rep@weights)
      }
    ),
    # Late MNIX
    MNIX_ideal_integrator = new_mvbu_family(
      name_short = "MNIX_ideal_integrator",
      name_long = "MNIX Ideal Adaptor with Late Ideal Integrator",
      stanmodel = list(
        stanmodel_name = "MNIX_ideal_adaptor_ideal_integrator",
        stanfit_class = "MNIX_IdealAdaptor_IdealIntegratorStanfit"
      ),
      name_details = list(
        CategoryRepresentation_name = character(0),
        CategoryTemplate_name = character(0),
        CognitiveModel_name = "MNIX_IdealAdaptor_IdealIntegrator",
        CognitiveModelList_name = character(0)
      ),
      density_function = function(x, rep, log = FALSE,
                                  noise_treatment = "no_noise",
                                  Sigma_noise = NULL) {
        .mnix_late_density_function(
          x, rep,
          weights = NULL, log = log,
          noise_treatment = noise_treatment, Sigma_noise = Sigma_noise
        )
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        .build_builtin_families()$MNIX@sampling_function(
          rep, n = n, with_replacement = with_replacement, ...
        )
      },
      expected_moments = function(rep) {
        list(mu = rep@m, Sigma = get_expected_sigma(rep), weights = rep@weights)
      },
      marginal_moments = function(rep) {
        list(mu = rep@m, Sigma = get_marginal_sigma(rep), weights = rep@weights)
      }
    ),
    MNIX_weighted_integrator = new_mvbu_family(
      name_short = "MNIX_weighted_integrator",
      name_long = "MNIX Ideal Adaptor with Late Weighted Integrator",
      stanmodel = list(
        stanmodel_name = "MNIX_ideal_adaptor_weighted_integrator",
        stanfit_class = "MNIX_IdealAdaptor_WeightedIntegratorStanfit"
      ),
      name_details = list(
        CategoryRepresentation_name = character(0),
        CategoryTemplate_name = character(0),
        CognitiveModel_name = "MNIX_IdealAdaptor_WeightedIntegrator",
        CognitiveModelList_name = character(0)
      ),
      density_function = function(x, rep, log = FALSE,
                                  noise_treatment = "no_noise",
                                  Sigma_noise = NULL) {
        .mnix_late_density_function(
          x, rep,
          weights = rep@weights, log = log,
          noise_treatment = noise_treatment, Sigma_noise = Sigma_noise
        )
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        .build_builtin_families()$MNIX@sampling_function(
          rep, n = n, with_replacement = with_replacement, ...
        )
      },
      expected_moments = function(rep) {
        list(mu = rep@m, Sigma = get_expected_sigma(rep), weights = rep@weights)
      },
      marginal_moments = function(rep) {
        list(mu = rep@m, Sigma = get_marginal_sigma(rep), weights = rep@weights)
      }
    ),
    MNIX_category_weighted_integrator = new_mvbu_family(
      name_short = "MNIX_category_weighted_integrator",
      name_long = "MNIX Ideal Adaptor with Late Category-Weighted Integrator",
      stanmodel = list(
        stanmodel_name = "MNIX_ideal_adaptor_category_weighted_integrator",
        stanfit_class = "MNIX_IdealAdaptor_CategoryWeightedIntegratorStanfit"
      ),
      name_details = list(
        CategoryRepresentation_name = character(0),
        CategoryTemplate_name = character(0),
        CognitiveModel_name = "MNIX_IdealAdaptor_CategoryWeightedIntegrator",
        CognitiveModelList_name = character(0)
      ),
      density_function = function(x, rep, log = FALSE,
                                  noise_treatment = "no_noise",
                                  Sigma_noise = NULL) {
        .mnix_late_density_function(
          x, rep,
          weights = rep@weights, log = log,
          noise_treatment = noise_treatment, Sigma_noise = Sigma_noise
        )
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        .build_builtin_families()$MNIX@sampling_function(
          rep, n = n, with_replacement = with_replacement, ...
        )
      },
      expected_moments = function(rep) {
        list(mu = rep@m, Sigma = get_expected_sigma(rep), weights = rep@weights)
      },
      marginal_moments = function(rep) {
        list(mu = rep@m, Sigma = get_marginal_sigma(rep), weights = rep@weights)
      }
    ),
    MNIX_group_weighted_integrator = new_mvbu_family(
      name_short = "MNIX_group_weighted_integrator",
      name_long = "MNIX Ideal Adaptor with Late Group-Weighted Integrator",
      stanmodel = list(
        stanmodel_name = "MNIX_ideal_adaptor_group_weighted_integrator",
        stanfit_class = "MNIX_IdealAdaptor_GroupWeightedIntegratorStanfit"
      ),
      name_details = list(
        CategoryRepresentation_name = character(0),
        CategoryTemplate_name = character(0),
        CognitiveModel_name = character(0),
        CognitiveModelList_name = character(0)
      ),
      density_function = function(x, rep, log = FALSE,
                                  noise_treatment = "no_noise",
                                  Sigma_noise = NULL) {
        .mnix_late_density_function(
          x, rep,
          weights = rep@weights, log = log,
          noise_treatment = noise_treatment, Sigma_noise = Sigma_noise
        )
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        .build_builtin_families()$MNIX@sampling_function(
          rep, n = n, with_replacement = with_replacement, ...
        )
      },
      expected_moments = function(rep) {
        list(mu = rep@m, Sigma = get_expected_sigma(rep), weights = rep@weights)
      },
      marginal_moments = function(rep) {
        list(mu = rep@m, Sigma = get_marginal_sigma(rep), weights = rep@weights)
      }
    ),
    MNIX_group_category_weighted_integrator = new_mvbu_family(
      name_short = "MNIX_group_category_weighted_integrator",
      name_long = paste(
        "MNIX Ideal Adaptor with Late",
        "Group-Category-Weighted Integrator"
      ),
      stanmodel = list(
        stanmodel_name =
          "MNIX_ideal_adaptor_group_category_weighted_integrator",
        stanfit_class =
          "MNIX_IdealAdaptor_GroupCategoryWeightedIntegratorStanfit"
      ),
      name_details = list(
        CategoryRepresentation_name = character(0),
        CategoryTemplate_name = character(0),
        CognitiveModel_name = character(0),
        CognitiveModelList_name = character(0)
      ),
      density_function = function(x, rep, log = FALSE,
                                  noise_treatment = "no_noise",
                                  Sigma_noise = NULL) {
        .mnix_late_density_function(
          x, rep,
          weights = rep@weights, log = log,
          noise_treatment = noise_treatment, Sigma_noise = Sigma_noise
        )
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        .build_builtin_families()$MNIX@sampling_function(
          rep, n = n, with_replacement = with_replacement, ...
        )
      },
      expected_moments = function(rep) {
        list(mu = rep@m, Sigma = get_expected_sigma(rep), weights = rep@weights)
      },
      marginal_moments = function(rep) {
        list(mu = rep@m, Sigma = get_marginal_sigma(rep), weights = rep@weights)
      }
    ),
    MNIX_informativity_weighted_integrator = new_mvbu_family(
      name_short = "MNIX_informativity_weighted_integrator",
      name_long = paste(
        "MNIX Ideal Adaptor with Late",
        "Informativity-Weighted Integrator"
      ),
      stanmodel = list(
        stanmodel_name =
          "MNIX_ideal_adaptor_informativity_weighted_integrator",
        stanfit_class =
          "MNIX_IdealAdaptor_InformativityWeightedIntegratorStanfit"
      ),
      name_details = list(
        CategoryRepresentation_name = character(0),
        CategoryTemplate_name = character(0),
        CognitiveModel_name =
          "MNIX_IdealAdaptor_InformativityWeightedIntegrator",
        CognitiveModelList_name = character(0)
      ),
      density_function = function(x, rep, log = FALSE,
                                  noise_treatment = "no_noise",
                                  Sigma_noise = NULL) {
        .mnix_late_density_function(
          x, rep,
          weights = rep@weights, log = log,
          noise_treatment = noise_treatment, Sigma_noise = Sigma_noise
        )
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        .build_builtin_families()$MNIX@sampling_function(
          rep, n = n, with_replacement = with_replacement, ...
        )
      },
      expected_moments = function(rep) {
        list(mu = rep@m, Sigma = get_expected_sigma(rep), weights = rep@weights)
      },
      marginal_moments = function(rep) {
        list(mu = rep@m, Sigma = get_marginal_sigma(rep), weights = rep@weights)
      }
    ),
    # Early MNIX
    ideal_integrator_MNIX = new_mvbu_family(
      name_short = "ideal_integrator_MNIX",
      name_long = "Early Ideal Integrator with MNIX Ideal Adaptor",
      stanmodel = list(
        stanmodel_name = "ideal_integrator_MNIX_ideal_adaptor",
        stanfit_class = "IdealIntegrator_MNIX_IdealAdaptorStanfit"
      ),
      name_details = list(
        CategoryRepresentation_name = character(0),
        CategoryTemplate_name = character(0),
        CognitiveModel_name = "IdealIntegrator_MNIX_IdealAdaptor",
        CognitiveModelList_name = character(0)
      ),
      density_function = function(x, rep, log = FALSE,
                                  noise_treatment = "no_noise",
                                  Sigma_noise = NULL) {
        .mnix_early_density_function(
          x, rep,
          log = log, noise_treatment = noise_treatment,
          Sigma_noise = Sigma_noise
        )
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        .build_builtin_families()$MNIX@sampling_function(
          rep, n = n, with_replacement = with_replacement, ...
        )
      },
      expected_moments = function(rep) {
        list(mu = rep@m, Sigma = get_expected_sigma(rep), weights = rep@weights)
      },
      marginal_moments = function(rep) {
        list(mu = rep@m, Sigma = get_marginal_sigma(rep), weights = rep@weights)
      }
    ),
    weighted_integrator_MNIX = new_mvbu_family(
      name_short = "weighted_integrator_MNIX",
      name_long = "Early Weighted Integrator with MNIX Ideal Adaptor",
      stanmodel = list(
        stanmodel_name = "weighted_integrator_MNIX_ideal_adaptor",
        stanfit_class = "WeightedIntegrator_MNIX_IdealAdaptorStanfit"
      ),
      name_details = list(
        CategoryRepresentation_name = character(0),
        CategoryTemplate_name = character(0),
        CognitiveModel_name = "WeightedIntegrator_MNIX_IdealAdaptor",
        CognitiveModelList_name = character(0)
      ),
      density_function = function(x, rep, log = FALSE,
                                  noise_treatment = "no_noise",
                                  Sigma_noise = NULL) {
        .mnix_early_density_function(
          x, rep,
          log = log, noise_treatment = noise_treatment,
          Sigma_noise = Sigma_noise
        )
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        .build_builtin_families()$MNIX@sampling_function(
          rep, n = n, with_replacement = with_replacement, ...
        )
      },
      expected_moments = function(rep) {
        list(mu = rep@m, Sigma = get_expected_sigma(rep), weights = rep@weights)
      },
      marginal_moments = function(rep) {
        list(mu = rep@m, Sigma = get_marginal_sigma(rep), weights = rep@weights)
      }
    ),
    category_weighted_integrator_MNIX = new_mvbu_family(
      name_short = "category_weighted_integrator_MNIX",
      name_long = "Early Category-Weighted Integrator with MNIX Ideal Adaptor",
      stanmodel = list(
        stanmodel_name = "category_weighted_integrator_MNIX_ideal_adaptor",
        stanfit_class = "CategoryWeightedIntegrator_MNIX_IdealAdaptorStanfit"
      ),
      name_details = list(
        CategoryRepresentation_name = character(0),
        CategoryTemplate_name = character(0),
        CognitiveModel_name = "CategoryWeightedIntegrator_MNIX_IdealAdaptor",
        CognitiveModelList_name = character(0)
      ),
      density_function = function(x, rep, log = FALSE,
                                  noise_treatment = "no_noise",
                                  Sigma_noise = NULL) {
        .mnix_early_density_function(
          x, rep,
          log = log, noise_treatment = noise_treatment,
          Sigma_noise = Sigma_noise
        )
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        .build_builtin_families()$MNIX@sampling_function(
          rep, n = n, with_replacement = with_replacement, ...
        )
      },
      expected_moments = function(rep) {
        list(mu = rep@m, Sigma = get_expected_sigma(rep), weights = rep@weights)
      },
      marginal_moments = function(rep) {
        list(mu = rep@m, Sigma = get_marginal_sigma(rep), weights = rep@weights)
      }
    ),
    group_weighted_integrator_MNIX = new_mvbu_family(
      name_short = "group_weighted_integrator_MNIX",
      name_long = "Early Group-Weighted Integrator with MNIX Ideal Adaptor",
      stanmodel = list(
        stanmodel_name = "group_weighted_integrator_MNIX_ideal_adaptor",
        stanfit_class = "GroupWeightedIntegrator_MNIX_IdealAdaptorStanfit"
      ),
      name_details = list(
        CategoryRepresentation_name = character(0),
        CategoryTemplate_name = character(0),
        CognitiveModel_name = character(0),
        CognitiveModelList_name = character(0)
      ),
      density_function = function(x, rep, log = FALSE,
                                  noise_treatment = "no_noise",
                                  Sigma_noise = NULL) {
        .mnix_early_density_function(
          x, rep,
          log = log, noise_treatment = noise_treatment,
          Sigma_noise = Sigma_noise
        )
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        .build_builtin_families()$MNIX@sampling_function(
          rep, n = n, with_replacement = with_replacement, ...
        )
      },
      expected_moments = function(rep) {
        list(mu = rep@m, Sigma = get_expected_sigma(rep), weights = rep@weights)
      },
      marginal_moments = function(rep) {
        list(mu = rep@m, Sigma = get_marginal_sigma(rep), weights = rep@weights)
      }
    ),
    group_category_weighted_integrator_MNIX = new_mvbu_family(
      name_short = "group_category_weighted_integrator_MNIX",
      name_long = paste(
        "Early Group-Category-Weighted Integrator with",
        "MNIX Ideal Adaptor"
      ),
      stanmodel = list(
        stanmodel_name =
          "group_category_weighted_integrator_MNIX_ideal_adaptor",
        stanfit_class =
          "GroupCategoryWeightedIntegrator_MNIX_IdealAdaptorStanfit"
      ),
      name_details = list(
        CategoryRepresentation_name = character(0),
        CategoryTemplate_name = character(0),
        CognitiveModel_name = character(0),
        CognitiveModelList_name = character(0)
      ),
      density_function = function(x, rep, log = FALSE,
                                  noise_treatment = "no_noise",
                                  Sigma_noise = NULL) {
        .mnix_early_density_function(
          x, rep,
          log = log, noise_treatment = noise_treatment,
          Sigma_noise = Sigma_noise
        )
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        .build_builtin_families()$MNIX@sampling_function(
          rep, n = n, with_replacement = with_replacement, ...
        )
      },
      expected_moments = function(rep) {
        list(mu = rep@m, Sigma = get_expected_sigma(rep), weights = rep@weights)
      },
      marginal_moments = function(rep) {
        list(mu = rep@m, Sigma = get_marginal_sigma(rep), weights = rep@weights)
      }
    ),
    informativity_weighted_integrator_MNIX = new_mvbu_family(
      name_short = "informativity_weighted_integrator_MNIX",
      name_long = paste(
        "Early Informativity-Weighted Integrator with",
        "MNIX Ideal Adaptor"
      ),
      stanmodel = list(
        stanmodel_name =
          "informativity_weighted_integrator_MNIX_ideal_adaptor",
        stanfit_class =
          "InformativityWeightedIntegrator_MNIX_IdealAdaptorStanfit"
      ),
      name_details = list(
        CategoryRepresentation_name = character(0),
        CategoryTemplate_name = character(0),
        CognitiveModel_name =
          "InformativityWeightedIntegrator_MNIX_IdealAdaptor",
        CognitiveModelList_name = character(0)
      ),
      density_function = function(x, rep, log = FALSE,
                                  noise_treatment = "no_noise",
                                  Sigma_noise = NULL) {
        .mnix_early_density_function(
          x, rep,
          log = log, noise_treatment = noise_treatment,
          Sigma_noise = Sigma_noise
        )
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        .build_builtin_families()$MNIX@sampling_function(
          rep, n = n, with_replacement = with_replacement, ...
        )
      },
      expected_moments = function(rep) {
        list(mu = rep@m, Sigma = get_expected_sigma(rep), weights = rep@weights)
      },
      marginal_moments = function(rep) {
        list(mu = rep@m, Sigma = get_marginal_sigma(rep), weights = rep@weights)
      }
    ),
    MVG = new_mvbu_family(
      name_short = "MVG",
      name_long = "Multivariate Gaussian Ideal Observer",
      stanmodel = NULL,
      name_details = list(
        CategoryRepresentation_name = "MVG_CategoryRepresentation",
        CategoryTemplate_name = "MVG_CategoryRepresentationTemplate",
        CognitiveModel_name = "MVG_IdealObserver",
        CognitiveModelList_name = "MVG_CognitiveModelList"
      ),
      density_function = function(x, rep, log = FALSE,
                                  noise_treatment = "no_noise",
                                  Sigma_noise = NULL) {
        d <- length(rep@mu)
        x <- .as_observation_matrix(x, d = d, arg_name = "x")
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(
            n = nrow(x),
            mean = rep(0, ncol(x)),
            sigma = Sigma_noise
          )
        }
        Sigma_eff <- if (!is.null(Sigma_noise) &&
          (identical(noise_treatment, "sample") ||
            identical(noise_treatment, "marginalize"))) {
          rep@Sigma + Sigma_noise
        } else {
          rep@Sigma
        }
        mvtnorm::dmvnorm(x, mean = rep@mu, sigma = Sigma_eff, log = log)
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        mvtnorm::rmvnorm(n, mean = rep@mu, sigma = rep@Sigma)
      },
      expected_moments = function(rep) {
        list(mu = rep@mu, Sigma = rep@Sigma, weights = NULL)
      }
    ),
    NIW = new_mvbu_family(
      name_short = "NIW",
      name_long = "Normal-Inverse-Wishart Ideal Adaptor",
      stanmodel = list(
        stanmodel_name = "NIW_ideal_adaptor",
        stanfit_class = "NIW_IdealAdaptorStanfit"
      ),
      name_details = list(
        CategoryRepresentation_name = "NIW_CategoryRepresentation",
        CategoryTemplate_name = "NIW_CategoryRepresentationTemplate",
        CognitiveModel_name = "NIW_IdealAdaptor",
        CognitiveModelList_name = "NIW_CognitiveModelList"
      ),
      density_function = function(x, rep, log = FALSE,
                                  noise_treatment = "no_noise",
                                  Sigma_noise = NULL) {
        d <- length(rep@m)
        x <- .as_observation_matrix(x, d = d, arg_name = "x")
        df <- rep@nu - d + 1
        Sigma_pred <- ((rep@kappa + 1) / (rep@kappa * df)) * rep@S
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(
            n = nrow(x),
            mean = rep(0, ncol(x)),
            sigma = Sigma_noise
          )
        }
        Sigma_eff <- if (!is.null(Sigma_noise) &&
          (identical(noise_treatment, "sample") ||
            identical(noise_treatment, "marginalize"))) {
          Sigma_pred + Sigma_noise
        } else {
          Sigma_pred
        }
        .dmvt_density(x, mean = rep@m, Sigma = Sigma_eff, df = df, log = log)
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        d <- length(rep@m)
        df <- rep@nu - d + 1
        Sigma_pred <- ((rep@kappa + 1) / (rep@kappa * df)) * rep@S
        .rmvt(n = n, delta = rep@m, sigma = Sigma_pred, df = df)
      },
      expected_moments = function(rep) {
        list(
          mu = get_expected_mu_from_m(rep@m),
          Sigma = get_expected_Sigma_from_S(rep@S, rep@nu),
          weights = NULL
        )
      },
      marginal_moments = function(rep) {
        list(
          mu = get_marginal_mu(rep), Sigma = get_marginal_sigma(rep),
          weights = NULL
        )
      }
    ),
    EXEMPLAR = new_mvbu_family(
      name_short = "EXEMPLAR",
      name_long = "Exemplar Model",
      stanmodel = NULL,
      name_details = list(
        CategoryRepresentation_name = "Exemplar_CategoryRepresentation",
        CategoryTemplate_name = "Exemplar_CategoryRepresentationTemplate",
        CognitiveModel_name = "Exemplar_Model",
        CognitiveModelList_name = "Exemplar_CognitiveModelList"
      ),
      density_function = function(x, rep, log = FALSE,
                                  noise_treatment = "no_noise",
                                  Sigma_noise = NULL) {
        ex0 <- rep@exemplars
        d0 <- ncol(ex0)
        n0 <- nrow(ex0)
        w0 <- as.numeric(rep@exemplar_weights)
        bandwidth_scale <- if (n0 > 1) {
          (4 / ((d0 + 2) * n0))^(2 / (d0 + 4))
        } else {
          1.0
        }
        Sigma0 <- if (n0 > 1) {
          stats::cov(ex0) * bandwidth_scale
        } else {
          diag(1, d0) * bandwidth_scale
        }
        if (!is.matrix(Sigma0) || any(!is.finite(Sigma0))) {
          Sigma0 <- diag(1, d0) * bandwidth_scale
        }
        Sigma0 <- Sigma0 + diag(MVBU_PROB_TOL, d0)

        x <- .as_observation_matrix(x, d = d0, arg_name = "x")
        if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
          x <- x + mvtnorm::rmvnorm(
            n = nrow(x),
            mean = rep(0, ncol(x)),
            sigma = Sigma_noise
          )
        }
        Sigma_eff <- if (!is.null(Sigma_noise) &&
          (identical(noise_treatment, "sample") ||
            identical(noise_treatment, "marginalize"))) {
          Sigma0 + Sigma_noise
        } else {
          Sigma0
        }

        chol_sigma <- tryCatch(chol(Sigma_eff), error = function(e) NULL)
        if (is.null(chol_sigma)) {
          chol_sigma <- chol(Sigma_eff + diag(MVBU_PROB_TOL, d0))
        }
        log_det <- 2 * sum(log(diag(chol_sigma)))

        x_std <- t(backsolve(chol_sigma, t(x), transpose = TRUE))
        ex_std <- t(backsolve(chol_sigma, t(ex0), transpose = TRUE))

        sq_x <- rowSums(x_std^2)
        sq_ex <- rowSums(ex_std^2)
        dist_mat <- outer(sq_x, sq_ex, "+") - 2 * tcrossprod(x_std, ex_std)
        dist_mat[dist_mat < 0] <- 0

        log_dens_by_exemplar <- -0.5 * (d0 * log(2 * pi) + log_det + dist_mat)
        weighted_logdens <- sweep(log_dens_by_exemplar, 2, log(w0), "+")
        log_mix <- .logsumexp_rows(weighted_logdens)
        if (isTRUE(log)) log_mix else exp(log_mix)
      },
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        n_avail <- nrow(rep@exemplars)
        if (isTRUE(with_replacement)) {
          idx <- sample.int(n_avail,
            size = n, replace = TRUE,
            prob = rep@exemplar_weights
          )
        } else {
          if (n > n_avail) {
            cat_lbl <- get_category_labels(rep)
            .stop(sprintf(
              paste0(
                "Cannot sample %d observations without replacement ",
                "from category '%s' which contains only %d exemplars."
              ),
              n, cat_lbl, n_avail
            ))
          }
          idx <- sample.int(n_avail,
            size = n, replace = FALSE,
            prob = rep@exemplar_weights
          )
        }
        rep@exemplars[idx, , drop = FALSE]
      },
      expected_moments = function(rep) {
        list(
          mu = as.numeric(rep@exemplar_weights %*% rep@exemplars),
          Sigma = stats::cov.wt(rep@exemplars, wt = rep@exemplar_weights)$cov,
          weights = NULL
        )
      }
    )
  )
}

#' Initialize the model family registry and Stan model bidirectional lookups
#'
#' @return Invisibly TRUE.
#' @keywords internal
#' @noRd
.init_family_registry <- function() {
  builtins <- .build_builtin_families()
  .mvbu_family_registry$families <- builtins
  .mvbu_family_registry$family_to_stanmodel <- list()
  .mvbu_family_registry$stanmodel_to_family <- list()

  for (fam_name in names(builtins)) {
    fam_obj <- builtins[[fam_name]]
    if (!is.null(fam_obj@stanmodel) &&
      !is.null(fam_obj@stanmodel$stanmodel_name) &&
      nzchar(fam_obj@stanmodel$stanmodel_name)) {
      sm_name <- fam_obj@stanmodel$stanmodel_name
      .mvbu_family_registry$family_to_stanmodel[[fam_name]] <- sm_name
      .mvbu_family_registry$stanmodel_to_family[[sm_name]] <- fam_name
    }
  }
  invisible(TRUE)
}

.init_family_registry()

# ------------------------------------------------------------------------------
# 3. Family Registration API
# ------------------------------------------------------------------------------

#' Register a model family in the MVBU family registry
#'
#' Registers a model family in the package registry. Accepts either an
#' [MVBU_Family] object directly (recommended) or a family name with associated
#' class names and Stan model configuration.
#'
#' @param family An [MVBU_Family] object, or character string specifying the
#'   model family short name.
#' @param category_representation_class Name of category-representation class.
#'   Used when \code{family} is a character string. Default: \code{NULL}.
#' @param cognitive_model_class Name of cognitive-model class.
#'   Used when \code{family} is a character string. Default: \code{NULL}.
#' @param stanmodel Optional list bundling `stanmodel_name` and `stanfit_class`
#'   for Stan model families, or \code{NULL}. Used when \code{family} is a
#'   character string. Default: \code{NULL}.
#' @param overwrite Logical indicating whether to overwrite an existing
#'   registration for this family. Default: \code{TRUE}.
#' @return Invisibly \code{TRUE}.
#' @seealso [MVBU_Family], [new_mvbu_family()], [get_model_family()],
#'   [list_model_families()]
#' @export
register_model_family <- function(
  family,
  category_representation_class = NULL,
  cognitive_model_class = NULL,
  stanmodel = NULL,
  overwrite = TRUE
) {
  if (S7::S7_inherits(family, MVBU_Family)) {
    family_obj <- family
  } else {
    .assert_character_scalar(family)
    family_name <- .normalize_family_name(family)
    .assert_character_scalar(category_representation_class)
    .assert_character_scalar(cognitive_model_class)
    if (!is.null(stanmodel)) {
      .assert_true(
        is.list(stanmodel),
        msg = "stanmodel must be a list or NULL."
      )
    }
    family_obj <- new_mvbu_family(
      name_short = family_name,
      name_details = list(
        CategoryRepresentation_name = category_representation_class,
        CognitiveModel_name = cognitive_model_class
      ),
      density_function = function(x, rep, log = FALSE) numeric(0),
      sampling_function = function(rep, n = 1L, with_replacement = TRUE, ...) {
        matrix(numeric(0), n, 0)
      },
      expected_moments = function(rep) {
        list(mu = numeric(0), Sigma = matrix(numeric(0), 0, 0), weights = NULL)
      },
      stanmodel = stanmodel
    )
  }
  .register_model_family(family_obj, overwrite = overwrite)
}

#' Internal registration of an MVBU_Family object
#'
#' @param family_obj An [MVBU_Family] object.
#' @param overwrite Logical indicating whether to overwrite
#'   an existing family.
#'   Default: \code{TRUE}.
#' @return Invisibly \code{TRUE}.
#' @keywords internal
#' @noRd
.register_model_family <- function(family_obj, overwrite = TRUE) {
  if (!S7::S7_inherits(family_obj, MVBU_Family)) {
    .stop("family_obj must be an MVBU_Family object.")
  }
  family <- family_obj@family_names$name_short
  family <- .normalize_family_name(family)
  if (!isTRUE(overwrite) &&
    !is.null(.mvbu_family_registry$families[[family]])) {
    .stop(sprintf(
      "Model family '%s' is already registered and overwrite = FALSE.",
      family
    ))
  }
  .mvbu_family_registry$families[[family]] <- family_obj
  if (!is.null(family_obj@stanmodel) &&
    !is.null(family_obj@stanmodel$stanmodel_name) &&
    nzchar(family_obj@stanmodel$stanmodel_name)) {
    sm_name <- family_obj@stanmodel$stanmodel_name
    .mvbu_family_registry$family_to_stanmodel[[family]] <- sm_name
    .mvbu_family_registry$stanmodel_to_family[[sm_name]] <- family
  }
  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# 4. Family Lookup & Introspection
# ------------------------------------------------------------------------------

#' List all registered model families
#'
#' Returns the short names of all currently registered model families in
#' `MVBeliefUpdatr`.
#'
#' @return Character vector of registered family short names.
#' @seealso [get_model_family()], [register_model_family()]
#' @export
list_model_families <- function() {
  sort(names(.mvbu_family_registry$families))
}

#' Get an MVBU_Family object from the registry
#'
#' @param family Model family short name (character string).
#' @return The [MVBU_Family] object.
#' @seealso [list_model_families()], [register_model_family()]
#' @export
get_model_family <- function(family) {
  .get_model_family(family)
}

#' Internal retrieval of an MVBU_Family object
#'
#' @param family Model family short name.
#' @return The [MVBU_Family] object.
#' @keywords internal
#' @noRd
.get_model_family <- function(family) {
  family <- .normalize_family_name(family)
  fam <- .mvbu_family_registry$families[[family]]
  if (is.null(fam)) {
    .stop(paste0("No model family registered for '", family, "'."))
  }
  fam
}

#' Normalize a model-family name to a canonical form
#'
#' Matches a candidate family name case-insensitively against currently
#' registered family names.
#'
#' @param family Model family name (character scalar).
#' @return Normalized character string representing the family name.
#' @keywords internal
#' @noRd
.normalize_family_name <- function(family) {
  .assert_character_scalar(family)
  family <- trimws(family)
  reg_names <- names(.mvbu_family_registry$families)
  idx <- which(tolower(reg_names) == tolower(family))
  if (length(idx) == 1L) {
    return(reg_names[idx])
  }
  if (!grepl("_", family)) {
    toupper(family)
  } else {
    family
  }
}

#' Infer canonical model family from an S7 class name
#'
#' @param class_name Character string specifying an S7 class name.
#' @return Character string specifying the canonical family short name.
#' @keywords internal
#' @noRd
.get_family_from_class_name <- function(class_name) {
  class_name <- sub("^.*::", "", class_name)
  families <- tryCatch(list_model_families(), error = function(e) character(0))
  for (fam in families) {
    reg <- tryCatch(
      .get_model_family_registration(fam),
      error = function(e) NULL
    )
    if (!is.null(reg) &&
      class_name %in% c(reg$category_representation, reg$cognitive_model)) {
      return(fam)
    }
  }
  clean_name <- gsub(
    paste0(
      "_(IdealObserver|IdealAdaptor|IdealAdaptorStanfit|",
      "CategoryRepresentation|CategoryRepresentationTemplate|Model|",
      "CognitiveModelList)$"
    ),
    "",
    class_name
  )
  for (fam in families) {
    if (tolower(fam) == tolower(clean_name)) {
      return(fam)
    }
  }
  if (grepl("MNIX", class_name)) {
    return("MNIX")
  }
  if (grepl("NIW", class_name)) {
    return("NIW")
  }
  if (grepl("NIX", class_name)) {
    return("NIX")
  }
  if (grepl("MUVG", class_name)) {
    return("MUVG")
  }
  if (grepl("UVG", class_name)) {
    return("UVG")
  }
  if (grepl("MVG", class_name)) {
    return("MVG")
  }
  if (toupper(clean_name) == "EXEMPLAR") {
    return("EXEMPLAR")
  }
  clean_name
}

#' Get the registration information for a model family
#'
#' @param family Model family name (character string).
#' @return A list with the registered class names and family object.
#' @seealso [get_model_family()], [list_model_families()]
#' @export
get_model_family_registration <- function(family) {
  .get_model_family_registration(family)
}

#' Internal retrieval of model family registration information
#'
#' @param family Model family name.
#' @return A list with registered class names and the family object.
#' @keywords internal
#' @noRd
.get_model_family_registration <- function(family) {
  fam <- .get_model_family(family)
  list(
    category_representation =
      fam@family_names$name_details$CategoryRepresentation_name,
    cognitive_model = fam@family_names$name_details$CognitiveModel_name,
    names = fam@family_names,
    family = fam
  )
}

# ------------------------------------------------------------------------------
# 5. Stan Model Lookups & Mapping
# ------------------------------------------------------------------------------

#' Get the Stan model name associated with a model family
#'
#' @param family Model family short name or [MVBU_Family] object.
#' @return Stan model name (character string) or \code{NULL} if the family
#'   has no Stan model.
#' @seealso [get_model_family()], [list_model_families()]
#' @export
get_stanmodel_name <- function(family) {
  if (S7::S7_inherits(family, MVBU_Family)) {
    return(if (!is.null(family@stanmodel)) {
      family@stanmodel$stanmodel_name
    } else {
      NULL
    })
  }
  fam <- .get_model_family(family)
  if (!is.null(fam@stanmodel)) fam@stanmodel$stanmodel_name else NULL
}

#' Resolve a model family name or object to its associated Stan model name
#'
#' @param family Character string or [MVBU_Family] object.
#' @return Character string specifying the Stan model name.
#' @keywords internal
#' @noRd
.to_stanmodel <- function(family) {
  if (S7::S7_inherits(family, MVBU_Family)) {
    sm_name <- if (!is.null(family@stanmodel)) {
      family@stanmodel$stanmodel_name
    } else {
      NULL
    }
    if (!is.null(sm_name) && nzchar(sm_name)) {
      return(sm_name)
    } else {
      .stop(sprintf(
        "Model family '%s' does not have an associated Stan model.",
        family@family_names$name_short
      ))
    }
  }
  .assert_character_scalar(family)
  family <- trimws(family)
  # If already a recognized Stan model name, return it directly
  if (!is.null(.mvbu_family_registry$stanmodel_to_family[[family]])) {
    return(family)
  }
  fam_norm <- tryCatch(
    .normalize_family_name(family),
    error = function(e) character(0)
  )
  if (length(fam_norm) > 0L) {
    sm <- .mvbu_family_registry$family_to_stanmodel[[fam_norm]]
    if (!is.null(sm)) {
      return(sm)
    }
    if (!is.null(.mvbu_family_registry$families[[fam_norm]])) {
      .stop(sprintf(
        "Model family '%s' does not have an associated Stan model.",
        fam_norm
      ))
    }
  }
  valid_fams <- sort(names(.mvbu_family_registry$family_to_stanmodel))
  valid_sms <- sort(names(.mvbu_family_registry$stanmodel_to_family))
  .stop(sprintf(
    paste0(
      "Unknown model family or Stan model '%s'. Available Stan model ",
      "families: %s (Stan models: %s)."
    ),
    family,
    paste(valid_fams, collapse = ", "),
    paste(valid_sms, collapse = ", ")
  ))
}

#' Resolve a Stan model name or model family to its canonical model family name
#'
#' @param stanmodel Character string or [MVBU_Family] object.
#' @return Character string specifying the model family short name.
#' @keywords internal
#' @noRd
.to_family <- function(stanmodel) {
  if (S7::S7_inherits(stanmodel, MVBU_Family)) {
    return(stanmodel@family_names$name_short)
  }
  .assert_character_scalar(stanmodel)
  stanmodel <- trimws(stanmodel)
  fam_norm <- tryCatch(
    .normalize_family_name(stanmodel),
    error = function(e) character(0)
  )
  if (length(fam_norm) > 0L &&
    !is.null(.mvbu_family_registry$families[[fam_norm]])) {
    return(fam_norm)
  }
  fam <- .mvbu_family_registry$stanmodel_to_family[[stanmodel]]
  if (!is.null(fam)) {
    return(fam)
  }
  valid_sms <- sort(names(.mvbu_family_registry$stanmodel_to_family))
  valid_fams <- sort(names(.mvbu_family_registry$families))
  .stop(sprintf(
    paste0(
      "Unknown Stan model or model family '%s'. Available Stan ",
      "models: %s (all families: %s)."
    ),
    stanmodel,
    paste(valid_sms, collapse = ", "),
    paste(valid_fams, collapse = ", ")
  ))
}

#' Check if a model family is a registered Stan model family
#'
#' @param family Character string or [MVBU_Family] object.
#' @return Logical indicating whether the family is registered and has an
#'   associated Stan model.
#' @keywords internal
#' @noRd
.is_registered_stan_family <- function(family) {
  if (S7::S7_inherits(family, MVBU_Family)) {
    sm_name <- if (!is.null(family@stanmodel)) {
      family@stanmodel$stanmodel_name
    } else {
      NULL
    }
    return(!is.null(sm_name) && nzchar(sm_name))
  }
  if (!.is_scalar_character(family)) {
    return(FALSE)
  }
  fam_norm <- tryCatch(
    .normalize_family_name(family),
    error = function(e) character(0)
  )
  if (length(fam_norm) == 0L) {
    return(FALSE)
  }
  fam_obj <- .mvbu_family_registry$families[[fam_norm]]
  if (is.null(fam_obj)) {
    return(FALSE)
  }
  sm_name <- if (!is.null(fam_obj@stanmodel)) {
    fam_obj@stanmodel$stanmodel_name
  } else {
    NULL
  }
  !is.null(sm_name) && nzchar(sm_name)
}

#' List all registered Stan model names
#'
#' @return Character vector of registered Stan model names.
#' @keywords internal
#' @noRd
.get_registered_stan_models <- function() {
  sort(names(.mvbu_family_registry$stanmodel_to_family))
}

#' Get the Stanfit class name associated with a model family
#'
#' @param family Model family short name, Stan model name, or [MVBU_Family]
#'   object.
#' @return Character string specifying the Stanfit class name, or \code{NULL}.
#' @keywords internal
#' @noRd
.get_stanfit_class <- function(family) {
  fam_obj <- if (S7::S7_inherits(family, MVBU_Family)) {
    family
  } else {
    fam_name <- .to_family(family)
    .get_model_family(fam_name)
  }
  if (!is.null(fam_obj@stanmodel)) {
    fam_obj@stanmodel$stanfit_class
  } else {
    NULL
  }
}

#' Resolve the concrete IdealAdaptorStanfit constructor for a family/Stan model
#'
#' @param family Model family short name, Stan model name, or [MVBU_Family]
#'   object.
#' @return Constructor function for the concrete IdealAdaptorStanfit class.
#' @keywords internal
#' @noRd
.get_stanfit_constructor <- function(family) {
  sfc <- .get_stanfit_class(family)
  if (is.null(sfc)) {
    default_ctor <- get0(
      "IdealAdaptorStanfit",
      envir = asNamespace("MVBeliefUpdatr")
    )
    return(if (!is.null(default_ctor)) default_ctor else NULL)
  }
  if (is.function(sfc) || inherits(sfc, "S7_class")) {
    return(sfc)
  }
  if (is.character(sfc) && length(sfc) == 1L && nzchar(sfc)) {
    ctor <- get0(sfc, envir = asNamespace("MVBeliefUpdatr"), inherits = TRUE)
    if (is.null(ctor)) {
      ctor <- get0(sfc, mode = "function")
    }
    if (!is.null(ctor)) {
      return(ctor)
    }
    fam_id <- if (S7::S7_inherits(family, MVBU_Family)) {
      family@family_names$name_short
    } else {
      family
    }
    .stop(sprintf(
      "Constructor class '%s' registered for family '%s' not found.",
      sfc,
      fam_id
    ))
  }
  default_ctor <- get0(
    "IdealAdaptorStanfit",
    envir = asNamespace("MVBeliefUpdatr")
  )
  if (!is.null(default_ctor)) default_ctor else NULL
}

# ------------------------------------------------------------------------------
# 6. Family Computational & Density Primitives
# ------------------------------------------------------------------------------

#' Generate likelihood closure from an MVBU_Family object
#'
#' @param rep Category representation object.
#' @param fam Optional [MVBU_Family] object. Default: \code{NULL}.
#' @return A likelihood evaluation function
#'   \code{function(x, log = FALSE, ...)}.
#' @keywords internal
#' @noRd
.make_family_likelihood_function <- function(rep, fam = NULL) {
  if (is.null(fam)) {
    fam <- if (!is.null(rep@metadata$family) &&
      S7::S7_inherits(rep@metadata$family, MVBU_Family)) {
      rep@metadata$family
    } else {
      fam_name <- .get_family_from_class_name(class(rep)[1])
      .get_model_family(fam_name)
    }
  }
  fn <- fam@density_function
  fn_formals <- names(formals(fn))
  has_noise <- "noise_treatment" %in% fn_formals || "..." %in% fn_formals

  function(x,
           log = FALSE,
           noise_treatment = "no_noise",
           Sigma_noise = NULL) {
    if (has_noise) {
      fn(
        x = x,
        rep = rep,
        log = log,
        noise_treatment = noise_treatment,
        Sigma_noise = Sigma_noise
      )
    } else {
      if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
        x <- as.matrix(x) + mvtnorm::rmvnorm(
          n = nrow(as.matrix(x)),
          mean = rep(0, ncol(as.matrix(x))),
          sigma = Sigma_noise
        )
      }
      fn(x = x, rep = rep, log = log)
    }
  }
}

#' Normalize input observations to a numeric matrix
#'
#' @param x Observations as vector, matrix, data.frame, or list.
#' @param d Expected number of dimensions (columns), or \code{NULL}.
#'   Default: \code{NULL}.
#' @param arg_name Name of argument for informative error messages.
#'   Default: \code{"x"}.
#' @return Numeric matrix of observations.
#' @keywords internal
#' @noRd
.as_observation_matrix <- function(x, d = NULL, arg_name = "x") {
  if (is.matrix(x) && is.numeric(x) && (is.null(d) || ncol(x) == d)) {
    return(x)
  }

  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (is.list(x) && !is.data.frame(x)) {
    if (length(x) == 0) {
      return(matrix(numeric(0), nrow = 0, ncol = if (!is.null(d)) d else 0))
    }
    valid_list <- vapply(x, function(el) {
      if (is.matrix(el) || is.data.frame(el)) {
        m <- as.matrix(el)
        return(is.numeric(m) && (is.null(d) || ncol(m) == d))
      }
      if (is.numeric(el) && is.vector(el)) {
        return(is.null(d) || length(el) == d)
      }
      FALSE
    }, logical(1))

    if (all(valid_list)) {
      mats <- lapply(x, function(el) {
        if (is.matrix(el)) {
          return(el)
        }
        if (is.data.frame(el)) {
          return(as.matrix(el))
        }
        ncol_val <- if (!is.null(d)) d else length(el)
        return(matrix(as.numeric(el), ncol = ncol_val, byrow = TRUE))
      })
      return(do.call(rbind, mats))
    }
  }

  if (is.matrix(x)) {
    if (!is.numeric(x)) {
      .stop(arg_name, " must be numeric.")
    }
    if (!is.null(d) && ncol(x) != d) {
      .stop(arg_name, " must have ", d, " column(s).")
    }
    return(x)
  }

  if (is.atomic(x) && is.numeric(x)) {
    if (is.null(d) || d == 1) {
      return(matrix(as.numeric(x), ncol = 1))
    }
    if (length(x) %% d != 0) {
      .stop(arg_name, " length must be a multiple of ", d, ".")
    }
    return(matrix(as.numeric(x), ncol = d, byrow = TRUE))
  }

  .stop(
    arg_name,
    paste0(
      " must be a numeric vector, matrix, data frame, ",
      "or list of numeric vectors."
    )
  )
}

#' Log-sum-exp helper across matrix rows
#'
#' @param matrix_log Matrix of log-scale values.
#' @return Numeric vector of log-sum-exp values per row.
#' @keywords internal
#' @noRd
.logsumexp_rows <- function(matrix_log) {
  max_val <- apply(matrix_log, 1, max)
  max_val + log(rowSums(exp(matrix_log - max_val)))
}

#' Multivariate Student-t density evaluation
#'
#' @param x Observations matrix.
#' @param mean Location vector.
#' @param Sigma Scale covariance matrix.
#' @param df Degrees of freedom.
#' @param log Logical indicating whether to return log density.
#'   Default: \code{FALSE}.
#' @return Numeric vector of densities.
#' @keywords internal
#' @noRd
.dmvt_density <- function(x, mean, Sigma, df, log = FALSE) {
  x <- .as_observation_matrix(x, d = length(mean), arg_name = "x")
  d <- ncol(x)
  diff <- sweep(x, 2, mean, "-")
  chol_sigma <- tryCatch(chol(Sigma), error = function(e) NULL)

  if (is.null(chol_sigma)) {
    ridge <- diag(MVBU_PROB_TOL, d)
    chol_sigma <- chol(Sigma + ridge)
  }

  log_det <- 2 * sum(log(diag(chol_sigma)))
  half_sol <- backsolve(chol_sigma, t(diff), transpose = TRUE)
  quad <- colSums(half_sol^2)

  log_const <- lgamma((df + d) / 2) - lgamma(df / 2) -
    0.5 * (d * log(df * pi) + log_det)
  log_d <- log_const - 0.5 * (df + d) * log1p(quad / df)

  if (isTRUE(log)) {
    log_d
  } else {
    exp(log_d)
  }
}

# Evaluate observation densities for early integration MNIX category
# representations by projecting observations onto weighted composite cue z = w'x,
# computing Student-t density with combined scale and degrees of freedom.
#
# Arguments:
#   x: Observation matrix or vector.
#   rep: MNIX category representation object (contains m, sigma2, kappa, nu, weights).
#   log: Logical; return log density if TRUE. Default: FALSE.
#   noise_treatment: Perceptual noise mode ("no_noise", "sample", "marginalize").
#   When noise_treatment != "no_noise", noise is added to the input cues.
#   Sigma_noise: Optional noise covariance matrix or variance vector.
# Returns:
#   Numeric vector of observation densities (or log densities).
.mnix_early_density_function <- function(
  x, rep,
  log = FALSE,
  noise_treatment = "no_noise",
  Sigma_noise = NULL
) {
  k <- length(rep@m)
  x <- .as_observation_matrix(x, d = k, arg_name = "x")
  if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
    x <- x + mvtnorm::rmvnorm(
      n = nrow(x),
      mean = rep(0, ncol(x)),
      sigma = Sigma_noise
    )
  }
  w <- rep@weights
  z <- as.numeric(x %*% w)
  m_z <- sum(w * rep@m)
  s2_z <- sum(w^2 * rep@sigma2)
  noise_var <- if (!is.null(Sigma_noise) &&
    (identical(noise_treatment, "sample") ||
      identical(noise_treatment, "marginalize"))) {
    if (is.matrix(Sigma_noise)) {
      as.numeric(t(w) %*% Sigma_noise %*% w)
    } else {
      sum(w^2 * as.numeric(Sigma_noise))
    }
  } else {
    0
  }
  scale <- sqrt((s2_z + noise_var) * (rep@kappa + 1) / rep@kappa)
  d <- stats::dt((z - m_z) / scale, df = rep@nu, log = TRUE) - log(scale)
  if (log) d else exp(d)
}

# Evaluate observation densities for late integration MNIX category
# representations by computing univariate Student-t densities along each cue
# dimension and combining them (weighted log-sum or product).
#
# Arguments:
#   x: Observation matrix or vector.
#   rep: MNIX category representation object.
#   weights: Optional cue dimension weights vector. Default: rep@weights.
#   log: Logical; return log density if TRUE. Default: FALSE.
#   noise_treatment: Perceptual noise mode ("no_noise", "sample", "marginalize").
#   When noise_treatment != "no_noise", noise is added to the input cues.
#   Sigma_noise: Optional noise covariance matrix or variance vector.
# Returns:
#   Numeric vector of observation densities (or log densities).
.mnix_late_density_function <- function(
  x, rep,
  weights = rep@weights,
  log = FALSE,
  noise_treatment = "no_noise",
  Sigma_noise = NULL
) {
  k <- length(rep@m)
  x <- .as_observation_matrix(x, d = k, arg_name = "x")
  if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
    x <- x + mvtnorm::rmvnorm(
      n = nrow(x),
      mean = rep(0, ncol(x)),
      sigma = Sigma_noise
    )
  }
  noise_diag <- if (!is.null(Sigma_noise) &&
    (identical(noise_treatment, "sample") ||
      identical(noise_treatment, "marginalize"))) {
    if (is.matrix(Sigma_noise)) {
      diag(Sigma_noise)
    } else {
      as.numeric(Sigma_noise)
    }
  } else {
    rep(0, k)
  }
  kappa_vec <- if (length(rep@kappa) == 1L) {
    rep(rep@kappa, k)
  } else {
    rep@kappa
  }
  nu_vec <- if (length(rep@nu) == 1L) {
    rep(rep@nu, k)
  } else {
    rep@nu
  }
  d_mat <- vapply(seq_len(k), function(i) {
    scale <- sqrt((rep@sigma2[i] + noise_diag[i]) *
      (kappa_vec[i] + 1) / kappa_vec[i])
    stats::dt(
      (x[, i] - rep@m[i]) / scale,
      df = nu_vec[i],
      log = TRUE
    ) - log(scale)
  }, numeric(nrow(x)))
  d_mat <- matrix(d_mat, nrow = nrow(x), ncol = k)

  log_dens <- if (is.null(weights)) {
    as.numeric(rowSums(d_mat))
  } else {
    as.numeric(d_mat %*% as.numeric(weights))
  }
  if (log) log_dens else exp(log_dens)
}
