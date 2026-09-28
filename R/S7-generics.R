#' @include S7-class.R
NULL

# -------------------------
# Core generics
# -------------------------

#' Get the model family from an object
#'
#' @param x An object.
#' @param ... Additional arguments passed to methods.
#' @return A character string representing the model family.
#' @export
get_model_family <- S7::new_generic("get_model_family", "x")

#' Get metadata from an object
#'
#' @param x An object.
#' @param ... Additional arguments passed to methods.
#' @return A list containing metadata associated with the object.
#' @export
get_metadata <- S7::new_generic("get_metadata", "x")

#' Extract parameters or expected parameters from representations, templates, models, or model distributions
#'
#' @title Extract parameters from models, templates, or representations
#' @param x A category representation, template, cognitive model, or stanfit object.
#' @param original_pars Logical scalar; if \code{TRUE}, return original
#'   (unmapped / top-level) Stan model parameter names from the underlying
#'   Stan program (as stored in \code{stanfit@model_pars}) rather than the
#'   expanded, indexed, and post-processed/renamed parameter draw names (as
#'   returned by \code{names(stanfit)}). Default: \code{FALSE}.
#' @param ... Additional arguments passed to methods.
#' @return A list, vector, matrix, or data frame of parameters or expected statistics.
#' @export
get_parameters <- S7::new_generic("get_parameters", "x", function(x, ...) S7::S7_dispatch())

#' @rdname get_parameters
#' @export
get_parameter_names <- S7::new_generic(
  "get_parameter_names",
  "x",
  function(x, original_pars = FALSE, ...) S7::S7_dispatch()
)

#' Aggregate representations, templates, or cognitive models
#'
#' @param x An S7 representation, template, cognitive model, or list of such objects.
#'   When `x` is not an S7 object or list of S7 objects, dispatch falls back to
#'   \code{\link[stats]{aggregate}}.
#' @param ... Additional objects of the same class (for variadic usage) or arguments
#'   passed to methods.
#' @param weights Optional numeric vector of positive weights with the same length as
#'   the number of objects being aggregated. If \code{NULL} (default), uniform weights
#'   \code{1/N} are used.
#' @return An aggregated S7 object of the same class.
#'
#' @rdname aggregate_models
#' @export
aggregate <- S7::new_generic("aggregate", "x", function(x, ...) {
  if (S7::S7_inherits(x, MVBU_Object) || (is.list(x) && length(x) > 0 && S7::S7_inherits(x[[1]], MVBU_Object))) {
    S7::S7_dispatch()
  } else {
    stats::aggregate(x, ...)
  }
})


get_category_representations <- S7::new_generic("get_category_representations", "x")
get_category_template <- S7::new_generic("get_category_template", "x")

#' Add a category representation to a template, model, or model list
#'
#' Adds a category representation to an existing template, cognitive model, or
#' each model in a model list. When adding to a cognitive model, category
#' priors and lapse biases for the new representation default to 0 (rescaling
#' existing categories proportionally), or can be specified for the new
#' category or across all categories.
#'
#' @param x An [MVBU_CategoryRepresentationTemplate], [MVBU_CognitiveModel],
#'   or [MVBU_ModelList] object.
#' @param representation An [MVBU_CategoryRepresentation] object to add.
#' @param name Optional character scalar specifying the name of the category.
#'   Default: `NULL`.
#' @param category_prior Optional numeric vector of category prior
#'   probabilities. Can be a single value for the new category (defaulting to 0,
#'   which rescales existing categories by
#'   `previous_value * (1 - category_prior)`), or a full vector of prior
#'   probabilities across all categories (including the new one) summing to 1.
#'   Default: `NULL`.
#' @param lapse_bias Optional numeric vector of lapse bias probabilities.
#'   Can be a single value for the new category (defaulting to 0, which rescales
#'   existing categories by `previous_value * (1 - lapse_bias)`), or a full
#'   vector of lapse biases across all categories (including the new one)
#'   summing to 1. Default: `NULL`.
#' @param ... Additional arguments passed to methods.
#' @return An updated object of the same class with the category representation
#'   added.
#' @seealso [get_category_labels()]
#' @export
add_category_representation <- S7::new_generic(
  "add_category_representation",
  "x",
  function(
    x,
    representation,
    name = NULL,
    category_prior = NULL,
    lapse_bias = NULL,
    ...
  ) S7::S7_dispatch()
)


#' Get the stored Stan fit from an S7 fit object
#'
#' @param x An S7 fit object.
#' @param ... Additional arguments passed to methods.
#' @return The stored [rstan::stanfit] object when present.
#' @export
get_stanfit <- S7::new_generic("get_stanfit", "x")

#' Set the stored Stan fit on an S7 fit object
#'
#' @param x An S7 fit object.
#' @param stanfit An [rstan::stanfit] object.
#' @return The updated S7 fit object.
#' @rdname get_stanfit
#' @aliases set_stanfit
#' @export
set_stanfit <- S7::new_generic("set_stanfit", c("x", "stanfit"))

#' Extract control parameters of the NUTS sampler
#'
#' @param x An S7 fit object or stanfit object.
#' @param pars Optional parameter names.
#' @param ... Additional arguments.
#' @return A named list of control parameters.
#' @export
control_params <- S7::new_generic(
  "control_params",
  "x",
  function(x, pars = NULL, ...) S7::S7_dispatch()
)

#' Extract diagnostic information from MVBU_Stanfit object
#'
#' @param x An S7 fit object or stanfit object.
#' @param ... Additional arguments.
#' @return Log posterior draws.
#' @rdname diagnostic-quantities
#' @export
log_posterior <- S7::new_generic("log_posterior", "x", function(x, ...) {
  if (inherits(x, "stanfit") || inherits(x, "CmdStanFit")) {
    bayesplot::log_posterior(x, ...)
  } else {
    S7::S7_dispatch()
  }
})

#' Extract NUTS sampler parameters
#'
#' @param x An S7 fit object or stanfit object.
#' @param pars Optional parameter names.
#' @param ... Additional arguments.
#' @return NUTS sampler parameters.
#' @rdname diagnostic-quantities
#' @export
nuts_params <- S7::new_generic(
  "nuts_params",
  "x",
  function(x, pars = NULL, ...) {
    if (inherits(x, "stanfit") || inherits(x, "CmdStanFit")) {
      bayesplot::nuts_params(x, pars = pars, ...)
    } else {
      S7::S7_dispatch()
    }
  }
)

#' Extract Rhat diagnostic values
#'
#' @param x An S7 fit object, draws object, or array.
#' @param pars Optional parameter names.
#' @param ... Additional arguments passed to methods.
#' @return Named numeric vector of Rhat values.
#' @rdname diagnostic-quantities
#' @export
rhat <- S7::new_generic(
  "rhat",
  "x",
  function(x, pars = NULL, ...) S7::S7_dispatch()
)

#' Extract effective sample size ratio
#'
#' @param x An S7 fit object or stanfit object.
#' @param pars Optional parameter names.
#' @param ... Additional arguments.
#' @return Named numeric vector of Neff ratios.
#' @rdname diagnostic-quantities
#' @export
neff_ratio <- S7::new_generic("neff_ratio", "x", function(x, pars = NULL, ...) {
  if (inherits(x, "stanfit") || inherits(x, "CmdStanFit")) {
    bayesplot::neff_ratio(x, pars = pars, ...)
  } else {
    S7::S7_dispatch()
  }
})


#' Get the stored Stan input from an S7 fit or fit-input object
#'
#' @param x An S7 fit or fit-input object.
#' @param ... Additional arguments passed to methods.
#' @return The stored S7 Stan input object.
#' @export
get_staninput <- S7::new_generic("get_staninput", "x")

#' Set the stored Stan input on an S7 fit or fit-input object
#'
#' @param x An S7 fit or fit-input object.
#' @param staninput An S7 Stan input object.
#' @param ... Additional arguments passed to methods.
#' @return The updated S7 object.
#' @export
set_staninput <- S7::new_generic("set_staninput", c("x", "staninput"))

#' Get transform metadata from an S7 fit or fit-input object
#'
#' @param x An S7 fit or fit-input object.
#' @param ... Additional arguments passed to methods.
#' @return A transform-information object.
#' @export
get_transform_information <- S7::new_generic("get_transform_information", "x")

#' Get category-prior values from a cognitive model
#'
#' Compatibility methods accept older list/data-frame shapes such as scalars,
#' named vectors, or table-like objects with a category column. These shims are
#' transitional and will be removed once S7-only representations are the only
#' supported interface.
#' @param x A cognitive model or legacy input object.
#' @param categories Optional category labels used to resolve values.
#' @return A numeric vector of category-prior values.
#' @export
get_category_prior <- S7::new_generic("get_category_prior", c("x", "categories"))

#' Get lapse-rate values from a cognitive model
#'
#' Compatibility methods handle older list/data-frame inputs while the S7 API
#' becomes the standard interface.
#' @param x A cognitive model or legacy input object.
#' @return A numeric lapse-rate value.
#' @export
get_lapse_rate <- S7::new_generic("get_lapse_rate", "x")

#' Get lapse-bias values from a cognitive model
#'
#' Compatibility methods handle older list/data-frame inputs while the S7 API
#' becomes the standard interface.
#' @param x A cognitive model or legacy input object.
#' @param categories Optional category labels used to resolve values.
#' @return A numeric vector of lapse-bias values.
#' @export
get_lapse_bias <- S7::new_generic("get_lapse_bias", c("x", "categories"))

#' Get perceptual noise covariance matrix from a cognitive model
#'
#' @param x A cognitive model or legacy input object.
#' @return A matrix or numeric value representing perceptual noise covariance (\eqn{\Sigma_{\text{noise}}}), or `NULL` if no noise is specified.
#' @export
get_noise <- S7::new_generic("get_noise", "x")

#' Get lapse treatment from a cognitive model
#'
#' @param x A cognitive model or legacy input object.
#' @return A character string representing the lapse treatment (`"no_lapses"`, `"sample"`, or `"marginalize"`).
#' @export
get_lapse_treatment <- S7::new_generic("get_lapse_treatment", "x")

#' Get perceptual noise treatment from a cognitive model
#'
#' @param x A cognitive model or legacy input object.
#' @return A character string representing the noise treatment (`"no_noise"`, `"sample"`, or `"marginalize"`).
#' @export
get_noise_treatment <- S7::new_generic("get_noise_treatment", "x")

#' Get cue labels from a representation, template, or model
#'
#' @param x A representation, representation template, or cognitive model.
#' @param indices An optional integer vector of indices to subset the returned labels.
#' @param ... Additional arguments passed to methods.
#' @return A character vector of cue labels.
#' @export
get_cue_labels <- S7::new_generic("get_cue_labels", "x", function(x, indices = NULL, ...) S7::S7_dispatch())

#' Get category labels from a representation, template, or model
#'
#' Compatibility methods also accept older list/data-frame inputs for
#' transitional support while the S7 API becomes the canonical interface.
#' @param x A representation, representation template, or cognitive model.
#' @param indices An optional integer vector of indices to subset the returned labels.
#' @param ... Additional arguments passed to methods.
#' @return A character vector of category labels.
#' @export
get_category_labels <- S7::new_generic("get_category_labels", "x", function(x, indices = NULL, ...) S7::S7_dispatch())

#' Get response category labels from an object
#'
#' @param x An object (e.g., \code{\link{IdealAdaptorStanfitInput}},
#'   \code{\link{MVBU_Stanfit}}, or \code{\link{MVBU_StanfitPosterior}}).
#' @param indices An optional integer vector of indices to subset the returned labels.
#' @param ... Additional arguments passed to methods.
#' @return A character vector of response category labels.
#' @seealso \code{\link{get_category_labels}}, \code{\link{get_group_labels}}, \code{\link{get_cue_labels}}
#' @export
get_response_category_labels <- S7::new_generic(
  "get_response_category_labels",
  "x",
  function(x, indices = NULL, ...) S7::S7_dispatch()
)

#' Get group labels from a model, stanfit, or stanfit-input object
#'
#' @param x A model, stanfit, or stanfit-input object.
#' @param indices An optional integer vector of indices to subset the returned labels.
#' @param include_prior Logical; whether to include the prior group label. Default: \code{FALSE}.
#'   Only applicable if `x` is an \code{\link{MVBU_Stanfit}} or \code{\link{MVBU_StanfitPosterior}}
#'   object.
#' @param ... Additional arguments passed to methods.
#' @return A character vector of group labels.
#' @seealso \code{\link{get_cue_labels}}, \code{\link{get_category_labels}}, \code{\link{get_labels}}
#' @export
get_group_labels <- S7::new_generic(
  "get_group_labels",
  "x",
  function(x, indices = NULL, include_prior = FALSE, ...) S7::S7_dispatch()
)

#' Get all labels from a representation, template, model, stanfit, or fit-input object
#'
#' @param x A representation, representation template, cognitive model, stanfit,
#'   or fit-input object.
#' @param ... Additional arguments passed to methods.
#' @return A named list of character vectors with elements \code{cue},
#'   \code{category}, \code{response_category}, and \code{group}.
#' @seealso \code{\link{set_labels}}, \code{\link{get_cue_labels}},
#'   \code{\link{get_category_labels}},
#'   \code{\link{get_response_category_labels}},
#'   \code{\link{get_group_labels}}
#' @export
get_labels <- S7::new_generic(
  "get_labels",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' Set label information on an MVBU object or metadata list
#'
#' @param x An [MVBU_Object] or a metadata list.
#' @param cue Character vector of cue labels. Defaults to \code{character(0)}.
#' @param category Character vector of category labels. Defaults to
#'   \code{character(0)}.
#' @param response_category Character vector of response category labels.
#'   Defaults to \code{category}.
#' @param group Character vector of group labels. Defaults to
#'   \code{character(0)}.
#' @param ... Additional arguments passed to methods.
#' @return The modified object or list with updated \code{label_information}.
#' @seealso \code{\link{get_labels}}, \code{\link{get_cue_labels}},
#'   \code{\link{get_category_labels}},
#'   \code{\link{get_response_category_labels}},
#'   \code{\link{get_group_labels}}
#' @export
set_labels <- S7::new_generic("set_labels", "x")


#' Extract the category-likelihood function from a representation, template, or model
#'
#' For NIX, MNIX, and NIW families the returned function evaluates the posterior
#' predictive, which is the category likelihood for those families.
#'
#' @param x A category representation, representation template, or cognitive model.
#' @param ... Additional arguments passed to methods.
#' @return A function of `(new_data, log, noise_treatment, Sigma_noise)` returning category likelihoods.
#' @export
get_category_likelihood_function <- S7::new_generic("get_category_likelihood_function", "x")

#' Extract a category-posterior function from a cognitive model
#'
#' @param x A cognitive model object.
#' @param noise_treatment Optional noise-handling mode for the returned function.
#' @param lapse_treatment Optional lapse-handling mode for the returned function.
#' @param ... Additional arguments passed to methods.
#' @return A function that computes category posterior probabilities for one or more observations.
#' @export
get_category_posterior_function <- S7::new_generic("get_category_posterior_function", c("x", "noise_treatment", "lapse_treatment"))

#' Compute category likelihoods for one or more observations
#'
#' For NIX, MNIX, and NIW families the likelihood is the posterior predictive of
#' the category, so this generic supersedes the family-specific posterior
#' predictive helpers.
#'
#' @param x A category representation, representation template, or cognitive model.
#' @param new_data A numeric matrix of observations, a vector, or a list of observations.
#' @param categories An optional subset of categories to restrict the result to.
#' @param ... Additional arguments passed to methods.
#' @return A matrix of category likelihoods with one row per observation and one column per category.
#' @export
likelihood <- S7::new_generic("likelihood", c("x", "new_data", "categories"))

#' Compute posterior category probabilities for one or more observations
#'
#' @param x A cognitive model object.
#' @param new_data A numeric matrix of observations or a list of matrices.
#' @param categories An optional subset of categories to restrict the posterior to.
#' @param ... Additional arguments passed to methods.
#' @return A matrix of posterior category probabilities with one row per observation.
#' @export
posterior <- S7::new_generic("posterior", c("x", "new_data", "categories"))

#' Categorize one or more observations under the model's decision rule
#'
#' Evaluates the model's response category decisions and response probabilities.
#'
#' @param x A cognitive model object, \code{\link{MVBU_Stanfit}}, or \code{\link{MVBU_StanfitPosterior}}.
#' @param new_data A numeric matrix of observations, data frame, or list of observation matrices.
#' @param decision_rule An optional override for the model's decision rule
#'   (\code{"criterion"}, \code{"sample"} / \code{"sampling"}, or \code{"proportional"}).
#'   Defaults to the model's configured decision rule.
#' @param simplify Logical scalar; if \code{TRUE}, returns a character vector of
#'   chosen \code{response_category} values. If \code{FALSE}, returns a data frame
#'   with columns \code{response_category} and \code{response_probability}.
#'   Defaults to \code{FALSE} when \code{decision_rule == "proportional"} and \code{TRUE} otherwise.
#' @param ... Additional arguments passed to methods.
#' @return If \code{simplify = TRUE}, a character vector of response categories.
#'   If \code{simplify = FALSE}, a data frame containing columns \code{response_category}
#'   and \code{response_probability} (which equals 1 for deterministic or sampled decisions,
#'   and equals the category posterior probability under proportional decision rules).
#' @export
categorize <- S7::new_generic(
  "categorize",
  c("x", "new_data", "decision_rule"),
  fun = function(
    x,
    new_data,
    decision_rule = "criterion",
    simplify = NULL,
    ...
  ) {
    S7::S7_dispatch()
  }
)


#' Update a model's category-representation template from observations
#'
#' Currently implemented for `NIW_IdealAdaptor`; methods for other model types
#' will be added as their update workflows are migrated to S7.
#' @param x A model object.
#' @param observations Observations used for updating.
#' @param ... Additional arguments passed to methods.
#' @return An updated model, or a history of updated models when requested by a method.
#' @export
update_template <- S7::new_generic(
  "update_template",
  c("x", "observations"),
  function(
    x,
    observations,
    updating = c("batch", "incremental"),
    keep_history = FALSE,
    lapse_treatment = "no_lapses",
    noise_treatment = "no_noise",
    update_method = "label-certain",
    ...
  ) {
    S7::S7_dispatch()
  }
)

#' Update a category representation from observations
#'
#' @param x A category representation object.
#' @param x_N Number of observations represented by the sufficient statistics.
#' @param x_mean Mean vector of the represented observations.
#' @param x_SS Centered sum-of-squares matrix of the represented observations.
#' @param ... Additional arguments passed to methods.
#' @return An updated category representation.
#' @export
update_category_representation <- S7::new_generic("update_category_representation", c("x", "x_N", "x_mean", "x_SS"))

get_expected_category <- S7::new_generic("get_expected_category", "x")

#' Get the model type of a cognitive model or model distribution
#'
#' @param x A model object.
#' @param ... Additional arguments passed to methods.
#' @return A character scalar identifying the model family/type.
#' @export
get_model_type <- S7::new_generic("get_model_type", "x")

#' Get the representation type of a category representation
#'
#' @param x A representation object.
#' @param ... Additional arguments passed to methods.
#' @return A character scalar identifying the representation family/type.
#' @export
get_representation_type <- S7::new_generic("get_representation_type", "x")

#' @rdname get_parameters
#' @export
get_expected_mu <- S7::new_generic("get_expected_mu", "x", function(x, ...) S7::S7_dispatch())

#' @rdname get_parameters
#' @export
get_expected_sigma <- S7::new_generic("get_expected_sigma", "x", function(x, ...) S7::S7_dispatch())

#' @rdname get_parameters
#' @export
get_marginal_mu <- S7::new_generic("get_marginal_mu", "x", function(x, ...) S7::S7_dispatch())

#' @rdname get_parameters
#' @export
get_marginal_sigma <- S7::new_generic("get_marginal_sigma", "x", function(x, ...) S7::S7_dispatch())

#' Extract expected category statistics from an S7 model, representation, or fit object
#'
#' Computes the expected value of category parameters such as mean vectors (\eqn{\mu}) or
#' covariance/scatter matrices (\eqn{\Sigma} / \eqn{S}) under the model's distribution.
#'
#' @param x An S7 object (cognitive model, category representation, template, or stanfit).
#' @param statistic Character scalar indicating the statistic to extract: \code{"mu"} (expected mean),
#'   \code{"Sigma"} or \code{"sigma"} (expected covariance), \code{"m"}, \code{"S"}, \code{"kappa"}, or \code{"nu"}.
#' @param categories Optional character vector of category names. Default: \code{NULL} (all categories).
#' @param groups Optional character vector of group names. Default: \code{NULL} (all groups).
#' @param ... Additional arguments passed to methods.
#' @return A numeric vector, matrix, or list of expected values.
#' @seealso \code{\link{get_expected_mu}}, \code{\link{get_expected_sigma}}, \code{\link{get_marginal_category_statistic}}, \code{\link{get_parameters}}
#' @export
get_expected_category_statistic <- S7::new_generic(
  "get_expected_category_statistic",
  "x",
  function(
    x,
    categories = NULL,
    groups = NULL,
    statistic = c("mu", "Sigma"),
    ...
  ) {
    S7::S7_dispatch()
  }
)

#' Extract marginal category statistics (predictive moments of observations) from an S7 object
#'
#' Computes the marginal predictive moments of observable stimulus tokens under the model's
#' predictive distribution, marginalizing over parameter uncertainty within the category belief state.
#'
#' @param x An S7 object (cognitive model, category representation, template, or stanfit).
#' @param statistic Character scalar indicating the statistic to extract: \code{"mu"} (marginal mean),
#'   \code{"Sigma"} or \code{"sigma"} (marginal predictive covariance).
#' @param categories Optional character vector of category names. Default: \code{NULL} (all categories).
#' @param groups Optional character vector of group names. Default: \code{NULL} (all groups).
#' @param ... Additional arguments passed to methods.
#' @return A numeric vector, matrix, or data.frame of marginal values.
#' @seealso \code{\link{get_expected_category_statistic}}, \code{\link{get_marginal_mu}}, \code{\link{get_marginal_sigma}}
#' @export
get_marginal_category_statistic <- S7::new_generic(
  "get_marginal_category_statistic",
  "x",
  function(
    x,
    categories = NULL,
    groups = NULL,
    statistic = c("mu", "Sigma"),
    ...
  ) {
    S7::S7_dispatch()
  }
)

#' Get MCMC prior or posterior draws from a stanfit object
#'
#' Get MCMC draws of all parameters from incremental Bayesian belief-updating (IBBU) as a tibble in long format.
#' By default all post-warmup draws are returned, but if \code{summarize = TRUE} then just the mean of each parameter is returned instead.
#'
#' By default, the category means and scatter matrices are nested, rather than each of their elements being
#' stored separately (\code{nest = TRUE}).
#'
#' @param fit An \code{\link{IdealAdaptorStanfit}} (or \code{\link{MVBU_Stanfit}}) object.
#' @param categories Character vector of category names for which draws should be returned. (default: all categories in model)
#' @param groups Character vector of group names for which draws should be returned. (default: all groups in model including prior)
#' @param which DEPRECATED. Use \code{groups} instead. Should parameters for the prior, posterior, or both be returned? (default: \code{"posterior"})
#' @param ndraws Number of random draws or \code{NULL} if all draws are to be returned. (default: \code{NULL})
#' @param summarize Should the mean of the draws be returned instead of all of the draws? (default: \code{FALSE})
#' @param nest Should the category mean vectors and scatter matrices be nested into one cell each, or should each element
#'   be stored in a separate column? (default: \code{TRUE})
#' @param seed Optional seed for reproducibility when sampling draws. (default: \code{NULL})
#' @param ... Additional arguments passed to methods.
#'
#' @return A tibble with MCMC draws of the model parameters.
#' @seealso \code{\link{get_parameters}}, \code{\link{get_stanfit}}
#' @export
get_draws <- S7::new_generic(
  "get_draws",
  "fit",
  function(
    fit,
    categories = NULL,
    groups = NULL,
    which = "posterior",
    ndraws = NULL,
    summarize = FALSE,
    nest = TRUE,
    seed = NULL,
    ...
  ) {
    S7::S7_dispatch()
  }
)

#' Get original variable names from an S7 model or stanfit input object
#'
#' Retrieves the original variable names that were mapped to standard roles
#' (\code{"group"}, \code{"group_unique"}, \code{"category"}, \code{"response_category"}, and \code{"cues"})
#' when the model or input object was created.
#'
#' @param x An S7 object containing original variable name metadata.
#' @param variable Character vector subsetting the variable types to retrieve.
#'   Must be one or more of \code{"group"}, \code{"group_unique"}, \code{"category"},
#'   \code{"response_category"}, or \code{"cues"}. Defaults to all of them.
#' @param ... Additional arguments passed to methods.
#' @return A named list of character vectors (or a character vector if a single variable is requested) containing the original variable names.
#' @seealso \code{\link{get_data}}, \code{\link{get_exposure_data}}, \code{\link{get_test_data}}
#' @export
get_original_variable_names <- S7::new_generic(
  "get_original_variable_names",
  "x",
  function(x, variable = c("group", "group_unique", "category", "response_category", "cues"), ...) S7::S7_dispatch()
)

#' Get data from an S7 stanfit or stanfit-input object
#'
#' Retrieves underlying experimental data (all observations across phases, exposure
#' observations only, or test observations only) stored in an S7 \code{\link{MVBU_Stanfit}}
#' or \code{\link{IdealAdaptorStanfitInput}} object.
#'
#' @param x An S7 stanfit or stanfit-input object.
#' @param groups Optional character vector of group labels to filter observations by.
#'   Defaults to all available groups.
#' @param categories Optional character vector of category labels to filter exposure
#'   observations by. Defaults to all available category levels.
#' @param response_categories Optional character vector of response category labels
#'   to filter test observations by. Defaults to all available response category levels.
#' @param n_samples Optional integer; subsample up to \code{n_samples} rows. Defaults to \code{NULL} (all observations).
#' @param original_names Logical; if \code{TRUE}, renames standardized columns back to
#'   original variable names. Defaults to \code{FALSE}.
#' @param ... Additional arguments passed to methods.
#' @return A \code{\link[tibble]{tibble}} containing the requested data.
#' @seealso \code{\link{get_original_variable_names}}
#' @rdname get_data
#' @export
get_data <- S7::new_generic(
  "get_data",
  "x",
  function(
    x,
    groups = NULL,
    categories = NULL,
    response_categories = NULL,
    n_samples = NULL,
    original_names = FALSE,
    ...
  ) {
    S7::S7_dispatch()
  }
)

#' @rdname get_data
#' @export
get_exposure_data <- S7::new_generic(
  "get_exposure_data",
  "x",
  function(
    x,
    groups = NULL,
    categories = NULL,
    n_samples = NULL,
    original_names = FALSE,
    ...
  ) {
    S7::S7_dispatch()
  }
)

#' @rdname get_data
#' @export
get_test_data <- S7::new_generic(
  "get_test_data",
  "x",
  function(
    x,
    groups = NULL,
    response_categories = NULL,
    n_samples = NULL,
    original_names = FALSE,
    ...
  ) {
    S7::S7_dispatch()
  }
)


#' Get exposure category statistics from an S7 stanfit or stanfit-input object
#'
#' @param x An S7 stanfit or stanfit-input object.
#' @param statistic Character scalar or vector specifying statistic(s) to extract:
#'   \code{"n"}, \code{"mean"}, \code{"cov"}, \code{"css"}, \code{"uss"}. Default: \code{"mean"}.
#' @param categories Optional character vector of category names. Default: \code{NULL}.
#' @param groups Optional character vector of group names. Default: \code{NULL}.
#' @param untransform_cues Logical; whether to return cues in original space. Default: \code{FALSE}.
#' @param ... Additional arguments passed to methods.
#' @return A tibble, vector, or matrix of exposure category statistics.
#' @seealso \code{\link{get_exposure_category_mean}}, \code{\link{get_exposure_category_cov}}
#' @export
get_exposure_category_statistic <- S7::new_generic(
  "get_exposure_category_statistic",
  "x",
  function(
    x,
    categories = NULL,
    groups = NULL,
    statistic = c("n", "mean", "css", "uss", "cov"),
    untransform_cues = FALSE,
    ...
  ) {
    S7::S7_dispatch()
  }
)

#' @rdname get_exposure_category_statistic
#' @export
get_exposure_category_mean <- S7::new_generic(
  "get_exposure_category_mean",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' @rdname get_exposure_category_statistic
#' @export
get_exposure_category_css <- S7::new_generic(
  "get_exposure_category_css",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' @rdname get_exposure_category_statistic
#' @export
get_exposure_category_uss <- S7::new_generic(
  "get_exposure_category_uss",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' @rdname get_exposure_category_statistic
#' @export
get_exposure_category_cov <- S7::new_generic(
  "get_exposure_category_cov",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' Get the transform/untransform function from an object
#'
#' @param x An S7 fit, fit-input, or transform-information object.
#' @param ... Additional arguments.
#' @return A function.
#' @rdname get_transform_function
#' @export
get_transform_function <- S7::new_generic("get_transform_function", "x")

#' @rdname get_transform_function
#' @export
get_untransform_function <- S7::new_generic("get_untransform_function", "x")

#' Get number of post-warmup MCMC samples from stanfit
#'
#' @param fit A \code{\link[rstan:stanfit-class]{stanfit}} object or \code{\link{MVBU_Stanfit}} object.
#' @param ... Additional arguments passed to methods.
#' @rdname get_number_of_draws
#' @export
get_number_of_draws <- S7::new_generic("get_number_of_draws", "fit")

#' Get indices for random MCMC draws from stanfit
#'
#' @param fit A \code{\link[rstan:stanfit-class]{stanfit}} object or \code{\link{MVBU_Stanfit}} object.
#' @param ndraws Number of indices to be returned.
#' @param ... Additional arguments passed to methods.
#' @rdname get_number_of_draws
#' @export
get_random_draw_indices <- S7::new_generic(
  "get_random_draw_indices",
  "fit",
  function(fit, ndraws = NULL, ...) S7::S7_dispatch()
)


#' Plot category representations and distributions
#'
#' Generic plotting function for visualizing category representations (densities,
#' contour lines, confidence ellipses, exemplar clouds, and 3D density surfaces/ellipsoids)
#' across single category representations, multi-category templates, cognitive decision models,
#' and fitted Stan models.
#'
#' @param x An MVBU object (\code{\link{MVBU_CategoryRepresentation}},
#'   \code{\link{MVBU_CategoryRepresentationTemplate}},
#'   \code{\link{MVBU_CognitiveModel}}, or \code{\link{MVBU_Stanfit}}).
#' @param cues Character vector of cue names to plot (1, 2, or 3 cues). Defaults
#'   to all cues of the model up to 3 (or the first 3 cues if the model has more
#'   than 3 cues, analytically marginalizing out all remaining dimensions). At most
#'   3 cue dimensions can be plotted simultaneously.
#' @param categories Character vector of category names to include. Defaults to
#'   \code{NULL} (all categories).
#' @param aes Plot aesthetic: \code{"contour"}, \code{"fill"},
#'   \code{"fill-gradient"}, \code{"fill-discrete"}, \code{"scatter"}, or
#'   combinations like \code{c("fill-gradient", "contour")}. For 1D and 2D
#'   categories, defaults to \code{c("fill-gradient", "contour")}. For 3D slices,
#'   defaults to \code{"contour"}. For 3D interactive parametric representations,
#'   defaults to \code{"fill-discrete"}. For exemplar representations, defaults to
#'   \code{"scatter"}.
#' @param levels Numeric vector or list specifying density / probability levels.
#'   Defaults to two-tailed central probabilities corresponding to 1, 2, 3 sigma
#'   (\code{2 * stats::pnorm(1:3) - 1}) for 1D/2D and sliced plots, and 2 sigma
#'   (\code{2 * stats::pnorm(2) - 1}) for 3D interactive ellipsoids. For a
#'   bivariate normal distribution \eqn{\mathcal{N}(\boldsymbol{\mu}, \boldsymbol{\Sigma})},
#'   the squared Mahalanobis distance follows \eqn{\chi^2_2}, meaning the ellipse
#'   enclosing mass \eqn{p \in (0, 1)} has radius
#'   \eqn{r = \sqrt{F_{\chi^2_2}^{-1}(p)} = \sqrt{-2 \ln(1 - p)}}.
#' @param limits Optional named list or numeric vector specifying axis / cue
#'   limits (e.g. \code{list(F1 = c(200, 800), F2 = c(800, 2500))}). Defaults to
#'   \code{NULL} (automatically computed from category distributions or
#'   exemplars).
#' @param resolution Integer specifying the grid resolution along continuous cue
#'   dimensions. Defaults to \code{100L} for 1D/2D category plots and \code{60L}
#'   for 3D sliced category plots.
#' @param n_exemplars Integer specifying max number of exemplars to randomly
#'   sample and plot for exemplar representations. Defaults to \code{0L} for 1D/2D
#' (no exemplars plotted) and \code{100L} for 3D interactive scatter plots.
#' @param uncertainty_treatment Character string specifying treatment of posterior parameter
#'   uncertainty in fitted Stan models: \code{"marginalize"} (default) shows the expected/marginal
#'   category density/decision function with posterior credible intervals, whereas \code{"sample"}
#'   overlays individual posterior sample draws.
#' @param groups Character vector of group names to plot for
#'   \code{\link{MVBU_Stanfit}} objects. Defaults to all exposure groups.
#' @param ndraws Integer specifying the number of posterior MCMC draws to use
#'   when extracting parameters or plotting from \code{\link{MVBU_Stanfit}}.
#'   Defaults to \code{20L} (\code{NULL} uses all draws).
#' @param show_exposure_data Logical; if \code{TRUE}, overlays observed exposure
#'   data points on \code{\link{MVBU_Stanfit}} category plots. Defaults to
#'   \code{FALSE}.
#' @param show_test_data Logical; if \code{TRUE}, overlays observed test data
#'   points on \code{\link{MVBU_Stanfit}} category plots. Defaults to
#'   \code{FALSE}.
#' @param parallel Logical; if \code{TRUE}, parallelizes grid evaluations over
#'   multiple CPU chunks using \pkg{parallel}. Defaults to \code{FALSE}.
#' @param n_cores Integer specifying the number of CPU cores to use when
#'   \code{parallel = TRUE}. Defaults to \code{NULL} (which automatically uses
#'   \code{max(1L, parallel::detectCores() - 1L)}).
#' @param noise_treatment Optional character string specifying treatment of
#'   perceptual noise (\code{"no_noise"}, \code{"sample"}, or \code{"marginalize"}).
#'   Defaults to \code{NULL} (uses the model's configured noise treatment).
#' @param ... Additional arguments passed to specific plotting methods (e.g.,
#'   \code{interactive = TRUE} for interactive 3D WebGL scenes via Plotly;
#'   \code{slices} and \code{slice_cue} for 3D sliced category plots).
#' @return A \code{ggplot} object (or a \code{plotly} htmlwidget if
#'   \code{interactive = TRUE}). Because ggplot output is a standard \code{ggplot}
#'   / \pkg{patchwork} object, it can be further customized with additional layers,
#'   scales, themes, or facets.
#' @details
#' In 3D interactive category plots (\code{plot_categories} with \code{interactive = TRUE}),
#' large diamond markers represent category means. For exemplar representations, the
#' mean is computed across all stored exemplars.
#' @seealso \code{\link{plot_categorization_functions}}, \code{\link{plot_parameters}},
#'   \code{\link{likelihood}}, \code{\link{posterior}}
#' @rdname plot_categories
#' @export
plot_categories <- S7::new_generic(
  "plot_categories",
  "x",
  function(
    x,
    cues = NULL,
    categories = NULL,
    groups = NULL,
    aes = NULL,
    levels = NULL,
    limits = NULL,
    resolution = 100,
    n_exemplars = 0L,
    ...
  ) {
    S7::S7_dispatch()
  }
)

#' Plot categorization decision functions and response boundaries
#'
#' Visualizes predicted categorization response probabilities and decision boundaries
#' across cue dimensions for cognitive models and fitted Stan models. Supports
#' proportional probability matching, deterministic criterion/MAP classification, and
#' stochastic sampling decision rules, fully accounting for perceptual noise and lapses.
#'
#' @inheritParams plot_categories
#' @param categories Character vector of category names to evaluate and plot.
#'   Defaults to the first category (\code{get_category_labels(x)[1L]}).
#' @param lapse_treatment Optional character string specifying treatment of
#'   lapses (\code{"no_lapses"}, \code{"sample"}, or \code{"marginalize"}).
#'   Defaults to \code{NULL} (uses the model's configured lapse treatment).
#' @param decision_rule Character string specifying the decision rule to use
#'   for categorization predictions (\code{"proportional"}, \code{"criterion"},
#'   or \code{"sampling"}). Defaults to \code{"proportional"}.
#' @param aes Plot aesthetic: \code{"contour"}, \code{"fill"}, \code{"fill-gradient"},
#'   or \code{"fill-discrete"}. Defaults to \code{"contour"} for 1D/2D static plots,
#'   and \code{"fill-discrete"} for 2D interactive plots (\code{interactive = TRUE}).
#'   When \code{interactive = TRUE} (or \code{aes = "interactive"}), renders a 3D
#'   response probability surface via Plotly.
#' @param levels Numeric vector specifying response probability contour levels.
#'   Defaults to \code{c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99)}.
#' @param resolution Integer specifying grid resolution along continuous cue dimensions.
#'   Defaults to \code{100L} for 1D and \code{60L} for 2D/3D.
#' @return A \code{ggplot} object (or a \code{plotly} htmlwidget if \code{interactive = TRUE}).
#' @seealso \code{\link{plot_categories}}, \code{\link{plot_parameters}},
#'   \code{\link{categorize}}, \code{\link{posterior}}
#' @rdname plot_categorization_functions
#' @export
plot_categorization_functions <- S7::new_generic(
  "plot_categorization_functions",
  "x",
  function(
    x,
    cues = NULL,
    categories = NULL,
    groups = NULL,
    aes = NULL,
    levels = NULL,
    limits = NULL,
    decision_rule = "proportional",
    resolution = 100,
    noise_treatment = NULL,
    lapse_treatment = NULL,
    ...
  ) {
    S7::S7_dispatch()
  }
)

#' Plot model parameter estimates, correlations, and pairwise posterior distributions
#'
#' Visualizes category parameter estimates (means, covariance ellipses, SDs, correlations),
#' parameter correlation matrices across MCMC draws, and pairwise posterior parameter
#' scatter/contour matrices for category representations, templates, cognitive models,
#' and fitted Stan models.
#'
#' @inheritParams plot_categories
#' @param pars Character vector of parameter names to filter in \code{plot_parameters},
#'   \code{plot_parameter_correlations}, and \code{plot_parameters_pairwise}. Defaults
#'   to \code{NULL} (all parameters).
#' @param groups Character vector of group names to plot for \code{\link{MVBU_Stanfit}}
#'   objects. For \code{plot_parameters}, defaults to all groups including the prior.
#'   For \code{plot_parameters_pairwise} and \code{plot_parameter_correlations}, defaults
#'   to \code{"prior"}.
#' @param combine_into_single_plot Logical; if \code{TRUE} (default), combines
#'   multi-panel parameter subplots into a single unified figure via \pkg{patchwork}.
#'   If \code{FALSE}, returns a named list of individual \code{ggplot} panel objects.
#' @param index_panels Logical; if \code{TRUE}, adds panel index labels (A, B, ...)
#'   to multi-panel figures. Defaults to \code{FALSE}.
#' @param ... Additional arguments passed to specific parameter plotting methods.
#' @return A \code{ggplot} object (or list of \code{ggplot} objects if
#'   \code{combine_into_single_plot = FALSE}).
#' @details
#' In parameter plots, scale is displayed as standard deviation
#' \eqn{\tau = \text{SD} = \sqrt{\text{diag}(\boldsymbol{\Sigma})}} (or the
#' expected standard deviation
#' \eqn{\sqrt{\text{diag}(\mathbf{S}) / (\nu - D - 1)}} for NIW/NIX models),
#' while off-diagonal covariances are plotted as dimensionless correlation
#' coefficients \eqn{\rho \in [-1, 1]}.
#'
#' \code{plot_parameter_correlations} displays a correlation heatmap of
#' Pearson correlation coefficients (\eqn{\rho}) computed across MCMC draws
#' for all parameters across categories and selected groups. Parameters for
#' distinct groups form separate rows and columns without averaging across
#' groups. Parameter names indicate category, cue, and group index in their
#' subscripts (with subscript \code{0} denoting the prior group). Gray bounding
#' boxes along the diagonal group parameters belonging to each category.
#'
#' \code{plot_parameters_pairwise} displays a matrix of pairwise parameter
#' relationships across categories and groups. The lower triangle shows
#' sample scatterplots with LOESS trendlines; the diagonal displays marginal
#' posterior auto-densities scaled to panel height; and the upper triangle
#' displays 2D posterior density contour lines (\code{geom_density_2d} with
#' normalized density contours). Points and curves are colored by group.
#' @aliases plot_parameter_correlations plot_parameters_pairwise
#' @seealso \code{\link{plot_categories}}, \code{\link{plot_categorization_functions}},
#'   \code{\link{get_parameters}}, \code{\link{get_draws}}
#' @rdname plot_parameters
#' @export
plot_parameters <- S7::new_generic(
  "plot_parameters",
  "x",
  function(
    x,
    pars = NULL,
    categories = NULL,
    groups = NULL,
    ndraws = 100,
    combine_into_single_plot = TRUE,
    index_panels = FALSE,
    ...
  ) {
    S7::S7_dispatch()
  }
)

#' @rdname plot_parameters
#' @export
plot_parameter_correlations <- S7::new_generic(
  "plot_parameter_correlations",
  "x",
  function(
    x,
    pars = NULL,
    categories = NULL,
    groups = NULL,
    cues = NULL,
    ndraws = 100,
    ...
  ) {
    S7::S7_dispatch()
  }
)

#' @rdname plot_parameters
#' @export
plot_parameters_pairwise <- S7::new_generic(
  "plot_parameters_pairwise",
  "x",
  function(
    x,
    pars = NULL,
    categories = NULL,
    groups = NULL,
    cues = NULL,
    ndraws = 100,
    ...
  ) {
    S7::S7_dispatch()
  }
)

#' Plot empirical sample distributions
#'
#' Visualizes the empirical distribution of cues across categories from observed exposure or test data
#' in a fitted model or stanfit-input object.
#'
#' @param sample Data source to plot: \code{"exposure"}, \code{"test"}, or both \code{c("exposure", "test")} (default).
#' @inheritParams plot_categories
#' @param aes Aesthetic mapping specification. Can be a character vector or a named list with elements
#'   \code{"exposure"} and \code{"test"}. Supported values include \code{"points"} (observed sample points),
#'   \code{"contour"} (empirical density contour lines), \code{"fill-gradient"} (continuous density fill),
#'   and \code{"fill-discrete"} (binned density levels). Defaults to \code{list(exposure = c("fill-gradient", "contour"), test = "points")}.
#' @param densities Density estimation method for contour and fill aesthetics: \code{"gaussian"} (default) evaluates
#'   parametric Gaussian distributions based on sample means and covariances, whereas \code{"kernel"} performs non-parametric
#'   kernel density estimation.
#' @param n_samples Optional integer; maximum number of points to plot when \code{"points"} is included in \code{aes}. Defaults to \code{NULL} (all observations).
#' @param ... Additional arguments passed to methods.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @seealso \code{\link{plot_exposure_sample}}, \code{\link{plot_test_sample}}, \code{\link{plot_categories}}
#' @rdname plot_sample
#' @export
plot_sample <- S7::new_generic(
  "plot_sample",
  "x",
  function(
    x,
    sample = c("exposure", "test"),
    cues = NULL,
    categories = NULL,
    groups = NULL,
    aes = NULL,
    densities = c("gaussian", "kernel"),
    n_samples = NULL,
    levels = NULL,
    limits = NULL,
    resolution = 100,
    ...
  ) {
    S7::S7_dispatch()
  }
)

#' @rdname plot_sample
#' @export
plot_exposure_sample <- function(x, sample = "exposure", ...) {
  plot_sample(x, sample = sample, ...)
}

#' @rdname plot_sample
#' @export
plot_test_sample <- function(x, sample = "test", ...) {
  plot_sample(x, sample = sample, ...)
}

#' Plot Stan MCMC diagnostics
#'
#' @inheritParams plot_categories
#' @rdname plot_diagnostics
#' @export
plot_diagnostics <- S7::new_generic("plot_diagnostics", "x")

#' Evaluate a cognitive model or fitted model
#'
#' Evaluates the predictions of a cognitive model or fitted Stan model against observed categorization responses.
#' When applied to fitted model objects with attached test data (such as \code{\link{MVBU_Stanfit}} or \code{\link{MVBU_StanfitPosterior}}),
#' leaving \code{x = NULL} and \code{response_category = NULL} automatically evaluates the model against the model's attached test data.
#'
#' Supported evaluation methods include:
#' \itemize{
#'   \item \code{"log_lik"}: Trial-by-trial categorical log-likelihood (default):
#'     \deqn{\log L_{\text{cat}} = \sum_{i=1}^N \log P(Y = y_i \mid X = x_i) = \sum_{u=1}^U \sum_{j=1}^K n_{u,j} \log p_{u,j}}
#'     where \eqn{u} indexes unique stimulus locations, \eqn{j} indexes categories, \eqn{n_{u,j}} is the count of responses
#'     for category \eqn{j} at stimulus \eqn{u}, and \eqn{p_{u,j}} is the model's predicted probability of category \eqn{j} at \eqn{u}.
#'     This represents the log probability of observing the exact sequence of individual responses \eqn{(y_1, \dots, y_N)}.
#'     This is the standard likelihood metric in GLMs, logistic regression, and Bayesian Stan models.
#'   \item \code{"accuracy"}: Proportion of correct / matching responses.
#'   \item \code{"likelihood-up-to-constant"}: Deprecated alias for \code{"log_lik"} (supported until version 0.1).
#' }
#'
#' Additionally, the following method is available to support calculation of order invariant log-likelihoods:
#' \itemize{
#'   \item \code{"log_lik_permutation_constant"}: The log multinomial combinatorial permutation factor:
#'     \deqn{C = \sum_{u=1}^U \left( \log(N_u!) - \sum_{j=1}^K \log(n_{u,j}!) \right) = \sum_{u=1}^U \log \binom{N_u}{n_{u,1}, \dots, n_{u,K}}}
#'     where \eqn{N_u = \sum_j n_{u,j}} is the total number of presentations at stimulus location \eqn{u}.
#'     This factor accounts for all permutations of response order at each stimulus location. Crucially, \eqn{C}
#'     depends purely on the observed data and is completely independent of the cognitive model.
#' }
#'
#' @param model A cognitive model object (e.g., \code{\link{MVG_IdealObserver}}, \code{\link{NIW_IdealAdaptor}}) or \code{\link{MVBU_Stanfit}}.
#' @param x Observations/stimuli (matrix, data frame, list of numeric vectors, or numeric vector).
#'   If \code{NULL} and \code{model} has attached test data (e.g., \code{\link{MVBU_Stanfit}} or \code{\link{MVBU_StanfitPosterior}}),
#'   stimuli \code{x} are automatically extracted from the model's attached test data.
#' @param response_category Vector of observed category responses for each observation.
#'   If \code{NULL} and \code{model} has attached test data (e.g., \code{\link{MVBU_Stanfit}} or \code{\link{MVBU_StanfitPosterior}}),
#'   responses are automatically extracted from the model's attached test data column \code{"response"}.
#' @param method Evaluation method(s): \code{"log_lik"} (default), \code{"log_lik_permutation_constant"},
#'   \code{"accuracy"}, or \code{"likelihood-up-to-constant"} (deprecated alias).
#'   Can be a character vector to return multiple metrics.
#' @param decision_rule Decision rule to apply: \code{"proportional"}, \code{"criterion"}, or \code{"sampling"}.
#'   Defaults to \code{"criterion"} when \code{method == "accuracy"}, and \code{"proportional"} otherwise.
#' @param return_by_x Logical; if `TRUE`, returns evaluation metrics grouped by unique stimulus values \eqn{x}.
#'   Defaults to `FALSE`.
#' @param ... Additional arguments passed to methods.
#' @return A numeric scalar, a tibble (if \code{return_by_x = TRUE}), or a named list of these if multiple methods were requested.
#' @rdname evaluate_model
#' @export
evaluate_model <- S7::new_generic("evaluate_model", "model", function(
  model,
  x = NULL,
  response_category = NULL,
  method = "log_lik",
  decision_rule = if (identical(method, "accuracy")) "criterion" else "proportional",
  return_by_x = FALSE,
  ...
) {
  S7::S7_dispatch()
})

#' Sample observations from category representations, templates, or cognitive models
#'
#' Generates simulated cue observations by sampling from an individual category representation,
#' a multi-category template, or a complete cognitive model.
#' \itemize{
#'   \item For single representations (\code{\link{MVBU_CategoryRepresentation}}), samples \eqn{n} observations from the representation's distribution.
#'   \item For templates (\code{\link{MVBU_CategoryRepresentationTemplate}}), samples \eqn{n} total observations distributed uniformly across categories.
#'   \item For cognitive models (\code{\link{MVBU_CognitiveModel}}), samples \eqn{n} total observations distributed across categories proportionally to the model's \code{category_prior}.
#' }
#'
#' @param x An \code{\link{MVBU_CategoryRepresentation}}, \code{\link{MVBU_CategoryRepresentationTemplate}}, or \code{\link{MVBU_CognitiveModel}} object.
#' @param n Non-negative integer specifying the total number of observations to sample. Defaults to \code{1L}.
#' @param with_replacement Logical; for exemplar representations, specifies whether to sample exemplars with replacement (\code{TRUE}, default) or without replacement (\code{FALSE}). Ignored for parametric families.
#' @param randomize_order Logical; whether to randomize the row order of the returned observations. Defaults to \code{TRUE}.
#' @param ... Additional arguments passed to methods (e.g., \code{Ns} or \code{randomize.order} for backwards compatibility).
#' @return A tibble with \eqn{n} rows containing a \code{category} factor column (with levels matching category labels) and one numeric column per cue.
#' @seealso \code{\link{plot_categories}}, \code{\link{likelihood}}, \code{\link{posterior}}
#' @rdname sample_observations
#' @export
sample_observations <- S7::new_generic(
  "sample_observations",
  "x",
  function(x, n = 1L, with_replacement = TRUE, randomize_order = TRUE, ...) {
    S7::S7_dispatch()
  }
)

#' Plot model belief updates across sequence steps or prior-to-posterior transitions
#'
#' Visualizes how category representations or decision surfaces evolve across update
#' steps (for sequential model lists / `update_template()` outputs) or across MCMC
#' prior-to-posterior transitions (for `MVBU_Stanfit` objects).
#'
#' @name plot_model_updates
#' @param x A list of S7 cognitive models/templates or an `MVBU_Stanfit` object.
#' @param what Character string indicating what to plot: `"categories"` (default)
#'   or `"categorization_functions"`.
#' @param step_labels Optional character vector of names or labels for each update step.
#' @param cues Optional character vector of cue names to plot.
#' @param categories Optional character vector of category names to plot.
#' @param aes Plot aesthetic (e.g. `"contour"`, `"fill-gradient"`).
#' @param levels Numeric vector of probability / density levels to plot.
#' @param limits Optional named list or numeric vector of plot axis limits.
#' @param resolution Integer grid resolution for continuous cue dimensions. Defaults to `100`.
#' @param ncol Integer number of facet columns. Defaults to `NULL`.
#' @param groups Optional character vector of group names to plot for `MVBU_Stanfit`.
#' @param ... Additional arguments passed to methods.
#' @return A `ggplot2` plot object faceted across update steps/stages.
#' @export
plot_model_updates <- S7::new_generic(
  "plot_model_updates",
  "x",
  function(
    x,
    what = c("categories", "categorization_functions"),
    groups = NULL,
    step_labels = NULL,
    cues = NULL,
    categories = NULL,
    aes = NULL,
    levels = NULL,
    limits = NULL,
    resolution = 100,
    ncol = NULL,
    ...
  ) {
    S7::S7_dispatch()
  }
)


#' Coerce a Stanfit model object to a lightweight StanfitPosterior representation
#'
#' Extracts posterior draws of category representation parameters and creates a
#' lightweight, self-contained [MVBU_StanfitPosterior] object suitable for rapid
#' categorization and visualization without retaining the full Stan sampling engine.
#'
#' @param x An [MVBU_Stanfit] object.
#' @param ... Additional arguments passed to methods.
#' @return An [MVBU_StanfitPosterior] object.
#' @export
as_MVBU_stanfit_posterior <- S7::new_generic(
  "as_MVBU_stanfit_posterior",
  "x",
  function(x, ...) S7::S7_dispatch()
)


#' Reconstruct exposure update history as a model list
#'
#' Reconstructs sequential cognitive model states or prior-to-posterior update checkpoints
#' into an [MVBU_ModelList].
#'
#' @param object An [MVBU_Stanfit] object or cognitive model.
#' @param uncertainty_treatment Character string specifying treatment of posterior parameter
#'   uncertainty: \code{"marginalize"} (default) reconstructs trajectories across posterior
#'   draws and summarizes them, while \code{"discard"} reconstructs a single trajectory
#'   from the prior point estimate.
#' @param ndraws Number of posterior draws to use when \code{uncertainty_treatment = "marginalize"}.
#'   Defaults to \code{20L}.
#' @param ... Additional arguments passed to methods.
#' @return An [MVBU_ModelList] object.
#' @export
reconstruct_update_history <- S7::new_generic(
  "reconstruct_update_history",
  "object",
  function(
    object,
    uncertainty_treatment = c("marginalize", "discard"),
    ndraws = 20L,
    step_size = 10L,
    groups = NULL,
    categories = NULL,
    parallel = FALSE,
    seed = 42L,
    ...
  ) {
    S7::S7_dispatch()
  }
)

#' Get sufficient category statistics from observations, models, or fit objects
#'
#' Computes or extracts summary exposure statistics (sample size, means, and centered
#' sums of squares / covariance matrices) across categories. Either an object (such as
#' a cognitive model, Stanfit input, Stanfit, or Stanfit posterior) or a data frame is
#' passed as `x`, which determines the method variant dispatched.
#'
#' @param x An object containing exposure data or category representations: a `data.frame`,
#'   an [MVBU_CognitiveModel], an `IdealAdaptorStanfitInput`, an [MVBU_Stanfit], or an [MVBU_StanfitPosterior].
#' @param untransform_cues Logical scalar; if \code{TRUE}, return cue
#'   statistics transformed back into original cue space. Default: \code{FALSE}.
#' @param ... Additional arguments passed to methods.
#' @return A [data.frame] of sufficient statistics across categories (and groups, if applicable).
#'   Columns include `category`, `x_N`, `x_mean`, `x_ss`, `x_css`, and optionally
#'   `group`, `x_uss`, and `x_cov`.
#' @export
get_sufficient_category_statistics <- S7::new_generic(
  "get_sufficient_category_statistics",
  "x",
  function(
    x,
    categories = NULL,
    groups = NULL,
    untransform_cues = FALSE,
    ...
  ) {
    S7::S7_dispatch()
  }
)
