#' @include S7-core-classes.R
NULL

# -------------------------
# Core generics
# -------------------------

construct_mvbu <- S7::new_generic("construct_mvbu", "x")
validate_mvbu <- S7::new_generic("validate_mvbu", "x")
plot_prep_mvbu <- S7::new_generic("plot_prep_mvbu", "x")
get_model_family <- S7::new_generic("get_model_family", "x")
get_metadata <- S7::new_generic("get_metadata", "x")

#' Extract parameters or expected parameters from representations, templates, models, or model distributions.
#'
#' @name get_parameters
#' @title Extract parameters or expected category statistics
#' @param x A category representation, template, cognitive model, or stanfit object.
#' @param statistic Character scalar indicating the statistic to extract (e.g. \code{"mu"} or \code{"Sigma"}).
#' @param ... Additional arguments passed to methods.
#' @return A list, vector, matrix, or data frame of parameters or expected statistics.
#' @export
get_parameters <- S7::new_generic("get_parameters", "x", function(x, ...) S7::S7_dispatch())

#' @rdname get_parameters
#' @export
get_parameter_names <- S7::new_generic("get_parameter_names", "x", function(x, ...) S7::S7_dispatch())


get_category_representations <- S7::new_generic("get_category_representations", "x")
get_category_template <- S7::new_generic("get_category_template", "x")


#' Get the stored Stan fit from an S7 fit object.
#'
#' @param x An S7 fit object.
#' @return The stored [rstan::stanfit] object when present.
#' @export
get_stanfit <- S7::new_generic("get_stanfit", "x")

#' Set the stored Stan fit on an S7 fit object.
#'
#' @param x An S7 fit object.
#' @param stanfit An [rstan::stanfit] object.
#' @return The updated S7 fit object.
#' @export
set_stanfit <- S7::new_generic("set_stanfit", c("x", "stanfit"))

#' Extract control parameters of the NUTS sampler
#'
#' @param x An S7 fit object or stanfit object.
#' @param pars Optional parameter names.
#' @param ... Additional arguments.
#' @return A named list of control parameters.
#' @export
control_params <- S7::new_generic("control_params", "x", function(x, ...) S7::S7_dispatch())

#' Extract log posterior draws
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
nuts_params <- S7::new_generic("nuts_params", "x", function(x, ...) {
  if (inherits(x, "stanfit") || inherits(x, "CmdStanFit")) {
    bayesplot::nuts_params(x, ...)
  } else {
    S7::S7_dispatch()
  }
})

#' Extract Rhat diagnostic values
#'
#' @param x An S7 fit object, draws object, or array.
#' @param pars Optional parameter names.
#' @param ... Additional arguments.
#' @return Named numeric vector of Rhat values.
#' @rdname diagnostic-quantities
#' @export
rhat <- S7::new_generic("rhat", "x", function(x, ...) {
  if (inherits(x, c("draws", "stanfit", "CmdStanFit")) || is.matrix(x) || is.array(x) || is.numeric(x)) {
    posterior::rhat(x, ...)
  } else {
    S7::S7_dispatch()
  }
})

#' Extract effective sample size ratio
#'
#' @param x An S7 fit object or stanfit object.
#' @param pars Optional parameter names.
#' @param ... Additional arguments.
#' @return Named numeric vector of Neff ratios.
#' @rdname diagnostic-quantities
#' @export
neff_ratio <- S7::new_generic("neff_ratio", "x", function(x, ...) {
  if (inherits(x, "stanfit") || inherits(x, "CmdStanFit")) {
    bayesplot::neff_ratio(x, ...)
  } else {
    S7::S7_dispatch()
  }
})


#' Get the stored Stan input from an S7 fit or fit-input object.
#'
#' @param x An S7 fit or fit-input object.
#' @return The stored S7 Stan input object.
#' @export
get_staninput <- S7::new_generic("get_staninput", "x")

#' Set the stored Stan input on an S7 fit or fit-input object.
#'
#' @param x An S7 fit or fit-input object.
#' @param staninput An S7 Stan input object.
#' @return The updated S7 object.
#' @export
set_staninput <- S7::new_generic("set_staninput", c("x", "staninput"))

#' Get transform metadata from an S7 fit or fit-input object.
#'
#' @param x An S7 fit or fit-input object.
#' @return A transform-information object.
#' @export
get_transform_information <- S7::new_generic("get_transform_information", "x")

#' Get category-prior values from a cognitive model or legacy input.
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

#' Get lapse-rate values from a cognitive model or legacy input.
#'
#' Compatibility methods handle older list/data-frame inputs while the S7 API
#' becomes the standard interface.
#' @param x A cognitive model or legacy input object.
#' @return A numeric lapse-rate value.
#' @export
get_lapse_rate <- S7::new_generic("get_lapse_rate", "x")

#' Get lapse-bias values from a cognitive model or legacy input.
#'
#' Compatibility methods handle older list/data-frame inputs while the S7 API
#' becomes the standard interface.
#' @param x A cognitive model or legacy input object.
#' @param categories Optional category labels used to resolve values.
#' @return A numeric vector of lapse-bias values.
#' @export
get_lapse_bias <- S7::new_generic("get_lapse_bias", c("x", "categories"))

#' Get perceptual noise covariance matrix from a cognitive model.
#'
#' @param x A cognitive model or legacy input object.
#' @return A matrix or numeric value representing perceptual noise covariance (\eqn{\Sigma_{\text{noise}}}), or `NULL` if no noise is specified.
#' @export
get_noise <- S7::new_generic("get_noise", "x")

#' Get cue labels from a representation, template, model, or model distribution.
#'
#' @param x A representation, representation template, or cognitive model.
#' @param indices An optional integer vector of indices to subset the returned labels.
#' @param ... Additional arguments passed to methods.
#' @return A character vector of cue labels.
#' @export
get_cue_labels <- S7::new_generic("get_cue_labels", "x", function(x, indices = NULL, ...) S7::S7_dispatch())

#' Get category labels from a representation, template, model, or model distribution.
#'
#' Compatibility methods also accept older list/data-frame inputs for
#' transitional support while the S7 API becomes the canonical interface.
#' @param x A representation, representation template, or cognitive model.
#' @param indices An optional integer vector of indices to subset the returned labels.
#' @param ... Additional arguments passed to methods.
#' @return A character vector of category labels.
#' @export
get_category_labels <- S7::new_generic("get_category_labels", "x", function(x, indices = NULL, ...) S7::S7_dispatch())

#' Get group labels from a model, stanfit, or fit-input object.
#'
#' @param x A model, stanfit, or fit-input object.
#' @param indices An optional integer vector of indices to subset the returned labels.
#' @param ... Additional arguments passed to methods (e.g., `include_prior = TRUE`).
#' @return A character vector of group labels.
#' @seealso \code{\link{get_cue_labels}}, \code{\link{get_category_labels}}, \code{\link{get_labels}}
#' @export
get_group_labels <- S7::new_generic("get_group_labels", "x", function(x, indices = NULL, ...) S7::S7_dispatch())

#' Get all labels from a representation, template, model, stanfit, or fit-input object.
#'
#' @param x A representation, representation template, cognitive model, stanfit, or fit-input object.
#' @param ... Additional arguments passed to methods.
#' @return A named list of character vectors with elements \code{cue}, \code{category}, and \code{group}.
#' @seealso \code{\link{get_cue_labels}}, \code{\link{get_category_labels}}, \code{\link{get_group_labels}}
#' @export
get_labels <- S7::new_generic("get_labels", "x", function(x, ...) S7::S7_dispatch())


#' Extract the category-likelihood function from a representation, template, or model.
#'
#' For NIX, MNIX, and NIW families the returned function evaluates the posterior
#' predictive, which is the category likelihood for those families.
#'
#' @param x A category representation, representation template, or cognitive model.
#' @return A function of `(new_data, log, noise_treatment, Sigma_noise)` returning category likelihoods.
#' @export
get_category_likelihood_function <- S7::new_generic("get_category_likelihood_function", "x")

#' Extract a category-posterior function from a cognitive model.
#'
#' @param x A cognitive model object.
#' @param noise_treatment Optional noise-handling mode for the returned function.
#' @param lapse_treatment Optional lapse-handling mode for the returned function.
#' @return A function that computes category posterior probabilities for one or more observations.
#' @export
get_category_posterior_function <- S7::new_generic("get_category_posterior_function", c("x", "noise_treatment", "lapse_treatment"))

#' Compute category likelihoods for one or more observations.
#'
#' For NIX, MNIX, and NIW families the likelihood is the posterior predictive of
#' the category, so this generic supersedes the family-specific posterior
#' predictive helpers.
#'
#' @param x A category representation, representation template, or cognitive model.
#' @param new_data A numeric matrix of observations, a vector, or a list of observations.
#' @param categories An optional subset of categories to restrict the result to.
#' @return A matrix of category likelihoods with one row per observation and one column per category.
#' @export
likelihood <- S7::new_generic("likelihood", c("x", "new_data", "categories"))

#' Compute posterior category probabilities for one or more observations.
#'
#' @param x A cognitive model object.
#' @param new_data A numeric matrix of observations or a list of matrices.
#' @param categories An optional subset of categories to restrict the posterior to.
#' @return A matrix of posterior category probabilities with one row per observation.
#' @export
posterior <- S7::new_generic("posterior", c("x", "new_data", "categories"))

#' Categorize one or more observations under the model's decision rule.
#'
#' @param x A cognitive model object.
#' @param new_data A numeric matrix of observations or a list of matrices.
#' @param decision_rule An optional override for the model's decision rule.
#' @return A data frame with one row per observation and the winning category plus its probability.
#' @export
categorize <- S7::new_generic("categorize", c("x", "new_data", "decision_rule"))

#' Update a model's category-representation template from observations.
#'
#' Currently implemented for `NIW_IdealAdaptor`; methods for other model types
#' will be added as their update workflows are migrated to S7.
#' @param x A model object.
#' @param observations Observations used for updating.
#' @return An updated model, or a history of updated models when requested by a method.
#' @export
update_template <- S7::new_generic("update_template", c("x", "observations"))

#' Update a category representation from observations.
#'
#' @param x A category representation object.
#' @param x_N Number of observations represented by the sufficient statistics.
#' @param x_mean Mean vector of the represented observations.
#' @param x_SS Centered sum-of-squares matrix of the represented observations.
#' @return An updated category representation.
#' @export
update_category_representation <- S7::new_generic("update_category_representation", c("x", "x_N", "x_mean", "x_SS"))

get_expected_category <- S7::new_generic("get_expected_category", "x")

#' Get the model type of a cognitive model or model distribution.
#'
#' @param x A model object.
#' @return A character scalar identifying the model family/type.
#' @export
get_model_type <- S7::new_generic("get_model_type", "x")

#' Get the representation type of a category representation.
#'
#' @param x A representation object.
#' @return A character scalar identifying the representation family/type.
#' @export
get_representation_type <- S7::new_generic("get_representation_type", "x")

#' @rdname get_parameters
#' @export
get_expected_mu <- S7::new_generic("get_expected_mu", "x", function(x, ...) S7::S7_dispatch())

#' @rdname get_parameters
#' @export
get_expected_sigma <- S7::new_generic("get_expected_sigma", "x", function(x, ...) S7::S7_dispatch())

#' Extract expected category statistics from an S7 model, representation, or fit object.
#'
#' Computes the expected value of category parameters such as mean vectors (\eqn{\mu}) or
#' covariance/scatter matrices (\eqn{\Sigma} / \eqn{S}) under the model's distribution.
#'
#' @param x An S7 object (cognitive model, category representation, template, or stanfit).
#' @param statistic Character scalar indicating the statistic to extract: \code{"mu"} (expected mean),
#'   \code{"Sigma"} or \code{"sigma"} (expected covariance), \code{"m"}, \code{"S"}, \code{"kappa"}, or \code{"nu"}.
#' @param ... Additional arguments passed to methods.
#' @return A numeric vector, matrix, or list of expected values.
#' @seealso \code{\link{get_expected_mu}}, \code{\link{get_expected_sigma}}, \code{\link{get_parameters}}
#' @export
get_expected_category_statistic <- S7::new_generic(
  "get_expected_category_statistic",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' Get MCMC prior or posterior draws from a stanfit object
#'
#' Get MCMC draws of all parameters from incremental Bayesian belief-updating (IBBU) as a tibble. Both wide
#' (\code{wide = TRUE}) or long format (\code{wide = FALSE}) can be chosen as output. By default all post-warmup draws are
#' returned, but if \code{summarize = TRUE} then just the mean of each parameter is returned instead.
#'
#' By default, the category means and scatter matrices are nested, rather than each of their elements being
#' stored separately (\code{nest = TRUE}).
#'
#' @param fit An \code{\link{IdealAdaptorStanfit}} (or \code{\link{MVBU_Stanfit}}) object.
#' @param categories Character vector of category names for which draws should be returned. (default: all categories in model)
#' @param groups Character vector of group names for which draws should be returned. (default: all groups in model including prior)
#' @param which DEPRECATED. Use \code{groups} instead. Should parameters for the prior, posterior, or both be returned? (default: \code{"posterior"})
#' @param ndraws Number of random draws or \code{NULL} if all draws are to be returned. (default: \code{NULL})
#' @param untransform_cues Should m_0 and S_0 be transformed back into the original cue space? (default: \code{FALSE})
#' @param summarize Should the mean of the draws be returned instead of all of the draws? (default: \code{FALSE})
#' @param wide Should all parameters be returned in one row? (default: \code{FALSE})
#' @param nest Should the category mean vectors and scatter matrices be nested into one cell each, or should each element
#'   be stored in a separate column? (default: \code{TRUE})
#' @param seed Optional seed for reproducibility when sampling draws. (default: \code{NULL})
#' @param ... Additional arguments passed to methods.
#'
#' @return A tibble with MCMC draws of the model parameters.
#' @seealso \code{\link{get_parameters}}, \code{\link{get_stanfit}}
#' @export
get_draws <- S7::new_generic("get_draws", "fit", function(fit, ...) S7::S7_dispatch())

#' Get the data from an S7 stanfit or stanfit-input object.
#'
#' @param x An S7 stanfit or stanfit-input object.
#' @param ... Additional arguments passed to methods.
#' @return A tibble containing the underlying data.
#' @seealso \code{\link{get_exposure_data}}, \code{\link{get_test_data}}
#' @export
get_data <- S7::new_generic("get_data", "x", function(x, ...) S7::S7_dispatch())

#' Get the exposure data from an S7 stanfit or stanfit-input object.
#'
#' @param x An S7 stanfit or stanfit-input object.
#' @param ... Additional arguments passed to methods.
#' @return A tibble containing the exposure data.
#' @seealso \code{\link{get_data}}, \code{\link{get_test_data}}
#' @export
get_exposure_data <- S7::new_generic(
  "get_exposure_data",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' Get the test data from an S7 stanfit or stanfit-input object.
#'
#' @param x An S7 stanfit or stanfit-input object.
#' @param groups Optional character vector of group labels to filter the test data by.
#' @param .recover_from_staninput Logical scalar; if \code{TRUE} and the \code{@data} property
#'   is missing or empty, attempts to reconstruct test data from \code{@staninput}. (default: \code{FALSE})
#' @param ... Additional arguments passed to methods.
#' @return A tibble containing the test data.
#' @seealso \code{\link{get_data}}, \code{\link{get_exposure_data}}
#' @export
get_test_data <- S7::new_generic(
  "get_test_data",
  "x",
  function(x, ...) S7::S7_dispatch()
)

#' Get exposure category statistics from an S7 stanfit or stanfit-input object.
#'
#' @param x An S7 stanfit or stanfit-input object.
#' @param ... Additional arguments passed to methods.
#' @return A tibble, vector, or matrix of exposure category statistics.
#' @seealso \code{\link{get_exposure_category_mean}}, \code{\link{get_exposure_category_cov}}
#' @export
get_exposure_category_statistic <- S7::new_generic(
  "get_exposure_category_statistic",
  "x",
  function(x, ...) S7::S7_dispatch()
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

#' Get the transform/untransform function from an object.
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
#' @param fit A \code{\link{stanfit}} object or \code{\link{MVBU_Stanfit}} object.
#' @rdname get_number_of_draws
#' @export
get_number_of_draws <- S7::new_generic("get_number_of_draws", "fit")

#' Get indices for random MCMC draws from stanfit
#'
#' @param fit A \code{\link{stanfit}} object or \code{\link{MVBU_Stanfit}} object.
#' @param ndraws Number of indices to be returned.
#' @rdname get_number_of_draws
#' @export
get_random_draw_indices <- S7::new_generic("get_random_draw_indices", "fit")


#' Plot functions for MVBU objects
#'
#' Generic plotting functions for MVBU objects, including cognitive models,
#' category representations, templates, stanfits, and staninput objects.
#'
#' @param x An MVBU object (\code{\link{MVBU_CategoryRepresentation}},
#'   \code{\link{MVBU_CategoryRepresentationTemplate}},
#'   \code{\link{MVBU_CognitiveModel}}, or \code{\link{MVBU_Stanfit}}).
#' @param cues Character vector of cue names to plot (1D or 2D).
#' @param categories Character vector of category names. For
#'   \code{plot_categories}, defaults to \code{NULL} (all categories). For
#'   \code{plot_categorization_function}, defaults to the first category
#'   (\code{get_category_labels(x)[1L]}).
#' @param aes Plot aesthetic: \code{"contour"}, \code{"fill-gradient"},
#'   \code{"fill-discrete"}, or combinations like
#'   \code{c("fill-gradient", "contour")}. Defaults to \code{"contour"} for 2D
#'   plots. \code{"fill"} defaults to \code{"fill-gradient"} in all cases.
#'   \code{"fill-discrete"} is only available for parametric representations.
#' @param levels Numeric vector or list specifying density / probability
#'   levels. For \code{plot_categories}, defaults to two-tailed central
#'   probabilities corresponding to 1, 2, 3, 4 sigma
#'   (\code{2 * stats::pnorm(1:4) - 1}). For
#'   \code{plot_categorization_function}, defaults to
#'   \code{c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99)}.
#' @param limits Optional named list or numeric vector specifying axis / cue
#'   limits (e.g. \code{list(F1 = c(200, 800), F2 = c(800, 2500))}). Defaults to
#'   \code{NULL} (computed from category statistics or exemplars).
#' @param groups Character vector of group names to plot for
#'   \code{\link{MVBU_Stanfit}} objects. For \code{plot_parameters_pairwise}
#'   and \code{plot_parameter_correlations}, defaults to \code{"prior"}
#'   (available groups include \code{"prior"} and all exposure groups).
#'   For \code{plot_parameters}, defaults to all groups, including the prior.
#' @param pars Character vector of parameter names to filter in
#'   \code{plot_parameters}, \code{plot_parameter_correlations}, and
#'   \code{plot_parameters_pairwise}. Defaults to \code{NULL} (all
#'   parameters).
#' @param combine_into_single_plot Logical; if \code{TRUE} (default), combines
#'   multi-panel parameter subplots into a single plot via \pkg{patchwork}.
#' @param n_exemplars Integer specifying max number of exemplars to randomly
#'   sample and plot for exemplar representations. Defaults to 100 (\code{NULL}
#'   plots all exemplars).
#' @param ndraws Number of posterior draws to use when plotting from
#'   \code{\link{MVBU_Stanfit}}. Defaults to 100 (\code{NULL} uses all draws).
#' @param ... Additional arguments passed to specific plotting methods.
#' @return A \code{ggplot} object (or list of \code{ggplot} objects if
#'   \code{combine_into_single_plot = FALSE}). Because output is a standard
#'   \code{ggplot} / \pkg{patchwork} object, it can be further customized with
#'   additional layers, scales, themes, or facets (e.g.
#'   \code{+ facet_wrap(~ Category)}).
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
#' @seealso \code{\link{categorize}}, \code{\link{posterior}},
#'   \code{\link{get_parameters}}, \code{\link{get_draws}}
#' @rdname plot_mvbu
#' @export
plot_categories <- S7::new_generic("plot_categories", "x")

#' @rdname plot_mvbu
#' @export
plot_categorization_function <- S7::new_generic(
  "plot_categorization_function",
  "x"
)

#' @rdname plot_mvbu
#' @export
plot_parameters <- S7::new_generic("plot_parameters", "x")

#' @rdname plot_mvbu
#' @export
plot_parameter_correlations <- S7::new_generic(
  "plot_parameter_correlations",
  "x"
)

#' @rdname plot_mvbu
#' @export
plot_parameters_pairwise <- S7::new_generic(
  "plot_parameters_pairwise",
  "x"
)

#' @rdname plot_mvbu
#' @export
plot_cues <- S7::new_generic("plot_cues", "x")

#' @rdname plot_mvbu
#' @export
plot_diagnostics <- S7::new_generic("plot_diagnostics", "x")

#' Evaluate a cognitive model or fitted model
#'
#' Evaluates the predictions of a cognitive model or fitted Stan model against observed categorization responses.
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
#' @param response_category Vector of observed category responses for each observation.
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
) S7::S7_dispatch())
