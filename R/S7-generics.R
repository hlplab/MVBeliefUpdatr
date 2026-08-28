#' @include S7-core-classes.R
NULL

# -------------------------
# Core generics
# -------------------------

construct_mvbu <- S7::new_generic("construct_mvbu", "x")
validate_mvbu <- S7::new_generic("validate_mvbu", "x")
summarize_mvbu <- S7::new_generic("summarize_mvbu", "x")
print_mvbu <- S7::new_generic("print_mvbu", "x")
plot_prep_mvbu <- S7::new_generic("plot_prep_mvbu", "x")

get_model_family <- S7::new_generic("get_model_family", "x")
get_metadata <- S7::new_generic("get_metadata", "x")
get_parameters <- S7::new_generic("get_parameters", "x")

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

#' Get cue labels from a representation, template, model, or model distribution.
#'
#' @param x A representation, representation template, or cognitive model.
#' @param indices An optional integer vector of indices to subset the returned labels.
#' @return A character vector of cue labels.
#' @export
get_cue_labels <- S7::new_generic("get_cue_labels", c("x", "indices"))

#' Get category labels from a representation, template, model, or model distribution.
#'
#' Compatibility methods also accept older list/data-frame inputs for
#' transitional support while the S7 API becomes the canonical interface.
#' @param x A representation, representation template, or cognitive model.
#' @param indices An optional integer vector of indices to subset the returned labels.
#' @return A character vector of category labels.
#' @export
get_category_labels <- S7::new_generic("get_category_labels", c("x", "indices"))
get_group_labels <- S7::new_generic("get_group_labels", c("x", "indices"))


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

plot_categories <- S7::new_generic("plot_categories", "x")
plot_parameters <- S7::new_generic("plot_parameters", "x")
plot_diagnostics <- S7::new_generic("plot_diagnostics", "x")
