#' @include S7-classes.R
NULL

# -------------------------
# Core generics
# -------------------------

construct_mvbu <- S7::new_generic("construct_mvbu", "x")
validate_mvbu <- S7::new_generic("validate_mvbu", "x")
summarize_mvbu <- S7::new_generic("summarize_mvbu", "x")
print_mvbu <- S7::new_generic("print_mvbu", "x")
categorize_mvbu <- S7::new_generic("categorize_mvbu", c("x", "new_data"))
predict_mvbu <- S7::new_generic("predict_mvbu", c("x", "new_data"))
posterior_mvbu <- S7::new_generic("posterior_mvbu", "x")
plot_prep_mvbu <- S7::new_generic("plot_prep_mvbu", "x")

get_model_family <- S7::new_generic("get_model_family", "x")
get_category_likelihood <- S7::new_generic("get_category_likelihood", "x")
get_parameters <- S7::new_generic("get_parameters", "x")
get_category_prior <- S7::new_generic("get_category_prior", "x")
get_cue_labels <- S7::new_generic("get_cue_labels", "x")
get_category_labels <- S7::new_generic("get_category_labels", "x")
get_group_labels <- S7::new_generic("get_group_labels", "x")

#' Get posterior category probabilities for one or more observations.
#'
#' @param x A cognitive model object.
#' @param new_data A numeric matrix of observations or a list of matrices.
#' @return A matrix of posterior category probabilities with one row per observation.
#' @export
get_category_posterior <- S7::new_generic("get_category_posterior", c("x", "new_data"))

#' Get the winning category and its probability under the model's decision rule.
#'
#' The meaning of "winning" depends on the model's decision rule or an explicit
#' override supplied by the caller.
#' @param x A cognitive model object.
#' @param new_data A numeric matrix of observations or a list of matrices.
#' @return A data frame with one row per observation and the winning category plus its probability.
#' @export
get_category <- S7::new_generic("get_category", c("x", "new_data"))

#' Get both posterior category probabilities and category predictions.
#'
#' @param x A cognitive model object.
#' @param new_data A numeric matrix of observations or a list of matrices.
#' @return A list containing the posterior matrix and the category prediction.
#' @export
get_category_posterior_prediction <- S7::new_generic("get_category_posterior_prediction", c("x", "new_data"))

update_model <- S7::new_generic("update_model", c("x", "data"))
update_category_likelihood <- S7::new_generic("update_category_likelihood", c("x", "data"))

get_posterior <- S7::new_generic("get_posterior", "x")
get_expected_category <- S7::new_generic("get_expected_category", "x")

plot_categories <- S7::new_generic("plot_categories", "x")
plot_parameters <- S7::new_generic("plot_parameters", "x")
plot_diagnostics <- S7::new_generic("plot_diagnostics", "x")