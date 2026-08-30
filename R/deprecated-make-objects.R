#' Make exemplar models from data.
#'
#' @description `r lifecycle::badge("deprecated")` Use \code{\link{new_exemplar_category_representation_template_from_data}} instead.
#' @param data The tibble or data.frame from which to construct the exemplar representations.
#' @param group No longer supported; must be `NULL`.
#' @param category Name of variable in \code{data} that contains the category information. (default: "category")
#' @param cues Name(s) of variables in \code{data} that contain the cue information.
#' @param sim_function No longer supported; custom similarity functions are ignored.
#' @param verbose If true provides more information. (default: `FALSE`)
#' @return An \code{MVBU_CategoryRepresentationTemplate} of Exemplar category representations.
#' @keywords internal
#' @export
make_exemplars_from_data <- function(data, group = NULL, category = "category", cues, sim_function = NULL, verbose = F) {
  lifecycle::deprecate_warn("0.0.9", "make_exemplars_from_data()", "new_exemplar_category_representation_template_from_data()")
  if (!is.null(group)) .stop("group is no longer supported. Call this function once per group instead.")
  if (!is.null(sim_function)) .warning("sim_function is no longer supported and will be ignored.")
  new_exemplar_category_representation_template_from_data(data = data, category = category, cues = cues, verbose = verbose)
}

#' Make exemplar model from data.
#'
#' @description `r lifecycle::badge("deprecated")` Use \code{\link{new_exemplar_model_from_data}} instead.
#' @inheritParams make_exemplars_from_data
#' @param prior Optional category-prior vector, forwarded as `category_prior`. (default: uniform over categories)
#' @param lapse_rate Optional lapse rate. (default: 0)
#' @param lapse_bias Optional lapse-bias vector. (default: same as `prior`)
#' @param Sigma_noise Optional perceptual-noise covariance matrix. (default: `NULL`)
#' @return An \code{Exemplar_Model} object.
#' @keywords internal
#' @export
make_exemplar_model_from_data <- function(data, group = NULL, category = "category", cues, sim_function = NULL, prior = NULL, lapse_rate = 0, lapse_bias = NULL, Sigma_noise = NULL, verbose = F) {
  lifecycle::deprecate_warn("0.0.9", "make_exemplar_model_from_data()", "new_exemplar_model_from_data()")
  if (!is.null(group)) .stop("group is no longer supported. Call this function once per group instead.")
  if (!is.null(sim_function)) .warning("sim_function is no longer supported and will be ignored.")
  new_exemplar_model_from_data(
    data = data, category = category, cues = cues,
    category_prior = prior, lapse_rate = lapse_rate, lapse_bias = lapse_bias, Sigma_noise = Sigma_noise,
    verbose = verbose)
}


#' Make multivariate Gaussian representation(s) from data.
#'
#' @description `r lifecycle::badge("deprecated")` Use \code{\link{new_mvg_category_representation_template_from_data}} instead.
#' @inheritParams make_exemplars_from_data
#' @return An \code{MVBU_CategoryRepresentationTemplate} of MVG category representations.
#' @keywords internal
#' @export
make_MVG_from_data <- function(data, group = NULL, category = "category", cues, verbose = F) {
  lifecycle::deprecate_warn("0.0.9", "make_MVG_from_data()", "new_mvg_category_representation_template_from_data()")
  if (!is.null(group)) .stop("group is no longer supported. Call this function once per group instead.")
  new_mvg_category_representation_template_from_data(data = data, category = category, cues = cues, verbose = verbose)
}

#' Make multivariate Gaussian ideal observer(s) from data.
#'
#' @description `r lifecycle::badge("deprecated")` Use \code{\link{new_mvg_ideal_observer_from_data}} instead.
#' @inheritParams make_exemplar_model_from_data
#' @return An \code{MVG_IdealObserver} object.
#' @keywords internal
#' @export
make_MVG_ideal_observer_from_data <- function(data, group = NULL, category = "category", cues, prior = NULL, lapse_rate = 0, lapse_bias = NULL, Sigma_noise = NULL, verbose = F) {
  lifecycle::deprecate_warn("0.0.9", "make_MVG_ideal_observer_from_data()", "new_mvg_ideal_observer_from_data()")
  if (!is.null(group)) .stop("group is no longer supported. Call this function once per group instead.")
  new_mvg_ideal_observer_from_data(
    data = data, category = category, cues = cues,
    category_prior = prior, lapse_rate = lapse_rate, lapse_bias = lapse_bias, Sigma_noise = Sigma_noise,
    verbose = verbose)
}


#' Make NIW belief from data.
#'
#' @description `r lifecycle::badge("deprecated")` Use \code{\link{new_niw_category_representation_template_from_data}} instead.
#' @inheritParams make_exemplars_from_data
#' @param kappa Strength of belief (pseudocount) about the category mean. (default: same as `nu`)
#' @param nu Strength of belief (pseudocount) about the category covariance matrix. (default: number of cues + 2)
#' @return An \code{MVBU_CategoryRepresentationTemplate} of NIW category representations.
#' @keywords internal
#' @export
make_NIW_belief_from_data <- function(data, group = NULL, category = "category", cues, kappa = nu, nu = length(cues) + 2, verbose = F) {
  lifecycle::deprecate_warn("0.0.9", "make_NIW_belief_from_data()", "new_niw_category_representation_template_from_data()")
  if (!is.null(group)) .stop("group is no longer supported. Call this function once per group instead.")
  new_niw_category_representation_template_from_data(data = data, category = category, cues = cues, kappa = kappa, nu = nu, verbose = verbose)
}

#' @rdname make_NIW_belief_from_data
#' @keywords internal
#' @export
make_NIW_prior_from_data <- make_NIW_belief_from_data

#' Make NIW ideal adaptor from data.
#'
#' @description `r lifecycle::badge("deprecated")` Use \code{\link{new_niw_ideal_adaptor_from_data}} instead.
#' @inheritParams make_NIW_belief_from_data
#' @param prior Optional category-prior vector, forwarded as `category_prior`. (default: uniform over categories)
#' @param lapse_rate Optional lapse rate. (default: 0)
#' @param lapse_bias Optional lapse-bias vector. (default: same as `prior`)
#' @param Sigma_noise Optional perceptual-noise covariance matrix. (default: `NULL`)
#' @return A \code{NIW_IdealAdaptor} object.
#' @keywords internal
#' @export
make_NIW_ideal_adaptor_from_data <- function(data, group = NULL, category = "category", cues, kappa = nu, nu = length(cues) + 2, prior = NULL, lapse_rate = 0, lapse_bias = NULL, Sigma_noise = NULL, verbose = F) {
  lifecycle::deprecate_warn("0.0.9", "make_NIW_ideal_adaptor_from_data()", "new_niw_ideal_adaptor_from_data()")
  if (!is.null(group)) .stop("group is no longer supported. Call this function once per group instead.")
  new_niw_ideal_adaptor_from_data(
    data = data, category = category, cues = cues, kappa = kappa, nu = nu,
    category_prior = prior, lapse_rate = lapse_rate, lapse_bias = lapse_bias, Sigma_noise = Sigma_noise,
    verbose = verbose)
}


#' Turn a representation (template) into a cognitive model.
#'
#' @description `r lifecycle::badge("deprecated")` Use \code{\link{new_mvg_ideal_observer}}, \code{\link{new_niw_ideal_adaptor}},
#'   or \code{\link{new_exemplar_model}} with an S7 \code{category_template} instead.
#' @param x An MVG, NIW_belief, or exemplars legacy tibble.
#' @param group No longer supported; must be `NULL`.
#' @param prior Optional category-prior vector. (default: uniform over categories)
#' @param lapse_rate Optional lapse rate. (default: 0)
#' @param lapse_bias Optional lapse-bias vector. (default: same as `prior`)
#' @param Sigma_noise Optional perceptual-noise covariance matrix. (default: `NULL`)
#' @param verbose Ignored.
#' @return An S7 cognitive-model object of the family matching `x`.
#' @keywords internal
#' @export
#' @rdname lift_likelihood_to_model
lift_likelihood_to_model <- function(x, group = NULL, prior = NULL, lapse_rate = 0, lapse_bias = NULL, Sigma_noise = NULL, verbose = F) {
  lifecycle::deprecate_warn(
    "0.0.9", "lift_likelihood_to_model()",
    details = "Use new_mvg_ideal_observer(), new_niw_ideal_adaptor(), or new_exemplar_model() with an S7 category_template instead."
  )
  if (!is.null(group)) .stop("group is no longer supported. Use purrr::map() to call this function once per group instead.")
  representations <-
    if (is.MVG(x)) as_s7_mvg_representations(x)
    else if (is.NIW_belief(x)) as_s7_niw_representations(x)
    else if (is.exemplars(x)) as_s7_exemplar_representations(x)
    else .stop("x must be an MVG, NIW_belief, or exemplars object.")
  constructor <-
    if (is.MVG(x)) new_mvg_ideal_observer
    else if (is.NIW_belief(x)) new_niw_ideal_adaptor
    else new_exemplar_model
  constructor(
    category_template = new_category_representation_template(representations),
    category_prior = prior, lapse_rate = lapse_rate, lapse_bias = lapse_bias, Sigma_noise = Sigma_noise)
}


#' @keywords internal
#' @export
#' @rdname lift_likelihood_to_model
lift_exemplars_to_exemplar_model <- function(x, group = NULL, prior = NULL, lapse_rate = 0, lapse_bias = NULL, Sigma_noise = NULL, verbose = F) {
  lifecycle::deprecate_warn("0.0.9", "lift_exemplars_to_exemplar_model()", "new_exemplar_model()")
  if (!is.null(group)) .stop("group is no longer supported. Use purrr::map() to call this function once per group instead.")
  new_exemplar_model(
    category_template = new_category_representation_template(as_s7_exemplar_representations(x)),
    category_prior = prior, lapse_rate = lapse_rate, lapse_bias = lapse_bias, Sigma_noise = Sigma_noise)
}


#' @keywords internal
#' @export
#' @rdname lift_likelihood_to_model
lift_MVG_to_MVG_ideal_observer <- function(x, group = NULL, prior = NULL, lapse_rate = 0, lapse_bias = NULL, Sigma_noise = NULL, verbose = F) {
  lifecycle::deprecate_warn("0.0.9", "lift_MVG_to_MVG_ideal_observer()", "new_mvg_ideal_observer()")
  if (!is.null(group)) .stop("group is no longer supported. Use purrr::map() to call this function once per group instead.")
  new_mvg_ideal_observer(
    category_template = new_category_representation_template(as_s7_mvg_representations(x)),
    category_prior = prior, lapse_rate = lapse_rate, lapse_bias = lapse_bias, Sigma_noise = Sigma_noise)
}

#' @keywords internal
#' @export
#' @rdname lift_likelihood_to_model
lift_NIW_belief_to_NIW_ideal_adaptor <- function(x, group = NULL, prior = NULL, lapse_rate = 0, lapse_bias = NULL, Sigma_noise = NULL, verbose = F) {
  lifecycle::deprecate_warn("0.0.9", "lift_NIW_belief_to_NIW_ideal_adaptor()", "new_niw_ideal_adaptor()")
  if (!is.null(group)) .stop("group is no longer supported. Use purrr::map() to call this function once per group instead.")
  new_niw_ideal_adaptor(
    category_template = new_category_representation_template(as_s7_niw_representations(x)),
    category_prior = prior, lapse_rate = lapse_rate, lapse_bias = lapse_bias, Sigma_noise = Sigma_noise)
}

#' @keywords internal
#' @export
#' @rdname lift_likelihood_to_model
lift_MVG_ideal_observer_to_NIW_ideal_adaptor <- function(x, group = NULL, kappa, nu, verbose = F) {
  lifecycle::deprecate_warn("0.0.9", "lift_MVG_ideal_observer_to_NIW_ideal_adaptor()", "as_niw_ideal_adaptor()")
  if (!is.null(group)) .stop("group is no longer supported. Use purrr::map() to call this function once per group instead.")
  model <- if (S7::S7_inherits(x, MVG_IdealObserver)) x else as_s7_mvg_ideal_observer(x)
  as_niw_ideal_adaptor(model, kappa = kappa, nu = nu)
}


#' Aggregate models by group structure (Deprecated)
#'
#' @description `r lifecycle::badge("deprecated")` Deprecated in MVBeliefUpdatr 0.0.9.
#'   Use \code{purrr::map()} to create lists of models and \code{\link{aggregate_models}} instead.
#'
#' @param x An MVG, MVG_ideal_observer, NIW_belief, or NIW_ideal_adaptor object.
#' @param group_structure The group structure that will be used for aggregation.
#'
#' @return The aggregated object.
#' @keywords internal
#' @export
aggregate_models_by_group_structure <- function(
  x,
  group_structure = NULL
) {
  lifecycle::deprecate_warn(
    "0.0.9",
    "aggregate_models_by_group_structure()",
    "aggregate_models()",
    details = "Use purrr::map() to create lists of models, and aggregate_models() to combine them."
  )
  .assert_that(all(is.character(group_structure)),
              msg = "Group structure must be a vector of characters.")
  .assert_that(all(group_structure %in% names(x)),
              msg = "All variables in group_structure must be contained in the x.")

  x_names <- setdiff(names(x), group_structure)
  while(length(group_structure) > 0) {
    group_structure <- group_structure[-1]
    x <- x %>%
      dplyr::group_by(!!! rlang::syms(group_structure), .data$category) %>%
      dplyr::summarise(
        dplyr::across(
          dplyr::intersect(names(x), c("kappa", "nu", "prior", "lapse_rate", "lapse_bias")),
          ~ mean(.x, na.rm = TRUE)),
        dplyr::across(
          dplyr::intersect(names(x), c("m", "mu", "S", "Sigma", "Sigma_noise")),
          ~ if (any(unlist(is.null(.x)))) { NULL } else { list(purrr::reduce(.x, `+`) / length(.x)) } ),
        .groups = "drop"
      )
  }

  x %>% dplyr::relocate(!!! rlang::syms(x_names))
}

