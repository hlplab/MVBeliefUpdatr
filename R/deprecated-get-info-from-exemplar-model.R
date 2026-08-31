#' @include S7-core-classes.R
#' @importFrom lifecycle deprecate_warn
#' @importFrom dplyr group_by mutate ungroup select arrange rename pull
#' @importFrom rlang :=
#' @importFrom stats rbinom runif
NULL

# deprecated ------------------------------------------------------------------

#' Deprecated: get_likelihood_from_exemplars
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_likelihood_from_exemplars()` is deprecated. Use \code{\link{likelihood}} instead.
#'
#' @param x Observations.
#' @param model Exemplar model.
#' @param noise_treatment Noise treatment.
#' @param log Logical; whether log likelihood is returned.
#' @param category Category column name.
#' @param category.label Category labels.
#' @return Likelihood data frame.
#' @rdname get_likelihood_from_exemplars
#' @export
get_likelihood_from_exemplars <- function(
  x,
  model,
  noise_treatment = "no_noise",
  log = TRUE,
  category = "category",
  category.label = NULL
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_likelihood_from_exemplars()",
    with = "likelihood()"
  )

  if (S7::S7_inherits(model, MVBU_CognitiveModel)) {
    d <- length(get_cue_labels(model))
    x_mat <- .as_observation_matrix(x, d = d, arg_name = "x")
    lik <- likelihood(model, x_mat, categories = category.label)
    if (log) lik <- log(lik)
    cats <- colnames(lik)
    val_name <- if (log) "log_likelihood" else "likelihood"
    return(tibble::tibble(
      !!rlang::sym(val_name) := as.vector(lik),
      !!rlang::sym(category) := rep(cats, each = nrow(lik))
    ))
  }

  tibble::tibble()
}

#' Deprecated: get_categorization_from_exemplar_model
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_categorization_from_exemplar_model()` is deprecated. Use \code{\link{categorize}} instead.
#'
#' @param x Observations.
#' @param model Model object.
#' @param decision_rule Decision rule.
#' @param noise_treatment Noise treatment.
#' @param lapse_treatment Lapse treatment.
#' @param simplify Logical; whether to simplify output.
#' @return Categorization data frame or vector.
#' @rdname get_categorization_from_model
#' @export
#' @keywords internal
get_categorization_from_exemplar_model <- function(
  x,
  model,
  decision_rule,
  noise_treatment = if (decision_rule == "sampling") "sample" else "no_noise",
  lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  simplify = FALSE
) {
  lifecycle::deprecate_warn(
    "0.0.3",
    "get_categorization_from_exemplar_model()",
    with = "categorize()"
  )

  if (S7::S7_inherits(model, MVBU_CognitiveModel)) {
    d.response <-
      .legacy_long_posterior(model, x, noise_treatment, lapse_treatment) %>%
      .legacy_apply_decision_rule(decision_rule)

    if (simplify) return(.legacy_simplify_categorization(d.response, decision_rule)) else return(d.response)
  }

  .assert_that(is.exemplar_model(model))
  .assert_that(decision_rule %in% c("criterion", "proportional", "sampling"),
    msg = "Decision rule must be one of: 'criterion', 'proportional', or 'sampling'."
  )
  .assert_that(any(lapse_treatment %in% c("no_lapses", "sample", "marginalize")),
    msg = "lapse_treatment must be one of 'no_lapses', 'sample' or 'marginalize'."
  )

  if (!is.list(x)) x <- list(x)

  posterior_probabilities <-
    suppressWarnings(get_likelihood_from_exemplars(x = x, model = model, log = FALSE, noise_treatment = noise_treatment)) %>%
    dplyr::group_by(.data$category) %>%
    dplyr::mutate(
      observationID = seq_along(x),
      x = x,
      lapse_rate = get_lapse_rate(model),
      lapse_bias = get_lapse_bias(model, categories = .data$category),
      prior = get_category_prior(model, categories = .data$category)
    ) %>%
    dplyr::group_by(.data$observationID) %>%
    dplyr::mutate(posterior_probability = (.data$likelihood * .data$prior) / sum(.data$likelihood * .data$prior))

  if (lapse_treatment == "sample") {
    posterior_probabilities %<>%
      dplyr::mutate(
        posterior_probability = ifelse(
          rep(stats::rbinom(1, 1, .data$lapse_rate), length(get_category_labels(model))),
          .data$lapse_bias,
          .data$posterior_probability
        )
      )
  } else if (lapse_treatment == "marginalize") {
    posterior_probabilities %<>%
      dplyr::mutate(posterior_probability = .data$lapse_rate * .data$lapse_bias + (1 - .data$lapse_rate) * .data$posterior_probability)
  }

  if (decision_rule == "criterion") {
    posterior_probabilities %<>%
      dplyr::mutate(
        posterior_probability = ifelse(
          rep(sum(.data$posterior_probability == max(.data$posterior_probability)) > 1, length(get_category_labels(model))),
          .data$posterior_probability + stats::runif(length(get_category_labels(model)), min = 0, max = 0),
          .data$posterior_probability
        ),
        response = ifelse(.data$posterior_probability == max(.data$posterior_probability), 1, 0)
      )
  } else if (decision_rule == "sampling") {
    posterior_probabilities %<>%
      dplyr::mutate(response = .rmultinom(1, 1, .data$posterior_probability) %>% as.vector())
  } else if (decision_rule == "proportional") {
    posterior_probabilities %<>%
      dplyr::mutate(response = .data$posterior_probability)
  }

  posterior_probabilities %<>%
    dplyr::ungroup() %>%
    dplyr::select(-dplyr::any_of(c("likelihood", "posterior_probability"))) %>%
    dplyr::select(dplyr::all_of(c("observationID", "x", "category", "response")))

  if (simplify) {
    .assert_that(decision_rule %in% c("criterion", "sampling"),
      msg = "For simplify = T, decision rule must be either criterion or sampling."
    )
    return(posterior_probabilities %>%
      dplyr::filter(.data$response == 1) %>%
      dplyr::select(dplyr::all_of(c("observationID", "category"))) %>%
      dplyr::arrange(.data$observationID) %>%
      dplyr::rename(response = "category") %>%
      dplyr::ungroup() %>%
      dplyr::pull(.data$response))
  } else {
    return(posterior_probabilities)
  }
}
