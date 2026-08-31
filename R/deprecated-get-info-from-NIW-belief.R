#' @include S7-core-classes.R
#' @importFrom lifecycle deprecate_warn
#' @importFrom foreach foreach %do%
#' @importFrom dplyr filter rename_with mutate
#' @importFrom tidyr as_tibble pivot_wider unnest
#' @importFrom purrr reduce
#' @importFrom rlang := sym
NULL

# deprecated ------------------------------------------------------------------

#' Deprecated: get_NIW_posterior_predictive.pmap
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_NIW_posterior_predictive.pmap()` is deprecated. Use \code{\link{get_NIW_posterior_predictive}} or \code{\link{likelihood}} instead.
#'
#' @inheritParams get_NIW_posterior_predictive
#' @param ... Additional arguments.
#' @export
get_NIW_posterior_predictive.pmap <- function(x, m, S, kappa, nu, ...) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_NIW_posterior_predictive.pmap()",
    with = "get_NIW_posterior_predictive()"
  )
  get_NIW_posterior_predictive(x = x, m = m, S = S, kappa = kappa, nu = nu, ...)
}

#' Deprecated: get_posterior_predictive_from_NIW_belief
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_posterior_predictive_from_NIW_belief()` is deprecated. Use \code{\link{likelihood}} instead.
#'
#' @param x Observations.
#' @param model An NIW belief or adaptor model.
#' @param noise_treatment Noise treatment.
#' @param log Logical; whether log likelihood is returned.
#' @param category Category column name.
#' @param category.label Category labels.
#' @param wide Logical; whether wide format is returned.
#' @return Posterior predictive data frame.
#' @rdname get_posterior_predictive_from_NIW_belief
#' @export
get_posterior_predictive_from_NIW_belief <- function(
  x,
  model,
  noise_treatment = "no_noise",
  log = TRUE,
  category = "category",
  category.label = NULL,
  wide = FALSE
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_posterior_predictive_from_NIW_belief()",
    with = "likelihood()"
  )

  if (S7::S7_inherits(model, MVBU_CognitiveModel)) {
    d <- length(get_cue_labels(model))
    x_mat <- .as_observation_matrix(x, d = d, arg_name = "x")
    lik <- likelihood(model, x_mat, categories = category.label)
    if (log) lik <- log(lik)
    cats <- colnames(lik)
    val_name <- if (log) "log_posterior_predictive" else "posterior_predictive"
    if (wide) {
      out <- tibble::as_tibble(as.data.frame(lik), .name_repair = "minimal")
      names(out) <- paste0(val_name, ".", cats)
      return(out)
    }
    out <- tibble::tibble(
      !!rlang::sym(val_name) := as.vector(lik),
      !!rlang::sym(category) := rep(cats, each = nrow(lik))
    )
    return(out)
  }

  .assert_that(is.NIW_belief(model))
  .assert_optional_character(category.label)

  if (is.null(category.label)) {
    model %<>% droplevels()
    category.label <- model %>%
      dplyr::pull(!!rlang::sym(category)) %>%
      unique()
  }

  posterior_predictive <- foreach::foreach(c = category.label) %do% {
    m <- model %>% dplyr::filter(!!rlang::sym(category) == c)
    get_NIW_posterior_predictive(
      x = x,
      m = m$m[[1]],
      S = m$S[[1]],
      kappa = m$kappa[[1]],
      nu = m$nu[[1]],
      log = log,
      noise_treatment = noise_treatment,
      Sigma_noise = if (noise_treatment == "no_noise") NULL else m$Sigma_noise[[1]]
    ) %>%
      tibble::as_tibble(.name_repair = "minimal") %>%
      dplyr::rename_with(~ if (log) "log_posterior_predictive" else "posterior_predictive") %>%
      dplyr::mutate(!!rlang::sym(category) := c)
  }

  posterior_predictive %<>% purrr::reduce(rbind)
  if (wide) {
    posterior_predictive %<>%
      tidyr::pivot_wider(
        values_from = if (log) "log_posterior_predictive" else "posterior_predictive",
        names_from = !!rlang::sym(category),
        names_prefix = if (log) "log_posterior_predictive." else "posterior_predictive."
      ) %>%
      tidyr::unnest()
  }

  return(posterior_predictive)
}

#' Deprecated: get_posterior_predictives_from_NIW_beliefs
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_posterior_predictives_from_NIW_beliefs()` is deprecated. Use \code{\link{likelihood}} instead.
#'
#' @param grouping.var Grouping variable name.
#' @rdname get_posterior_predictive_from_NIW_belief
#' @export
get_posterior_predictives_from_NIW_beliefs <- function(
  x,
  model,
  noise_treatment = "no_noise",
  log = TRUE,
  category = "category",
  category.label = NULL,
  grouping.var = NULL,
  wide = FALSE
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_posterior_predictives_from_NIW_beliefs()",
    with = "likelihood()"
  )

  if (is.null(grouping.var)) {
    return(suppressWarnings(get_posterior_predictive_from_NIW_belief(
      x,
      model,
      log = log,
      category = category,
      category.label = category.label,
      wide = wide
    )))
  } else {
    .assert_that(grouping.var %in% names(x),
      msg = "Grouping variable not found in the NIW belief object."
    )

    foreach::foreach(i = unique(x[[grouping.var]])) %do% {
      suppressWarnings(get_posterior_predictive_from_NIW_belief(
        x,
        model %>% dplyr::filter(!!rlang::sym(grouping.var) == i),
        log = log,
        category = category,
        category.label = category.label,
        wide = wide
      )) %>%
        dplyr::mutate(!!rlang::sym(grouping.var) := i)
    } %>%
      purrr::reduce(rbind)
  }
}
