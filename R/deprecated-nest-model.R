#' @include S7-core-classes.R
#' @importFrom lifecycle deprecate_warn
#' @importFrom dplyr group_by across any_of arrange summarise relocate select mutate ungroup all_of starts_with
#' @importFrom tidyr unnest pivot_longer
#' @importFrom rlang sym syms .data :=
NULL

# deprecated ------------------------------------------------------------------

# deprecated ------------------------------------------------------------------

#' Deprecated: make_named_vector
#'
#' @description `r lifecycle::badge("deprecated")`
#' `make_named_vector()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0.
#'
#' @param x A vector.
#' @param names Character vector of names.
#' @return A named vector.
#' @seealso [stats::setNames()]
#' @keywords internal
#' @export
make_named_vector <- function(x, names) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "make_named_vector()"
  )
  x <- as.vector(x)
  names(x) <- names
  return(x)
}

#' Deprecated: make_named_square_matrix
#'
#' @description `r lifecycle::badge("deprecated")`
#' `make_named_square_matrix()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0.
#'
#' @param x Values for the matrix.
#' @param names Character vector of row and column names.
#' @return A named square matrix.
#' @keywords internal
#' @export
make_named_square_matrix <- function(x, names) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "make_named_square_matrix()"
  )
  x <- matrix(x, nrow = sqrt(length(x)), dimnames = list(names, names))
  return(x)
}

#' Deprecated: nest_cue_information_in_model
#'
#' @description `r lifecycle::badge("deprecated")`
#' `nest_cue_information_in_model()` and `unnest_cue_information_in_model()` were deprecated
#' in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0. S7 cognitive models store native
#' parameters directly within category representations.
#'
#' @param model A model object or tibble.
#'
#' @keywords internal
#' @rdname nest_model
#' @export
nest_cue_information_in_model <- function(model) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "nest_cue_information_in_model()"
  )

  if (all(c("mu", "Sigma") %in% names(model))) {
    m <- "mu"
    S <- "Sigma"
  } else if (all(c("m", "S") %in% names(model))) {
    m <- "m"
    S <- "S"
  } else if (is.MVG(model) || is.MVG_ideal_observer(model)) {
    m <- "mu"
    S <- "Sigma"
  } else if (is.NIW_belief(model) || is.NIW_ideal_adaptor(model) || is.ideal_adaptor_stanfit(model)) {
    m <- "m"
    S <- "S"
  } else {
    stop("Object not recognized.")
  }

  .assert_that(all(c("cue", "cue2") %in% names(model)),
    msg = "cue and cue2 columns not found. There is nothing to nest."
  )
  model %>%
    dplyr::group_by(dplyr::across(-dplyr::any_of(c("cue", "cue2", m, S)))) %>%
    dplyr::arrange(.data$cue, .data$cue2, .by_group = TRUE) %>%
    dplyr::summarise(
      !!rlang::sym(m) := list(suppressWarnings(make_named_vector(unique(!!rlang::sym(m)), unique(.data$cue)))),
      !!rlang::sym(S) := list(suppressWarnings(make_named_square_matrix(!!rlang::sym(S), unique(.data$cue)))),
      .groups = "drop"
    ) %>%
    dplyr::relocate(dplyr::starts_with(c("lapse_", "prior")), .after = !!rlang::sym(S))
}

#' Deprecated: unnest_cue_information_in_model
#'
#' @description `r lifecycle::badge("deprecated")`
#' `unnest_cue_information_in_model()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0.
#'
#' @keywords internal
#' @rdname nest_model
#' @export
unnest_cue_information_in_model <- function(model) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "unnest_cue_information_in_model()"
  )

  cue <- cue2 <- NULL

  if (all(c("mu", "Sigma") %in% names(model))) {
    m <- "mu"
    S <- "Sigma"
  } else if (all(c("m", "S") %in% names(model))) {
    m <- "m"
    S <- "S"
  } else if (is.MVG(model) || is.MVG_ideal_observer(model)) {
    m <- "mu"
    S <- "Sigma"
  } else if (is.NIW_belief(model) || is.NIW_ideal_adaptor(model) || is.ideal_adaptor_stanfit(model)) {
    m <- "m"
    S <- "S"
  } else {
    stop("Object not recognized.")
  }

  .assert_that(all(c("cue", "cue2") %nin% names(model)),
    msg = "Cannot create cue and cue2 columns since they already exist in the model."
  )

  cue.labels <- if (is.data.frame(model) && m %in% names(model) && length(model[[m]]) > 0) {
    names(model[[m]][[1]])
  } else {
    NULL
  }
  if (is.null(cue.labels) && is.data.frame(model) && S %in% names(model) && length(model[[S]]) > 0) {
    cue.labels <- colnames(model[[S]][[1]])
  }
  if (is.null(cue.labels)) {
    cue.labels <- get_cue_labels_from_model(model)
  }

  model <- model %>%
    tidyr::unnest(c(!!rlang::sym(m), !!rlang::sym(S))) %>%
    dplyr::group_by(dplyr::across(-dplyr::any_of(c(m, S)))) %>%
    dplyr::mutate(cue = cue.labels)

  for (i in seq_along(cue.labels)) {
    model <- model %>% dplyr::mutate(!!rlang::sym(cue.labels[i]) := (!!rlang::sym(S))[, i])
  }

  model %>%
    dplyr::select(-dplyr::all_of(S)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(cue.labels), values_to = S, names_to = "cue2") %>%
    dplyr::ungroup() %>%
    dplyr::relocate(dplyr::all_of(c("cue", "cue2")), .after = dplyr::any_of("nu"))
}

