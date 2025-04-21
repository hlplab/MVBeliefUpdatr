#' @importFrom dplyr near tally distinct pull
NULL

#' @export
get_class <- function(x) {
  class <-
    case_when(
      is.NIW_ideal_adaptor(x) ~ "NIW_ideal_adaptor",
      is.NIW_belief(x) ~ "NIW_belief",
      is.MVG_ideal_observer(x) ~ "MVG_ideal_observer",
      is.MVG(x) ~ "MVG",
      is.exemplar_model(x) ~ "exemplar_model",
      is.exemplars(x) ~ "exemplars",
      is.ideal_adaptor_stanfit(x) ~ "ideal_adaptor_stanfit",
      T ~ "Unrecognized class")

  return(class)
}

#' @export
is.Sigma <- function(x) {
  if (is.null(x)) stop2("Expected a covariance matrix, but got NULL.")
  if (is.matrix(x)) {
    if (all(x == 0) | is.positive.definite(x)) return(T) else return(F)
  } else {
    if (is_scalar_double(x)) return(T) else return(F)
  }
}


get_expected_columns_for_model <- function() c("prior", "lapse_rate", "lapse_bias", "Sigma_noise")

#' Is this an MVBeliefUpdatr representation?
#'
#' Check whether \code{x} is recognized as an MVBeliefUpdatr category representation.
#'
#' @param x Object to be checked.
#' @param group Name of one or more group variables, each unique combination of which describes a model. (default: NULL)
#' @param verbose Should verbose output be provided? (default: `TRUE`)
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#'
#' @importFrom dplyr group_map
#' @export
is.MVBU_representation <- function(x, group = NULL, verbose = F, tolerance = 1e-5) {
  name_of_x <- deparse(substitute(x))

  if (!is.null(group)) {
    if (verbose) message("Checking whether ", name_of_x, " is a model within each unique combination of group values.")
    x %<>% group_by(!!! syms(group))
  }

  if (any(unlist(group_map(x, .f = ~ !(is.exemplars(.x) || is.MVG(.x) || is.NIW_belief(.x)))))) {
    return(FALSE)
  }

  return(TRUE)
}

#' Is this an MVBeliefUpdatr model?
#'
#' Check whether \code{x} is a model with lapse rates, biases, priors, and perceptual noise.
#'
#' @param x Object to be checked.
#' @param group Name of one or more group variables, each unique combination of which describes a model. (default: NULL)
#' @param verbose Should verbose output be provided? (default: `TRUE`)
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#'
#' @importFrom purrr map_lgl
#' @export
is.MVBU_model <- function(x, group = NULL, verbose = F, tolerance = 1e-5) {
  name_of_x <- deparse(substitute(x))

  if (!is_tibble(x)) {
    if (verbose) message("Object is not a tibble. All MVBeliefUpdatr models are stored in tibbles.")
    return(FALSE)
  }

  if (verbose) message("Checking whether ", name_of_x, " has all the column names required for a model.")
  if (!all(c("prior", "lapse_rate", "lapse_bias", "Sigma_noise") %in% names(x))) return(FALSE)

  if (!is.null(group)) {
    if (verbose) message("Checking whether ", name_of_x, " is a model within each unique combination of group values.")
    x %<>% group_by(!!! syms(group))
  }

  # Check that all entries of noise are either NULL or a matrix
  if (!all(map_lgl(x$Sigma_noise, ~ if (is.null(.x) | is.matrix(.x)) { TRUE } else { FALSE }))) {
    if (verbose) message("If not NULL, Sigma_noise must be a matrix.")
    return(FALSE)
  }

  # Check that noise is constant across categories (within each group)
  if (any(x %>% distinct(Sigma_noise) %>% tally(name = "n.sigma_noise") %>% filter(n.sigma_noise != 1) %>% nrow() != 0)) {
    if (verbose) message(paste("Noise covariance matrix Sigma_noise in", name_of_x, "is not constant across categories."))
    return(FALSE)
  }

  # Since noise is tested to be constant, it is sufficient to test the details of the first Sigma_noise
  # If noise is not NULL, check its dimensionality and dimension names
  if (!is.null(first(x$Sigma_noise))) {
    Sigma_noise <- first(x$Sigma_noise)
    assert_that(all(dim(Sigma_noise) == rep(get_cue_dimensionality_from_model(x), 2)),
                msg = paste("If not NULL, Sigma_noise must match the dimensionality of other parameters in the model (here: a",
                            get_cue_dimensionality_from_model(x), "x", get_cue_dimensionality_from_model(x), " matrix)."))
    assert_that(!is.null(dimnames(Sigma_noise)),
                msg = "If not NULL, Sigma_noise must have non-NULL dimnames.")
    assert_that(map(dimnames(Sigma_noise), ~ .x == get_cue_labels_from_model(x)) %>% reduce(all),
                msg = "If not NULL, the dimnames of Sigma_noise must match the cue names used in the model.")
  }

  # Check that the prior probabilities are all between 0 and 1
  if (any(!between(x$prior, 0, 1))) {
    if (verbose) message(paste("Prior probabilities in", name_of_x, "are not all between 0 and 1: ", paste(x$prior, collapse = ",")))
    return(FALSE)
  }

  # Check that the prior probabilities add up to 1
  if (any(!near(x %>% summarise(sum_prior = sum(prior)) %>% pull(sum_prior), 1, tol = tolerance))) {
    if (verbose) message(paste("Prior probabilities in", name_of_x, "do not add up to 1: ", sum(x$prior)))
    return(FALSE)
  }

  # Check that the lapse rates are all between 0 and 1
  if (any(!between(x$lapse_rate, 0, 1))) {
    if (verbose) message(paste("Lapse rates in", name_of_x, "are not all between 0 and 1: ", paste(x$lapse_rate, collapse = ",")))
    return(FALSE)
  }

  # Check that the lapse rate is constant across categories
  if (any(x %>% distinct(lapse_rate) %>% tally(name = "n.lapse_rate") %>% filter(n.lapse_rate != 1) %>% nrow() != 0)) {
    if (verbose)
      message(
        paste(
          "Lapse rates in",
          name_of_x,
          "are not constant across categories: ",
          paste(x$lapse_rate, collapse = ", ")))
    return(FALSE)
  }

  # Check that the lapse bias probabilities are all between 0 and 1
  if (any(!between(x$lapse_bias, 0, 1))) {
    if (verbose) message(paste("Lapse bias probabilities in", name_of_x, "are not all between 0 and 1: ", paste(x$lapse_bias, collapse = ",")))
    return(FALSE)
  }

  # Check that the lapse bias probabilities add up to 1
  if (any(!near(x %>% summarise(sum_lapse_bias = sum(lapse_bias)) %>% pull(sum_lapse_bias), 1, tol = tolerance))) {
    if (verbose) message(paste("Lapse bias probabilities in", name_of_x, "do not add up to 1: ", sum(x$lapse_bias)))
    return(FALSE)
  }

  return(TRUE)
}

#' Print MVBeliefUpdatr model
#'
#' Specifies reasonable defaults for the parameters to be summarized for the MVBeliefUpdatr_model object.
#'
#' @param x An \code{\link{MVBeliefUpdatr_model}} object.
#'
#' @export
print.MVBU_model <- function(x, ...) {
  assertthat(is.MVBU_model(x))

  if (get_cue_dimensionality_from_model(x) == 1)
    x %<>%
    mutate(across(is.numeric, ~ unlist(.x)))

  print(x)
}
