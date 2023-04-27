#' @importFrom dplyr near

#' @export
get_class <- function(x) {
  class <-
    case_when(
      is.NIW_ideal_adaptor(x) ~ "NIW_ideal_adaptor",
      is.NIW_belief(x) ~ "NIW_belief",
      is.MVG_ideal_observer(x) ~ "MVG_ideal_observer",
      is.MVG(x) ~ "MVG",
      is.NIW_ideal_adaptor_stanfit(x) ~ new_stanfit_class_name,
      T ~ "Unrecognized class")

  return(class)
}

#' @export
is.Sigma <- function(x) {
  if (is.matrix(x)) {
    if (all(x == 0) | is.positive.definite(x)) return(T) else return(F)
  } else {
    if (is_scalar_double(x)) return(T) else return(F)
  }
}


get_expected_columns_for_model <- function() c("prior", "lapse_rate", "lapse_bias", "Sigma_noise")

#' Is this a model?
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
#' @examples
#' TBD
#' @export
is.model <- function(x, group = NULL, verbose = F, tolerance = 1e-5) {
  name_of_x <- deparse(substitute(x))

  if (verbose) message("Checking whether ", name_of_x, " has all the column names required for a model.")
  if (!all(c("prior", "lapse_rate", "lapse_bias", "Sigma_noise") %in% names(x))) return(FALSE)

  if (!is.null(group)) {
    if (verbose) message("Checking whether ", name_of_x, " is a model within each unique combination of group values.")
    x %<>% group_by(!!! syms(group))
  }

  # Check that noise is constant across categories (within each group)
  if (any(x %>% distinct(Sigma_noise) %>% tally(name = "n.sigma_noise") %>% filter(n.sigma_noise != 1) %>% nrow() != 0)) {
    if (verbose) message(paste("Noise covariance matrix Sigma_noise in", name_of_x, "is not constant across categories."))
    return(FALSE)
  }

  # Check that all entries of noise are either NULL or a matrix
  if (!all(map(x$Sigma_noise, ~ if (is.null(.x) | is.matrix(.x)) { TRUE } else { FALSE }) %>% unlist())) {
    if (verbose) message("If not NULL, Sigma_noise must be a matrix.")
    return(FALSE)
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
