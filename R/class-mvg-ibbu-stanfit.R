new_stanfit_class_name = "NIW_ideal_adaptor_stanfit"

#' An S4 class for stanfit objects that use one of the NIW_ideal_adaptor stan programs.
#'
#' @slot input_data list containing the data handed to rstan through \code{compose_data} function.
#' @slot transform_functions list containing elements transform.function and untransform.function.
#' @slot labels list
#' @export
NIW_ideal_adaptor_stanfit <-
  setClass(
    new_stanfit_class_name,
    slots = c(input_data = "list", transform_functions = "list", labels = "list"),
    contains = "stanfit",
    package = "MVBeliefUpdatr")

# Call class constructor function
NIW_ideal_adaptor_stanfit


#' Coerce stanfit into NIW_ideal_adaptor_stanfit
#'
#' Combines a \code{\link[rstan]{stanfit}} object and its input into an \code{\link{NIW_ideal_adaptor_stanfit}} object.
#'
#' @param stanfit stanfit object
#' @param input Input for a call to rstan prepared for one of the acceptable
#' stan programs. See \code{\link{compose_data}}.
#' @param transform_functions Optionally, a list of transform (and corresponding untransform) functions of the
#' type returned by \code{\link[transform_cues]{transform_cues}}.
#'
#' @return NIW_ideal_adaptor_stanfit object
#'
#' @seealso \code{\link[rstan]{stanfit}}
#' @examples
#' TBD
#' @export
as.NIW_ideal_adaptor_stanfit <- function(stanfit, input_data, transform_functions = NULL) {
  assert_that(class(stanfit) %in% c("stanfit", "NIW_ideal_adaptor_stanfit"),
              msg = paste0("Only stanfit and NIW_ideal_adaptor_stanfit objects can be converted into ", new_stanfit_class_name, " objects."))
  assert_that(stanfit@model_name %in% names(MVBeliefUpdatr:::stanmodels),
              msg = paste0("stanfit object was not created by one of the accepted stancodes:\n\t",
                           paste(names(MVBeliefUpdatr:::stanmodels), collapse = "\n\t"),
                           "\n(you can get the name of your model from your_stanfit@model_name)."))

  class(stanfit) <- new_stanfit_class_name
  stanfit %<>% attach_stanfit_input_data(input_data)

  if (!is.null(transform_functions)) stanfit %<>% attach_stanfit_transform(transform_functions)

  return(stanfit)
}

#' Is this an NIW IBBU stanfit?
#'
#' Check whether \code{x} is of class \code{NIW_ideal_adaptor_stanfit}.
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
is.NIW_ideal_adaptor_stanfit <- function(x, verbose = F) {
  if (all(class(x) %in% c("stanfit", new_stanfit_class_name)))
    return(TRUE) else return(FALSE)
}


#' Is this a tibble of MCMC draws of an NIW ideal adaptor?
#'
#' Check whether \code{x} is a tibble of post-warmup draws of parameters obtained from incremental
#' conjugate Bayesian belief-updating (IBBU) over a Normal-Inverse-Wishart (NIW) prior.
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
is.NIW_ideal_adaptor_MCMC <- function(x, is.nested = T, is.long = T, with.lapse = if (with.lapse_bias) T else F, with.lapse_bias = F) {
  if(
    all(
       is.NIW_ideal_adaptor(x, is.long = is.long, category = "category", with.lapse = with.lapse, with.lapse_bias = with.lapse_bias),
       all(c(".chain", ".iteration", ".draw",
             "group", "lapse_rate") %in% names(x)),
       xor(is.nested, all(c("cue", "cue2") %in% names(x)))
    )
  ) return(T) else return(F)
}


#' Is this a list of NIW IBBU inputs?
#'
#' Check whether \code{x} is of class \code{NIW_ideal_adaptor_input}.
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
is.NIW_ideal_adaptor_input <- function(x) {
  message("Test of NIW_ideal_adaptor_input class not yet implemented. Always returning T.")

  # Proposed names for slides in input object (at least internally / not necessarily handed to stan like this:
  #
  #   exposure_N (N)
  #   exposure_cue_mean (x_mean)
  #   exposure_cue_ss (x_ss)
  #   test_N (N_test)
  #   test_cue (x_test)
  #   test_response (z_test_counts)
  #   test_group (y_test)

  return(TRUE)
}

