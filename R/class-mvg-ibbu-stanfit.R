#' @import rstan

new_stanfit_class_name = "mvg_ibbu_stanfit"

#' An S4 class for stanfit objects that use one of the mvg_ibbu stan programs.
#'
#' @slot input_data A list
#' @slot labels A list
mvg_ibbu_stanfit <- setClass(new_stanfit_class_name,
         slots = c(input_data = "list", labels = "list"),
         contains = "stanfit",
         package = "MVBeliefUpdatr")

#' Make an mvg_ibbu_stanfit
#'
#' Combines a stanfit object and its input into an mvg_ibbu_stanfit object.
#'
#' @param stanfit stanfit object
#' @param input Input for a call to rstan prepared for one of the acceptable
#' stan programs. See \code{\link{compose_data}}.
#'
#' @return An mvg_ibbu_stanfit object
#'
#' @seealso \code{\link[rstan]{stanfit}}
#' @examples
#' TBD
#' @export
#'
as.mvg_ibbu_stanfit = function(stanfit, input) {
  assert_that(class(stanfit) == "stanfit",
              msg = paste0("Only stanfit objects can be converted into ", class_name, " objects."))
  assert_that(stanfit@stanmodel@model_name %in% names(MVBeliefUpdatr:::stanmodels),
              msg = paste0("stanfit object was not created by one of the accepted stancodes: code is ",
                           stanfit@stanmodel@model_name, " but would have to be one of ", names(MVBeliefUpdatr:::stanmodels)))

  class(stanfit) <- new_stanfit_class_name
  stanfit %<>%
    attach_stanfit_input_data(input)

  return(stanfit)
}

#' Is this an MV IBBU stanfit?
#'
#' Check whether \code{x} is of class \code{mvg_ibbu_stanfit}.
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
is.mvg_ibbu_stanfit = function(x) {
  if (all(c(new_stanfit_class_name, "stanfit")  %in% class(x)))
    return(TRUE) else return(FALSE)
}


#' Is this a tibble of MCMC draws of NIW beliefs?
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
is.NIW_belief_MCMC = function(x, is.nested = T, is.long = T, with.lapse = if (with.bias) T else F, with.bias = F) {
  if(
    all(
       is.NIW_belief(x, is.long = is.long, category = "category", with.lapse = with.lapse, with.bias = with.bias),
       all(c(".chain", ".iteration", ".draw",
             "group", "lapse_rate") %in% names(x)),
       xor(is.nested, all(c("cue", "cue2") %in% names(x)))
    )
  ) return(T) else return(F)
}

#' @describeIn is.NIW_belief_MCMC Also checks whether x has a lapse term.
#' @export
is.NIW_belief_w_lapse_MCMC = function(x, is.nested = T, is.long = T, with.bias = F) {
  is.NIW_belief_MCMC(x, is.nested = is.nested, is.long = is.long, with.lapse = T, with.bias = with.bias)
}

#' @describeIn is.NIW_belief_MCMC Also checks whether x has a bias term.
#' @export
is.NIW_belief_w_bias_MCMC = function(x, is.nested = T, is.long = T) {
  is.NIW_belief_MCMC(x, is.nested = is.nested, is.long = is.long, with.lapse = T, with.bias = T)
}


#' Is this a list of MV IBBU inputs?
#'
#' Check whether \code{x} is of class \code{mvg_ibbu_input}.
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
is.mvg_ibbu_input = function(x) {
  message("Test of mvg_ibbu_input class not yet implemented. Always returning T.")

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

