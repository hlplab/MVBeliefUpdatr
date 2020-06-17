#' @import assertthat
#' @import rstan

class_name = "mv_ibbu_stanfit"

setClass(class_name,
        contains = "stanfit")

#' Is this an MV IBBU stanfit?
#'
#' Check whether \code{x} is of class \code{mv_ibbu_stanfit}.
is.mv_ibbu_stanfit = function(x) {
  if ("stanfit" %in% class(x)) message("Accepting stanfit as valid class. This might change in the future.")
  if (class_name %in% class(x) | "stanfit" %in% class(x))
    return(TRUE) else return(FALSE)
}


#' Is this a tibble of NIW beliefs?
#'
#' Check whether \code{x} is a tibble of NIW beliefs.
is.NIW_belief = function(x, is.long = T, category = "category") {
  if (
    any(
      # Currently only IBBU MCMC tibbles in long format can be recognized.
      !is.flag(is.long), !is.long == T,
      !is_tibble(x),
      !all(c(category, "kappa", "nu", "M", "S") %in% names(x))
    )
  ) return(FALSE)

  # Check whether category is a factor only after everything else is checked.
  if (
    any(
      !is.factor(get(category, x))
    )
  ) return(FALSE)

  return(TRUE)
}


#' Is this a tibble of MV IBBU draws?
#'
#' Check whether \code{x} is a tibble of post-warmup draws of parameters obtained from incremental
#' Bayesian belief-updating (IBBU).
is.mv_ibbu_MCMC = function(x, is.nested = T, is.long = T) {
  if(
    any(
       !is.NIW_belief(x, category = "category", is.long = is.long),
       !all(c(".chain", ".iteration", ".draw",
             "group", "lapse_rate") %in% names(x)),
       !(!is.nested | all(c("cue", "cue2") %in% names(x)))
    )
  ) return(F) else return(T)
}



#' Is this a list of MV IBBU inputs?
#'
#' Check whether \code{x} is of class \code{mv_ibbu_input}.
is.mv_ibbu_input = function(x) {
  warning("Test of mv_ibbu_input class not yet implemented. Always returning T.")

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

