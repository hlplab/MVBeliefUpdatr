#' @import assertthat
#' @import rstan

class_name = "mv_ibbu_stanfit"

setClass(class_name,
        contains = "stanfit")

#' Is this an MV IBBU stanfit?
#'
is.mv_ibbu_stanfit = function(x) {
  if ("stanfit" %in% class(x)) message("Accepting stanfit as valid class. This might change in the future.")
  if (class_name %in% class(x) | "stanfit" %in% class(x))
    return(TRUE) else return(FALSE)
}

#' Is this a tibble of IBBU draws?
#'
#' Check whether \code{x} is a tibble of post-warmup draws of parameters obtained from incremental
#' Bayesian belief-updating (IBBU).
is.mv_ibbu_MCMC = function(x, nested = T, long = T) {
  assert_that(is.flag(long) & long == T,
              msg = "Currently only IBBU MCMC tibbles in long format can be recognized.")
  flag = all(
    is_tibble(x),
    c(".chain", ".iteration", ".draw",
      "group", "category",
      "kappa", "nu", "M", "S", "lapse_rate") %in% names(x))

  if (nested) return(flag) else return(all(flag, c("cue", "cue2") %in% names(x)))
}
