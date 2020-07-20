#' @import rstan

class_name = "mv_ibbu_stanfit"

setClass(class_name,
        contains = "stanfit")

#' Is this an MV IBBU stanfit?
#'
#' Check whether \code{x} is of class \code{mv_ibbu_stanfit}.
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
#'
is.mv_ibbu_stanfit = function(x) {
  if ("stanfit" %in% class(x)) message("Accepting stanfit as valid class. This might change in the future.")
  if (class_name %in% class(x) | "stanfit" %in% class(x))
    return(TRUE) else return(FALSE)
}


#' Is this a tibble of NIW beliefs?
#'
#' Check whether \code{x} is a tibble of NIW beliefs. Optionally, one can also check whether a lapse rate
#' and bias is part of the belief.
#'
#' @param x Object to be checked.
#' @param category Name of the category variable. (default: "category")
#' @param is.long Is this check assessing whether the belief is in long format (`TRUE`) or wide format (`FALSE`)?
#' (default: `TRUE`)
#' @param w
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @rdname is_NIW_belief
#' @export
#'
is.NIW_belief = function(x, category = "category", is.long = T, with.lapse = if (with.bias) T else F, with.bias = F) {
  assert_that(all(is.flag(is.long), is.flag(with.lapse), is.flag(with.bias)))

  if (
    any(
      !is.long == T,
      !is_tibble(x),
      !all(c("kappa", "nu", "M", "S") %in% names(x))
    )
  ) {
    warning("Currently only NIW beliefs in long format can be recognized.")
    return(FALSE)
  }

  if (!(category %in% names(x))) {
    warning("Column category not found. Did you use another name for this column? You can use the category
            argument to specify the name of that column.")
    return(FALSE)
  }

  if (
    any(
      with.lapse & "lapse" %nin% names(x),
      with.bias & "bias" %nin% names(x)
    )
  ) return(FALSE)

  # Check that category is a factor only after everything else is checked.
  if (any(!is.factor(get(category, x)))) return(FALSE)

  # Check that M and S contain the cue names and that those cue names match.
  names_M = names(x$M[[1]])
  names_S = dimnames(x$S[[1]])
  if (!all(
    names_S[[1]] == names_S[[2]],
    names_S[[1]] == names_M)) return(FALSE)

  return(TRUE)
}

#' @rdname is_NIW_belief
#' @export
is.NIW_belief_w_lapse = function(x, category = "category", is.long = T, with.bias = F) {
  is.NIW_belief(x, category = category, is.long = is.long, with.lapse = T, with.bias = with.bias)
}

#' @rdname is_NIW_belief
#' @export
is.NIW_belief_w_bias = function(x, category = "category", is.long = T) {
  is.NIW_belief(x, category = category, is.long = is.long, with.lapse = T, with.bias = T)
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
#' @rdname is.NIW_belief_MCMC
#' @export
#'
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

#' @rdname is.NIW_belief_MCMC
#' @export
is.NIW_belief_w_lapse_MCMC = function(x, is.nested = T, is.long = T, with.bias = F) {
  is.NIW_belief_MCMC(x, is.nested = is.nested, is.long = is.long, with.lapse = T, with.bias = with.bias)
}

#' @rdname is.NIW_belief_MCMC
#' @export
is.NIW_belief_w_bias_MCMC = function(x, is.nested = T, is.long = T) {
  is.NIW_belief_MCMC(x, is.nested = is.nested, is.long = is.long, with.lapse = T, with.bias = T)
}


#' Is this a list of MV IBBU inputs?
#'
#' Check whether \code{x} is of class \code{mv_ibbu_input}.
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
#'
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

