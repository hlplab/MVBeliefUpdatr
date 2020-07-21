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
      all(!is_tibble(x), is.data.frame(x))
    )
  ) {
    message("Currently only NIW beliefs in long format can be recognized.")
    return(FALSE)
  }

  if (category %nin% names(x)) {
    message("Column category not found. Did you use another name for this column? You can use the category
            argument to specify the name of that column.")
    return(FALSE)
  }

  if (
    any(
      any(c("kappa", "nu", "M", "S") %nin% names(x)),
      with.lapse & "lapse_rate" %nin% names(x),
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
