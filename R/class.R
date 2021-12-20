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
    if (is.positive.definite(x)) return(T) else return(F)
  } else {
    if (is_scalar_double(x)) return(T) else return(F)
  }
}
