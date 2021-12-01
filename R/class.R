#' @export
get_class <- function(x) {
  class <-
    case_when(
      is.NIW_ideal_adaptor(x) ~ "NIW_ideal_adaptor",
      is.NIW_belief(x) ~ "NIW_belief",
      is.MVG_ideal_observer(x) ~ "MVG_ideal_observer",
      is.MVG(x) ~ "MVG",
      T ~ "Unrecognized class")

  return(class)
}
