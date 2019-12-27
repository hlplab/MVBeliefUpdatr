#' @import assertthat
#' @import rstan

class_name = "mv_ibbu_stanfit"

setClass(class_name,
        contains = "stanfit")

is.mv_ibbu_stanfit = function(x) {
  if (class(x) == "stanfit") message("Accepting stanfit as valid class. This might change in the future.")
  if (class(x) %in% c(class_name, "stanfit"))
    return(TRUE) else return(FALSE)
}
