#' @import assertthat
#' @import rstan

class_name = "mv_ibbu_stanfit"

setClass(class_name,
        contains = "stanfit")

is.mv_ibbu_stanfit = function(x) {
  if (class(x) == class_name)
    return(TRUE) else return(FALSE)
}
