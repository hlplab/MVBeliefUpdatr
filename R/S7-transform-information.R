#' @include S7-core-classes.R

#' @name MVBU-TransformInformation-class
#' @title S7 base class for transformation information
#' @docType class
#'
#' @slot transform.parameters A list of sufficient parameters for the transform.
#' @slot transform.function Function to transform parameters.
#' @slot untransform.function Function to untransform parameters.
#' @export
MVBU_TransformInformation <- S7::new_class(
  "MVBU_TransformInformation",
  package = NULL,
  parent = MVBU_Object,
  properties = list(
    transform.parameters = S7::class_list,
    transform.function = S7::class_function,
    untransform.function = S7::class_function
  ),
  constructor = function(transform.parameters = list(), transform.function = function(x) x, untransform.function = function(x) x) {
    S7::new_object(
      MVBU_Object(),
      transform.parameters = as.list(transform.parameters),
      transform.function = transform.function,
      untransform.function = untransform.function
    )
  },
  validator = function(self) {
    if (!is.list(self@transform.parameters)) {
      return("`transform.parameters` must be a list")
    }
    if (!is.function(self@transform.function)) {
      return("`transform.function` must be a function")
    }
    if (!is.function(self@untransform.function)) {
      return("`untransform.function` must be a function")
    }
    NULL
  }
)

#' Construct transform-information from an affine-transform object.
#'
#' @param transform An affine-transform object produced by get_affine_transform().
#' @return An S7 object of class MVBU_TransformInformation.
#' @keywords internal
#' @noRd
new_transform_information <- function(transform) {
  MVBU_TransformInformation(
    transform.parameters = transform$transform.parameters,
    transform.function = transform$transform.function,
    untransform.function = transform$untransform.function
  )
}
