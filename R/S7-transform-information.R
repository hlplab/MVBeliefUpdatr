#' @include S7-class.R
NULL

#' Helper to construct affine transformation function on data
#'
#' @param cues Character vector of cue labels.
#' @param shift Numeric shift vector.
#' @param SCALE Numeric scale matrix or scalar.
#'
#' @return A function \code{f(data, return_type)}.
#' @keywords internal
#' @noRd
.make_affine_transform_function <- function(cues, shift, SCALE) {
  force(cues)
  force(shift)
  force(SCALE)
  n_cues <- length(cues)

  function(data, return_type = c("replace", "add", "cues only")[1]) {
    return_type <- match.arg(
      return_type,
      c("replace", "add", "cues only", "only cues")
    )
    .assert_data_contains_cols(data, cues)
    newdata <- as.matrix(data[, cues, drop = FALSE])

    if (n_cues == 1L) {
      trans_cues <- (newdata[, 1] + shift) * as.numeric(SCALE)
      newdata <- matrix(trans_cues, ncol = 1L)
      colnames(newdata) <- cues
    } else {
      newdata <- sweep(newdata, 2, -shift) %*% t(SCALE)
      colnames(newdata) <- cues
    }

    newdata <- as.data.frame(newdata)
    if (return_type == "replace") {
      data[, cues] <- newdata
    } else if (return_type == "add") {
      names(newdata) <- paste0(cues, "_transformed")
      data <- cbind(data, newdata)
    } else if (return_type %in% c("cues only", "only cues")) {
      names(newdata) <- cues
      data <- newdata
    }
    data
  }
}

#' S7 class for transformation information
#'
#' Represents an invertible affine transformation \eqn{f(x) = A(x + b)},
#' where \eqn{b} is the shift vector and \eqn{A} is the scale/rotation matrix.
#'
#' @param cue.labels Character vector of cue variable names.
#' @param type Transformation type (e.g., \code{"identity"}, \code{"center"},
#'   \code{"standardize"}, \code{"PCA whiten"}, \code{"ZCA whiten"}).
#' @param shift Numeric shift vector added before scaling.
#' @param SCALE Numeric square invertible matrix or scalar scale factor.
#' @param transform.function Function to apply transformation to data.
#'
#' @docType class
#' @name MVBU-TransformInformation-class
#' @rdname MVBU-TransformInformation-class
#' @export
MVBU_TransformInformation <- S7::new_class(
  "MVBU_TransformInformation",
  package = NULL,
  parent = MVBU_Object,
  properties = list(
    cue.labels = S7::class_character,
    type = S7::class_character,
    shift = S7::class_numeric,
    SCALE = S7::new_property(
      validator = function(value) {
        if (!is.numeric(value)) {
          return("`SCALE` must be numeric")
        }
        NULL
      }
    ),
    transform.function = S7::class_function
  ),
  constructor = function(
    cue.labels = character(),
    type = "identity",
    shift = numeric(),
    SCALE = 1,
    transform.function = NULL
  ) {
    t_fn <- if (!is.null(transform.function) &&
      is.function(transform.function)) {
      transform.function
    } else if (length(cue.labels) > 0L && length(shift) > 0L) {
      .make_affine_transform_function(cue.labels, shift, SCALE)
    } else {
      function(x, ...) x
    }

    S7::new_object(
      MVBU_Object(),
      cue.labels = as.character(cue.labels),
      type = as.character(type),
      shift = as.numeric(shift),
      SCALE = SCALE,
      transform.function = t_fn
    )
  },
  validator = function(self) {
    if (length(self@cue.labels) != length(self@shift)) {
      return("Length of `cue.labels` must match length of `shift`.")
    }
    NULL
  }
)

#' Invert an affine transformation
#'
#' Computes the exact algebraic inverse of an affine transformation
#' \eqn{f(x) = A(x + b)}. The inverse transformation is itself an affine
#' transformation \eqn{f^{-1}(y) = A^{-1}(y - A b) = A_{inv}(y + b_{inv})},
#' where \eqn{A_{inv} = A^{-1}} and \eqn{b_{inv} = -A b}.
#'
#' @param transform An \code{\link{MVBU_TransformInformation}} object.
#' @return An inverted \code{\link{MVBU_TransformInformation}} object.
#'
#' @export
invert_transform <- S7::new_generic("invert_transform", "transform")

#' @rdname invert_transform
#' @export
S7::method(invert_transform, MVBU_TransformInformation) <- function(transform) {
  k <- length(transform@shift)
  scale_val <- transform@SCALE
  shift_val <- transform@shift

  if (k == 1L) {
    inv_scale <- 1 / as.numeric(scale_val)
    inv_shift <- -as.numeric(scale_val) * shift_val
  } else {
    inv_scale <- solve(scale_val)
    inv_shift <- -as.vector(scale_val %*% shift_val)
  }

  inv_type <- if (grepl("^inverse_", transform@type)) {
    sub("^inverse_", "", transform@type)
  } else {
    paste0("inverse_", transform@type)
  }

  MVBU_TransformInformation(
    cue.labels = transform@cue.labels,
    type = inv_type,
    shift = inv_shift,
    SCALE = inv_scale
  )
}
