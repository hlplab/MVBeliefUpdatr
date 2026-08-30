#' @include S7-make-objects.R
NULL

#' Example S7 objects built from the bundled Chodroff-Wilson data.
#'
#' `n_cues` selects VOT, VOT plus f0, or all three bundled cues. UVG and NIX
#' examples require `n_cues = 1`.
#' @name example-s7-objects
NULL

.example_data <- function(n_cues, categories = NULL) {
  .assert_true(.is_non_NA_scalar_count(n_cues), msg = "n_cues must be a positive whole number.")
  .assert_true(n_cues %in% 1:3, msg = "n_cues must be one of 1, 2, or 3.")
  data("ChodroffWilson2018", package = "MVBeliefUpdatr", envir = environment())
  if (is.null(categories)) categories <- levels(ChodroffWilson2018$category)
  ChodroffWilson2018[ChodroffWilson2018$category %in% categories, c("category", c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)])]
}

#' @rdname example-s7-objects
#' @export
example_uvg_category_representation <- function(n_cues = 1, category = "/b/") {
  .assert_true(n_cues == 1, msg = "UVG examples require n_cues = 1.")
  new_uvg_category_representation_from_data(.example_data(n_cues, category), cues = "VOT")
}

#' @rdname example-s7-objects
#' @export
example_nix_category_representation <- function(n_cues = 1, category = "/b/", kappa = 10, nu = 30) {
  .assert_true(n_cues == 1, msg = "NIX examples require n_cues = 1.")
  new_nix_category_representation_from_data(.example_data(n_cues, category), cues = "VOT", kappa = kappa, nu = nu)
}

#' @rdname example-s7-objects
#' @export
example_muvg_category_representation <- function(n_cues = 2, category = "/b/") {
  .assert_true(n_cues >= 2, msg = "MUVG examples require n_cues >= 2.")
  new_muvg_category_representation_from_data(.example_data(n_cues, category), cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)])
}

#' @rdname example-s7-objects
#' @export
example_mnix_category_representation <- function(n_cues = 2, category = "/b/", kappa = 10, nu = 30) {
  .assert_true(n_cues >= 2, msg = "MNIX examples require n_cues >= 2.")
  new_mnix_category_representation_from_data(.example_data(n_cues, category), cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)], kappa = kappa, nu = nu)
}

#' @rdname example-s7-objects
#' @export
example_mvg_category_representation <- function(n_cues = 1, category = "/b/") {
  new_mvg_category_representation_from_data(.example_data(n_cues, category), cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)])
}

#' @rdname example-s7-objects
#' @export
example_niw_category_representation <- function(n_cues = 1, category = "/b/", kappa = 10, nu = NULL) {
  if (is.null(nu)) nu <- n_cues + 2
  new_niw_category_representation_from_data(.example_data(n_cues, category), cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)], kappa = kappa, nu = nu)
}

#' @rdname example-s7-objects
#' @export
example_exemplar_category_representation <- function(n_cues = 1, category = "/b/") {
  new_exemplar_category_representation_from_data(.example_data(n_cues, category), cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)])
}

#' @rdname example-s7-objects
#' @export
example_category_representation <- function(type, n_cues = NULL, category = "/b/", ...) {
  type <- toupper(type)
  if (is.null(n_cues)) {
    n_cues <- if (type %in% c("MUVG", "MNIX")) 2 else 1
  }
  constructor <- switch(
    type,
    UVG = example_uvg_category_representation,
    NIX = example_nix_category_representation,
    MUVG = example_muvg_category_representation,
    MNIX = example_mnix_category_representation,
    MVG = example_mvg_category_representation,
    NIW = example_niw_category_representation,
    EXEMPLAR = example_exemplar_category_representation,
    .stop("type must be one of UVG, NIX, MUVG, MNIX, MVG, NIW, or EXEMPLAR.")
  )
  constructor(n_cues = n_cues, category = category, ...)
}

#' @rdname example-s7-objects
#' @export
example_uvg_category_representation_template <- function(n_cues = 1, categories = NULL) {
  .assert_true(n_cues == 1, msg = "UVG examples require n_cues = 1.")
  new_uvg_category_representation_template_from_data(.example_data(n_cues, categories), cues = "VOT")
}

#' @rdname example-s7-objects
#' @export
example_nix_category_representation_template <- function(n_cues = 1, categories = NULL, kappa = 10, nu = 30) {
  .assert_true(n_cues == 1, msg = "NIX examples require n_cues = 1.")
  new_nix_category_representation_template_from_data(.example_data(n_cues, categories), cues = "VOT", kappa = kappa, nu = nu)
}

#' @rdname example-s7-objects
#' @export
example_muvg_category_representation_template <- function(n_cues = 2, categories = NULL) {
  .assert_true(n_cues >= 2, msg = "MUVG examples require n_cues >= 2.")
  new_muvg_category_representation_template_from_data(.example_data(n_cues, categories), cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)])
}

#' @rdname example-s7-objects
#' @export
example_mnix_category_representation_template <- function(n_cues = 2, categories = NULL, kappa = 10, nu = 30) {
  .assert_true(n_cues >= 2, msg = "MNIX examples require n_cues >= 2.")
  new_mnix_category_representation_template_from_data(.example_data(n_cues, categories), cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)], kappa = kappa, nu = nu)
}

#' @rdname example-s7-objects
#' @export
example_mvg_category_representation_template <- function(n_cues = 1, categories = NULL) {
  new_mvg_category_representation_template_from_data(.example_data(n_cues, categories), cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)])
}

#' @rdname example-s7-objects
#' @export
example_niw_category_representation_template <- function(n_cues = 1, categories = NULL, kappa = 10, nu = NULL) {
  if (is.null(nu)) nu <- n_cues + 2
  new_niw_category_representation_template_from_data(.example_data(n_cues, categories), cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)], kappa = kappa, nu = nu)
}

#' @rdname example-s7-objects
#' @export
example_exemplar_category_representation_template <- function(n_cues = 1, categories = NULL) {
  new_exemplar_category_representation_template_from_data(.example_data(n_cues, categories), cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)])
}

#' @rdname example-s7-objects
#' @export
example_category_representation_template <- function(type, n_cues = NULL, categories = NULL, ...) {
  type <- toupper(type)
  if (is.null(n_cues)) {
    n_cues <- if (type %in% c("MUVG", "MNIX")) 2 else 1
  }
  constructor <- switch(
    type,
    UVG = example_uvg_category_representation_template,
    NIX = example_nix_category_representation_template,
    MUVG = example_muvg_category_representation_template,
    MNIX = example_mnix_category_representation_template,
    MVG = example_mvg_category_representation_template,
    NIW = example_niw_category_representation_template,
    EXEMPLAR = example_exemplar_category_representation_template,
    .stop("type must be one of UVG, NIX, MUVG, MNIX, MVG, NIW, or EXEMPLAR.")
  )
  constructor(n_cues = n_cues, categories = categories, ...)
}

#' @rdname example-s7-objects
#' @export
example_uvg_ideal_observer <- function(n_cues = 1, categories = NULL, ...) {
  .assert_true(n_cues == 1, msg = "UVG examples require n_cues = 1.")
  new_uvg_ideal_observer_from_data(.example_data(n_cues, categories), cues = "VOT", ...)
}

#' @rdname example-s7-objects
#' @export
example_nix_ideal_adaptor <- function(n_cues = 1, categories = NULL, kappa = 10, nu = 30, ...) {
  .assert_true(n_cues == 1, msg = "NIX examples require n_cues = 1.")
  new_nix_ideal_adaptor_from_data(.example_data(n_cues, categories), cues = "VOT", kappa = kappa, nu = nu, ...)
}

#' @rdname example-s7-objects
#' @export
example_muvg_ideal_observer <- function(n_cues = 2, categories = NULL, ...) {
  .assert_true(n_cues >= 2, msg = "MUVG examples require n_cues >= 2.")
  new_muvg_ideal_observer_from_data(.example_data(n_cues, categories), cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)], ...)
}

#' @rdname example-s7-objects
#' @export
example_mnix_ideal_adaptor <- function(n_cues = 2, categories = NULL, kappa = 10, nu = 30, ...) {
  .assert_true(n_cues >= 2, msg = "MNIX examples require n_cues >= 2.")
  new_mnix_ideal_adaptor_from_data(.example_data(n_cues, categories), cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)], kappa = kappa, nu = nu, ...)
}

#' @rdname example-s7-objects
#' @export
example_mvg_ideal_observer <- function(n_cues = 1, categories = NULL, ...) {
  new_mvg_ideal_observer_from_data(.example_data(n_cues, categories), cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)], ...)
}

#' @rdname example-s7-objects
#' @export
example_niw_ideal_adaptor <- function(n_cues = 1, categories = NULL, kappa = 10, nu = NULL, ...) {
  if (is.null(nu)) nu <- n_cues + 2
  new_niw_ideal_adaptor_from_data(.example_data(n_cues, categories), cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)], kappa = kappa, nu = nu, ...)
}

#' @rdname example-s7-objects
#' @export
example_exemplar_model <- function(n_cues = 1, categories = NULL, ...) {
  new_exemplar_model_from_data(.example_data(n_cues, categories), cues = c("VOT", "f0_semitones", "vowel_duration")[seq_len(n_cues)], ...)
}

#' @rdname example-s7-objects
#' @export
example_model <- function(type, n_cues = NULL, categories = NULL, ...) {
  type <- toupper(type)
  if (is.null(n_cues)) {
    n_cues <- if (type %in% c("MUVG", "MNIX")) 2 else 1
  }
  constructor <- switch(
    type,
    UVG = example_uvg_ideal_observer,
    NIX = example_nix_ideal_adaptor,
    MUVG = example_muvg_ideal_observer,
    MNIX = example_mnix_ideal_adaptor,
    MVG = example_mvg_ideal_observer,
    NIW = example_niw_ideal_adaptor,
    EXEMPLAR = example_exemplar_model,
    .stop("type must be one of UVG, NIX, MUVG, MNIX, MVG, NIW, or EXEMPLAR.")
  )
  constructor(n_cues = n_cues, categories = categories, ...)
}
