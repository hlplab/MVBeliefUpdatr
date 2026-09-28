# =============================================================================
# Factory & Generalized Dispatch Helpers for S7 Cognitive Models
# =============================================================================

#' @include S7-generics.R
#' @include S7-class.R
#' @include S7-class-nix.R
#' @include S7-class-mnix.R
#' @include S7-class-niw.R
NULL

#' Internal helper: Create an ideal adaptor category representation dynamically
#'
#' Factory dispatcher that constructs a category representation (NIX, MNIX, or NIW)
#' without hardcoding `if/else` checks across calling functions.
#'
#' @param model_family Character string indicating model type ("NIX", "MNIX", "NIW", etc.).
#' @param category_labels Character vector of category labels.
#' @param cue_labels Character vector of cue labels.
#' @param m Mean vector or scalar.
#' @param s Scale matrix, variance vector, or scalar.
#' @param kappa Prior weight parameter.
#' @param nu Degrees of freedom parameter.
#' @return An S7 category representation object.
#' @keywords internal
.create_ideal_adaptor_representation <- function(
  model_family,
  category_labels,
  cue_labels,
  m,
  s,
  kappa,
  nu
) {
  fam <- toupper(as.character(model_family))

  if (grepl("NIX", fam) && !grepl("MNIX", fam)) {
    new_nix_category_representation(
      category_labels = category_labels,
      cue_labels = cue_labels,
      m = as.numeric(m)[1],
      sigma2 = as.numeric(s)[1] / as.numeric(nu),
      kappa = as.numeric(kappa),
      nu = as.numeric(nu)
    )
  } else if (grepl("MNIX", fam)) {
    k_comp <- length(as.vector(m))
    new_mnix_category_representation(
      category_labels = category_labels,
      cue_labels = cue_labels,
      m = as.vector(m),
      sigma2 = as.vector(s) / as.numeric(nu),
      kappa = rep_len(as.numeric(kappa), k_comp),
      nu = rep_len(as.numeric(nu), k_comp)
    )
  } else if (grepl("NIW", fam)) {
    new_niw_category_representation(
      category_labels = category_labels,
      cue_labels = cue_labels,
      m = as.vector(m),
      S = as.matrix(s),
      kappa = as.numeric(kappa),
      nu = as.numeric(nu)
    )
  } else {
    .stop(sprintf("Unsupported ideal adaptor model family for dynamic representation creation: '%s'.", model_family))
  }
}


.create_ideal_adaptor_model <- function(model_family, category_template, ...) {
  fam <- toupper(as.character(model_family))
  if (grepl("NIX", fam) && !grepl("MNIX", fam)) {
    new_nix_ideal_adaptor(category_template = category_template, ...)
  } else if (grepl("MNIX", fam)) {
    new_mnix_ideal_adaptor(category_template = category_template, ...)
  } else if (grepl("NIW", fam)) {
    new_niw_ideal_adaptor(category_template = category_template, ...)
  } else {
    .stop(sprintf("Unsupported ideal adaptor model family for dynamic model creation: '%s'.", model_family))
  }
}
