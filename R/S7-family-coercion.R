#' @include S7-core-classes.R
#' @include S7-core-uvg-classes.R
#' @include S7-core-nix-classes.R
#' @include S7-core-muvg-classes.R
#' @include S7-core-mnix-classes.R
#' @include S7-core-mvg-classes.R
#' @include S7-core-niw-classes.R
#' @include S7-core-exemplar-classes.R
NULL

# Cross-family coercion between S7 category representations, templates, and cognitive models.
#
# Supported coercions are: same-family (identity), matched point/uncertain pairs
# (UVG<->NIX, MUVG<->MNIX, MVG<->NIW), and any family <-> EXEMPLAR (by summarizing exemplars
# into point/uncertain parameters, or by sampling exemplars from a parametric representation).

.mvbu_family_of_representation <- function(x) {
  if (S7::S7_inherits(x, UVG_CategoryRepresentation)) return("UVG")
  if (S7::S7_inherits(x, NIX_CategoryRepresentation)) return("NIX")
  if (S7::S7_inherits(x, MUVG_CategoryRepresentation)) return("MUVG")
  if (S7::S7_inherits(x, MNIX_CategoryRepresentation)) return("MNIX")
  if (S7::S7_inherits(x, MVG_CategoryRepresentation)) return("MVG")
  if (S7::S7_inherits(x, NIW_CategoryRepresentation)) return("NIW")
  if (S7::S7_inherits(x, Exemplar_CategoryRepresentation)) return("EXEMPLAR")
  .stop("x must be an MVBU_CategoryRepresentation of a supported family.")
}

# Draw n cue observations (a matrix with one column per cue) from a single category representation.
.mvbu_sample_from_representation <- function(x, from, n, with_replacement = TRUE) {
  switch(
    from,
    UVG = matrix(stats::rnorm(n, mean = x@mu, sd = sqrt(x@sigma2)), ncol = 1),
    NIX = {
      scale <- sqrt(x@sigma2 * (x@kappa + 1) / x@kappa)
      matrix(x@m + scale * stats::rt(n, df = x@nu), ncol = 1)
    },
    MVG = mvtnorm::rmvnorm(n, mean = x@mu, sigma = x@Sigma),
    NIW = {
      d <- length(x@m)
      df <- x@nu - d + 1
      Sigma_pred <- ((x@kappa + 1) / (x@kappa * df)) * x@S
      .rmvt(n = n, delta = x@m, sigma = Sigma_pred, df = df)
    },
    MUVG = {
      k <- length(x@component_mu)
      out <- vapply(seq_len(k), function(i) stats::rnorm(n, mean = x@component_mu[i], sd = sqrt(x@component_sigma2[i])), numeric(n))
      matrix(out, nrow = n, ncol = k)
    },
    MNIX = {
      k <- length(x@component_m)
      out <- vapply(seq_len(k), function(i) {
        scale <- sqrt(x@component_sigma2[i] * (x@component_kappa[i] + 1) / x@component_kappa[i])
        x@component_m[i] + scale * stats::rt(n, df = x@component_nu[i])
      }, numeric(n))
      matrix(out, nrow = n, ncol = k)
    },
    EXEMPLAR = {
      n_avail <- nrow(x@exemplars)
      if (isTRUE(with_replacement)) {
        idx <- sample.int(n_avail, size = n, replace = TRUE, prob = x@exemplar_weights)
      } else {
        if (n > n_avail) {
          cat_lbl <- .mvbu_extract_label_metadata(x)$category
          .stop(sprintf(
            "Cannot sample %d observations without replacement from category '%s' which contains only %d exemplars.",
            n, cat_lbl, n_avail
          ))
        }
        idx <- sample.int(n_avail, size = n, replace = FALSE, prob = x@exemplar_weights)
      }
      x@exemplars[idx, , drop = FALSE]
    },
    .stop("Unsupported source family for exemplar sampling.")
  )
}

# Canonical (mu, Sigma, component weights) summary of a non-exemplar representation, used as the
# common intermediate for all coercions between parametric families.
.mvbu_point_summary_from_representation <- function(x, from) {
  switch(
    from,
    UVG = list(mu = x@mu, Sigma = matrix(x@sigma2, 1, 1), weights = 1),
    NIX = list(mu = x@m, Sigma = matrix(x@sigma2 / (x@nu - 2), 1, 1), weights = 1),
    MUVG = list(mu = x@component_mu, Sigma = diag(x@component_sigma2, nrow = length(x@component_sigma2)), weights = x@component_weights),
    MNIX = list(
      mu = x@component_m,
      Sigma = diag(x@component_sigma2 / (x@component_nu - 2), nrow = length(x@component_sigma2)),
      weights = x@component_weights
    ),
    MVG = list(mu = x@mu, Sigma = x@Sigma, weights = NULL),
    NIW = list(mu = get_expected_mu_from_m(x@m), Sigma = get_expected_Sigma_from_S(x@S, x@nu), weights = NULL),
    .stop("Unsupported source family.")
  )
}

.mvbu_coerce_category_representation <- function(x, to, kappa = NULL, nu = NULL, n = NULL, component_weights = NULL) {
  allowed <- c("UVG", "NIX", "MUVG", "MNIX", "MVG", "NIW", "EXEMPLAR")
  to <- toupper(to)
  .assert_true(to %in% allowed, msg = paste0("to must be one of: ", paste(allowed, collapse = ", ")))
  .assert_true(S7::S7_inherits(x, MVBU_CategoryRepresentation), msg = "x must be an MVBU_CategoryRepresentation.")

  from <- .mvbu_family_of_representation(x)
  if (from == to) return(x)

  matched_pairs <- c(UVG = "NIX", NIX = "UVG", MUVG = "MNIX", MNIX = "MUVG", MVG = "NIW", NIW = "MVG")
  is_matched_pair <- identical(unname(matched_pairs[from]), to)
  involves_exemplar <- (from == "EXEMPLAR") || (to == "EXEMPLAR")

  if (!is_matched_pair && !involves_exemplar) {
    .stop(
      paste0(
        "Coercion from ", from, " to ", to, " category representations is not supported. ",
        "Supported coercions are: same-family, matched point/uncertain pairs (UVG<->NIX, MUVG<->MNIX, MVG<->NIW), ",
        "and any family <-> EXEMPLAR."
      )
    )
  }

  labels <- .mvbu_extract_label_metadata(x)
  category_labels <- labels$category
  cue_labels <- labels$cue

  if (to == "EXEMPLAR") {
    .assert_true(.is_non_NA_scalar_count(n) && n >= 1, msg = "n must be a non-NA, positive whole number giving the number of exemplars to sample.")
    exemplars <- .mvbu_sample_from_representation(x, from = from, n = n)
    return(new_exemplar_category_representation(category_labels = category_labels, cue_labels = cue_labels, exemplars = exemplars))
  }

  if (from == "EXEMPLAR") {
    mu <- as.numeric(x@exemplar_weights %*% x@exemplars)
    Sigma <- stats::cov.wt(x@exemplars, wt = x@exemplar_weights)$cov
  } else {
    point <- .mvbu_point_summary_from_representation(x, from = from)
    mu <- point$mu
    Sigma <- point$Sigma
    if (is.null(component_weights)) component_weights <- point$weights
  }

  if (to %in% c("UVG", "NIX") && length(mu) != 1) {
    .stop(paste0("Coercion to ", to, " requires a single-cue representation (found ", length(mu), " cues)."))
  }

  switch(
    to,
    UVG = new_uvg_category_representation(category_labels, cue_labels, mu = mu[1], sigma2 = Sigma[1, 1]),
    MVG = new_mvg_category_representation(category_labels, cue_labels, mu = mu, Sigma = Sigma),
    MUVG = new_muvg_category_representation(
      category_labels, cue_labels,
      component_mu = mu, component_sigma2 = diag(as.matrix(Sigma)),
      component_weights = component_weights
    ),
    NIX = {
      .assert_non_NA_scalar_numeric(kappa, msg = "kappa must be a non-NA scalar numeric value.")
      .assert_non_NA_scalar_numeric(nu, msg = "nu must be a non-NA scalar numeric value greater than 2.")
      .assert_true(nu > 2, msg = "nu must be greater than 2 for a univariate NIX representation.")
      new_nix_category_representation(category_labels, cue_labels, m = mu[1], kappa = kappa, nu = nu, sigma2 = Sigma[1, 1] * (nu - 2))
    },
    NIW = {
      .assert_non_NA_scalar_numeric(kappa, msg = "kappa must be a non-NA scalar numeric value.")
      .assert_non_NA_scalar_numeric(nu, msg = "nu must be a non-NA scalar numeric value.")
      .assert_true(nu > length(cue_labels) + 1, msg = paste0("nu must be larger than dimensionality of cues + 1 (>", length(cue_labels) + 1, ")."))
      new_niw_category_representation(category_labels, cue_labels, m = mu, kappa = kappa, nu = nu, S = get_S_from_expected_Sigma(Sigma, nu))
    },
    MNIX = {
      k <- length(mu)
      kappa <- rep_len(kappa, k)
      nu <- rep_len(nu, k)
      .assert_true(all(!is.na(kappa)) && is.numeric(kappa), msg = "kappa must be numeric.")
      .assert_true(all(!is.na(nu)) && is.numeric(nu) && all(nu > 2), msg = "nu must be numeric and greater than 2 for each component.")
      component_sigma2 <- diag(as.matrix(Sigma))
      new_mnix_category_representation(
        category_labels, cue_labels,
        component_m = mu, component_kappa = kappa, component_nu = nu,
        component_sigma2 = component_sigma2 * (nu - 2),
        component_weights = component_weights
      )
    },
    .stop("Unsupported target family.")
  )
}

.mvbu_coerce_category_representation_template <- function(x, to, kappa = NULL, nu = NULL, n = NULL, component_weights = NULL) {
  .assert_true(S7::S7_inherits(x, MVBU_CategoryRepresentationTemplate), msg = "x must be an MVBU_CategoryRepresentationTemplate.")
  representations <- lapply(
    x@representations,
    .mvbu_coerce_category_representation,
    to = to, kappa = kappa, nu = nu, n = n, component_weights = component_weights
  )
  new_category_representation_template(representations)
}

.mvbu_coerce_cognitive_model <- function(x, to, kappa = NULL, nu = NULL, n = NULL, component_weights = NULL) {
  .assert_true(S7::S7_inherits(x, MVBU_CognitiveModel), msg = "x must be an MVBU_CognitiveModel.")
  template <- .mvbu_coerce_category_representation_template(
    x@category_template,
    to = to, kappa = kappa, nu = nu, n = n, component_weights = component_weights
  )

  constructor <- switch(
    to,
    UVG = new_uvg_ideal_observer, NIX = new_nix_ideal_adaptor,
    MUVG = new_muvg_ideal_observer, MNIX = new_mnix_ideal_adaptor,
    MVG = new_mvg_ideal_observer, NIW = new_niw_ideal_adaptor,
    EXEMPLAR = new_exemplar_model,
    .stop("Unsupported target family.")
  )

  constructor(
    category_template = template,
    decision_rule = x@decision_rule,
    category_prior = x@category_prior,
    lapse_rate = get_lapse_rate(x),
    lapse_bias = get_lapse_bias(x),
    Sigma_noise = get_noise(x),
    noise_treatment = get_noise_treatment(x),
    lapse_treatment = get_lapse_treatment(x),
    metadata = x@metadata
  )
}

#' Coerce category representations, templates, and cognitive models between families
#'
#' These functions coerce between the point-estimate families (UVG, MUVG, MVG), the corresponding
#' uncertainty families (NIX, MNIX, NIW), and the non-parametric Exemplar family. Supported coercions
#' are: same-family (returned unchanged), matched pairs (UVG<->NIX, MUVG<->MNIX, MVG<->NIW), and any
#' family <-> EXEMPLAR. Coercions into NIX/MNIX/NIW require `kappa` and `nu`; coercions into EXEMPLAR
#' require `n`, the number of exemplars to sample.
#'
#' @param x An MVBU_CategoryRepresentation, MVBU_CategoryRepresentationTemplate, or MVBU_CognitiveModel object.
#' @param kappa Strength of belief (pseudocount) about the category mean.
#' @param nu Strength of belief (pseudocount) about the category covariance/variance.
#' @param n Number of exemplars to sample when coercing into the EXEMPLAR family.
#' @param component_weights Optional per-component/per-cue weights for MUVG/MNIX targets. (default: preserved
#'   from `x` if available, otherwise the family constructor's default.)
#' @return An object of the same kind as `x` (representation, template, or model), coerced to the target family.
#' @name as_family_coercion
NULL

#' @rdname as_family_coercion
#' @export
as_uvg_category_representation <- function(x) .mvbu_coerce_category_representation(x, to = "UVG")
#' @rdname as_family_coercion
#' @export
as_nix_category_representation <- function(x, kappa, nu) .mvbu_coerce_category_representation(x, to = "NIX", kappa = kappa, nu = nu)
#' @rdname as_family_coercion
#' @export
as_muvg_category_representation <- function(x, component_weights = NULL) .mvbu_coerce_category_representation(x, to = "MUVG", component_weights = component_weights)
#' @rdname as_family_coercion
#' @export
as_mnix_category_representation <- function(x, kappa, nu, component_weights = NULL) .mvbu_coerce_category_representation(x, to = "MNIX", kappa = kappa, nu = nu, component_weights = component_weights)
#' @rdname as_family_coercion
#' @export
as_mvg_category_representation <- function(x) .mvbu_coerce_category_representation(x, to = "MVG")
#' @rdname as_family_coercion
#' @export
as_niw_category_representation <- function(x, kappa, nu) .mvbu_coerce_category_representation(x, to = "NIW", kappa = kappa, nu = nu)
#' @rdname as_family_coercion
#' @export
as_exemplar_category_representation <- function(x, n) .mvbu_coerce_category_representation(x, to = "EXEMPLAR", n = n)

#' @rdname as_family_coercion
#' @export
as_uvg_category_representation_template <- function(x) .mvbu_coerce_category_representation_template(x, to = "UVG")
#' @rdname as_family_coercion
#' @export
as_nix_category_representation_template <- function(x, kappa, nu) .mvbu_coerce_category_representation_template(x, to = "NIX", kappa = kappa, nu = nu)
#' @rdname as_family_coercion
#' @export
as_muvg_category_representation_template <- function(x, component_weights = NULL) .mvbu_coerce_category_representation_template(x, to = "MUVG", component_weights = component_weights)
#' @rdname as_family_coercion
#' @export
as_mnix_category_representation_template <- function(x, kappa, nu, component_weights = NULL) .mvbu_coerce_category_representation_template(x, to = "MNIX", kappa = kappa, nu = nu, component_weights = component_weights)
#' @rdname as_family_coercion
#' @export
as_mvg_category_representation_template <- function(x) .mvbu_coerce_category_representation_template(x, to = "MVG")
#' @rdname as_family_coercion
#' @export
as_niw_category_representation_template <- function(x, kappa, nu) .mvbu_coerce_category_representation_template(x, to = "NIW", kappa = kappa, nu = nu)
#' @rdname as_family_coercion
#' @export
as_exemplar_category_representation_template <- function(x, n) .mvbu_coerce_category_representation_template(x, to = "EXEMPLAR", n = n)

#' @rdname as_family_coercion
#' @export
as_uvg_ideal_observer <- function(x) .mvbu_coerce_cognitive_model(x, to = "UVG")
#' @rdname as_family_coercion
#' @export
as_nix_ideal_adaptor <- function(x, kappa, nu) .mvbu_coerce_cognitive_model(x, to = "NIX", kappa = kappa, nu = nu)
#' @rdname as_family_coercion
#' @export
as_muvg_ideal_observer <- function(x, component_weights = NULL) .mvbu_coerce_cognitive_model(x, to = "MUVG", component_weights = component_weights)
#' @rdname as_family_coercion
#' @export
as_mnix_ideal_adaptor <- function(x, kappa, nu, component_weights = NULL) .mvbu_coerce_cognitive_model(x, to = "MNIX", kappa = kappa, nu = nu, component_weights = component_weights)
#' @rdname as_family_coercion
#' @export
as_mvg_ideal_observer <- function(x) .mvbu_coerce_cognitive_model(x, to = "MVG")
#' @rdname as_family_coercion
#' @export
as_niw_ideal_adaptor <- function(x, kappa, nu) .mvbu_coerce_cognitive_model(x, to = "NIW", kappa = kappa, nu = nu)
#' @rdname as_family_coercion
#' @export
as_exemplar_model <- function(x, n) .mvbu_coerce_cognitive_model(x, to = "EXEMPLAR", n = n)

