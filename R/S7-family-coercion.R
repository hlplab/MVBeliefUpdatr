#' @include S7-family.R
#' @include S7-class.R
#' @include S7-class-uvg.R
#' @include S7-class-nix.R
#' @include S7-class-muvg.R
#' @include S7-class-mnix.R
#' @include S7-class-mvg.R
#' @include S7-class-niw.R
#' @include S7-class-exemplar.R
NULL

# Cross-family coercion between S7 category representations, templates, and cognitive models.
#
# Supported coercions are: same-family (identity), matched point/uncertain pairs
# (UVG<->NIX, MUVG<->MNIX, MVG<->NIW), and any family <-> EXEMPLAR (by summarizing exemplars
# into point/uncertain parameters, or by sampling exemplars from a parametric representation).

.mvbu_coerce_category_representation <- function(x, to, kappa = NULL, nu = NULL, n = NULL, weights = NULL) {
  allowed <- c("UVG", "NIX", "MUVG", "MNIX", "MVG", "NIW", "EXEMPLAR")
  to <- toupper(to)
  .assert_true(to %in% allowed, msg = paste0("to must be one of: ", paste(allowed, collapse = ", ")))
  .assert_true(S7::S7_inherits(x, MVBU_CategoryRepresentation), msg = "x must be an MVBU_CategoryRepresentation.")

  from <- get_family_name(x)
  if (from == to) {
    return(x)
  }

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

  labels <- get_labels(x)
  category_labels <- labels$category
  cue_labels <- labels$cue

  fam <- get_family(x)

  if (to == "EXEMPLAR") {
    .assert_true(
      .is_scalar_count(n) && n >= 1,
      msg = "n must be a non-NA, positive whole number giving the number of exemplars to sample."
    )
    exemplars <- fam@sampling_function(x, n = n)
    return(
      new_exemplar_category_representation(
        category_labels = category_labels,
        cue_labels = cue_labels,
        exemplars = exemplars
      )
    )
  }

  point <- fam@expected_moments(x)
  mu <- point$mu
  Sigma <- as.matrix(point$Sigma)
  if (is.null(weights)) weights <- point$weights

  if (to %in% c("UVG", "NIX") && length(mu) != 1) {
    .stop(paste0("Coercion to ", to, " requires a single-cue representation (found ", length(mu), " cues)."))
  }

  switch(to,
    UVG = new_uvg_category_representation(category_labels, cue_labels, mu = mu[1], sigma2 = Sigma[1, 1]),
    MVG = new_mvg_category_representation(category_labels, cue_labels, mu = mu, Sigma = Sigma),
    MUVG = new_muvg_category_representation(
      category_labels, cue_labels,
      mu = mu, sigma2 = diag(as.matrix(Sigma)),
      weights = weights
    ),
    NIX = {
      .assert_numeric_scalar(kappa, msg = "kappa must be a non-NA scalar numeric value.")
      .assert_numeric_scalar(nu, msg = "nu must be a non-NA scalar numeric value greater than 2.")
      .assert_true(nu > 2, msg = "nu must be greater than 2 for a univariate NIX representation.")
      new_nix_category_representation(
        category_labels,
        cue_labels,
        m = mu[1],
        kappa = kappa,
        nu = nu,
        sigma2 = Sigma[1, 1] * (nu - 2) / nu
      )
    },
    MNIX = {
      k <- length(mu)
      kappa <- rep_len(kappa, k)
      nu <- rep_len(nu, k)
      .assert_true(all(!is.na(kappa)) && is.numeric(kappa), msg = "kappa must be numeric.")
      .assert_true(all(!is.na(nu)) && is.numeric(nu) && all(nu > 2), msg = "nu must be numeric and greater than 2 for each component.")
      sigma2 <- diag(as.matrix(Sigma))
      new_mnix_category_representation(
        category_labels,
        cue_labels,
        m = mu,
        kappa = kappa,
        nu = nu,
        sigma2 = sigma2 * (nu - 2) / nu,
        weights = weights
      )
    },
    NIW = {
      .assert_numeric_scalar(kappa, msg = "kappa must be a non-NA scalar numeric value.")
      .assert_numeric_scalar(nu, msg = "nu must be a non-NA scalar numeric value.")
      .assert_true(nu > length(cue_labels) + 1, msg = paste0("nu must be larger than dimensionality of cues + 1 (>", length(cue_labels) + 1, ")."))
      new_niw_category_representation(
        category_labels, cue_labels,
        m = mu,
        kappa = kappa,
        nu = nu,
        S = get_S_from_expected_Sigma(Sigma, nu)
      )
    },
    .stop("Unsupported target family.")
  )
}

.mvbu_coerce_category_representation_template <- function(x, to, kappa = NULL, nu = NULL, n = NULL, weights = NULL) {
  .assert_true(S7::S7_inherits(x, MVBU_CategoryRepresentationTemplate), msg = "x must be an MVBU_CategoryRepresentationTemplate.")
  representations <- lapply(
    x@representations,
    .mvbu_coerce_category_representation,
    to = to, kappa = kappa, nu = nu, n = n, weights = weights
  )
  new_category_representation_template(representations)
}

.mvbu_coerce_cognitive_model <- function(x, to, kappa = NULL, nu = NULL, n = NULL, weights = NULL) {
  .assert_true(S7::S7_inherits(x, MVBU_CognitiveModel), msg = "x must be an MVBU_CognitiveModel.")
  template <- .mvbu_coerce_category_representation_template(
    x@category_template,
    to = to, kappa = kappa, nu = nu, n = n, weights = weights
  )

  constructor <- switch(to,
    UVG = new_uvg_ideal_observer,
    NIX = new_nix_ideal_adaptor,
    MUVG = new_muvg_ideal_observer,
    MNIX = new_mnix_ideal_adaptor,
    MVG = new_mvg_ideal_observer,
    NIW = new_niw_ideal_adaptor,
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

#' Coerce representations, templates, and cognitive models between families
#'
#' Coerces between point-estimate families (UVG, MUVG, MVG), uncertainty
#' families (NIX, MNIX, NIW), and the non-parametric Exemplar family.
#' Supported coercions are: same-family (returned unchanged), matched
#' pairs (UVG<->NIX, MUVG<->MNIX, MVG<->NIW), and any family <-> EXEMPLAR.
#' Coercions into NIX/MNIX/NIW require `kappa` and `nu`; coercions into
#' EXEMPLAR require `n`, the number of exemplars to sample.
#'
#' Also provides generic lifters `as_ideal_adaptor()` and
#' `as_ideal_observer()` to lift representations or templates into full
#' cognitive models.
#'
#' @param x An `MVBU_CategoryRepresentation`,
#'   `MVBU_CategoryRepresentationTemplate`, or `MVBU_CognitiveModel` object.
#' @param kappa Strength of belief (pseudocount) about the category mean
#'   \eqn{\mu}.
#' @param nu Degrees of freedom (pseudocount) for covariance/variance beliefs
#'   \eqn{\Sigma}.
#' @param n Number of exemplars to sample when coercing into the EXEMPLAR family.
#' @param weights Optional per-cue weights for MUVG/MNIX targets.
#'   (default: preserved from `x` if available, otherwise the family
#'   constructor's default.)
#' @param model_family Optional character string specifying target model family
#'   ("NIX", "MNIX", "NIW" for adaptors; "UVG", "MUVG", "MVG" for observers).
#' @param ... Additional arguments passed to constructors.
#' @return An object of the same kind as `x` (representation, template, or
#'   model), coerced to the target family.
#' @examples
#' # Coerce between UVG and NIX representations
#' uvg_rep <- new_uvg_category_representation(
#'   category_labels = "A",
#'   cue_labels = "cue1",
#'   mu = 0,
#'   sigma2 = 1
#' )
#' nix_rep <- as_nix_category_representation(uvg_rep, kappa = 10, nu = 15)
#' nix_rep
#'
#' # Coerce back to UVG
#' as_uvg_category_representation(nix_rep)
#'
#' # Sample exemplars from a parametric representation
#' ex_rep <- as_exemplar_category_representation(uvg_rep, n = 20)
#' ex_rep
#'
#' # Lift representation to an ideal observer model
#' obs <- as_ideal_observer(uvg_rep)
#' obs
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
as_muvg_category_representation <- function(x, weights = NULL) .mvbu_coerce_category_representation(x, to = "MUVG", weights = weights)
#' @rdname as_family_coercion
#' @export
as_mnix_category_representation <- function(x, kappa, nu, weights = NULL) .mvbu_coerce_category_representation(x, to = "MNIX", kappa = kappa, nu = nu, weights = weights)
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
as_muvg_category_representation_template <- function(x, weights = NULL) .mvbu_coerce_category_representation_template(x, to = "MUVG", weights = weights)
#' @rdname as_family_coercion
#' @export
as_mnix_category_representation_template <- function(x, kappa, nu, weights = NULL) .mvbu_coerce_category_representation_template(x, to = "MNIX", kappa = kappa, nu = nu, weights = weights)
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
as_muvg_ideal_observer <- function(x, weights = NULL) .mvbu_coerce_cognitive_model(x, to = "MUVG", weights = weights)
#' @rdname as_family_coercion
#' @export
as_mnix_ideal_adaptor <- function(x, kappa, nu, weights = NULL) .mvbu_coerce_cognitive_model(x, to = "MNIX", kappa = kappa, nu = nu, weights = weights)
#' @rdname as_family_coercion
#' @export
as_mvg_ideal_observer <- function(x) .mvbu_coerce_cognitive_model(x, to = "MVG")
#' @rdname as_family_coercion
#' @export
as_niw_ideal_adaptor <- function(x, kappa, nu) .mvbu_coerce_cognitive_model(x, to = "NIW", kappa = kappa, nu = nu)
#' @rdname as_family_coercion
#' @export
as_exemplar_model <- function(x, n) {
  .mvbu_coerce_cognitive_model(x, to = "EXEMPLAR", n = n)
}

# Internal helper: instantiate an ideal observer model given family and template
.create_ideal_observer_model <- function(model_family,
                                         category_template,
                                         ...) {
  fam <- toupper(as.character(model_family))
  if (grepl("UVG", fam) && !grepl("MUVG", fam)) {
    new_uvg_ideal_observer(category_template = category_template, ...)
  } else if (grepl("MUVG", fam)) {
    new_muvg_ideal_observer(category_template = category_template, ...)
  } else if (grepl("MVG", fam)) {
    new_mvg_ideal_observer(category_template = category_template, ...)
  } else {
    .stop(sprintf(
      paste0(
        "Unsupported ideal observer model family for dynamic model creation: ",
        "'%s'."
      ),
      model_family
    ))
  }
}

# Internal helper: instantiate an ideal adaptor model given family and template
.create_ideal_adaptor_model <- function(model_family,
                                        category_template,
                                        ...) {
  fam <- toupper(as.character(model_family))
  if (grepl("NIX", fam) && !grepl("MNIX", fam)) {
    new_nix_ideal_adaptor(category_template = category_template, ...)
  } else if (grepl("MNIX", fam)) {
    new_mnix_ideal_adaptor(category_template = category_template, ...)
  } else if (grepl("NIW", fam)) {
    new_niw_ideal_adaptor(category_template = category_template, ...)
  } else {
    .stop(sprintf(
      paste0(
        "Unsupported ideal adaptor model family for dynamic model creation: ",
        "'%s'."
      ),
      model_family
    ))
  }
}

#' @rdname as_family_coercion
#' @export
as_ideal_adaptor <- function(x, model_family = NULL, ...) {
  .assert_true(
    !missing(x) && !is.null(x),
    msg = "x must be provided to as_ideal_adaptor()."
  )

  if (S7::S7_inherits(x, MVBU_CategoryRepresentation)) {
    x <- new_category_representation_template(list(x))
  }
  if (S7::S7_inherits(x, MVBU_CategoryRepresentationTemplate)) {
    if (is.null(model_family)) {
      model_family <- get_family_name(x)
    }
    return(.create_ideal_adaptor_model(model_family, category_template = x, ...))
  }
  if (S7::S7_inherits(x, MVBU_CognitiveModel)) {
    if (!is.null(model_family)) {
      return(.mvbu_coerce_cognitive_model(x, to = model_family, ...))
    }
    return(x)
  }
  .stop("x must be an MVBU_CategoryRepresentationTemplate, MVBU_CategoryRepresentation, or MVBU_CognitiveModel.")
}

#' @rdname as_family_coercion
#' @export
as_ideal_observer <- function(x, model_family = NULL, ...) {
  .assert_true(
    !missing(x) && !is.null(x),
    msg = "x must be provided to as_ideal_observer()."
  )

  if (S7::S7_inherits(x, MVBU_CategoryRepresentation)) {
    x <- new_category_representation_template(list(x))
  }
  if (S7::S7_inherits(x, MVBU_CategoryRepresentationTemplate)) {
    if (is.null(model_family)) {
      model_family <- get_family_name(x)
    }
    return(.create_ideal_observer_model(model_family, category_template = x, ...))
  }
  if (S7::S7_inherits(x, MVBU_CognitiveModel)) {
    if (!is.null(model_family)) {
      return(.mvbu_coerce_cognitive_model(x, to = model_family, ...))
    }
    target_fam <- switch(get_family_name(x),
      NIX = "UVG",
      MNIX = "MUVG",
      NIW = "MVG",
      get_family_name(x)
    )
    return(.mvbu_coerce_cognitive_model(x, to = target_fam, ...))
  }
  .stop("x must be an MVBU_CategoryRepresentationTemplate, MVBU_CategoryRepresentation, or MVBU_CognitiveModel.")
}
