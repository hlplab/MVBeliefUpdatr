#' @include asserts.R
#' @include S7-class.R
#' @include S7-generics.R
#' @importFrom S7 new_class new_object props class_character S7_inherits method
#' @importFrom tibble tibble as_tibble
NULL

# -------------------------------------------------------------------------
# S7 Prior Class and Constructor
# -------------------------------------------------------------------------

#' S7 class for prior distribution specifications
#'
#' Represents prior distribution specifications in a dual-space (original vs.
#' transformed cue space) framework, supporting flexible brms-style
#' declarations across parameter classes.
#'
#' @param declaration Character vector of prior distribution strings
#'   (e.g., \code{"cauchy(0, 5)"}, \code{"lkj(1)"}, \code{"normal(0, 1)"}).
#' @param par_class Character vector of parameter classes (e.g., \code{"tau_0"},
#'   \code{"cor_0"}, \code{"m_0"}, \code{"weight"}, \code{"lapse_rate"},
#'   \code{"lapse_bias"}, \code{"category_prior"}, \code{"kappa_0"},
#'   \code{"nu_0"}).
#' @param category Character vector specifying category names (default:
#'   \code{""} for class-wide or category-invariant priors).
#' @param cue Character vector specifying cue names (default: \code{""} for
#'   all cues).
#' @param par Character vector specifying parameter or cue names (default:
#'   \code{""} for class-wide priors).
#' @param cue_space Character vector specifying the cue space for the prior:
#'   \code{"transformed"} (default) or \code{"original"}.
#'
#' @docType class
#' @name MVBU_Prior
#' @rdname MVBU_Prior
#' @export
#' @examples
#' # Specify class-wide priors
#' p1 <- parameter_prior("cauchy(0, 5)", class = "tau_0")
#' p2 <- parameter_prior("lkj(1)", class = "cor_0")
#' p3 <- parameter_prior("normal(0, 1)", class = "m_0")
#'
#' # Combine priors into a list
#' c(p1, p2, p3)
MVBU_Prior <- S7::new_class(
  "MVBU_Prior",
  package = NULL,
  parent = MVBU_Object,
  properties = list(
    declaration = S7::class_character,
    par_class = S7::class_character,
    category = S7::class_character,
    cue = S7::class_character,
    par = S7::class_character,
    cue_space = S7::class_character
  ),
  constructor = function(
    declaration = character(),
    par_class = character(),
    category = character(),
    cue = character(),
    par = character(),
    cue_space = character()
  ) {
    if (length(declaration) == 0L) {
      return(
        S7::new_object(
          MVBU_Object(),
          declaration = character(),
          par_class = character(),
          category = character(),
          cue = character(),
          par = character(),
          cue_space = character()
        )
      )
    }
    n <- length(declaration)
    if (length(par_class) == 1L && n > 1L) par_class <- rep(par_class, n)
    if (length(category) == 0L) category <- rep("", n)
    if (length(category) == 1L && n > 1L) category <- rep(category, n)
    if (length(cue) == 0L) cue <- rep("", n)
    if (length(cue) == 1L && n > 1L) cue <- rep(cue, n)
    if (length(par) == 0L) {
      par <- ifelse(nzchar(category), category, cue)
    }
    if (length(par) == 1L && n > 1L) par <- rep(par, n)
    if (length(cue_space) == 0L) cue_space <- rep("transformed", n)
    if (length(cue_space) == 1L && n > 1L) cue_space <- rep(cue_space, n)

    S7::new_object(
      MVBU_Object(),
      declaration = as.character(declaration),
      par_class = as.character(par_class),
      category = as.character(category),
      cue = as.character(cue),
      par = as.character(par),
      cue_space = as.character(cue_space)
    )
  },
  validator = function(self) {
    n <- length(self@declaration)
    if (length(self@par_class) != n || length(self@category) != n ||
        length(self@cue) != n || length(self@par) != n ||
        length(self@cue_space) != n) {
      return("All properties of MVBU_Prior must have the same length.")
    }
    allowed_classes <- c(
      "tau_0", "cor_0", "m_0", "weight", "weights", "lapse_rate", "lapse_bias",
      "category_prior", "kappa_0", "nu_0"
    )
    invalid_classes <- setdiff(self@par_class, allowed_classes)
    if (length(invalid_classes) > 0L) {
      return(
        sprintf(
          "Invalid prior class: '%s'. Allowed classes are: %s.",
          invalid_classes[1],
          paste(allowed_classes, collapse = ", ")
        )
      )
    }
    allowed_spaces <- c("transformed", "original")
    invalid_spaces <- setdiff(self@cue_space, allowed_spaces)
    if (length(invalid_spaces) > 0L) {
      return(
        sprintf(
          "Invalid cue_space: '%s'. Must be 'transformed' or 'original'.",
          invalid_spaces[1]
        )
      )
    }
    NULL
  }
)

#' Declare a prior distribution
#'
#' Constructs a prior specification in the style of \pkg{brms}.
#'
#' @param declaration Character string or numeric scalar specifying the prior
#'   distribution (e.g., \code{"cauchy(0, 5)"}, \code{"lkj(1)"},
#'   \code{"normal(0, 1)"}, or numeric \code{5}).
#' @param class Character string specifying the parameter class. Supported
#'   classes include:
#'   \itemize{
#'     \item \code{"tau_0"}: Scale / standard deviation for category covariance.
#'     \item \code{"cor_0"}: Correlation matrix / Cholesky factor
#'       \code{L_omega}.
#'     \item \code{"m_0"}: Category prior mean \code{m_0}.
#'     \item \code{"weight"}: Cue integration weights.
#'     \item \code{"lapse_rate"}: Lapse probability.
#'     \item \code{"lapse_bias"}: Lapse response distribution simplex.
#'     \item \code{"category_prior"}: Category prior probability simplex.
#'     \item \code{"kappa_0"}: Prior pseudocount for mean.
#'     \item \code{"nu_0"}: Prior pseudocount for covariance.
#'   }
#' @param category Character string specifying category name (default:
#'   \code{""} or \code{NULL} for class-wide priors).
#' @param cue Character string specifying cue name (default: \code{""} or
#'   \code{NULL} for all cues).
#' @param par Character string specifying an optional coefficient, cue, or
#'   category name (default: \code{NULL} for class-wide priors).
#' @param cue_space Character string specifying whether the prior is formulated
#'   in transformed cue space (\code{"transformed"}, the default) or original
#'   physical cue space (\code{"original"}).
#'
#' @return An object of class \code{\link{MVBU_Prior}}.
#' @rdname MVBU_Prior
#' @seealso \code{\link{fit_ideal_adaptor}}, \code{\link{fix_parameter}}
#' @export
parameter_prior <- function(
  declaration,
  class = c(
    "tau_0", "cor_0", "m_0", "weight", "weights", "lapse_rate", "lapse_bias",
    "category_prior", "kappa_0", "nu_0"
  ),
  category = NULL,
  cue = NULL,
  par = NULL,
  cue_space = c("transformed", "original")
) {
  .assert_true(
    !missing(declaration),
    msg = "Argument `declaration` must be provided to `parameter_prior()`."
  )
  class <- match.arg(class)
  if (class == "weights") class <- "weight"
  cue_space <- match.arg(cue_space)

  if (is.null(category)) category <- ""
  if (is.null(cue)) cue <- ""

  # Handle legacy par argument
  if (!is.null(par) && nzchar(trimws(par))) {
    par_val <- trimws(par)
    if (class %in% c("tau_0", "m_0")) {
      if (!nzchar(cue)) cue <- par_val
    } else if (class %in% c("category_prior", "lapse_bias")) {
      if (!nzchar(category)) category <- par_val
    }
  }

  par_str <- ifelse(nzchar(category), category, cue)

  if (is.numeric(declaration)) {
    if (class == "tau_0") {
      declaration <- sprintf("cauchy(0, %s)", as.character(declaration))
    } else if (class == "cor_0") {
      declaration <- sprintf("lkj(%s)", as.character(declaration))
    } else {
      declaration <- as.character(declaration)
    }
  }

  .assert_true(
    is.character(declaration) && length(declaration) == 1L &&
      nzchar(trimws(declaration)),
    msg = "`declaration` must be a non-empty character string or numeric value."
  )
  .assert_true(
    is.character(category) && length(category) == 1L,
    msg = "`category` must be a character string of length 1."
  )
  .assert_true(
    is.character(cue) && length(cue) == 1L,
    msg = "`cue` must be a character string of length 1."
  )

  MVBU_Prior(
    declaration = trimws(declaration),
    par_class = class,
    category = trimws(category),
    cue = trimws(cue),
    par = trimws(par_str),
    cue_space = cue_space
  )
}

# -------------------------------------------------------------------------
# S3 / S7 Methods for MVBU_Prior
# -------------------------------------------------------------------------

#' Combine prior specifications
#'
#' @param ... \code{\link{MVBU_Prior}} objects or lists of them.
#' @param recursive Logical (default: FALSE).
#'
#' @return A single combined \code{\link{MVBU_Prior}} object.
#' @rdname MVBU_Prior
#' @export
c.MVBU_Prior <- function(..., recursive = FALSE) {
  dots <- list(...)
  decl <- character()
  cls <- character()
  cats <- character()
  cues <- character()
  pr <- character()
  spc <- character()

  for (item in dots) {
    if (is.null(item)) next
    if (S7::S7_inherits(item, MVBU_Prior)) {
      decl <- c(decl, item@declaration)
      cls <- c(cls, item@par_class)
      cats <- c(cats, item@category)
      cues <- c(cues, item@cue)
      pr <- c(pr, item@par)
      spc <- c(spc, item@cue_space)
    } else if (is.list(item)) {
      for (sub_item in item) {
        if (S7::S7_inherits(sub_item, MVBU_Prior)) {
          decl <- c(decl, sub_item@declaration)
          cls <- c(cls, sub_item@par_class)
          cats <- c(cats, sub_item@category)
          cues <- c(cues, sub_item@cue)
          pr <- c(pr, sub_item@par)
          spc <- c(spc, sub_item@cue_space)
        }
      }
    } else {
      .stop(
        "All arguments to `c()` for MVBU_Prior must be MVBU_Prior objects."
      )
    }
  }

  MVBU_Prior(
    declaration = decl,
    par_class = cls,
    category = cats,
    cue = cues,
    par = pr,
    cue_space = spc
  )
}

#' Subset prior specifications
#'
#' @param x An \code{\link{MVBU_Prior}} object.
#' @param i Indices specifying elements to extract.
#' @param ... Additional arguments (ignored).
#'
#' @return A subsetted \code{\link{MVBU_Prior}} object.
#' @rdname MVBU_Prior
#' @export
`[.MVBU_Prior` <- function(x, i, ...) {
  MVBU_Prior(
    declaration = x@declaration[i],
    par_class = x@par_class[i],
    category = x@category[i],
    cue = x@cue[i],
    par = x@par[i],
    cue_space = x@cue_space[i]
  )
}

#' Number of prior specifications
#'
#' @param x An \code{\link{MVBU_Prior}} object.
#'
#' @return An integer length.
#' @rdname MVBU_Prior
#' @export
length.MVBU_Prior <- function(x) {
  length(x@declaration)
}

#' @export
names.MVBU_Prior <- function(x) {
  unique(x@par_class)
}

#' Convert prior specification to data.frame
#'
#' @param x An \code{\link{MVBU_Prior}} object.
#' @param row.names Optional row names (default: NULL).
#' @param optional Logical (default: FALSE).
#' @param ... Additional arguments passed to data.frame.
#'
#' @return A \code{data.frame} of prior specifications.
#' @rdname MVBU_Prior
#' @export
as.data.frame.MVBU_Prior <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  ...
) {
  data.frame(
    declaration = x@declaration,
    class = x@par_class,
    category = ifelse(nzchar(x@category), x@category, "(all)"),
    cue = ifelse(nzchar(x@cue), x@cue, "(all)"),
    par = x@par,
    cue_space = x@cue_space,
    stringsAsFactors = FALSE
  )
}

#' Convert prior specification to tibble
#'
#' @param x An \code{\link{MVBU_Prior}} object.
#' @param ... Additional arguments passed to \code{as_tibble}.
#'
#' @return A \code{\link[tibble]{tibble}} of prior specifications.
#' @rdname MVBU_Prior
#' @export
as_tibble.MVBU_Prior <- function(x, ...) {
  tibble::tibble(
    declaration = x@declaration,
    class = x@par_class,
    category = ifelse(nzchar(x@category), x@category, "(all)"),
    cue = ifelse(nzchar(x@cue), x@cue, "(all)"),
    par = x@par,
    cue_space = x@cue_space
  )
}

#' Print prior specifications
#'
#' @param x An \code{\link{MVBU_Prior}} object.
#' @param ... Additional arguments passed to print.
#'
#' @return The object \code{x}, invisibly.
#' @rdname MVBU_Prior
#' @export
print.MVBU_Prior <- function(x, ...) {
  cat("MVBeliefUpdatr Prior Specifications:\n")
  if (length(x) == 0L) {
    cat("  (empty prior specification)\n")
    return(invisible(x))
  }
  df <- as.data.frame(x)
  print(df, row.names = FALSE, ...)
  invisible(x)
}

S7::`method<-`(print, MVBU_Prior, value = print.MVBU_Prior)
S7::`method<-`(as.data.frame, MVBU_Prior, value = as.data.frame.MVBU_Prior)
S7::`method<-`(c, MVBU_Prior, value = function(x, ...) c.MVBU_Prior(x, ...))
S7::`method<-`(`[`, MVBU_Prior, value = function(x, i, ...) {
  `[.MVBU_Prior`(x, i, ...)
})
S7::`method<-`(length, MVBU_Prior, value = length.MVBU_Prior)
S7::`method<-`(names, MVBU_Prior, value = names.MVBU_Prior)

# -------------------------------------------------------------------------
# Internal Helper Functions for Prior Resolution
# -------------------------------------------------------------------------

#' Extract scale parameter from a tau prior declaration
#'
#' @param declaration Character string.
#' @return Numeric scale value.
#' @keywords internal
#' @noRd
.extract_tau_scale <- function(declaration) {
  if (grepl("cauchy\\s*\\(", declaration, ignore.case = TRUE)) {
    val_str <- sub(
      "^.*cauchy\\s*\\([^,]+,\\s*([^\\)]+)\\).*$",
      "\\1",
      declaration,
      ignore.case = TRUE
    )
    val <- suppressWarnings(as.numeric(trimws(val_str)))
  } else if (grepl("^[0-9.]+$", declaration)) {
    val <- suppressWarnings(as.numeric(trimws(declaration)))
  } else {
    val <- NA_real_
  }
  .assert_true(
    !is.na(val) && val > 0,
    msg = sprintf(
      paste0(
        "Could not extract a positive scale parameter from tau prior ",
        "declaration '%s'. Expected format: 'cauchy(0, <scale>)' or a ",
        "positive number."
      ),
      declaration
    )
  )
  val
}

#' Extract eta parameter from an LKJ correlation prior declaration
#'
#' @param declaration Character string.
#' @return Numeric eta value.
#' @keywords internal
#' @noRd
.extract_lkj_eta <- function(declaration) {
  pattern <- "(?:lkj|lkj_corr_cholesky)\\s*\\("
  if (grepl(pattern, declaration, ignore.case = TRUE)) {
    val_str <- sub(
      "^.*(?:lkj|lkj_corr_cholesky)\\s*\\(\\s*([^\\)]+)\\).*$",
      "\\1",
      declaration,
      ignore.case = TRUE
    )
    val <- suppressWarnings(as.numeric(trimws(val_str)))
  } else if (grepl("^[0-9.]+$", declaration)) {
    val <- suppressWarnings(as.numeric(trimws(declaration)))
  } else {
    val <- NA_real_
  }
  .assert_true(
    !is.na(val) && val > 0,
    msg = sprintf(
      paste0(
        "Could not extract a positive eta parameter from cor prior ",
        "declaration '%s'. Expected format: 'lkj(<eta>)' or a positive number."
      ),
      declaration
    )
  )
  val
}

#' Extract location and scale from a Cauchy prior declaration
#' @keywords internal
#' @noRd
.extract_cauchy_params <- function(
  declaration,
  default_loc = 0,
  default_scale = 5
) {
  declaration <- trimws(declaration)
  if (grepl("cauchy\\s*\\(", declaration, ignore.case = TRUE)) {
    inner <- sub(
      "^.*cauchy\\s*\\((.*)\\).*$",
      "\\1",
      declaration,
      ignore.case = TRUE
    )
    parts <- trimws(strsplit(inner, ",")[[1]])
    if (length(parts) == 1L) {
      loc <- default_loc
      scale <- as.numeric(parts[1])
    } else {
      loc <- as.numeric(parts[1])
      scale <- as.numeric(parts[2])
    }
  } else if (grepl("^[0-9.]+$", declaration)) {
    loc <- default_loc
    scale <- as.numeric(declaration)
  } else {
    loc <- default_loc
    scale <- default_scale
  }
  list(loc = loc, scale = scale)
}

#' Extract mean and sd from a Normal prior declaration
#' @keywords internal
#' @noRd
.extract_normal_params <- function(
  declaration,
  default_mean = 0,
  default_sd = 1
) {
  declaration <- trimws(declaration)
  if (grepl("normal\\s*\\(", declaration, ignore.case = TRUE)) {
    inner <- sub(
      "^.*normal\\s*\\((.*)\\).*$",
      "\\1",
      declaration,
      ignore.case = TRUE
    )
    parts <- trimws(strsplit(inner, ",")[[1]])
    mean_val <- as.numeric(parts[1])
    sd_val <- if (length(parts) >= 2L) as.numeric(parts[2]) else default_sd
  } else if (grepl("^[0-9.-]+$", declaration)) {
    mean_val <- as.numeric(declaration)
    sd_val <- default_sd
  } else {
    mean_val <- default_mean
    sd_val <- default_sd
  }
  list(mean = mean_val, sd = sd_val)
}

#' Extract alpha and beta from a Beta prior declaration
#' @keywords internal
#' @noRd
.extract_beta_params <- function(
  declaration,
  default_alpha = 1,
  default_beta = 19
) {
  declaration <- trimws(declaration)
  if (grepl("beta\\s*\\(", declaration, ignore.case = TRUE)) {
    inner <- sub(
      "^.*beta\\s*\\((.*)\\).*$",
      "\\1",
      declaration,
      ignore.case = TRUE
    )
    parts <- trimws(strsplit(inner, ",")[[1]])
    alpha <- as.numeric(parts[1])
    beta <- as.numeric(parts[2])
  } else {
    alpha <- default_alpha
    beta <- default_beta
  }
  list(alpha = alpha, beta = beta)
}

#' Extract Dirichlet concentration vector from a prior declaration
#' @keywords internal
#' @noRd
.extract_dirichlet_alphas <- function(declaration, n, default_alpha = 1) {
  declaration <- trimws(declaration)
  if (grepl("dirichlet\\s*\\(", declaration, ignore.case = TRUE)) {
    inner <- sub(
      "^.*dirichlet\\s*\\((.*)\\).*$",
      "\\1",
      declaration,
      ignore.case = TRUE
    )
    vals <- tryCatch({
      eval(parse(text = inner), envir = baseenv())
    }, error = function(e) {
      suppressWarnings(as.numeric(trimws(strsplit(inner, "[, ]+")[[1]])))
    })
    if (length(vals) == 1L) {
      rep(as.numeric(vals), n)
    } else {
      as.numeric(vals)
    }
  } else {
    rep(default_alpha, n)
  }
}

#' Resolve full dual-space prior specifications and Stan hyperparameters
#'
#' Resolves user-declared priors against model defaults, translating between
#' transformed space (used directly in Stan) and original physical cue space
#' (used for reporting and transparency).
#'
#' @param priors User-provided prior specification (\code{NULL}, list, or
#'   \code{\link{MVBU_Prior}}).
#' @param stanmodel Character string identifying the Stan model.
#' @param cues Character vector of cue variable names.
#' @param category_labels Character vector of category names.
#' @param transform Affine transform object of class
#'   \code{\link{MVBU_TransformInformation}}.
#' @param max_n_exposure Maximum observation count across exposure groups.
#'
#' @return A named list with components:
#'   \item{original}{\code{\link{MVBU_Prior}} in original cue space.}
#'   \item{transformed}{\code{\link{MVBU_Prior}} in transformed space.}
#'   \item{stan}{Named list of resolved prior hyperparameters for Stan.}
#'
#' @keywords internal
#' @noRd
.resolve_priors <- function(
  priors,
  stanmodel,
  cues,
  category_labels,
  transform,
  max_n_exposure = NULL
) {
  n_cues <- length(cues)
  n_categories <- length(category_labels)

  # Determine cue scale multipliers: how 1 transformed unit maps to original
  if (length(transform@shift) == 1L) {
    orig_scale_mult <- rep(
      abs(1 / as.numeric(transform@SCALE)),
      n_cues
    )
  } else if (is.matrix(transform@SCALE)) {
    inv_scale_mat <- solve(transform@SCALE)
    orig_scale_mult <- abs(diag(inv_scale_mat))
  } else {
    orig_scale_mult <- rep(1, n_cues)
  }
  names(orig_scale_mult) <- cues

  # Convert user input to MVBU_Prior
  if (is.null(priors) || length(priors) == 0L) {
    user_prior <- MVBU_Prior()
  } else if (S7::S7_inherits(priors, MVBU_Prior)) {
    user_prior <- priors
  } else if (is.list(priors) && all(vapply(priors, function(p) {
    S7::S7_inherits(p, MVBU_Prior)
  }, logical(1)))) {
    user_prior <- do.call(c, priors)
  } else {
    .stop("Argument `priors` must be an MVBU_Prior object or NULL.")
  }

  # Validate par parameters
  if (length(user_prior) > 0L) {
    for (i in seq_along(user_prior@declaration)) {
      cls <- user_prior@par_class[i]
      pr <- user_prior@par[i]
      if (nzchar(pr)) {
        if (cls %in% c("tau_0", "m_0") && !pr %in% cues) {
          .stop(
            sprintf(
              paste0(
                "Parameter par = '%s' is not recognized. For class '%s', ",
                "par must be empty or one of the cues: %s."
              ),
              pr,
              cls,
              paste(cues, collapse = ", ")
            )
          )
        }
        if (cls %in% c("category_prior", "lapse_bias") &&
            !pr %in% category_labels) {
          .stop(
            sprintf(
              paste0(
                "Parameter par = '%s' is not recognized. For class '%s', ",
                "par must be empty or one of the categories: %s."
              ),
              pr,
              cls,
              paste(category_labels, collapse = ", ")
            )
          )
        }
      }
    }
  }

  # Resolve tau (default in transformed space: 5)
  stan_tau_scale <- rep(5, n_cues)
  names(stan_tau_scale) <- cues
  stan_tau_loc <- rep(0, n_cues)

  tau_idx <- which(user_prior@par_class == "tau_0")
  for (idx in tau_idx) {
    pr <- user_prior@par[idx]
    spc <- user_prior@cue_space[idx]
    sc <- .extract_tau_scale(user_prior@declaration[idx])
    if (nzchar(pr)) {
      if (spc == "original") {
        stan_tau_scale[pr] <- sc / orig_scale_mult[pr]
      } else {
        stan_tau_scale[pr] <- sc
      }
    } else {
      if (spc == "original") {
        stan_tau_scale <- rep(sc, n_cues) / orig_scale_mult
      } else {
        stan_tau_scale <- rep(sc, n_cues)
      }
    }
  }

  # Resolve cor (L_omega_eta, default: 1)
  stan_L_omega_eta <- 1
  cor_idx <- which(user_prior@par_class == "cor_0")
  if (length(cor_idx) > 0L) {
    last_cor_idx <- cor_idx[length(cor_idx)]
    stan_L_omega_eta <- .extract_lkj_eta(user_prior@declaration[last_cor_idx])
  }

  # Resolve pseudocounts kappa_0 and nu_0
  max_n <- if (!is.null(max_n_exposure) && length(max_n_exposure) > 0L) {
    max(max_n_exposure, na.rm = TRUE)
  } else {
    10
  }
  sigma_kappanu <- max(max_n * 4, 10)

  prior_loc_of_kappa <- rep(0, n_categories)
  prior_scale_of_kappa <- rep(sigma_kappanu, n_categories)
  prior_loc_of_nu <- rep(0, n_categories)
  prior_scale_of_nu <- rep(sigma_kappanu, n_categories)

  kappa_idx <- which(user_prior@par_class == "kappa_0")
  for (idx in kappa_idx) {
    res <- .extract_cauchy_params(
      user_prior@declaration[idx],
      default_loc = 0,
      default_scale = sigma_kappanu
    )
    cat_match <- user_prior@category[idx]
    if (nzchar(cat_match) && cat_match %in% category_labels) {
      c_i <- match(cat_match, category_labels)
      prior_loc_of_kappa[c_i] <- res$loc
      prior_scale_of_kappa[c_i] <- res$scale
    } else {
      prior_loc_of_kappa[] <- res$loc
      prior_scale_of_kappa[] <- res$scale
    }
  }

  nu_idx <- which(user_prior@par_class == "nu_0")
  for (idx in nu_idx) {
    res <- .extract_cauchy_params(
      user_prior@declaration[idx],
      default_loc = 0,
      default_scale = sigma_kappanu
    )
    cat_match <- user_prior@category[idx]
    if (nzchar(cat_match) && cat_match %in% category_labels) {
      c_i <- match(cat_match, category_labels)
      prior_loc_of_nu[c_i] <- res$loc
      prior_scale_of_nu[c_i] <- res$scale
    } else {
      prior_loc_of_nu[] <- res$loc
      prior_scale_of_nu[] <- res$scale
    }
  }

  # Resolve lapse_rate (default: beta(1, 19))
  prior_alpha_of_lapse_rate <- 1
  prior_beta_of_lapse_rate <- 19
  lr_idx <- which(user_prior@par_class == "lapse_rate")
  if (length(lr_idx) > 0L) {
    res <- .extract_beta_params(
      user_prior@declaration[lr_idx[length(lr_idx)]],
      default_alpha = 1,
      default_beta = 19
    )
    prior_alpha_of_lapse_rate <- res$alpha
    prior_beta_of_lapse_rate <- res$beta
  }

  # Resolve lapse_bias (default: dirichlet(rep(1, M)))
  prior_alphas_of_lapse_bias <- rep(1, n_categories)
  lb_idx <- which(user_prior@par_class == "lapse_bias")
  if (length(lb_idx) > 0L) {
    prior_alphas_of_lapse_bias <- .extract_dirichlet_alphas(
      user_prior@declaration[lb_idx[length(lb_idx)]],
      n = n_categories,
      default_alpha = 1
    )
  }

  # Resolve category_prior (default: dirichlet(rep(1, M)))
  prior_alphas_of_category_prior <- rep(1, n_categories)
  cp_idx <- which(user_prior@par_class == "category_prior")
  if (length(cp_idx) > 0L) {
    prior_alphas_of_category_prior <- .extract_dirichlet_alphas(
      user_prior@declaration[cp_idx[length(cp_idx)]],
      n = n_categories,
      default_alpha = 1
    )
  }

  # Resolve prior_mu_of_m (default: normal(0, tau) in transformed space)
  prior_mu_of_m <- matrix(0, nrow = n_categories, ncol = n_cues)
  m_idx <- which(user_prior@par_class == "m_0")
  for (idx in m_idx) {
    res <- .extract_normal_params(user_prior@declaration[idx])
    cat_match <- user_prior@category[idx]
    cue_match <- user_prior@cue[idx]
    c_inds <- if (nzchar(cat_match) && cat_match %in% category_labels) {
      match(cat_match, category_labels)
    } else {
      seq_len(n_categories)
    }
    k_inds <- if (nzchar(cue_match) && cue_match %in% cues) {
      match(cue_match, cues)
    } else {
      seq_len(n_cues)
    }
    prior_mu_of_m[c_inds, k_inds] <- res$mean
  }

  # Construct dual-space representations
  is_weighted <- grepl("weighted", stanmodel)
  has_cor <- (n_cues >= 2L) && !startsWith(stanmodel, "NIX") &&
    !startsWith(stanmodel, "MNIX")

  # Transformed space
  trans_decl <- character()
  trans_cls <- character()
  trans_par <- character()
  trans_spc <- character()

  all_same_tau <- length(unique(stan_tau_scale)) == 1L
  if (all_same_tau) {
    trans_decl <- c(
      trans_decl,
      sprintf("cauchy(0, %g)", round(stan_tau_scale[1], 4))
    )
    trans_cls <- c(trans_cls, "tau_0")
    trans_par <- c(trans_par, "")
    trans_spc <- c(trans_spc, "transformed")
  } else {
    for (k in seq_along(cues)) {
      trans_decl <- c(
        trans_decl,
        sprintf("cauchy(0, %g)", round(stan_tau_scale[k], 4))
      )
      trans_cls <- c(trans_cls, "tau_0")
      trans_par <- c(trans_par, cues[k])
      trans_spc <- c(trans_spc, "transformed")
    }
  }

  if (has_cor) {
    trans_decl <- c(
      trans_decl,
      sprintf("lkj(%g)", round(stan_L_omega_eta, 4))
    )
    trans_cls <- c(trans_cls, "cor_0")
    trans_par <- c(trans_par, "")
    trans_spc <- c(trans_spc, "transformed")
  }

  trans_decl <- c(trans_decl, "normal(0, tau_0)")
  trans_cls <- c(trans_cls, "m_0")
  trans_par <- c(trans_par, "")
  trans_spc <- c(trans_spc, "transformed")

  if (is_weighted) {
    trans_decl <- c(trans_decl, "dirichlet(1)")
    trans_cls <- c(trans_cls, "weight")
    trans_par <- c(trans_par, "")
    trans_spc <- c(trans_spc, "transformed")
  }

  trans_decl <- c(trans_decl, "dirichlet(1)")
  trans_cls <- c(trans_cls, "category_prior")
  trans_par <- c(trans_par, "")
  trans_spc <- c(trans_spc, "transformed")

  trans_decl <- c(trans_decl, "dirichlet(1)")
  trans_cls <- c(trans_cls, "lapse_bias")
  trans_par <- c(trans_par, "")
  trans_spc <- c(trans_spc, "transformed")

  # Original space
  orig_decl <- character()
  orig_cls <- character()
  orig_par <- character()
  orig_spc <- character()

  orig_tau_scales <- stan_tau_scale * orig_scale_mult
  all_same_orig_tau <- length(unique(round(orig_tau_scales, 4))) == 1L
  if (all_same_orig_tau && n_cues > 1L) {
    orig_decl <- c(
      orig_decl,
      sprintf("cauchy(0, %g)", round(orig_tau_scales[1], 4))
    )
    orig_cls <- c(orig_cls, "tau_0")
    orig_par <- c(orig_par, "")
    orig_spc <- c(orig_spc, "original")
  } else {
    for (k in seq_along(cues)) {
      orig_decl <- c(
        orig_decl,
        sprintf("cauchy(0, %g)", round(orig_tau_scales[k], 4))
      )
      orig_cls <- c(orig_cls, "tau_0")
      orig_par <- c(orig_par, cues[k])
      orig_spc <- c(orig_spc, "original")
    }
  }

  if (has_cor) {
    orig_decl <- c(
      orig_decl,
      sprintf("lkj(%g)", round(stan_L_omega_eta, 4))
    )
    orig_cls <- c(orig_cls, "cor_0")
    orig_par <- c(orig_par, "")
    orig_spc <- c(orig_spc, "original")
  }

  shift_vals <- if (length(transform@shift) == 1L) {
    rep(transform@shift, n_cues)
  } else {
    transform@shift
  }
  for (k in seq_along(cues)) {
    mean_center <- -shift_vals[k]
    orig_decl <- c(
      orig_decl,
      sprintf("normal(%g, tau_0)", round(mean_center, 4))
    )
    orig_cls <- c(orig_cls, "m_0")
    orig_par <- c(orig_par, cues[k])
    orig_spc <- c(orig_spc, "original")
  }

  if (is_weighted) {
    orig_decl <- c(orig_decl, "dirichlet(1)")
    orig_cls <- c(orig_cls, "weight")
    orig_par <- c(orig_par, "")
    orig_spc <- c(orig_spc, "original")
  }

  orig_decl <- c(orig_decl, "dirichlet(1)")
  orig_cls <- c(orig_cls, "category_prior")
  orig_par <- c(orig_par, "")
  orig_spc <- c(orig_spc, "original")

  orig_decl <- c(orig_decl, "dirichlet(1)")
  orig_cls <- c(orig_cls, "lapse_bias")
  orig_par <- c(orig_par, "")
  orig_spc <- c(orig_spc, "original")

  list(
    original = MVBU_Prior(
      declaration = orig_decl,
      par_class = orig_cls,
      par = orig_par,
      cue_space = orig_spc
    ),
    transformed = MVBU_Prior(
      declaration = trans_decl,
      par_class = trans_cls,
      par = trans_par,
      cue_space = trans_spc
    ),
    stan = list(
      prior_alphas_of_category_prior = as.array(
        as.numeric(prior_alphas_of_category_prior)
      ),
      prior_alpha_of_lapse_rate = as.numeric(prior_alpha_of_lapse_rate),
      prior_beta_of_lapse_rate = as.numeric(prior_beta_of_lapse_rate),
      prior_alphas_of_lapse_bias = as.array(
        as.numeric(prior_alphas_of_lapse_bias)
      ),
      prior_mu_of_m = prior_mu_of_m,
      prior_loc_of_tau = as.array(as.numeric(stan_tau_loc)),
      prior_scale_of_tau = as.array(as.numeric(stan_tau_scale)),
      prior_eta_of_cor = as.numeric(stan_L_omega_eta),
      prior_loc_of_kappa = as.array(as.numeric(prior_loc_of_kappa)),
      prior_scale_of_kappa = as.array(as.numeric(prior_scale_of_kappa)),
      prior_loc_of_nu = as.array(as.numeric(prior_loc_of_nu)),
      prior_scale_of_nu = as.array(as.numeric(prior_scale_of_nu))
    )
  )
}

# -------------------------------------------------------------------------
# S7 Method for get_stan_priors on MVBU_Prior
# -------------------------------------------------------------------------

#' @rdname get_parameters
#' @export
S7::method(get_stan_priors, MVBU_Prior) <- function(
  x,
  space = c("original", "transformed"),
  ...
) {
  space <- match.arg(space)
  idx <- which(x@cue_space == space)
  if (length(idx) == 0L) {
    MVBU_Prior()
  } else {
    x[idx]
  }
}
