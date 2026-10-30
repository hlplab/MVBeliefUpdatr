#' @include asserts.R
#' @include S7-class.R
#' @include S7-generics.R
#' @importFrom S7 new_class new_object props class_character class_list
#' @importFrom S7 S7_inherits method
#' @importFrom tibble tibble as_tibble
NULL

# -------------------------------------------------------------------------
# S7 Fixed Parameter Class and Constructor
# -------------------------------------------------------------------------

#' S7 class for fixed cognitive parameter specifications
#'
#' Represents fixed parameter specifications in a dual-space (original vs.
#' transformed cue space) framework, supporting declarative parameter fixing
#' across global, category-specific, and cue-specific cognitive parameters.
#'
#' @param par_class Character vector of parameter classes (e.g., \code{"m_0"},
#'   \code{"S_0"}, \code{"kappa_0"}, \code{"nu_0"}, \code{"expected_mu_0"},
#'   \code{"expected_Sigma_0"}, \code{"mu_0"}, \code{"Sigma_0"},
#'   \code{"lapse_rate"}, \code{"lapse_bias"}, \code{"category_prior"},
#'   \code{"weights"}).
#' @param value List of fixed parameter values (scalars, vectors, or matrices).
#' @param category Character vector specifying category names (default:
#'   \code{""} for class-wide or category-invariant parameters).
#' @param cue Character vector specifying cue names (default: \code{""} for
#'   all cues).
#' @param cue_space Character vector specifying the cue space:
#'   \code{"original"} (default) or \code{"transformed"}.
#'
#' @docType class
#' @name MVBU_FixedParameter
#' @rdname MVBU_FixedParameter
#' @export
MVBU_FixedParameter <- S7::new_class(
  "MVBU_FixedParameter",
  package = NULL,
  parent = MVBU_Object,
  properties = list(
    par_class = S7::class_character,
    value = S7::class_list,
    category = S7::class_character,
    cue = S7::class_character,
    cue_space = S7::class_character
  ),
  constructor = function(
    par_class = character(),
    value = list(),
    category = character(),
    cue = character(),
    cue_space = character()
  ) {
    if (length(par_class) == 0L) {
      return(
        S7::new_object(
          MVBU_Object(),
          par_class = character(),
          value = list(),
          category = character(),
          cue = character(),
          cue_space = character()
        )
      )
    }
    n <- length(par_class)
    if (length(value) == 1L && n > 1L) value <- rep(value, n)
    if (length(category) == 0L) category <- rep("", n)
    if (length(category) == 1L && n > 1L) category <- rep(category, n)
    if (length(cue) == 0L) cue <- rep("", n)
    if (length(cue) == 1L && n > 1L) cue <- rep(cue, n)
    if (length(cue_space) == 0L) cue_space <- rep("original", n)
    if (length(cue_space) == 1L && n > 1L) cue_space <- rep(cue_space, n)

    S7::new_object(
      MVBU_Object(),
      par_class = as.character(par_class),
      value = as.list(value),
      category = as.character(category),
      cue = as.character(cue),
      cue_space = as.character(cue_space)
    )
  },
  validator = function(self) {
    n <- length(self@par_class)
    if (length(self@value) != n || length(self@category) != n ||
      length(self@cue) != n || length(self@cue_space) != n) {
      return("All properties of MVBU_FixedParameter must have the same length.")
    }
    allowed_classes <- c(
      "m_0", "S_0", "kappa_0", "nu_0",
      "expected_mu_0", "expected_Sigma_0", "mu_0", "Sigma_0",
      "lapse_rate", "lapse_bias", "category_prior", "weights"
    )
    invalid_classes <- setdiff(self@par_class, allowed_classes)
    if (length(invalid_classes) > 0L) {
      return(
        sprintf(
          "Invalid fixed parameter class: '%s'. Allowed classes are: %s.",
          invalid_classes[1],
          paste(allowed_classes, collapse = ", ")
        )
      )
    }
    allowed_spaces <- c("original", "transformed")
    invalid_spaces <- setdiff(self@cue_space, allowed_spaces)
    if (length(invalid_spaces) > 0L) {
      return(
        sprintf(
          "Invalid cue_space: '%s'. Must be 'original' or 'transformed'.",
          invalid_spaces[1]
        )
      )
    }
    NULL
  }
)

#' Fix a model parameter to a known value
#'
#' Constructs a fixed parameter specification for
#' \code{\link{fit_ideal_adaptor}}.
#' Supports fixing parameters globally, per-category, or per-cue, in either
#' original or transformed cue space.
#'
#' @param value Numeric value, vector, matrix, or character assignment
#'   (e.g., \code{"lapse_rate = 0.05"}).
#' @param class Character string specifying parameter class:
#'   \code{"m_0"}, \code{"S_0"}, \code{"kappa_0"}, \code{"nu_0"},
#'   \code{"expected_mu_0"}, \code{"expected_Sigma_0"}, \code{"mu_0"},
#'   \code{"Sigma_0"}, \code{"lapse_rate"}, \code{"lapse_bias"},
#'   \code{"category_prior"}, or \code{"weights"}.
#' @param category Character string specifying category name (default:
#'   \code{""} or \code{NULL} for all categories).
#' @param cue Character string specifying cue name (default: \code{""} or
#'   \code{NULL} for all cues).
#' @param cue_space Character string specifying cue space:
#'   \code{"original"} (default) or \code{"transformed"}.
#'
#' @section Category Moment Parameter Interpretations:
#' When fixing category representation parameters, note the distinction
#' between prior hyperparameters, expected prior moments, and known category
#' parameters:
#' \describe{
#'   \item{\code{"m_0"} and \code{"S_0"}}{
#'     Prior hyperparameters of the conjugate prior distribution
#'     (Normal-Inverse-Wishart or Normal-Inverse-Chi-squared). \code{"m_0"} is
#'     the prior category mean vector (center of the normal prior on
#'     \eqn{\mu_0}), and \code{"S_0"} is the prior scatter matrix.
#'   }
#'   \item{\code{"expected_mu_0"} and \code{"expected_Sigma_0"}}{
#'     Expected moments of category representations under the prior. For
#'     category means, \eqn{E[\mu_0] = m_0}, so \code{"expected_mu_0"} is
#'     canonicalized directly to \code{"m_0"}. For covariance matrices,
#'     \eqn{E[\Sigma_0] = S_0 / (\nu_0 - K - 1)} for NIW (or
#'     \eqn{E[\Sigma_0] = S_0 / (\nu_0 - 2)} for NIX/MNIX). Specifying
#'     \code{"expected_Sigma_0"} allows users to fix the intuitive expected
#'     covariance. Stan input preparation automatically rescales this value to
#'     the conjugate scatter matrix \code{"S_0"} given \nu_0.
#'   }
#'   \item{\code{"mu_0"} and \code{"Sigma_0"}}{
#'     Fixed, known category representation parameters entering the likelihood
#'     directly without prior belief updating or uncertainty of that parameter.
#'   }
#' }
#'
#' @return An object of class \code{\link{MVBU_FixedParameter}}.
#' @rdname MVBU_FixedParameter
#' @seealso \code{\link{prior}}, \code{\link{fit_ideal_adaptor}}
#' @export
fix_parameter <- function(
  value,
  class = c(
    "m_0", "S_0", "kappa_0", "nu_0",
    "expected_mu_0", "expected_Sigma_0", "mu_0", "Sigma_0",
    "lapse_rate", "lapse_bias", "category_prior", "weights"
  ),
  category = NULL,
  cue = NULL,
  cue_space = c("original", "transformed")
) {
  .assert_true(
    !missing(value),
    msg = "Argument `value` must be provided to `fix_parameter()`."
  )
  if (is.null(category)) category <- ""
  if (is.null(cue)) cue <- ""
  cue_space <- match.arg(cue_space)

  all_classes <- c(
    "m_0", "S_0", "kappa_0", "nu_0",
    "expected_mu_0", "expected_Sigma_0", "mu_0", "Sigma_0",
    "lapse_rate", "lapse_bias", "category_prior", "weights"
  )
  if (is.character(value) && length(value) == 1L && value %in% all_classes &&
    !missing(class) && !is.character(class)) {
    tmp_val <- class
    class <- value
    value <- tmp_val
  }

  # Support "name = value" string specification
  if (is.character(value) && length(value) == 1L && grepl("=", value)) {
    lhs <- trimws(sub("=.*$", "", value))
    rhs_str <- trimws(sub("^[^=]*=", "", value))
    parsed_val <- eval(parse(text = rhs_str), envir = parent.frame())
    class <- lhs
    value <- parsed_val
  } else {
    class <- match.arg(class)
  }

  if (class %in% c("mu_0", "Sigma_0")) {
    .stop(
      sprintf(
        "Fixing known '%s' is not yet implemented (deferred to Stage 6).",
        class
      )
    )
  }

  # Canonicalize parameter class names
  if (class == "expected_mu_0") {
    class <- "m_0"
  }

  .assert_true(
    is.character(category) && length(category) == 1L,
    msg = "`category` must be a character string of length 1."
  )
  .assert_true(
    is.character(cue) && length(cue) == 1L,
    msg = "`cue` must be a character string of length 1."
  )

  category <- trimws(category)
  cue <- trimws(cue)

  if (nzchar(category) &&
    class %in% c("lapse_rate", "category_prior", "lapse_bias")) {
    .stop(
      sprintf(
        "Parameter '%s' is a global parameter and cannot be category-specific.",
        class
      )
    )
  }

  if (nzchar(cue) &&
    !class %in% c(
      "m_0", "expected_mu_0", "mu_0",
      "S_0", "expected_Sigma_0", "Sigma_0",
      "weights"
    )) {
    .stop(
      sprintf(
        "Parameter '%s' cannot be cue-specific.",
        class
      )
    )
  }

  MVBU_FixedParameter(
    par_class = class,
    value = list(value),
    category = category,
    cue = cue,
    cue_space = cue_space
  )
}

# -------------------------------------------------------------------------
# S3 / S7 Methods for MVBU_FixedParameter
# -------------------------------------------------------------------------

#' @export
c.MVBU_FixedParameter <- function(..., recursive = FALSE) {
  dots <- list(...)
  dots <- dots[!vapply(dots, is.null, logical(1))]
  if (length(dots) == 0L) {
    return(MVBU_FixedParameter())
  }
  for (item in dots) {
    .assert_true(
      S7::S7_inherits(item, MVBU_FixedParameter),
      msg = "All arguments to `c()` must be MVBU_FixedParameter objects."
    )
  }
  all_classes <- unlist(lapply(dots, function(x) x@par_class))
  all_values <- unlist(lapply(dots, function(x) x@value), recursive = FALSE)
  all_cats <- unlist(lapply(dots, function(x) x@category))
  all_cues <- unlist(lapply(dots, function(x) x@cue))
  all_spaces <- unlist(lapply(dots, function(x) x@cue_space))

  MVBU_FixedParameter(
    par_class = all_classes,
    value = all_values,
    category = all_cats,
    cue = all_cues,
    cue_space = all_spaces
  )
}

#' @export
length.MVBU_FixedParameter <- function(x) {
  length(x@par_class)
}

#' @export
names.MVBU_FixedParameter <- function(x) {
  unique(x@par_class)
}

#' @export
print.MVBU_FixedParameter <- function(x, ...) {
  n <- length(x)
  if (n == 0L) {
    cat("empty MVBU_FixedParameter()\n")
    return(invisible(x))
  }
  cat(sprintf("MVBU_FixedParameter with %d specification(s):\n", n))
  print(tibble::as_tibble(x), ...)
  invisible(x)
}

#' @importFrom tibble as_tibble
#' @export
as_tibble.MVBU_FixedParameter <- function(x, ...) {
  n <- length(x)
  if (n == 0L) {
    return(
      tibble::tibble(
        class = character(),
        category = character(),
        cue = character(),
        cue_space = character(),
        value = character()
      )
    )
  }
  formatted_vals <- vapply(x@value, function(v) {
    if (is.matrix(v)) {
      sprintf("[%dx%d matrix]", nrow(v), ncol(v))
    } else if (length(v) > 1L) {
      paste0("c(", paste(round(v, 3), collapse = ", "), ")")
    } else {
      as.character(round(v, 4))
    }
  }, character(1))

  tibble::tibble(
    class = x@par_class,
    category = ifelse(nzchar(x@category), x@category, "(all)"),
    cue = ifelse(nzchar(x@cue), x@cue, "(all)"),
    cue_space = x@cue_space,
    value = formatted_vals
  )
}

#' Subset fixed parameter specifications
#'
#' @param x An \code{\link{MVBU_FixedParameter}} object.
#' @param i Indices specifying elements to extract.
#' @param ... Additional arguments (ignored).
#'
#' @return A subsetted \code{\link{MVBU_FixedParameter}} object.
#' @rdname MVBU_FixedParameter
#' @export
`[.MVBU_FixedParameter` <- function(x, i, ...) {
  MVBU_FixedParameter(
    par_class = x@par_class[i],
    value = x@value[i],
    category = x@category[i],
    cue = x@cue[i],
    cue_space = x@cue_space[i]
  )
}

#' @export
as.data.frame.MVBU_FixedParameter <- function(x, ...) {
  as.data.frame(tibble::as_tibble(x, ...))
}

S7::method(print, MVBU_FixedParameter) <- print.MVBU_FixedParameter
S7::method(as.data.frame, MVBU_FixedParameter) <-
  as.data.frame.MVBU_FixedParameter
S7::method(c, MVBU_FixedParameter) <- function(x, ...) {
  c.MVBU_FixedParameter(x, ...)
}
S7::method(`[`, MVBU_FixedParameter) <- function(x, i, ...) {
  `[.MVBU_FixedParameter`(x, i, ...)
}
S7::method(length, MVBU_FixedParameter) <- function(x) {
  length.MVBU_FixedParameter(x)
}
S7::method(names, MVBU_FixedParameter) <- function(x) {
  names.MVBU_FixedParameter(x)
}

# -------------------------------------------------------------------------
# Internal Conversion & Resolution Helpers
# -------------------------------------------------------------------------

#' Convert list or existing MVBU_FixedParameter into MVBU_FixedParameter
#' @keywords internal
#' @noRd
.as_fixed_parameters <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(MVBU_FixedParameter())
  }
  if (S7::S7_inherits(x, MVBU_FixedParameter)) {
    return(x)
  }
  if (is.list(x) && all(vapply(x, function(i) {
    S7::S7_inherits(i, MVBU_FixedParameter)
  }, logical(1)))) {
    return(do.call(c, x))
  }
  .stop(
    "Argument `fixed_parameters` must be an MVBU_FixedParameter object ",
    "or NULL."
  )
}

#' Resolve fixed parameter specifications into structured Stan inputs
#'
#' @param fixed An \code{\link{MVBU_FixedParameter}} object or list.
#' @param categories Character vector of category names.
#' @param cues Character vector of cue names.
#' @param transform Affine transform object of class
#'   \code{\link{MVBU_TransformInformation}}.
#'
#' @return A list with structured partitioned components for Stan.
#' @keywords internal
#' @noRd
.resolve_fixed_parameters <- function(
  fixed,
  categories,
  cues,
  transform
) {
  fixed <- .as_fixed_parameters(fixed)
  n_categories <- length(categories)
  n_cues <- length(cues)

  # 1. Category prior
  cp_idx <- which(fixed@par_class == "category_prior")
  category_prior <- if (length(cp_idx) > 0L) {
    val <- fixed@value[[cp_idx[length(cp_idx)]]]
    .assert_true(
      length(val) == n_categories,
      msg = sprintf(
        "Fixed category_prior must be a simplex of length %d.",
        n_categories
      )
    )
    list(
      known = 1L,
      data = matrix(as.numeric(val), nrow = 1, ncol = n_categories)
    )
  } else {
    list(
      known = 0L,
      data = matrix(0, nrow = 0, ncol = n_categories)
    )
  }

  # 2. Lapse rate
  lr_idx <- which(fixed@par_class == "lapse_rate")
  lapse_rate <- if (length(lr_idx) > 0L) {
    val <- as.numeric(fixed@value[[lr_idx[length(lr_idx)]]][1])
    .assert_true(
      val >= 0 && val <= 1,
      msg = "Fixed lapse_rate must be between 0 and 1."
    )
    list(
      known = 1L,
      data = array(val, dim = 1)
    )
  } else {
    list(
      known = 0L,
      data = numeric(0)
    )
  }

  # 3. Lapse bias
  lb_idx <- which(fixed@par_class == "lapse_bias")
  lapse_bias <- if (length(lb_idx) > 0L) {
    val <- fixed@value[[lb_idx[length(lb_idx)]]]
    .assert_true(
      length(val) == n_categories,
      msg = sprintf(
        "Fixed lapse_bias must be a simplex of length %d.",
        n_categories
      )
    )
    list(
      known = 1L,
      data = matrix(as.numeric(val), nrow = 1, ncol = n_categories)
    )
  } else {
    list(
      known = 0L,
      data = matrix(0, nrow = 0, ncol = n_categories)
    )
  }

  # 4. m_0
  m_idx <- which(fixed@par_class == "m_0")
  m_known_list <- list()
  m_known_indices <- integer()

  for (idx in m_idx) {
    val <- fixed@value[[idx]]
    cat_spec <- fixed@category[idx]
    spc <- fixed@cue_space[idx]

    if (nzchar(cat_spec)) {
      .assert_true(
        cat_spec %in% categories,
        msg = sprintf("Category '%s' not found in model categories.", cat_spec)
      )
      c_i <- match(cat_spec, categories)
      val_vec <- as.numeric(val)
      .assert_true(
        length(val_vec) == n_cues,
        msg = sprintf("Fixed m_0 must have length %d.", n_cues)
      )
      if (spc == "original") {
        if (length(transform@shift) == 1L) {
          val_trans <- as.numeric(
            transform@SCALE * (val_vec - transform@shift)
          )
        } else {
          val_trans <- as.numeric(
            transform@SCALE %*% (val_vec - transform@shift)
          )
        }
      } else {
        val_trans <- val_vec
      }
      m_known_indices <- c(m_known_indices, c_i)
      m_known_list[[length(m_known_list) + 1L]] <- val_trans
    } else {
      if (is.list(val)) {
        for (c_nm in names(val)) {
          .assert_true(
            c_nm %in% categories,
            msg = sprintf("Category '%s' not found.", c_nm)
          )
          c_i <- match(c_nm, categories)
          val_vec <- as.numeric(val[[c_nm]])
          .assert_true(
            length(val_vec) == n_cues,
            msg = sprintf("Fixed mu_0 must have length %d.", n_cues)
          )
          if (spc == "original") {
            if (length(transform@shift) == 1L) {
              val_trans <- as.numeric(
                transform@SCALE * (val_vec - transform@shift)
              )
            } else {
              val_trans <- as.numeric(
                transform@SCALE %*% (val_vec - transform@shift)
              )
            }
          } else {
            val_trans <- val_vec
          }
          m_known_indices <- c(m_known_indices, c_i)
          m_known_list[[length(m_known_list) + 1L]] <- val_trans
        }
      } else if (n_categories == 1L && (is.vector(val) || is.numeric(val))) {
        val_vec <- as.numeric(val)
        .assert_true(
          length(val_vec) == n_cues,
          msg = sprintf("Fixed m_0 must have length %d.", n_cues)
        )
        if (spc == "original") {
          if (length(transform@shift) == 1L) {
            val_trans <- as.numeric(
              transform@SCALE * (val_vec - transform@shift)
            )
          } else {
            val_trans <- as.numeric(
              transform@SCALE %*% (val_vec - transform@shift)
            )
          }
        } else {
          val_trans <- val_vec
        }
        m_known_indices <- c(m_known_indices, 1L)
        m_known_list[[length(m_known_list) + 1L]] <- val_trans
      }
    }
  }

  if (length(m_known_indices) > 0L) {
    unique_inds <- unique(m_known_indices)
    final_m_list <- list()
    for (u_i in unique_inds) {
      pos <- which(m_known_indices == u_i)
      final_m_list[[length(final_m_list) + 1L]] <-
        m_known_list[[pos[length(pos)]]]
    }
    M_m_0_known <- length(unique_inds)
    idx_m_0_known <- as.array(as.integer(unique_inds))
    m_0_known_data <- matrix(
      unlist(final_m_list),
      nrow = M_m_0_known,
      ncol = n_cues,
      byrow = TRUE
    )
  } else {
    M_m_0_known <- 0L
    idx_m_0_known <- as.array(integer(0))
    m_0_known_data <- matrix(0, nrow = 0, ncol = n_cues)
  }

  m_0 <- list(
    M_known = M_m_0_known,
    idx_known = idx_m_0_known,
    data = m_0_known_data
  )

  # 5. S_0 / expected_Sigma_0
  s_idx <- which(fixed@par_class %in% c("S_0", "expected_Sigma_0"))
  s_known_list <- list()
  s_known_indices <- integer()
  s_known_modes <- integer()

  for (idx in s_idx) {
    cls <- fixed@par_class[idx]
    mode_val <- if (cls == "S_0") 1L else 2L
    val <- fixed@value[[idx]]
    cat_spec <- fixed@category[idx]
    spc <- fixed@cue_space[idx]

    if (nzchar(cat_spec)) {
      .assert_true(
        cat_spec %in% categories,
        msg = sprintf("Category '%s' not found.", cat_spec)
      )
      c_i <- match(cat_spec, categories)
      val_mat <- as.matrix(val)
      .assert_true(
        nrow(val_mat) == n_cues && ncol(val_mat) == n_cues,
        msg = sprintf("Fixed %s must be %d x %d.", cls, n_cues, n_cues)
      )
      if (spc == "original") {
        if (length(transform@shift) == 1L) {
          sc <- as.numeric(transform@SCALE)
          val_trans <- (sc^2) * val_mat
        } else {
          val_trans <- transform@SCALE %*% val_mat %*% t(transform@SCALE)
        }
      } else {
        val_trans <- val_mat
      }
      s_known_indices <- c(s_known_indices, c_i)
      s_known_modes <- c(s_known_modes, mode_val)
      s_known_list[[length(s_known_list) + 1L]] <- val_trans
    } else {
      if (is.list(val)) {
        for (c_nm in names(val)) {
          .assert_true(
            c_nm %in% categories,
            msg = sprintf("Category '%s' not found.", c_nm)
          )
          c_i <- match(c_nm, categories)
          val_mat <- as.matrix(val[[c_nm]])
          .assert_true(
            nrow(val_mat) == n_cues && ncol(val_mat) == n_cues,
            msg = sprintf("Fixed %s must be %d x %d.", cls, n_cues, n_cues)
          )
          if (spc == "original") {
            if (length(transform@shift) == 1L) {
              sc <- as.numeric(transform@SCALE)
              val_trans <- (sc^2) * val_mat
            } else {
              val_trans <- transform@SCALE %*% val_mat %*%
                t(transform@SCALE)
            }
          } else {
            val_trans <- val_mat
          }
          s_known_indices <- c(s_known_indices, c_i)
          s_known_modes <- c(s_known_modes, mode_val)
          s_known_list[[length(s_known_list) + 1L]] <- val_trans
        }
      } else if (n_categories == 1L && (is.matrix(val) || is.array(val))) {
        val_mat <- as.matrix(val)
        .assert_true(
          nrow(val_mat) == n_cues && ncol(val_mat) == n_cues,
          msg = sprintf("Fixed %s must be %d x %d.", cls, n_cues, n_cues)
        )
        if (spc == "original") {
          if (length(transform@shift) == 1L) {
            sc <- as.numeric(transform@SCALE)
            val_trans <- (sc^2) * val_mat
          } else {
            val_trans <- transform@SCALE %*% val_mat %*%
              t(transform@SCALE)
          }
        } else {
          val_trans <- val_mat
        }
        s_known_indices <- c(s_known_indices, 1L)
        s_known_modes <- c(s_known_modes, mode_val)
        s_known_list[[length(s_known_list) + 1L]] <- val_trans
      }
    }
  }

  if (length(s_known_indices) > 0L) {
    unique_inds <- unique(s_known_indices)
    final_s_list <- list()
    final_modes <- integer(length(unique_inds))
    for (u_pos in seq_along(unique_inds)) {
      u_i <- unique_inds[u_pos]
      pos <- which(s_known_indices == u_i)
      last_pos <- pos[length(pos)]
      final_s_list[[u_pos]] <- s_known_list[[last_pos]]
      final_modes[u_pos] <- s_known_modes[last_pos]
    }
    M_S_0_known <- length(unique_inds)
    idx_S_0_known <- as.array(as.integer(unique_inds))
    mode_S_0_known <- as.array(as.integer(final_modes))
    S_0_known_data <- array(0, dim = c(M_S_0_known, n_cues, n_cues))
    for (k in seq_len(M_S_0_known)) {
      S_0_known_data[k, , ] <- final_s_list[[k]]
    }
  } else {
    M_S_0_known <- 0L
    idx_S_0_known <- as.array(integer(0))
    mode_S_0_known <- as.array(integer(0))
    S_0_known_data <- array(0, dim = c(0, n_cues, n_cues))
  }

  S_0 <- list(
    M_known = M_S_0_known,
    idx_known = idx_S_0_known,
    mode_known = mode_S_0_known,
    data = S_0_known_data
  )

  # 6. kappa_0
  k_idx <- which(fixed@par_class == "kappa_0")
  k_known_list <- list()
  k_known_indices <- integer()
  is_cat_k <- FALSE

  for (idx in k_idx) {
    val <- as.numeric(fixed@value[[idx]])
    cat_spec <- fixed@category[idx]
    if (nzchar(cat_spec)) {
      is_cat_k <- TRUE
      .assert_true(
        cat_spec %in% categories,
        msg = sprintf("Category '%s' not found.", cat_spec)
      )
      c_i <- match(cat_spec, categories)
      k_known_indices <- c(k_known_indices, c_i)
      k_known_list[[length(k_known_list) + 1L]] <- val[1]
    } else {
      if (length(val) == n_categories && !is.null(names(val))) {
        is_cat_k <- TRUE
        for (c_nm in names(val)) {
          c_i <- match(c_nm, categories)
          k_known_indices <- c(k_known_indices, c_i)
          k_known_list[[length(k_known_list) + 1L]] <- as.numeric(val[c_nm])
        }
      } else {
        k_known_indices <- c(k_known_indices, 1L)
        k_known_list[[length(k_known_list) + 1L]] <- val[1]
      }
    }
  }

  if (length(k_known_indices) > 0L) {
    unique_inds <- unique(k_known_indices)
    final_k <- numeric(length(unique_inds))
    for (u_pos in seq_along(unique_inds)) {
      pos <- which(k_known_indices == unique_inds[u_pos])
      final_k[u_pos] <- k_known_list[[pos[length(pos)]]]
    }
    M_kappa_known <- length(unique_inds)
    idx_kappa_known <- as.array(as.integer(unique_inds))
    kappa_data <- as.array(as.numeric(final_k))
  } else {
    M_kappa_known <- 0L
    idx_kappa_known <- as.array(integer(0))
    kappa_data <- as.array(numeric(0))
  }

  kappa_0 <- list(
    is_category_specific = as.integer(is_cat_k),
    M_known = M_kappa_known,
    idx_known = idx_kappa_known,
    data = kappa_data
  )

  # 7. nu_0
  nu_idx <- which(fixed@par_class == "nu_0")
  nu_known_list <- list()
  nu_known_indices <- integer()
  is_cat_nu <- FALSE

  for (idx in nu_idx) {
    val <- as.numeric(fixed@value[[idx]])
    cat_spec <- fixed@category[idx]
    if (nzchar(cat_spec)) {
      is_cat_nu <- TRUE
      .assert_true(
        cat_spec %in% categories,
        msg = sprintf("Category '%s' not found.", cat_spec)
      )
      c_i <- match(cat_spec, categories)
      nu_known_indices <- c(nu_known_indices, c_i)
      nu_known_list[[length(nu_known_list) + 1L]] <- val[1]
    } else {
      if (length(val) == n_categories && !is.null(names(val))) {
        is_cat_nu <- TRUE
        for (c_nm in names(val)) {
          c_i <- match(c_nm, categories)
          nu_known_indices <- c(nu_known_indices, c_i)
          nu_known_list[[length(nu_known_list) + 1L]] <- as.numeric(val[c_nm])
        }
      } else {
        nu_known_indices <- c(nu_known_indices, 1L)
        nu_known_list[[length(nu_known_list) + 1L]] <- val[1]
      }
    }
  }

  if (length(nu_known_indices) > 0L) {
    unique_inds <- unique(nu_known_indices)
    final_nu <- numeric(length(unique_inds))
    for (u_pos in seq_along(unique_inds)) {
      pos <- which(nu_known_indices == unique_inds[u_pos])
      final_nu[u_pos] <- nu_known_list[[pos[length(pos)]]]
    }
    M_nu_known <- length(unique_inds)
    idx_nu_known <- as.array(as.integer(unique_inds))
    nu_data <- as.array(as.numeric(final_nu))
  } else {
    M_nu_known <- 0L
    idx_nu_known <- as.array(integer(0))
    nu_data <- as.array(numeric(0))
  }

  nu_0 <- list(
    is_category_specific = as.integer(is_cat_nu),
    M_known = M_nu_known,
    idx_known = idx_nu_known,
    data = nu_data
  )

  # 8. weights
  w_idx <- which(fixed@par_class == "weights")
  weights <- if (length(w_idx) > 0L) {
    val <- fixed@value[[w_idx[length(w_idx)]]]
    list(
      known = 1L,
      data = val
    )
  } else {
    list(
      known = 0L,
      data = NULL
    )
  }

  list(
    category_prior = category_prior,
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    m_0 = m_0,
    S_0 = S_0,
    kappa_0 = kappa_0,
    nu_0 = nu_0,
    weights = weights
  )
}
