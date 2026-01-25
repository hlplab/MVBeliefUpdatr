#' @importFrom S7 new_class new_object new_generic method class_name props
NULL

#' S7 class for processed Stan input used by ideal_adaptor_stanfit
#'
#' @name ideal_adaptor_staninput-class
#' @docType class
#' @property transformed A list containing transformed stan input.
#' @property untransformed A list containing untransformed stan input.
#' @export
ideal_adaptor_staninput <- S7::new_class(
  "ideal_adaptor_staninput",
  properties = list(
    transformed = "list",
    untransformed = "list"
  ),
  validator = function(self) {
    stopifnot(is.list(self@transformed))
    stopifnot(is.list(self@untransformed))
    TRUE
  }
)

# Constructor helper for staninput
#' @export
ideal_adaptor_staninput <- function(transformed = list(), untransformed = list()) {
  S7::new_object(
    "ideal_adaptor_staninput",
    transformed = transformed,
    untransformed = untransformed
  )
}

#' S7 class for transformation information for Stan models
#'
#' @name transform_information-class
#' @docType class
#'
#' @property transform.parameters A list of sufficient parameters for the transform.
#' @property transform.function Function to transform parameters.
#' @property untransform.function Function to untransform parameters.
#' @export
transform_information <- S7::new_class(
  "transform_information",
  properties = list(
    transform.parameters = "list",
    transform.function = "function",
    untransform.function = "function"
  ),
  validator = function(self) {
    stopifnot(is.list(self@transform.parameters))
    stopifnot(is.function(self@transform.function))
    stopifnot(is.function(self@untransform.function))
    TRUE
  }
)

# Constructor helper
#' @export
transform_information <- function(transform.parameters, transform.function, untransform.function) {
  S7::new_object(
    "transform_information",
    transform.parameters = transform.parameters,
    transform.function = transform.function,
    untransform.function = untransform.function
  )
}


#' An S7 class for ideal_adaptor stanfit objects that use one of the ideal_adaptor Stan programs.
#'
#' @name ideal_adaptor_stanfit-class
#' @aliases NIW_ideal_adaptor_stanfit
#' @docType class
#'
#' @details
#' See \code{methods(class = "ideal_adaptor_stanfit")} for an overview of available methods.
#'
#' @property data A \code{data.frame} containing the data used to fit the model.
#' @property staninput A named object of class \code{ideal_adaptor_staninput} containing the
#'   data handed to rstan through \code{\link{make_staninput}}. The staninput object
#'   contains at least two components: \code{transformed} and \code{untransformed}.
#' @property stanvars A \code{\link{stanvars}} object or \code{NULL}.
#' @property backend The name of the backend used to fit the model (character).
#' @property save_pars Optional storage for saved parameter names.
#' @property stan_args Named list of additional control arguments that were passed
#'   to the Stan backend directly. NOT YET USED
#' @property stanfit An object of class \code{\link[rstan:stanfit-class]{stanfit}}
#'   containing the posterior draws.
#' @property basis An object that contains a small subset of the Stan data
#'   created at fitting time, which is needed to process new data correctly. NOT YET USED
#' @property transform_information An object of type \code{\link{transform_information}}.
#' @property criteria An empty \code{list} for adding model fit criteria
#'   after estimation of the model. NOT YET USED
#' @property file Optional name of a file in which the model object was stored in
#'   or loaded from.
#' @property version The versions of \pkg{MVBeliefUpdatr} and \pkg{rstan} with
#'   which the model was fitted.
#' @property labels List of labels
#'
#' @rawNamespace if (getRversion() < "4.3.0") importFrom(S7, "@")
#' @importFrom purrr map_lgl map_chr
#' @importFrom utils packageVersion
#' @export
ideal_adaptor_stanfit <- S7::new_class(
  "ideal_adaptor_stanfit",
  properties = list(
    data = "data.frame",
    staninput = "ideal_adaptor_staninput",
    stanvars = "ANY",
    backend = "character",
    save_pars = "ANY",
    stan_args = "list",
    stanfit = "stanfit",
    basis = "ANY",
    transform_information = "transform_information",
    criteria = "list",
    file = "character",
    version = "ANY",
    labels = "list"
  ),
  validator = function(self) {
    stopifnot(is.data.frame(self@data))
    stopifnot(S7::is_object(self@staninput) && S7::class_name(self@staninput) == "ideal_adaptor_staninput")
    stopifnot(is.character(self@backend))
    stopifnot(is.list(self@stan_args))
    if (!is.null(stanfit)) {
      stopifnot(inherits(stanfit, "stanfit"))
      stopifnot(stanfit@model_name %in% names(MVBeliefUpdatr:::stanmodels))
    }
    stopifnot(is.list(self@criteria))
    stopifnot(is.character(self@file))
    stopifnot(S7::is_object(self@transform_information) && S7::class_name(self@transform_information) == "transform_information")
    TRUE
  }
)

# Constructor (modernized S7)
#' Constructor for ideal_adaptor_stanfit (S7)
#'
#' @param data data.frame
#' @param staninput ideal_adaptor_staninput or list; if list, will be coerced
#' @param stanvars backend stanvars or NULL
#' @param backend character backend name
#' @param save_pars optional
#' @param stan_args list of stan args
#' @param stanfit rstan stanfit or NULL
#' @param basis optional
#' @param transform_information transform_information
#' @param criteria list
#' @param file character path or NULL/"".
#' @param version package version info (from get_current_versions())
#' @param labels list
#' @export
ideal_adaptor_stanfit <- function(
    data = data.frame(),
    staninput = ideal_adaptor_staninput(),
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
    basis = NULL,
    transform_information = NULL,
    criteria = list(),
    file = NULL,
    version = get_current_versions(),
    labels = list()
) {
  # Coerce staninput if it's NULL or a list
  if (is.list(staninput) && !S7::is_object(staninput)) {
      staninput <- ideal_adaptor_staninput(
        transformed = staninput@transformed %||% list(),
        untransformed = staninput@untransformed %||% list()
      )
  }

  # Coerce transform_information if it's NULL or a list
  if (!S7::is_object(transform_information)) {
    transform_information <- transform_information(
      transform.parameters = transform_information$transform.parameters %||% list(),
      transform.function = transform_information$transform.function %||% (function(x) x),
      untransform.function = transform_information$untransform.function %||% (function(x) x)
    )
  }

  S7::new_object(
    "ideal_adaptor_stanfit",
    data = as.data.frame(data),
    staninput = staninput,
    stanvars = stanvars,
    backend = as.character(backend),
    save_pars = save_pars,
    stan_args = stan_args,
    stanfit = stanfit,
    basis = basis,
    transform_information = transform_information,
    criteria = criteria,
    file = as.character(file),
    version = version,
    labels = labels
  )
}

# -------------------------
# Utilities ported from original file
# -------------------------

make_parnames <- function(prefix, ...) {
  combinations <- expand.grid(..., stringsAsFactors = FALSE)
  paste0(prefix, "[", apply(combinations, 1, paste0, collapse = ","), "]")
}

rename_pars <- function(x, include_original_pars = FALSE) {
  stopifnot(S7::is_object(x) && S7::class_name(x) == "ideal_adaptor_stanfit")
  stanfit <- x@stanfit

  chains <- length(stanfit@sim$samples)

  .rename <- function(parname) {
    parname <- gsub("(t_scale)\\[", "\\1_transformed\\[", parname)
    parname <- gsub("(m|S|tau)(_(0|n))\\[", "\\1\\2_transformed\\[", parname)
    parname <- gsub("(m|S|tau)(_(0|n))_original\\[", "\\1\\2\\[", parname)
    parname <- gsub("p_cat\\[", "p_category\\[", parname)
    parname
  }

  if (include_original_pars) stanfit@model_pars <- .rename(stanfit@model_pars)
  stanfit@sim$fnames_oi <- vapply(stanfit@sim$fnames_oi, .rename, FUN.VALUE = character(1))

  for (i in seq_len(chains)) names(stanfit@sim$samples[[i]])  <- vapply(names(stanfit@sim$samples[[i]]), .rename, FUN.VALUE = character(1))

  x@stanfit <- stanfit
  x
}

#' Is this an NIW ideal adaptor stanfit?
#'
#' Check whether \code{x} is of class \code{\link{ideal_adaptor_stanfit}}.
#'
#' @param x Object to be checked.
#' @param verbose Currently being ignored.
#' @return A logical.
#' @export
is.ideal_adaptor_stanfit <- function(x, verbose = FALSE) {
  inherits(x, "ideal_adaptor_stanfit") || (S7::is_object(x) && S7::class_name(x) == "ideal_adaptor_stanfit")
}

is.stanfit <- function(x) {
  inherits(x, "stanfit")
}


is.ideal_adaptor_staninput <- function(x) {
  # Accept either the S7 staninput object or an equivalent list structure
  if (!S7::is_object(x) && !S7::class_name(x) == "ideal_adaptor_staninput") return(FALSE)

  return(TRUE)
}

#' Is this a list of NIW ideal adaptor stanfit inputs?
#' @export
is.ideal_adaptor_stanfit_input <- function(x, verbose = FALSE) {
  if (!is.list(x)) {
    if (verbose) message("Object x is not a list.")
    return(FALSE)
  }
  if (!all(c("staninput", "data", "transform_information") %in% names(x))) {
    if (verbose) message("Object x is missing one of the required components: staninput, data, transform_information.")
    return(FALSE)
  }
  if (!is.ideal_adaptor_staninput(x@staninput)) {
    if (verbose) message("Component staninput in object x is not a ideal_adaptor_staninput.")
    return(FALSE)
  }
  if (!is.data.frame(x@data)) {
    if (verbose) message("Component data in object x is not a data.frame.")
    return(FALSE)
  }
  if (!is.transform_information(x@transform_information)) {
    if (verbose) message("Component transform_information in object x is not a list.")
    return(FALSE)
  }

  TRUE
}

# --- contains_draws generics and methods (S7) ---
contains_draws <- S7::new_generic("contains_draws", function(x, ...) standardGeneric("contains_draws"))

# Method for stanfit S4 objects (delegated to S4 internals)
S7::method(contains_draws, "stanfit", function(x, ...) {
  if (!(length(x@sim))) return(FALSE)
  TRUE
})

# Method for ideal_adaptor_stanfit
S7::method(contains_draws, "ideal_adaptor_stanfit", function(x, ...) {
  stanfit <- x@stanfit
  contains_draws(stanfit)
})

# --- file helpers ---
check_stanfit_file <- function(file) {
  file <- as_one_character(file)
  file_ending <- tolower(get_matches("\\.[^\\.]+$", file))
  if (!isTRUE(file_ending == ".rds")) {
    file <- paste0(file, ".rds")
  }
  file
}

file_refit_options <- function() {
  c("never", "always", "on_change")
}

#' Check if cached \code{ideal_adaptor_stanfit} can be used.
#' @keywords internal
stanfit_needs_refit <- function(
    x,
    current_version = get_current_versions(),
    data = NULL, staninput = NULL,
    silent = FALSE, verbose = FALSE
) {
  stopifnot(is.ideal_adaptor_stanfit(x))
  silent <- as_one_logical(silent)
  verbose <- as_one_logical(verbose)

  # Backward compatibility: if object is an old S4 object where 'version' slot didn't exist
  if (is.null(x@version)) {
    if (!silent) {
      message("The model ", deparse1(substitute(x)), " was fit with an old version of MVBeliefUpdater (< 0.0.1.0015).")
    }
    return(TRUE)
  }
  if (!isTRUE(all.equal(x@version, current_version))) {
    if (!silent) {
      message("Version of MVBeliefUpdatr or rstan has changed (current version is", paste(purrr::map_chr(current_version, ~ paste(.x, collapse = ", ")), collapse = "; "), ").")
      if (verbose) {
        print(x@version)
      }
    }
    return(TRUE)
  }

  if (!is.null(staninput)) {
    stopifnot(is.list(staninput))
    cached_staninput <- get_staninput(x, which = "both")
  }
  if (!is.null(data)) {
    stopifnot(is.data.frame(data))
    cached_data <- x@data
  }

  refit <- FALSE

  if (!is.null(staninput)) {
    staninput_equality <- all.equal(staninput, cached_staninput, check.attributes = FALSE, use.names = TRUE)
    if (!isTRUE(staninput_equality)) {
      if (!silent) {
        message("The processed input for Stan has changed.")
        if (verbose) print(staninput_equality)
      }
      refit <- TRUE
    }
  }
  if (!is.null(data)) {
    factor_level_message <- FALSE
    for (var in names(cached_data)) {
      if (is_like_factor(cached_data[[var]])) {
        cached_levels <- levels(factor(cached_data[[var]]))
        new_levels <- levels(factor(data[[var]]))
        if (!is_equal(cached_levels, new_levels)) {
          if (!silent) {
            factor_level_message <- TRUE
            if (verbose) {
              cat(paste0(
                "Names of factor levels in data have changed for variable '", var, "' ",
                "with cached levels (", collapse_comma(cached_levels), ") ",
                "but new levels (", collapse_comma(new_levels), ").\n"
              ))
            }
          }
          refit <- TRUE
          if (!verbose) break
        }
      }
    }
    if (factor_level_message) message("Names of factor levels in data have changed.")
  }

  if (!silent && refit) message("Model needs to be refit.")
  refit
}

# read/write functions
read_ideal_adaptor_stanfit <- function(file) {
  file <- check_stanfit_file(file)
  dir <- dirname(file)
  if (!dir.exists(dir)) {
    stop2(
      "The directory '", dir, "' does not exist. Please choose an ",
      "existing directory where the model can be saved after fitting."
    )
  }

  x <- suppressWarnings(try(readRDS(file), silent = TRUE))
  if (!is_try_error(x)) {
    if (!is.ideal_adaptor_stanfit(x)) {
      stop2("Object loaded from 'file' is not recognized as class 'ideal_adaptor_stanfit'. This usually indicates that you need to refit the model.")
    }
    x@file <- file
  } else {
    x <- NULL
  }
  x
}

write_ideal_adaptor_stanfit <- function(x, file, compress = TRUE) {
  stopifnot(is.ideal_adaptor_stanfit(x))
  file <- check_stanfit_file(file)
  x@file <- file
  saveRDS(x, file = file, compress = compress)
  invisible(x)
}
