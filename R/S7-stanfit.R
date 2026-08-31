#' @include asserts.R
#' @include S7-core-classes.R
#' @include S7-transform-information.R
#' @include S7-staninput.R
#' @include S7-stanfit-input.R
#' @importFrom S7 new_class new_object new_generic method class_name props
NULL

#' An S7 base class for Stan fit objects
#'
#' @name MVBU_Stanfit
#' @rdname MVBU_Stanfit
#' @title MVBeliefUpdatr Stanfit Classes and Constructors
#' @docType class
#'
#' @slot data A \code{data.frame} containing the data used to fit the model.
#' @slot staninput A Stan input object containing the data handed to rstan through
#'   \code{\link{make_staninput}}. The staninput object contains at
#'   least two components: \code{transformed} and \code{untransformed}.
#' @slot stanvars A \code{\link{stanvars}} object or \code{NULL}.
#' @slot backend The name of the backend used to fit the model (character).
#' @slot save_pars Optional storage for saved parameter names.
#' @slot stan_args Named list of additional control arguments that were passed
#'   to the Stan backend directly. NOT YET USED
#' @slot stanfit An object of class
#'   \code{\link[rstan:stanfit-class]{stanfit}} containing posterior draws.
#' @slot basis An object that contains a small subset of the Stan data
#'   created at fitting time, needed to process new data. NOT YET USED
#' @slot transform_information An object of type \code{\link{MVBU_TransformInformation}}.
#' @slot criteria An empty \code{list} for adding model fit criteria
#'   after estimation of the model. NOT YET USED
#' @slot file Optional name of a file in which the model object was stored in
#'   or loaded from.
#' @slot version The versions of \pkg{MVBeliefUpdatr} and \pkg{rstan} with
#'   which the model was fitted.
#' @slot metadata List containing auxiliary information including
#'   \code{label_information}.
#'
#' @rawNamespace if (getRversion() < "4.3.0") importFrom(S7, "@")
#' @importFrom purrr map_lgl map_chr
#' @importFrom utils packageVersion
#' @export
MVBU_Stanfit <- S7::new_class(
  "MVBU_Stanfit",
  package = NULL,
  parent = MVBU_Object,
  properties = list(
    data = S7::new_S3_class("data.frame"),
    staninput = S7::class_any,
    stanvars = S7::class_any,
    backend = S7::class_character,
    save_pars = S7::class_any,
    stan_args = S7::class_list,
    stanfit = S7::class_any,
    basis = S7::class_any,
    transform_information = MVBU_TransformInformation,
    criteria = S7::class_list,
    file = S7::class_character,
    version = S7::class_any,
    metadata = S7::class_list
  ),
  constructor = function(
    data = data.frame(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
    basis = NULL,
    transform_information = NULL,
    criteria = list(),
    file = NULL,
    version = NULL,
    metadata = list()
  ) {
    if (is.null(transform_information)) {
      transform_information <- MVBU_TransformInformation()
    }

    if (is.null(staninput)) {
      staninput <- MVBU_Staninput(values = list())
    }

    if (is.null(version)) {
      version <- get_current_versions()
    }
    if (is.null(file)) {
      file <- character(0)
    }

    if (!is.list(metadata)) {
      metadata <- list()
    }

    label_info <- if (!is.null(metadata$label_information) &&
      is.list(metadata$label_information)) {
      metadata$label_information
    } else {
      list()
    }

    metadata$label_information <- list(
      cue = if (!is.null(label_info$cue)) as.character(label_info$cue) else character(0),
      category = if (!is.null(label_info$category)) as.character(label_info$category) else character(0),
      group = if (!is.null(label_info$group)) as.character(label_info$group) else character(0)
    )

    S7::new_object(
      MVBU_Object(),
      data = as.data.frame(data),
      staninput = staninput,
      stanvars = stanvars,
      backend = as.character(backend),
      save_pars = save_pars,
      stan_args = as.list(stan_args),
      stanfit = stanfit,
      basis = basis,
      transform_information = transform_information,
      criteria = as.list(criteria),
      file = as.character(file),
      version = version,
      metadata = metadata
    )
  },
  validator = function(self) {
    if (!is.data.frame(self@data)) {
      return("`data` must be a data.frame")
    }
    if (!is.null(self@staninput) && !S7::S7_inherits(self@staninput, MVBU_Staninput)) {
      return("`staninput` must be NULL or an S7 object")
    }
    if (!is.character(self@backend)) {
      return("`backend` must be a character")
    }
    if (!is.list(self@stan_args)) {
      return("`stan_args` must be a list")
    }
    if (!is.null(self@stanfit)) {
      if (!inherits(self@stanfit, "stanfit")) {
        return("`stanfit` must inherit from stanfit when provided")
      }
      if (!(self@stanfit@model_name %in% names(MVBeliefUpdatr:::stanmodels))) {
        return(
          paste0(
            "`stanfit` model_name is not recognized. `stanfit` has to be ",
            "created by one of the accepted stanmodels:\n\t",
            paste(names(MVBeliefUpdatr:::stanmodels), collapse = "\n\t"),
            "\n(you can get the name of your model from ",
            "your_stanfit@model_name)."
          )
        )
      }
    }
    if (!is.list(self@criteria)) {
      return("`criteria` must be a list")
    }
    if (!is.null(self@file) && !is.character(self@file)) {
      return("`file` must be NULL or a character")
    }
    if (!is.null(self@transform_information) &&
        !S7::S7_inherits(self@transform_information,
                         MVBU_TransformInformation)) {
      return("`transform_information` must inherit from MVBU_TransformInformation")
    }
    if (!is.list(self@metadata)) {
      return("`metadata` must be a list")
    }
    NULL
  }
)

#' @rdname MVBU_Stanfit
#' @export
IdealAdaptorStanfit <- S7::new_class(
  "IdealAdaptorStanfit",
  package = NULL,
  parent = MVBU_Stanfit,
  constructor = function(
    data = data.frame(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
    basis = NULL,
    transform_information = NULL,
    criteria = list(),
    file = NULL,
    version = NULL,
    metadata = list()
  ) {
    ti <- if (is.null(transform_information)) {
      MVBU_TransformInformation()
    } else {
      transform_information
    }
    S7::new_object(
      MVBU_Stanfit(
        data = data,
        staninput = staninput,
        stanvars = stanvars,
        backend = as.character(backend),
        save_pars = save_pars,
        stan_args = as.list(stan_args),
        stanfit = stanfit,
        basis = basis,
        transform_information = ti,
        criteria = as.list(criteria),
        file = as.character(file),
        version = if (is.null(version)) get_current_versions() else version,
        metadata = as.list(metadata)
      )
    )
  }
)

NIX_IdealAdaptorStanfit <- S7::new_class(
  "NIX_IdealAdaptorStanfit",
  package = NULL,
  parent = IdealAdaptorStanfit,
  constructor = function(
    data = data.frame(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
    basis = NULL,
    transform_information = NULL,
    criteria = list(),
    file = NULL,
    version = NULL,
    metadata = list()
  ) {
    ti <- if (is.null(transform_information)) {
      MVBU_TransformInformation()
    } else {
      transform_information
    }
    S7::new_object(
      IdealAdaptorStanfit(
        data = data,
        staninput = staninput,
        stanvars = stanvars,
        backend = as.character(backend),
        save_pars = save_pars,
        stan_args = as.list(stan_args),
        stanfit = stanfit,
        basis = basis,
        transform_information = ti,
        criteria = as.list(criteria),
        file = as.character(file),
        version = if (is.null(version)) get_current_versions() else version,
        metadata = as.list(metadata)
      )
    )
  }
)

MNIX_IdealAdaptorStanfit <- S7::new_class(
  "MNIX_IdealAdaptorStanfit",
  package = NULL,
  parent = IdealAdaptorStanfit,
  constructor = function(
    data = data.frame(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
    basis = NULL,
    transform_information = NULL,
    criteria = list(),
    file = NULL,
    version = NULL,
    metadata = list()
  ) {
    ti <- if (is.null(transform_information)) {
      MVBU_TransformInformation()
    } else {
      transform_information
    }
    S7::new_object(
      IdealAdaptorStanfit(
        data = data,
        staninput = staninput,
        stanvars = stanvars,
        backend = as.character(backend),
        save_pars = save_pars,
        stan_args = as.list(stan_args),
        stanfit = stanfit,
        basis = basis,
        transform_information = ti,
        criteria = as.list(criteria),
        file = as.character(file),
        version = if (is.null(version)) get_current_versions() else version,
        metadata = as.list(metadata)
      )
    )
  }
)

NIW_IdealAdaptorStanfit <- S7::new_class(
  "NIW_IdealAdaptorStanfit",
  package = NULL,
  parent = IdealAdaptorStanfit,
  constructor = function(
    data = data.frame(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
    basis = NULL,
    transform_information = NULL,
    criteria = list(),
    file = NULL,
    version = NULL,
    metadata = list()
  ) {
    ti <- if (is.null(transform_information)) {
      MVBU_TransformInformation()
    } else {
      transform_information
    }
    S7::new_object(
      IdealAdaptorStanfit(
        data = data,
        staninput = staninput,
        stanvars = stanvars,
        backend = as.character(backend),
        save_pars = save_pars,
        stan_args = as.list(stan_args),
        stanfit = stanfit,
        basis = basis,
        transform_information = ti,
        criteria = as.list(criteria),
        file = as.character(file),
        version = if (is.null(version)) get_current_versions() else version,
        metadata = as.list(metadata)
      )
    )
  }
)

#' @export
ideal_adaptor_stanfit <- function(
  data = data.frame(),
  staninput = NULL,
  stanvars = NULL,
  backend = "rstan",
  save_pars = NULL,
  stan_args = list(),
  stanfit = NULL,
  basis = NULL,
  transform_information = NULL,
  criteria = list(),
  file = NULL,
  version = NULL,
  metadata = list()
) {
  constructor <- get_ideal_adaptor_stanfit_constructor(staninput)

  constructor(
    data = data,
    staninput = staninput,
    stanvars = stanvars,
    backend = backend,
    save_pars = save_pars,
    stan_args = stan_args,
    stanfit = stanfit,
    basis = basis,
    transform_information = transform_information,
    criteria = criteria,
    file = file,
    version = version,
    metadata = metadata
  )
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
  inherits(x, "ideal_adaptor_stanfit") ||
    S7::S7_inherits(x, IdealAdaptorStanfit)
}

# -------------------------
# Utilities ported from original file
# -------------------------

#' Build parameter names for a rectangular grid of indices
#'
#' @param prefix Character prefix for the parameter names.
#' @param ... Arguments passed to \code{expand.grid}.
#' @return A character vector of parameter names.
#' @keywords internal
#' @noRd
.make_parnames <- function(prefix, ...) {
  combinations <- expand.grid(..., stringsAsFactors = FALSE)
  paste0(prefix, "[", apply(combinations, 1, paste0, collapse = ","), "]")
}

#' Rename fitted parameter names to their transformed equivalents
#'
#' @param x An ideal adaptor Stanfit object.
#' @param include_original_pars Whether to preserve original parameter names.
#' @return The updated Stanfit object.
#' @keywords internal
#' @noRd
.rename_pars <- function(x, include_original_pars = FALSE) {
  assert_IdealAdaptorStanfit(x)
  stanfit <- get_stanfit(x)

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

  for (i in seq_len(chains)) names(stanfit@sim$samples[[i]]) <- vapply(names(stanfit@sim$samples[[i]]), .rename, FUN.VALUE = character(1))

  set_stanfit(x, stanfit)
}

# --- contains_draws temporary dispatch shim ---
# NOTE: legacy S7/S4 mixed generic syntax is temporarily replaced to keep
# package loadable during S7 migration and roxygen generation.
#' Check whether a Stanfit object contains posterior draws
#'
#' @param x Object to inspect.
#' @param ... Additional arguments (currently unused).
#' @return A logical scalar.
#' @keywords internal
#' @noRd
.contains_draws <- function(x, ...) {
  if (S7::S7_inherits(x, IdealAdaptorStanfit)) {
    return(.contains_draws(x@stanfit))
  }

  if (inherits(x, "stanfit")) {
    return(length(x@sim) > 0)
  }

  FALSE
}

# --- file helpers ---
#' Normalize a fit file path to include the .rds suffix
#'
#' @param file File path or base name.
#' @return A character scalar with the normalized file path.
#' @keywords internal
#' @noRd
.check_stanfit_file <- function(file) {
  file <- .as_one_character(file)
  file_ending <- tolower(.get_matches("\\.[^\\.]+$", file))
  if (!isTRUE(file_ending == ".rds")) {
    file <- paste0(file, ".rds")
  }
  file
}

#' Return the supported refit-policy options for cached Stanfits
#'
#' @return A character vector of supported options.
#' @keywords internal
#' @noRd
.file_refit_options <- function() {
  c("never", "always", "on_change")
}

#' Check if cached \code{ideal_adaptor_stanfit} can be used
#'
#' Checks whether a given cached fit can be used without refitting when
#' \code{file_refit = "on_change"} is used.
#'
#' @param x Old \code{ideal_adaptor_stanfit} object (e.g., loaded from file).
#' @param current_version Current version of relevant packages. (default: will be automatically
#'  obtained from current packages).
#' @param data New data to check consistency of factor level names. (default: \code{NULL}))
#' @param staninput New Stan data (result of a call to \code{\link{make_staninput}}).
#'   Pass \code{NULL} to avoid this data check. (default: \code{NULL}))
#' @param silent Logical. If \code{TRUE}, no messages will be given. (default: \code{FALSE}))
#' @param verbose Logical. If \code{TRUE} detailed report of the differences
#'   is printed to the console. (default: \code{FALSE}))
#' @return A boolean indicating whether a refit is needed.
#'
#' @details
#' fit differs from the given data and code.
#'
#' @keywords internal
#' @noRd
.stanfit_needs_refit <- function(
  x,
  current_version = get_current_versions(),
  data = NULL, staninput = NULL,
  silent = FALSE, verbose = FALSE
) {
  assert_IdealAdaptorStanfit(x)
  silent <- .as_one_logical(silent)
  verbose <- .as_one_logical(verbose)

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
    assert_IdealAdaptorStaninput(staninput)
    cached_staninput <- get_staninput(x)
    staninput <- staninput@values
  }
  if (!is.null(data)) {
    .assert_data_frame_like(data)
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
      if (.is_like_factor(cached_data[[var]])) {
        cached_levels <- levels(factor(cached_data[[var]]))
        new_levels <- levels(factor(data[[var]]))
        if (!.is_equal(cached_levels, new_levels)) {
          if (!silent) {
            factor_level_message <- TRUE
            if (verbose) {
              cat(paste0(
                "Names of factor levels in data have changed for variable '", var, "' ",
                "with cached levels (", paste(as.character(cached_levels), collapse = ", "), ") ",
                "but new levels (", paste(as.character(new_levels), collapse = ", "), ").\n"
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
#' Read a cached ideal adaptor Stanfit from disk
#'
#' @param file File path to the cached fit.
#' @return A cached fit object or \code{NULL} if none is available.
#' @keywords internal
#' @noRd
.read_ideal_adaptor_stanfit <- function(file) {
  file <- .check_stanfit_file(file)
  if (!file.exists(file)) {
    return(NULL)
  }
  x <- suppressWarnings(try(readRDS(file), silent = TRUE))
  if (.is_try_error(x)) {
    return(NULL)
  }

  .assert_true(
    S7::S7_inherits(x, IdealAdaptorStanfit),
    msg = "Object loaded from 'file' is not an IdealAdaptorStanfit object. This indicates that it was fit with an outdated version of MVBeliefUpdatr and needs to be refit."
  )

  label_info <- try(x@metadata$label_information, silent = TRUE)
  .assert_true(
    !is.null(label_info) && is.list(label_info) &&
      !is.null(label_info$cue) && !is.null(label_info$category) && !is.null(label_info$group),
    msg = paste0(
      "Object loaded from 'file' is missing label_information in metadata. ",
      "Please refit the model."
    )
  )

  x@file <- file
  x
}

#' Write a fitted ideal adaptor Stanfit object to disk
#'
#' @param x The fitted object to save.
#' @param file File path for the saved object.
#' @param compress Compression level passed to \code{saveRDS}.
#' @return The saved object, invisibly.
#' @keywords internal
#' @noRd
.write_ideal_adaptor_stanfit <- function(x, file, compress = TRUE) {
  assert_IdealAdaptorStanfit(x)
  file <- .check_stanfit_file(file)
  x@file <- file
  saveRDS(x, file = file, compress = compress)
  invisible(x)
}
