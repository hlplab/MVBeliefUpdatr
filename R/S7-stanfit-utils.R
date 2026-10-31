#' @include asserts.R
#' @include internal-globals.R
#' @include S7-class.R
#' @include S7-stanfit.R
#' @include S7-staninput.R
#' @include S7-stanfit-input.R
#' @include S7-prior.R
#' @include S7-fixed.R
#' @importFrom S7 S7_inherits
#' @importFrom purrr map_chr
NULL

# -------------------------------------------------------------------------
# Path and Option Helpers
# -------------------------------------------------------------------------

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

# -------------------------------------------------------------------------
# Refit Detection
# -------------------------------------------------------------------------

#' Check if cached Stanfit can be used without refitting
#'
#' Checks whether a given cached fit can be used without refitting when
#' \code{file_refit = "on_change"} is used.
#'
#' @param x Old \code{\link{MVBU_Stanfit}} object (e.g., loaded from file).
#' @param current_version Current version of relevant packages (default: will
#'   be automatically obtained via \code{\link{get_current_versions}}).
#' @param data New data to check consistency of factor levels (default:
#'   \code{NULL}).
#' @param staninput New Stan data (result of a call to
#'   \code{\link{IdealAdaptorStaninput}}). Pass \code{NULL} to avoid this check
#' @param silent Logical. If \code{TRUE}, no messages will be given
#'   (default: \code{FALSE}).
#' @param verbose Logical. If \code{TRUE}, detailed report of differences
#'   is printed to the console (default: \code{FALSE}).
#' @return A boolean indicating whether a refit is needed.
#'
#' @keywords internal
#' @noRd
.stanfit_needs_refit <- function(
  x,
  current_version = get_current_versions(),
  data = NULL, staninput = NULL,
  silent = FALSE, verbose = FALSE
) {
  .assert_true(
    S7::S7_inherits(x, MVBU_Stanfit),
    msg = "x must inherit from MVBU_Stanfit."
  )
  silent <- .as_one_logical(silent)
  verbose <- .as_one_logical(verbose)

  core_pkgs <- intersect(c("MVBeliefUpdatr", "rstan"), names(current_version))
  if (!isTRUE(all.equal(x@version[core_pkgs], current_version[core_pkgs]))) {
    if (!silent) {
      message(
        "Version of MVBeliefUpdatr or rstan has changed (current version is ",
        paste(
          purrr::map_chr(
            current_version[core_pkgs],
            ~ paste(.x, collapse = ", ")
          ),
          collapse = "; "
        ),
        ")."
      )
      if (verbose) {
        print(x@version)
      }
    }
    return(TRUE)
  }

  has_orig_names <- !is.null(x@metadata$original_variable_names)
  required_cols <- c(
    "category", "response_category", "group", "group_unique", "Phase"
  )
  orig_data <- get_data(x, space = "original")
  has_std_cols <- !is.null(orig_data) &&
    all(required_cols %in% names(orig_data))
  if (!has_orig_names || !has_std_cols) {
    if (!silent) {
      message(
        "Cached model does not have standardized data columns or ",
        "original variable names metadata."
      )
    }
    return(TRUE)
  }

  if (!is.null(staninput)) {
    if (S7::S7_inherits(staninput, IdealAdaptorStaninput)) {
      staninput <- staninput@values
    }
    cached_staninput <- get_staninput(x)
    if (S7::S7_inherits(cached_staninput, IdealAdaptorStaninput)) {
      cached_staninput <- cached_staninput@values
    }
  }
  if (!is.null(data)) {
    .assert_data_frame_like(data)
    cached_data <- get_data(x, space = "original")
  }

  refit <- FALSE

  if (!is.null(staninput)) {
    si_new <- if (is.list(staninput) && !is.null(names(staninput))) {
      staninput[sort(names(staninput))]
    } else {
      staninput
    }
    si_cached <- if (
      is.list(cached_staninput) && !is.null(names(cached_staninput))
    ) {
      cached_staninput[sort(names(cached_staninput))]
    } else {
      cached_staninput
    }
    staninput_equality <- all.equal(
      si_new,
      si_cached,
      check.attributes = FALSE,
      use.names = TRUE
    )
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

# -------------------------------------------------------------------------
# Read and Write Functions
# -------------------------------------------------------------------------

#' Read and write cached Stanfit objects to disk
#'
#' Reads or saves a previously fitted \code{\link{MVBU_Stanfit}} object (or
#' subclass such as \code{\link{IdealAdaptorStanfit}}) from or to an RDS file.
#'
#' @param file File path to the cached fit. If the file extension is omitted,
#'   \code{.rds} is appended automatically.
#' @param x An object inheriting from \code{\link{MVBU_Stanfit}}.
#' @param compress Logical or character string specifying compression passed
#'   to \code{\link{saveRDS}} (default: \code{TRUE}).
#'
#' @return For \code{read_stanfit}, the loaded \code{\link{MVBU_Stanfit}} object
#'   with its \code{@file} slot updated, or \code{NULL} if the file does not
#'   exist or cannot be read. For \code{write_stanfit}, the saved object
#'   \code{x} invisibly with its \code{@file} slot updated.
#'
#' @seealso \code{\link{fit_ideal_adaptor}}
#' @examples
#' \donttest{
#' fit <- example_ideal_adaptor_stanfit()
#' tmp <- tempfile(fileext = ".rds")
#' write_stanfit(fit, tmp)
#' fit_loaded <- read_stanfit(tmp)
#' is(fit_loaded, "MVBU_Stanfit")
#' unlink(tmp)
#' }
#' @name read_stanfit
#' @rdname read_stanfit
#' @export
read_stanfit <- function(file) {
  file <- .check_stanfit_file(file)
  if (!file.exists(file)) {
    return(NULL)
  }
  x <- suppressWarnings(try(readRDS(file), silent = TRUE))
  if (.is_try_error(x)) {
    return(NULL)
  }

  .assert_true(
    S7::S7_inherits(x, MVBU_Stanfit),
    msg = paste0(
      "Object loaded from '", file, "' is not an MVBU_Stanfit object. ",
      "This indicates that it was fit with an outdated version of ",
      "MVBeliefUpdatr and needs to be refit."
    )
  )

  if (S7::S7_inherits(x, IdealAdaptorStanfit)) {
    label_info <- try(x@metadata$label_information, silent = TRUE)
    .assert_true(
      !is.null(label_info) && is.list(label_info) &&
        !is.null(label_info$cue) && !is.null(label_info$category) &&
        !is.null(label_info$group),
      msg = paste0(
        "Object loaded from '", file, "' is missing label_information in ",
        "metadata. Please refit the model."
      )
    )
  }

  x@file <- file
  x
}

#' @rdname read_stanfit
#' @export
write_stanfit <- function(x, file, compress = TRUE) {
  .assert_true(
    S7::S7_inherits(x, MVBU_Stanfit),
    msg = "x must inherit from MVBU_Stanfit."
  )
  file <- .check_stanfit_file(file)
  x@file <- file
  saveRDS(x, file = file, compress = compress)
  invisible(x)
}

# -------------------------------------------------------------------------
# Stanfit Parameter Extraction Helpers
# -------------------------------------------------------------------------

#' Extract fitted parameter names from Stanfit object
#'
#' Identifies canonical parameter names sampled in the Stan model draws.
#'
#' @param stanfit A \code{stanfit} object.
#'
#' @return A character vector of unique fitted parameter names
#'   (e.g., \code{c("kappa", "nu", "m", "S")}).
#'
#' @keywords internal
#' @noRd
.extract_fitted_parameter_names_from_stanfit <- function(stanfit) {
  if (is.null(stanfit)) {
    return(character(0))
  }
  fit_pars <- if (is.character(stanfit)) {
    stanfit
  } else if (inherits(stanfit, "stanfit")) {
    names(stanfit)
  } else {
    return(character(0))
  }
  matched <- grep(
    "^((kappa|nu|m|S|weights)_|weights\\[|lapse_rate|lapse_bias|prior_p|Sigma_noise)",
    fit_pars,
    value = TRUE
  )
  unique(
    gsub(
      "^(kappa|nu|m|S|weights|lapse_rate|lapse_bias|prior_p|Sigma_noise).*$",
      "\\1",
      matched
    )
  )
}

# -------------------------------------------------------------------------
# Stan Compilation and Threading Helpers
# -------------------------------------------------------------------------

#' Validate thread specification
#'
#' @param threads An integer count or list specifying thread options.
#' @return A validated list with element \code{threads}, or \code{NULL}.
#' @keywords internal
#' @noRd
.validate_threads <- function(threads) {
  if (is.null(threads)) {
    return(NULL)
  }
  if (is.list(threads) && !is.null(threads$threads)) {
    th <- .as_one_integer(threads$threads)
  } else if (is.numeric(threads)) {
    th <- .as_one_integer(threads)
  } else {
    .stop("Argument 'threads' must be a positive integer or a list containing 'threads'.")
  }
  if (th < 1L) {
    .stop("Number of threads must be a positive integer.")
  }
  list(threads = th)
}

#' Check whether threading should be used
#'
#' @param threads A validated or raw thread specification.
#' @param force Logical; if \code{TRUE}, ignore any force flag. Default: \code{FALSE}.
#' @return Logical scalar indicating whether threading is active.
#' @keywords internal
#' @noRd
.use_threading <- function(threads, force = FALSE) {
  threads <- .validate_threads(threads)
  if (is.null(threads)) {
    return(FALSE)
  }
  isTRUE(threads$threads > 1L)
}

#' Resolve file path for a Stan model
#'
#' @param stanmodel Character string naming the Stan model.
#' @return Character path to the .stan file, or error if not found.
#' @keywords internal
#' @noRd
.get_stan_file <- function(stanmodel) {
  stanmodel <- .as_one_character(stanmodel)
  if (!grepl("\\.stan$", stanmodel, ignore.case = TRUE)) {
    stan_filename <- paste0(stanmodel, ".stan")
  } else {
    stan_filename <- stanmodel
  }

  # First check installed package
  installed_file <- system.file("stan", stan_filename, package = "MVBeliefUpdatr")
  if (nzchar(installed_file) && file.exists(installed_file)) {
    return(installed_file)
  }

  # Fallback to source directory during development
  src_file <- file.path("inst", "stan", stan_filename)
  if (file.exists(src_file)) {
    return(normalizePath(src_file))
  }

  # Check direct path
  if (file.exists(stanmodel)) {
    return(normalizePath(stanmodel))
  }

  .stop("Could not locate Stan model file for '", stanmodel, "'.")
}

#' Get include paths for Stan compilation
#'
#' @return Character vector of directory paths containing Stan include headers.
#' @keywords internal
#' @noRd
.get_stan_include_paths <- function() {
  paths <- character(0)
  pkg_stan <- system.file("stan", package = "MVBeliefUpdatr")
  if (nzchar(pkg_stan) && dir.exists(pkg_stan)) {
    paths <- c(paths, pkg_stan)
  }
  dev_stan <- file.path("inst", "stan")
  if (dir.exists(dev_stan)) {
    paths <- c(paths, normalizePath(dev_stan))
  }
  unique(paths)
}

#' Normalize Stan code
#'
#' Normalizes Stan code to avoid unnecessary recompilation after whitespace
#' or comment changes.
#'
#' @param x A string containing the Stan code.
#' @return A normalized character string.
#' @keywords internal
#' @noRd
.normalize_stancode <- function(x) {
  x <- .as_one_character(x)
  # Remove single-line comments
  x <- gsub("//[^\n\r]*[\n\r]", " ", x)
  x <- gsub("//[^\n\r]*$", " ", x)
  # Remove multi-line comments
  x <- gsub("/\\*([^*]*(\\*[^/])?)*\\*/", " ", x)
  # Standardize whitespace (including newlines)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

#' Compile a Stan model with the requested backend
#'
#' @param model Stan model code (character string) or path to a \code{.stan} file.
#' @param backend Backend name to use for compilation: \code{"rstan"} or \code{"cmdstanr"}.
#' @param ... Additional arguments passed to the backend-specific compiler.
#' @return A compiled Stan model object (e.g., class \code{stanmodel} for rstan,
#'   or \code{CmdStanModel} for cmdstanr).
#' @keywords internal
#' @noRd
.compile_model <- function(model, backend = "rstan", ...) {
  compiler <- get(paste0(".compile_model_", backend), mode = "function")
  compiler(model, ...)
}

#' Compile a Stan model with rstan
#'
#' @param model Stan model code or path to a \code{.stan} file.
#' @param threads Optional thread specification (integer count or list with
#'   \code{threads}). Default: \code{NULL}.
#' @param silent Verbosity level for compilation output between 0 and 2.
#'   Default: \code{1}.
#' @param ... Additional arguments passed to \code{\link[rstan]{stan_model}}.
#' @return A compiled \code{stanmodel} object from \pkg{rstan}.
#' @keywords internal
#' @noRd
.compile_model_rstan <- function(
  model,
  threads = NULL,
  silent = 1,
  ...
) {
  args <- list(...)
  model_char <- .as_one_character(model)
  if (file.exists(model_char)) {
    args$file <- model_char
  } else {
    args$model_code <- model_char
  }

  inc_paths <- .get_stan_include_paths()
  if (length(inc_paths) > 0 && is.null(args$isystem)) {
    args$isystem <- inc_paths
  }

  if (silent < 2) {
    message("Compiling Stan program with rstan...")
  }

  threads_spec <- .validate_threads(threads)
  if (isTRUE(threads_spec$threads > 1L)) {
    if (utils::packageVersion("rstan") >= "2.26") {
      threads_per_chain_def <- rstan::rstan_options("threads_per_chain")
      on.exit(
        rstan::rstan_options(threads_per_chain = threads_per_chain_def),
        add = TRUE
      )
      rstan::rstan_options(threads_per_chain = threads_spec$threads)
    } else {
      .stop(
        "Threading is not supported by backend 'rstan' version ",
        utils::packageVersion("rstan"), "."
      )
    }
  }

  .eval_silent(
    .do_call(rstan::stan_model, args),
    type = "message",
    try = TRUE,
    silent = silent >= 2
  )
}

#' Compile a Stan model with cmdstanr
#'
#' @param model Stan model code or path to a \code{.stan} file.
#' @param threads Optional thread specification (integer count or list with \code{threads}).
#'   Default: \code{NULL}.
#' @param silent Verbosity level for compilation output between 0 and 2. Default: \code{1}.
#' @param ... Additional arguments passed to \code{cmdstanr::cmdstan_model()}.
#' @return A compiled \code{CmdStanModel} object from \pkg{cmdstanr}.
#' @keywords internal
#' @noRd
.compile_model_cmdstanr <- function(
  model,
  threads = NULL,
  silent = 1,
  ...
) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    .stop("Package 'cmdstanr' is required to use backend 'cmdstanr'.")
  }

  args <- list(...)
  model_char <- .as_one_character(model)
  if (file.exists(model_char)) {
    args$stan_file <- model_char
  } else {
    args$stan_file <- cmdstanr::write_stan_file(model_char)
  }

  inc_paths <- .get_stan_include_paths()
  if (length(inc_paths) > 0) {
    args$include_paths <- unique(c(args$include_paths, inc_paths))
  }

  if (.use_threading(threads)) {
    if (is.null(args$cpp_options)) {
      args$cpp_options <- list()
    }
    args$cpp_options$stan_threads <- TRUE
  }

  if (silent < 2) {
    message("Compiling Stan program with cmdstanr...")
  }

  .eval_silent(
    .do_call(cmdstanr::cmdstan_model, args),
    type = "message",
    try = TRUE,
    silent = silent >= 2
  )
}

#' Read CmdStan CSV files as an rstan stanfit object
#'
#' @param files Character vector of paths to CmdStan output CSV files.
#' @param exclude_pars Character vector of parameter names to drop.
#' @param model_name Optional model name to assign to the stanfit object.
#'   Default: \code{""}.
#' @return An object of class \code{\link[rstan:stanfit-class]{stanfit}}.
#' @keywords internal
#' @noRd
.read_cmdstan_csv_as_stanfit <- function(files, exclude_pars = character(0), model_name = "") {
  norm_files <- vapply(files, function(f) {
    lines <- readLines(f)
    c_idx <- which(startsWith(lines, "#"))
    if (length(c_idx) > 0) {
      lines[c_idx] <- gsub("save_warmup[[:space:]]*=[[:space:]]*false", "save_warmup = 0", lines[c_idx])
      lines[c_idx] <- gsub("save_warmup[[:space:]]*=[[:space:]]*true", "save_warmup = 1", lines[c_idx])
      lines[c_idx] <- gsub("opencl_device_id[[:space:]]*=[[:space:]]*false", "opencl_device_id = 0", lines[c_idx])
      lines[c_idx] <- gsub("opencl_platform_id[[:space:]]*=[[:space:]]*false", "opencl_platform_id = 0", lines[c_idx])
    }
    if (length(exclude_pars) > 0) {
      non_comment <- which(!startsWith(lines, "#"))
      if (length(non_comment) > 0) {
        h_idx <- min(non_comment)
        header_cols <- strsplit(lines[h_idx], ",")[[1]]
        col_bases <- sub("\\[.*", "", header_cols)
        drop_cols <- which(col_bases %in% exclude_pars)
        if (length(drop_cols) > 0) {
          keep_cols <- setdiff(seq_along(header_cols), drop_cols)
          data_indices <- h_idx:length(lines)
          valid_indices <- data_indices[!startsWith(lines[data_indices], "#") & nzchar(trimws(lines[data_indices]))]
          for (idx in valid_indices) {
            row_vals <- strsplit(lines[idx], ",")[[1]]
            lines[idx] <- paste(row_vals[keep_cols], collapse = ",")
          }
        }
      }
    }
    tmp <- tempfile(fileext = ".csv")
    writeLines(lines, tmp)
    tmp
  }, character(1))
  on.exit(unlink(norm_files), add = TRUE)
  res <- rstan::read_stan_csv(norm_files)
  if (nzchar(model_name)) {
    res@model_name <- model_name
  } else {
    res@model_name <- sub("_model$", "", res@model_name)
  }
  res
}

# -------------------------------------------------------------------------
# Model Specification Validation and Transformation
# -------------------------------------------------------------------------

#' Validate and transform fixed category moments
#'
#' Internal helper that checks the dimensionality of supplied category moments
#' and applies the affine transform to the category means and covariance
#' matrices.
#'
#' @param value Category representation moment values (mean vector or
#'   covariance matrix).
#' @param categories Character vector of category labels.
#' @param n_cues Number of cue dimensions expected for the supplied
#'   moment values.
#' @param transform Affine transform object.
#' @param which Whether the moment is for `mu_0` or `Sigma_0`.
#' @keywords internal
#' @noRd
.validate_and_transform_fixed_moments <- function(
  value,
  categories,
  n_cues,
  transform,
  which = c("mu_0", "Sigma_0")
) {
  .assert_numeric_scalar(n_cues)
  n_cues <- as.integer(n_cues)
  which <- match.arg(which)

  if (length(categories) == 1L) {
    if (which == "mu_0") {
      .assert_true(
        is.vector(value),
        msg = "mu_0 must be a vector when there is only one category."
      )
      value <- list(as.numeric(value))
    } else if (which == "Sigma_0") {
      .assert_true(
        is.array(value) || is.matrix(value),
        msg = "Sigma_0 must be a matrix when there is only one category."
      )
      value <- list(as.matrix(value))
    }
  } else {
    .assert_list(value)
    .assert_true(
      length(value) == length(categories),
      msg = sprintf(
        paste0(
          "%s must be a named list with names matching the categories ",
          "present in the exposure and test data."
        ),
        which
      )
    )
    value_names <- names(value)
    .assert_true(
      !is.null(value_names) && setequal(value_names, categories),
      msg = sprintf(
        paste0(
          "%s must be a named list with names matching the categories ",
          "present in the exposure and test data."
        ),
        which
      )
    )
  }

  if (which == "mu_0") {
    .assert_true(
      all(vapply(value, function(x) {
        is.numeric(x) && length(as.numeric(x)) == n_cues
      }, logical(1))),
      msg = sprintf(
        paste0(
          "mu_0 must be a named list with each entry a numeric vector ",
          "of length %s"
        ),
        n_cues
      )
    )
    value <- lapply(value, function(x) {
      transform_mean(as.numeric(x), transform)
    })
  } else {
    .assert_true(
      all(vapply(value, function(x) {
        x <- as.matrix(x)
        is.numeric(x) && nrow(x) == n_cues && ncol(x) == n_cues
      }, logical(1))),
      msg = sprintf(
        paste0(
          "Sigma_0 must be a named list with each entry a square matrix ",
          "of dimension %s x %s"
        ),
        n_cues,
        n_cues
      )
    )
    value <- lapply(value, function(x) {
      transform_cov(as.matrix(x), transform)
    })
  }
  value
}

#' Validate and transform model specifications for Stan fitting
#'
#' Internal helper that validates model compatibility, cue constraints,
#' and fixed parameters, transforming prior distributions and parameters
#' from original cue space to transformed space.
#'
#' @param stanmodel Character string naming the Stan model.
#' @param cues Character vector of cue names.
#' @param categories Character vector of category names.
#' @param groups Character vector of group names.
#' @param fixed_parameters Named list of fixed parameters.
#' @param priors List of prior specifications.
#' @param transform Affine transform object of class
#'   \code{\link{MVBU_TransformInformation}}.
#' @param control Optional list of control parameters.
#'
#' @return A list with components:
#'   \item{fixed_untransformed}{Original fixed parameters.}
#'   \item{fixed_transformed}{Transformed fixed parameters ready for Stan.}
#'   \item{priors_untransformed}{Priors in original space.}
#'   \item{priors_transformed}{Priors in transformed space.}
#'
#' @keywords internal
#' @noRd
.validate_and_transform_model_spec <- function(
  stanmodel,
  cues,
  categories,
  groups,
  fixed_parameters = NULL,
  priors = NULL,
  transform,
  control = list()
) {
  allowed_stanmodels <- c(
    "NIW_ideal_adaptor",
    "NIX_ideal_adaptor",
    # Late MNIX
    "MNIX_ideal_adaptor_ideal_integrator",
    "MNIX_ideal_adaptor_weighted_integrator",
    "MNIX_ideal_adaptor_category_weighted_integrator",
    "MNIX_ideal_adaptor_group_weighted_integrator",
    "MNIX_ideal_adaptor_group_category_weighted_integrator",
    "MNIX_ideal_adaptor_informativity_weighted_integrator",
    # Early MNIX
    "ideal_integrator_MNIX_ideal_adaptor",
    "weighted_integrator_MNIX_ideal_adaptor",
    "category_weighted_integrator_MNIX_ideal_adaptor",
    "group_weighted_integrator_MNIX_ideal_adaptor",
    "group_category_weighted_integrator_MNIX_ideal_adaptor",
    "informativity_weighted_integrator_MNIX_ideal_adaptor"
  )
  stanmodel <- match.arg(stanmodel, allowed_stanmodels)

  n_cues <- length(cues)
  n_categories <- length(categories)
  n_groups <- length(groups)

  .assert_any(
    stanmodel != "NIX_ideal_adaptor",
    n_cues == 1L,
    msg = "NIX_ideal_adaptor requires exactly one cue."
  )
  .assert_any(
    !startsWith(stanmodel, "MNIX_ideal_adaptor"),
    n_cues >= 2L,
    msg = sprintf("%s requires at least two cues.", stanmodel)
  )
  .assert_any(
    !grepl("MNIX_ideal_adaptor$", stanmodel),
    n_cues >= 2L,
    msg = sprintf("%s requires at least two cues.", stanmodel)
  )

  if (is.null(fixed_parameters)) {
    fixed_parameters <- MVBU_FixedParameter()
  } else if (!S7::S7_inherits(fixed_parameters, MVBU_FixedParameter)) {
    .stop(
      "Argument `fixed_parameters` must be an MVBU_FixedParameter ",
      "object or NULL."
    )
  }

  if (is.null(priors)) {
    priors <- MVBU_Prior()
  } else if (!S7::S7_inherits(priors, MVBU_Prior)) {
    .stop(
      "Argument `priors` must be an MVBU_Prior object or NULL."
    )
  }

  fixed_structured <- .resolve_fixed_parameters(
    fixed = fixed_parameters,
    categories = categories,
    cues = cues,
    transform = transform
  )

  # Extract bridge variables for NIX and MNIX models
  lapse_rate <- if (fixed_structured$lapse_rate$known == 1L) {
    as.numeric(fixed_structured$lapse_rate$data)
  } else {
    NULL
  }

  lapse_bias <- if (fixed_structured$lapse_bias$known == 1L) {
    as.numeric(fixed_structured$lapse_bias$data)
  } else {
    NULL
  }

  prior_p <- if (fixed_structured$category_prior$known == 1L) {
    as.numeric(fixed_structured$category_prior$data)
  } else {
    NULL
  }

  is_global_weighted <- stanmodel %in% c(
    "MNIX_ideal_adaptor_weighted_integrator",
    "weighted_integrator_MNIX_ideal_adaptor"
  )
  is_cat_weighted <- stanmodel %in% c(
    "MNIX_ideal_adaptor_category_weighted_integrator",
    "category_weighted_integrator_MNIX_ideal_adaptor"
  )
  is_grp_weighted <- stanmodel %in% c(
    "MNIX_ideal_adaptor_group_weighted_integrator",
    "group_weighted_integrator_MNIX_ideal_adaptor"
  )
  is_grpc_weighted <- stanmodel %in% c(
    "MNIX_ideal_adaptor_group_category_weighted_integrator",
    "group_category_weighted_integrator_MNIX_ideal_adaptor"
  )
  is_any_weighted <- is_global_weighted || is_cat_weighted ||
    is_grp_weighted || is_grpc_weighted

  w_idx <- which(fixed_parameters@par_class == "weights")
  weights <- if (length(w_idx) > 0L) {
    if (!is_any_weighted) {
      .stop(sprintf("Model '%s' does not use fixed weights.", stanmodel))
    }
    w_val <- fixed_parameters@value[[w_idx[length(w_idx)]]]
    if (is_global_weighted) {
      .assert_all(
        is.numeric(w_val),
        length(w_val) == n_cues,
        all(w_val >= 0),
        isTRUE(all.equal(sum(w_val), 1, tolerance = MVBU_PROB_TOL)),
        msg = sprintf(
          paste0(
            "weights must be a non-negative numeric vector of length %d ",
            "summing to 1."
          ),
          n_cues
        )
      )
    }
    w_val
  } else {
    NULL
  }

  # For NIX and MNIX (temporary bridge before Stages 2 and 3)
  mu_0_transformed <- if (fixed_structured$m_0$M_known == n_categories) {
    res_m <- vector("list", n_categories)
    names(res_m) <- categories
    for (i in seq_along(fixed_structured$m_0$idx_known)) {
      c_idx <- fixed_structured$m_0$idx_known[i]
      res_m[[c_idx]] <- as.numeric(fixed_structured$m_0$data[i, ])
    }
    res_m
  } else {
    NULL
  }

  Sigma_0_transformed <- if (fixed_structured$S_0$M_known == n_categories) {
    res_s <- vector("list", n_categories)
    names(res_s) <- categories
    for (i in seq_along(fixed_structured$S_0$idx_known)) {
      c_idx <- fixed_structured$S_0$idx_known[i]
      res_s[[c_idx]] <- matrix(
        fixed_structured$S_0$data[i, , ],
        nrow = n_cues,
        ncol = n_cues
      )
    }
    res_s
  } else {
    NULL
  }

  fixed_transformed <- list(
    lapse_rate = lapse_rate,
    lapse_bias = lapse_bias,
    prior_p = prior_p,
    mu_0 = mu_0_transformed,
    Sigma_0 = Sigma_0_transformed,
    weights = weights,
    structured = fixed_structured
  )

  resolved_priors <- .resolve_priors(
    priors = priors,
    stanmodel = stanmodel,
    cues = cues,
    category_labels = categories,
    transform = transform
  )

  list(
    fixed_untransformed = fixed_parameters,
    fixed_transformed = fixed_transformed,
    priors_untransformed = resolved_priors$original,
    priors_transformed = list(
      table = resolved_priors$transformed,
      stan = resolved_priors$stan
    )
  )
}
