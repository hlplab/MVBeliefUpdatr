#' @include asserts.R
#' @include S7-core-classes.R
#' @include S7-stanfit.R
#' @include S7-staninput.R
#' @include S7-stanfit-input.R
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
#' @param current_version Current version of relevant packages (default: will be automatically
#'   obtained from current packages via \code{\link{get_current_versions}}).
#' @param data New data to check consistency of factor level names (default: \code{NULL}).
#' @param staninput New Stan data (result of a call to \code{\link{make_staninput}} or
#'   \code{\link{IdealAdaptorStaninput}}). Pass \code{NULL} to avoid this check (default: \code{NULL}).
#' @param silent Logical. If \code{TRUE}, no messages will be given (default: \code{FALSE}).
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

  if (!isTRUE(all.equal(x@version, current_version))) {
    if (!silent) {
      message(
        "Version of MVBeliefUpdatr or rstan has changed (current version is ",
        paste(purrr::map_chr(current_version, ~ paste(.x, collapse = ", ")), collapse = "; "),
        ")."
      )
      if (verbose) {
        print(x@version)
      }
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

# -------------------------------------------------------------------------
# Read and Write Functions
# -------------------------------------------------------------------------

#' Read a cached Stanfit object from disk
#'
#' Reads a previously saved \code{\link{MVBU_Stanfit}} object (or subclass such as
#' \code{\link{IdealAdaptorStanfit}}) from an RDS file.
#'
#' @param file File path to the cached fit. If the file extension is omitted,
#'   \code{.rds} is appended automatically.
#'
#' @return The loaded \code{\link{MVBU_Stanfit}} object with its \code{@file}
#'   slot updated, or \code{NULL} if the file does not exist or cannot be read.
#'
#' @seealso \code{\link{write_stanfit}}, \code{\link{fit_ideal_adaptor}}
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
      "This indicates that it was fit with an outdated version of MVBeliefUpdatr and needs to be refit."
    )
  )

  if (S7::S7_inherits(x, IdealAdaptorStanfit)) {
    label_info <- try(x@metadata$label_information, silent = TRUE)
    .assert_true(
      !is.null(label_info) && is.list(label_info) &&
        !is.null(label_info$cue) && !is.null(label_info$category) && !is.null(label_info$group),
      msg = paste0(
        "Object loaded from '", file, "' is missing label_information in metadata. ",
        "Please refit the model."
      )
    )
  }

  x@file <- file
  x
}

#' Write a fitted Stanfit object to disk
#'
#' Saves an \code{\link{MVBU_Stanfit}} object (or subclass) to disk via \code{\link{saveRDS}},
#' and records the target file path in the object's \code{@file} slot.
#'
#' @param x An object inheriting from \code{\link{MVBU_Stanfit}}.
#' @param file File path for saving the object. If the file extension is omitted,
#'   \code{.rds} is appended automatically.
#' @param compress Logical or character string specifying compression passed
#'   to \code{\link{saveRDS}} (default: \code{TRUE}).
#'
#' @return The saved object \code{x} invisibly, with \code{@file} updated.
#'
#' @seealso \code{\link{read_stanfit}}, \code{\link{fit_ideal_adaptor}}
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
