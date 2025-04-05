#' An S4 class for NIW_ideal_adaptor stanfit objects that use one of the NIW_ideal_adaptor stan programs.
#'
#' @name NIW_ideal_adaptor_stanfit-class
#' @aliases NIW_ideal_adaptor_stanfit
#' @docType class
#'
#' @details
#' See \code{methods(class = "NIW_ideal_adaptor_stanfit")} for an overview of available methods.
#'
#' @slot data A \code{data.frame} containing the data used to fit the model.
#' @slot staninput A named \code{list} containing the data handed to rstan through
#'   \code{\link{make_staninput}}.
#' @slot stanvars A \code{\link{stanvars}} object.
#' @slot backend The name of the backend used to fit the model.
#' @slot stan_args Named list of additional control arguments that were passed
#'   to the Stan backend directly. NOT YET USED
#' @slot stanfit An object of class \code{\link[rstan:stanfit-class]{stanfit}}
#'   among others containing the posterior draws.
#' @slot basis An object that contains a small subset of the Stan data
#'   created at fitting time, which is needed to process new data correctly. NOT YET USED
#' @slot transform_information list containing elements transform.parameters, transform.function, and
#'   untransform.function.
#' @slot criteria An empty \code{list} for adding model fit criteria
#'   after estimation of the model. NOT YET USED
#' @slot file Optional name of a file in which the model object was stored in
#'   or loaded from. NOT YET USED
#' @slot version The versions of \pkg{MVBeliefUpdatr} and \pkg{rstan} with
#'   which the model was fitted.
#' @slot labels list
#'
#' @importFrom rstan nlist
#' @importFrom utils packageVersion
NULL

# NIW_ideal_adaptor_stanfit class
NIW_ideal_adaptor_stanfit <- function(
    data = data.frame(),
    staninput = list(),
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
    basis = NULL,
    transform_information = NULL,
    criteria = list(),
    file = NULL
) {
  if (!is.null(stanfit)) {
    assert_that(is.stanfit(stanfit))
    assert_that(
      stanfit@model_name %in% names(MVBeliefUpdatr:::stanmodels),
      msg = paste0("stanfit object was not created by one of the accepted stancodes:\n\t",
                   paste(names(MVBeliefUpdatr:::stanmodels), collapse = "\n\t"),
                   "\n(you can get the name of your model from your_stanfit@model_name)."))
  }
  assert_that(
    is.NIW_ideal_adaptor_staninput(staninput),
    msg = paste("staninput is not an acceptable input for NIW_ideal_adaptor_stanfit stan program."))

  version <- get_current_versions()

  x <-
    nlist(
      data,
      staninput,
      stanvars,
      save_pars,
      backend,
      stan_args,
      stanfit,
      basis,
      transform_information,
      criteria,
      file,
      version)

  # setClass(
  #   "NIW_ideal_adaptor_stanfit",
  #   slots = c(input_data = "list", transform_information = "list", labels = "list"),
  #   contains = "stanfit",
  #   package = "MVBeliefUpdatr")
  class(x) <- "NIW_ideal_adaptor_stanfit"

  x
}

# Ultimately, this could be turned into a method, but care would have to be taken to still make tidybayes::recover_types()
# work since that applies to objects of class "stanfit" (but is not defined as a method, I think)
#' @export
recover_types.NIW_ideal_adaptor_stanfit <- function(fit, staninput = NULL) {
  stanfit <- get_stanfit(fit)
  if (is.null(staninput)) staninput <- get_staninput(fit)

  # the levels information recovered below should probably should be stored in a more systematic way, either as attributes
  # to the model or as some part of a list
  stanfit %<>%
    recover_types(
      crossing(
        category = factor(colnames(staninput$z_test_counts), levels = colnames(staninput$z_test_counts)),
        group = factor(attr(staninput$y_test, "levels"), levels = attr(staninput$y_test, "levels")),
        cue = factor(dimnames(staninput$x_test)[[2]], levels = dimnames(staninput$x_test)[[2]]),
        cue2 = cue))

  fit <- set_stanfit(fit, stanfit)
  return(fit)
}

# Example usage
# make_parnames("lapse_rate_param", 1),
# make_parnames("m_0_param", 1:get_staninput(fit)$M, 1:get_staninput(fit)$K),
# make_parnames("tau_0_param", 1:get_staninput(fit)$M, 1:get_staninput(fit)$K),
# make_parnames("L_omega_0_param", 1:get_staninput(fit)$M, 1:get_staninput(fit)$K, 1:get_staninput(fit)$K),
# make_parnames("m_0_tau", 1:get_staninput(fit)$K),
# make_parnames("m_0_tau_param", 1:get_staninput(fit)$K),
# make_parnames("m_0_L_omega", 1:get_staninput(fit)$K, 1:get_staninput(fit)$K),
# make_parnames("m_0_L_omega_param", 1:get_staninput(fit)$K, 1:get_staninput(fit)$K),
# make_parnames("L_S_0", 1:get_staninput(fit)$M, 1:get_staninput(fit)$K, 1:get_staninput(fit)$K),
# make_parnames("L_S_n", 1:get_staninput(fit)$M, 1:get_staninput(fit)$L, 1:get_staninput(fit)$K, 1:get_staninput(fit)$K),
# make_parnames("L_t_scale", 1:get_staninput(fit)$M, 1:get_staninput(fit)$L, 1:get_staninput(fit)$K, 1:get_staninput(fit)$K),
# make_parnames("p_test_conj", 1:get_staninput(fit)$N_test, 1:get_staninput(fit)$M),
# make_parnames("log_p_test_conj", 1:get_staninput(fit)$N_test, 1:get_staninput(fit)$M))
make_parnames <- function(prefix, ...) {
  # Create all unique combinations of the values of ..., each of which is a list
  combinations <- expand.grid(..., stringsAsFactors = FALSE)
  paste0(prefix, "[", apply(combinations, 1, paste0, collapse = ","), "]")
}


#' Is this an NIW ideal adaptor stanfit?
#'
#' Check whether \code{x} is of class \code{\link{NIW_ideal_adaptor_stanfit}}.
#'
#' @param x Object to be checked.
#' @param verbose Currently being ignored.
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#' @export
is.NIW_ideal_adaptor_stanfit <- function(x, verbose = F) {
  inherits(x, "NIW_ideal_adaptor_stanfit")
}

is.stanfit <- function(x) {
  inherits(x, "stanfit")
}

#' Is this a tibble of MCMC draws of an NIW ideal adaptor?
#'
#' Check whether \code{x} is a tibble of post-warmup draws of parameters obtained from incremental
#' conjugate Bayesian belief-updating (IBBU) over a Normal-Inverse-Wishart (NIW) prior.
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#' @export
is.NIW_ideal_adaptor_MCMC <- function(x, is.nested = T, is.long = T, with.prior = F, with.lapse = if (with.lapse_bias) T else F, with.lapse_bias = F) {
  if(
    all(
       is.NIW_ideal_adaptor(x, is.long = is.long, category = "category", with.prior = with.prior, with.lapse = with.lapse, with.lapse_bias = with.lapse_bias),
       all(c(".chain", ".iteration", ".draw",
             "group") %in% names(x)),
       xor(is.nested, all(c("cue", "cue2") %in% names(x)))
    )
  ) return(T) else return(F)
}


#' Is this a list of NIW ideal adaptor stanfit inputs?
#'
#' Check whether \code{x} is of class \code{\link{NIW_ideal_adaptor_stanfit}}.
#'
#' @return A logical.
#'
#' @seealso TBD
#' @keywords TBD
#' @export
is.NIW_ideal_adaptor_staninput <- function(x) {
  # Test of NIW_ideal_adaptor_staninput class not yet implemented. Always returning T.

  # Proposed names for slots in input object (at least internally / not necessarily handed to stan like this:
  #
  #   exposure_N (N)
  #   exposure_category_mean (x_mean)
  #   exposure_cue_ss (x_ss)
  #   test_N (N_test)
  #   test_cue (x_test)
  #   test_response (z_test_counts)
  #   test_group (y_test)

  return(TRUE)
}



# from brms
# check validity of file name to store a stanfit object in
check_stanfit_file <- function(file) {
  file <- as_one_character(file)
  file_ending <- tolower(get_matches("\\.[^\\.]+$", file))
  if (!isTRUE(file_ending == ".rds")) {
    file <- paste0(file, ".rds")
  }
  file
}

# possible options for argument 'file_refit'
file_refit_options <- function() {
  c("never", "always", "on_change")
}

# Modified from brms
#' Check if cached \code{NIW_ideal_adaptor_stanfit} can be used.
#'
#' Checks whether a given cached fit can be used without refitting when
#' \code{file_refit = "on_change"} is used.
#'
#' @param x Old \code{NIW_ideal_adaptor_stanfit} object (e.g., loaded from file).
#' @param current_version Current version of relevant packages. (default: will be automatically
#'  obtained from current packages).
#' @param data New data to check consistency of factor level names. (default: \code{NULL}))
#' @param staninput New Stan data (result of a call to \code{\link[standata.default]{standata}}).
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
stanfit_needs_refit <- function(
    x,
    current_version = get_current_versions(),
    data = NULL, staninput = NULL,
    silent = FALSE, verbose = FALSE
) {
  stopifnot(is.NIW_ideal_adaptor_stanfit(x))
  silent <- as_one_logical(silent)
  verbose <- as_one_logical(verbose)

  # Check if the stanfit is older than version 0.0.1.0015 when we introduced the version slot
  if (!hasSlot(x, "version")) {
    if (!silent) {
      message("The model ", deparse1(substitute(x)), " was fit with an old version of MVBeliefUpdater (< 0.0.1.0015).")
    }
    refit <- TRUE
  }
  if (!all.equal(x$version, current_version)) {
    if (!silent) {
      message("Version of MVBeliefUpdatr or rstan has changed (current version is", paste(current_version, collapse = ", "), ").")
      if (verbose) {
        print(x$version)
      }
    }
    refit <- TRUE
  }


  # if (!is.null(scode)) {
  #   scode <- as_one_character(scode)
  #   cached_scode <- stancode(x)
  # }
  if (!is.null(staninput)) {
    stopifnot(is.list(staninput))
    cached_staninput <- get_staninput(x)
  }
  if (!is.null(data)) {
    stopifnot(is.data.frame(data))
    cached_data <- x$data
  }

  refit <- FALSE
  # if (!is.null(scode)) {
  #   if (normalize_stancode(scode) != normalize_stancode(cached_scode)) {
  #     if (!silent) {
  #       message("Stan code has changed beyond whitespace/comments.")
  #       if (verbose) {
  #         require_package("diffobj")
  #         print(diffobj::diffChr(scode, cached_scode, format = "ansi8"))
  #       }
  #     }
  #     refit <- TRUE
  #   }
  # }
  if (!is.null(staninput)) {
    staninput_equality <- all.equal(staninput, cached_staninput, check.attributes = FALSE)
    if (!isTRUE(staninput_equality)) {
      if (!silent) {
        message("The processed data for Stan has changed.")
        if (verbose) {
          print(staninput_equality)
        }
      }
      refit <- TRUE
    }
  }
  if (!is.null(data)) {
    # check consistency of factor names
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
                "Names of factor levels have changed for variable '", var, "' ",
                "with cached levels (", collapse_comma(cached_levels), ") ",
                "but new levels (", collapse_comma(new_levels), ").\n"
              ))
            }
          }
          refit <- TRUE
          if (!verbose) {
            # no need to check all variables if we trigger a refit anyway
            break
          }
        }
      }
    }
    if (factor_level_message) {
      message("Names of factor levels have changed.")
    }
  }
  refit
}

# modified from brms
# read a NIW_ideal_adaptor_stanfit object from a file
# @param file path to an rds file
# @return a NIW_ideal_adaptor_stanfit object or NULL
read_NIW_ideal_adaptor_stanfit <- function(file) {
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
    if (!is.NIW_ideal_adaptor_stanfit(x)) {
      stop2("Object loaded via 'file' is not of class 'NIW_ideal_adaptor_stanfit'.")
    }
    x$file <- file
  } else {
    x <- NULL
  }
  x
}

# modified from brms
# write an NIW_ideal_adaptor_stanfit object to a file.
# @param x an NIW_ideal_adaptor_stanfit object
# @param file path to an rds file
# @param compress compression format supported by saveRDS
# @return NULL
write_NIW_ideal_adaptor_stanfit <- function(x, file, compress = TRUE) {
  stopifnot(is.NIW_ideal_adaptor_stanfit(x))
  file <- check_stanfit_file(file)
  x$file <- file
  saveRDS(x, file = file, compress = compress)
  invisible(x)
}


