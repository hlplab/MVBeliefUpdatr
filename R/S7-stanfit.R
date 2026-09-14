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
    metadata = S7::class_list,
    cache = S7::class_list
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
    metadata = list(),
    cache = list()
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
    if (!is.list(cache)) {
      cache <- list()
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
      metadata = metadata,
      cache = as.list(cache)
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
      !S7::S7_inherits(
        self@transform_information,
        MVBU_TransformInformation
      )) {
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
    metadata = list(),
    cache = list()
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
        metadata = as.list(metadata),
        cache = as.list(cache)
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
    metadata = list(),
    cache = list()
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
        metadata = as.list(metadata),
        cache = as.list(cache)
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
    metadata = list(),
    cache = list()
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
        metadata = as.list(metadata),
        cache = as.list(cache)
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
    metadata = list(),
    cache = list()
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
        metadata = as.list(metadata),
        cache = as.list(cache)
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

#' Check whether an IdealAdaptorStanfit or stanfit object contains posterior draws
#'
#' Internal utility for verifying that a fitted Stanfit object has samples
#' available for downstream operations (e.g., draws extraction, diagnostics).
#'
#' @param x An \code{IdealAdaptorStanfit} or S4 \code{stanfit} object.
#' @return A logical scalar: `TRUE` if draws are present.
#' @keywords internal
#' @noRd
.contains_draws <- function(x) {
  sf <- if (S7::S7_inherits(x, IdealAdaptorStanfit)) x@stanfit else x
  if (inherits(sf, "stanfit")) {
    length(sf@sim) > 0 && length(sf@sim$samples) > 0
  } else {
    FALSE
  }
}
