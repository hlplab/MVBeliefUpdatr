#' @include asserts.R
#' @include S7-class.R
#' @include S7-transform-information.R
#' @include S7-staninput.R
#' @include S7-stanfit-input.R
#' @importFrom S7 new_class new_object new_generic method class_name props
NULL

#' MVBeliefUpdatr Stanfit Classes and Constructors
#'
#' An S7 base class for Stan fit objects.
#'
#' @param data A named \code{list} containing \code{$original} and
#'   \code{$transformed} data frames. Default: \code{list()}.
#' @param staninput A Stan input object containing the data handed to rstan
#'   (an \code{\link{IdealAdaptorStaninput}} object). Default: \code{NULL}.
#' @param stanvars A \code{stanvars} object or \code{NULL}. Default: \code{NULL}.
#' @param backend The name of the backend used to fit the model (character).
#'   Default: \code{"rstan"}.
#' @param save_pars Optional storage for saved parameter names. Default:
#'   \code{NULL}.
#' @param stan_args Named list of additional control arguments that were passed
#'   to the Stan backend directly. Default: \code{list()}.
#' @param stanfit An object of class \code{\link[rstan:stanfit-class]{stanfit}}
#'   containing posterior draws. Default: \code{NULL}.
#' @param transform_information An object of type
#'   \code{\link{MVBU_TransformInformation}} or \code{NULL}. Default: \code{NULL}.
#' @param criteria An empty \code{list} for adding model fit criteria after
#'   estimation of the model. Default: \code{list()}.
#' @param file Optional name of a file in which the model object was stored in
#'   or loaded from. Default: \code{NULL}.
#' @param version The versions of \pkg{MVBeliefUpdatr} and \pkg{rstan} with
#'   which the model was fitted. Default: \code{NULL}.
#' @param metadata List containing auxiliary information including
#'   \code{label_information}. Default: \code{list()}.
#' @param cache Internal cache list storing pre-extracted \code{$draws} and
#'   \code{$summary} objects populated by \code{\link{add_posterior_latents}}.
#'   Default: \code{list()}.
#' @seealso \code{\link{add_posterior_latents}}, \code{\link{fit_ideal_adaptor}},
#'   \code{vignette("fitting-and-working-with-stanfit-models")}
#'
#' @examples
#' # Instantiate an empty Stanfit container object
#' fit <- MVBU_Stanfit()
#' is(fit, "MVBU_Stanfit")
#'
#' \donttest{
#' # Load a pre-compiled example Stanfit object
#' ex_fit <- example_ideal_adaptor_stanfit()
#' class(ex_fit)
#' }
#'
#' @name MVBU_Stanfit
#' @rdname MVBU_Stanfit
#' @docType class
#' @rawNamespace if (getRversion() < "4.3.0") importFrom(S7, "@")
#' @importFrom purrr map_lgl map_chr
#' @importFrom utils packageVersion
#' @export
MVBU_Stanfit <- S7::new_class(
  "MVBU_Stanfit",
  package = NULL,
  parent = MVBU_Object,
  properties = list(
    data = S7::class_list,
    staninput = S7::class_any,
    stanvars = S7::class_any,
    backend = S7::class_character,
    save_pars = S7::class_any,
    stan_args = S7::class_list,
    stanfit = S7::class_any,
    transform_information = MVBU_TransformInformation,
    criteria = S7::class_list,
    file = S7::class_character,
    version = S7::class_any,
    metadata = S7::class_list,
    cache = S7::class_list
  ),
  constructor = function(
    data = list(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
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

    if (is.data.frame(data)) {
      data <- list(original = data, transformed = data)
    } else if (is.list(data)) {
      if (is.null(data$original) && is.null(data$transformed) &&
        length(data) == 0L) {
        data <- list(
          original = tibble::tibble(),
          transformed = tibble::tibble()
        )
      } else if (is.null(data$original) && !is.null(data$transformed)) {
        data$original <- data$transformed
      } else if (!is.null(data$original) && is.null(data$transformed)) {
        data$transformed <- data$original
      }
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
      response_category = if (!is.null(label_info$response_category)) as.character(label_info$response_category) else character(0),
      group = if (!is.null(label_info$group)) as.character(label_info$group) else character(0)
    )

    S7::new_object(
      MVBU_Object(),
      data = data,
      staninput = staninput,
      stanvars = stanvars,
      backend = as.character(backend),
      save_pars = save_pars,
      stan_args = as.list(stan_args),
      stanfit = stanfit,
      transform_information = transform_information,
      criteria = as.list(criteria),
      file = as.character(file),
      version = version,
      metadata = metadata,
      cache = as.list(cache)
    )
  },
  validator = function(self) {
    if (!is.list(self@data) || is.data.frame(self@data)) {
      return("`data` must be a list containing $original and $transformed.")
    }
    if (!inherits(self@data$original, "data.frame")) {
      return("`data$original` must inherit from data.frame.")
    }
    if (!inherits(self@data$transformed, "data.frame")) {
      return("`data$transformed` must inherit from data.frame.")
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
      if (!(self@stanfit@model_name %in% names(stanmodels))) {
        return(
          paste0(
            "`stanfit` model_name is not recognized. `stanfit` has to be ",
            "created by one of the accepted stanmodels:\n\t",
            paste(names(stanmodels), collapse = "\n\t"),
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

#' @docType class
#' @rdname MVBU_Stanfit
#' @export
IdealAdaptorStanfit <- S7::new_class(
  "IdealAdaptorStanfit",
  package = NULL,
  parent = MVBU_Stanfit,
  constructor = function(
    data = list(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
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
    data = list(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
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

# --- Base MNIX Stanfit Class ---

#' @docType class
#' @rdname MVBU_Stanfit
#' @export
MNIX_IdealAdaptorStanfit <- S7::new_class(
  "MNIX_IdealAdaptorStanfit",
  package = NULL,
  parent = IdealAdaptorStanfit,
  constructor = function(
    data = list(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
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

# --- Early Integrator Stanfit Classes ---

#' @docType class
#' @rdname MVBU_Stanfit
#' @export
IdealIntegrator_MNIX_IdealAdaptorStanfit <- S7::new_class(
  "IdealIntegrator_MNIX_IdealAdaptorStanfit",
  package = NULL,
  parent = MNIX_IdealAdaptorStanfit,
  constructor = function(
    data = list(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
    transform_information = NULL,
    criteria = list(),
    file = NULL,
    version = NULL,
    metadata = list(),
    cache = list()
  ) {
    S7::new_object(
      MNIX_IdealAdaptorStanfit(
        data = data,
        staninput = staninput,
        stanvars = stanvars,
        backend = backend,
        save_pars = save_pars,
        stan_args = stan_args,
        stanfit = stanfit,
        transform_information = transform_information,
        criteria = criteria,
        file = file,
        version = version,
        metadata = metadata,
        cache = cache
      )
    )
  }
)

#' @docType class
#' @rdname MVBU_Stanfit
#' @export
WeightedIntegrator_MNIX_IdealAdaptorStanfit <- S7::new_class(
  "WeightedIntegrator_MNIX_IdealAdaptorStanfit",
  package = NULL,
  parent = MNIX_IdealAdaptorStanfit,
  constructor = function(
    data = list(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
    transform_information = NULL,
    criteria = list(),
    file = NULL,
    version = NULL,
    metadata = list(),
    cache = list()
  ) {
    S7::new_object(
      MNIX_IdealAdaptorStanfit(
        data = data,
        staninput = staninput,
        stanvars = stanvars,
        backend = backend,
        save_pars = save_pars,
        stan_args = stan_args,
        stanfit = stanfit,
        transform_information = transform_information,
        criteria = criteria,
        file = file,
        version = version,
        metadata = metadata,
        cache = cache
      )
    )
  }
)

#' @docType class
#' @rdname MVBU_Stanfit
#' @export
CategoryWeightedIntegrator_MNIX_IdealAdaptorStanfit <- S7::new_class(
  "CategoryWeightedIntegrator_MNIX_IdealAdaptorStanfit",
  package = NULL,
  parent = MNIX_IdealAdaptorStanfit,
  constructor = function(
    data = list(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
    transform_information = NULL,
    criteria = list(),
    file = NULL,
    version = NULL,
    metadata = list(),
    cache = list()
  ) {
    S7::new_object(
      MNIX_IdealAdaptorStanfit(
        data = data,
        staninput = staninput,
        stanvars = stanvars,
        backend = backend,
        save_pars = save_pars,
        stan_args = stan_args,
        stanfit = stanfit,
        transform_information = transform_information,
        criteria = criteria,
        file = file,
        version = version,
        metadata = metadata,
        cache = cache
      )
    )
  }
)

#' @docType class
#' @rdname MVBU_Stanfit
#' @export
GroupWeightedIntegrator_MNIX_IdealAdaptorStanfit <- S7::new_class(
  "GroupWeightedIntegrator_MNIX_IdealAdaptorStanfit",
  package = NULL,
  parent = MNIX_IdealAdaptorStanfit,
  constructor = function(
    data = list(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
    transform_information = NULL,
    criteria = list(),
    file = NULL,
    version = NULL,
    metadata = list(),
    cache = list()
  ) {
    S7::new_object(
      MNIX_IdealAdaptorStanfit(
        data = data,
        staninput = staninput,
        stanvars = stanvars,
        backend = backend,
        save_pars = save_pars,
        stan_args = stan_args,
        stanfit = stanfit,
        transform_information = transform_information,
        criteria = criteria,
        file = file,
        version = version,
        metadata = metadata,
        cache = cache
      )
    )
  }
)

#' @docType class
#' @rdname MVBU_Stanfit
#' @export
GroupCategoryWeightedIntegrator_MNIX_IdealAdaptorStanfit <- S7::new_class(
  "GroupCategoryWeightedIntegrator_MNIX_IdealAdaptorStanfit",
  package = NULL,
  parent = MNIX_IdealAdaptorStanfit,
  constructor = function(
    data = list(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
    transform_information = NULL,
    criteria = list(),
    file = NULL,
    version = NULL,
    metadata = list(),
    cache = list()
  ) {
    S7::new_object(
      MNIX_IdealAdaptorStanfit(
        data = data,
        staninput = staninput,
        stanvars = stanvars,
        backend = backend,
        save_pars = save_pars,
        stan_args = stan_args,
        stanfit = stanfit,
        transform_information = transform_information,
        criteria = criteria,
        file = file,
        version = version,
        metadata = metadata,
        cache = cache
      )
    )
  }
)

#' @docType class
#' @rdname MVBU_Stanfit
#' @export
InformativityWeightedIntegrator_MNIX_IdealAdaptorStanfit <- S7::new_class(
  "InformativityWeightedIntegrator_MNIX_IdealAdaptorStanfit",
  package = NULL,
  parent = MNIX_IdealAdaptorStanfit,
  constructor = function(
    data = list(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
    transform_information = NULL,
    criteria = list(),
    file = NULL,
    version = NULL,
    metadata = list(),
    cache = list()
  ) {
    S7::new_object(
      MNIX_IdealAdaptorStanfit(
        data = data,
        staninput = staninput,
        stanvars = stanvars,
        backend = backend,
        save_pars = save_pars,
        stan_args = stan_args,
        stanfit = stanfit,
        transform_information = transform_information,
        criteria = criteria,
        file = file,
        version = version,
        metadata = metadata,
        cache = cache
      )
    )
  }
)

# --- Late Integrator Stanfit Classes ---

#' @docType class
#' @rdname MVBU_Stanfit
#' @export
MNIX_IdealAdaptor_IdealIntegratorStanfit <- S7::new_class(
  "MNIX_IdealAdaptor_IdealIntegratorStanfit",
  package = NULL,
  parent = MNIX_IdealAdaptorStanfit,
  constructor = function(
    data = list(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
    transform_information = NULL,
    criteria = list(),
    file = NULL,
    version = NULL,
    metadata = list(),
    cache = list()
  ) {
    S7::new_object(
      MNIX_IdealAdaptorStanfit(
        data = data,
        staninput = staninput,
        stanvars = stanvars,
        backend = backend,
        save_pars = save_pars,
        stan_args = stan_args,
        stanfit = stanfit,
        transform_information = transform_information,
        criteria = criteria,
        file = file,
        version = version,
        metadata = metadata,
        cache = cache
      )
    )
  }
)

#' @docType class
#' @rdname MVBU_Stanfit
#' @export
MNIX_IdealAdaptor_WeightedIntegratorStanfit <- S7::new_class(
  "MNIX_IdealAdaptor_WeightedIntegratorStanfit",
  package = NULL,
  parent = MNIX_IdealAdaptorStanfit,
  constructor = function(
    data = list(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
    transform_information = NULL,
    criteria = list(),
    file = NULL,
    version = NULL,
    metadata = list(),
    cache = list()
  ) {
    S7::new_object(
      MNIX_IdealAdaptorStanfit(
        data = data,
        staninput = staninput,
        stanvars = stanvars,
        backend = backend,
        save_pars = save_pars,
        stan_args = stan_args,
        stanfit = stanfit,
        transform_information = transform_information,
        criteria = criteria,
        file = file,
        version = version,
        metadata = metadata,
        cache = cache
      )
    )
  }
)

#' @docType class
#' @rdname MVBU_Stanfit
#' @export
MNIX_IdealAdaptor_CategoryWeightedIntegratorStanfit <- S7::new_class(
  "MNIX_IdealAdaptor_CategoryWeightedIntegratorStanfit",
  package = NULL,
  parent = MNIX_IdealAdaptorStanfit,
  constructor = function(
    data = list(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
    transform_information = NULL,
    criteria = list(),
    file = NULL,
    version = NULL,
    metadata = list(),
    cache = list()
  ) {
    S7::new_object(
      MNIX_IdealAdaptorStanfit(
        data = data,
        staninput = staninput,
        stanvars = stanvars,
        backend = backend,
        save_pars = save_pars,
        stan_args = stan_args,
        stanfit = stanfit,
        transform_information = transform_information,
        criteria = criteria,
        file = file,
        version = version,
        metadata = metadata,
        cache = cache
      )
    )
  }
)

#' @docType class
#' @rdname MVBU_Stanfit
#' @export
MNIX_IdealAdaptor_GroupWeightedIntegratorStanfit <- S7::new_class(
  "MNIX_IdealAdaptor_GroupWeightedIntegratorStanfit",
  package = NULL,
  parent = MNIX_IdealAdaptorStanfit,
  constructor = function(
    data = list(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
    transform_information = NULL,
    criteria = list(),
    file = NULL,
    version = NULL,
    metadata = list(),
    cache = list()
  ) {
    S7::new_object(
      MNIX_IdealAdaptorStanfit(
        data = data,
        staninput = staninput,
        stanvars = stanvars,
        backend = backend,
        save_pars = save_pars,
        stan_args = stan_args,
        stanfit = stanfit,
        transform_information = transform_information,
        criteria = criteria,
        file = file,
        version = version,
        metadata = metadata,
        cache = cache
      )
    )
  }
)

#' @docType class
#' @rdname MVBU_Stanfit
#' @export
MNIX_IdealAdaptor_GroupCategoryWeightedIntegratorStanfit <- S7::new_class(
  "MNIX_IdealAdaptor_GroupCategoryWeightedIntegratorStanfit",
  package = NULL,
  parent = MNIX_IdealAdaptorStanfit,
  constructor = function(
    data = list(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
    transform_information = NULL,
    criteria = list(),
    file = NULL,
    version = NULL,
    metadata = list(),
    cache = list()
  ) {
    S7::new_object(
      MNIX_IdealAdaptorStanfit(
        data = data,
        staninput = staninput,
        stanvars = stanvars,
        backend = backend,
        save_pars = save_pars,
        stan_args = stan_args,
        stanfit = stanfit,
        transform_information = transform_information,
        criteria = criteria,
        file = file,
        version = version,
        metadata = metadata,
        cache = cache
      )
    )
  }
)

#' @docType class
#' @rdname MVBU_Stanfit
#' @export
MNIX_IdealAdaptor_InformativityWeightedIntegratorStanfit <- S7::new_class(
  "MNIX_IdealAdaptor_InformativityWeightedIntegratorStanfit",
  package = NULL,
  parent = MNIX_IdealAdaptorStanfit,
  constructor = function(
    data = list(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
    transform_information = NULL,
    criteria = list(),
    file = NULL,
    version = NULL,
    metadata = list(),
    cache = list()
  ) {
    S7::new_object(
      MNIX_IdealAdaptorStanfit(
        data = data,
        staninput = staninput,
        stanvars = stanvars,
        backend = backend,
        save_pars = save_pars,
        stan_args = stan_args,
        stanfit = stanfit,
        transform_information = transform_information,
        criteria = criteria,
        file = file,
        version = version,
        metadata = metadata,
        cache = cache
      )
    )
  }
)

NIW_IdealAdaptorStanfit <- S7::new_class(
  "NIW_IdealAdaptorStanfit",
  package = NULL,
  parent = IdealAdaptorStanfit,
  constructor = function(
    data = list(),
    staninput = NULL,
    stanvars = NULL,
    backend = "rstan",
    save_pars = NULL,
    stan_args = list(),
    stanfit = NULL,
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
