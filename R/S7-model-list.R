# =============================================================================
# S7 Container Class for Lists of Cognitive Models / Posteriors
# =============================================================================

#' @include S7-generics.R
#' @include S7-class.R
#' @include S7-stanfit-posterior.R
#' @importFrom S7 new_class new_object method class_character class_list S7_inherits
NULL

#' S7 class for ordered collections of cognitive models
#'
#' Represents an ordered collection or list of cognitive model objects (which may
#' represent sequential update histories, parameter variations, or distinct model variants).
#' @param models A list of S7 cognitive model objects or `MVBU_StanfitPosterior` objects. Default: \code{list()}.
#' @param model_labels Character vector of labels corresponding to each model. Default: \code{character(0)}.
#' @param metadata Optional list of metadata attributes. Default: \code{list()}.
#'
#' @name MVBU_ModelList
#' @rdname MVBU_ModelList
#' @docType class
#' @export
MVBU_ModelList <- S7::new_class(
  "MVBU_ModelList",
  package = NULL,
  parent = MVBU_Object,
  properties = list(
    models = S7::class_list,
    model_labels = S7::class_character,
    metadata = S7::class_list
  ),
  validator = function(self) {
    if (length(self@models) == 0L) {
      return("`models` must contain at least one model object.")
    }
    for (i in seq_along(self@models)) {
      m <- self@models[[i]]
      if (!S7::S7_inherits(m, MVBU_Object)) {
        return(sprintf("Element %d in `models` is not an S7 MVBU_Object.", i))
      }
    }
    NULL
  }
)

#' Helper constructor to create an MVBU_ModelList
#'
#' @param models A list of S7 cognitive model objects or `MVBU_StanfitPosterior` objects.
#' @param model_labels Optional character vector of names or step labels for each model.
#' @param metadata Optional named list of general attributes.
#' @return An [MVBU_ModelList] object.
#' @export
as_model_list <- function(models, model_labels = NULL, metadata = list()) {
  if (S7::S7_inherits(models, MVBU_ModelList)) {
    return(models)
  }
  if (!is.list(models)) {
    models <- list(models)
  }
  if (is.null(model_labels)) {
    model_labels <- names(models)
    if (is.null(model_labels) || any(!nzchar(model_labels))) {
      model_labels <- paste0("Step_", seq_along(models) - 1L)
    }
  } else {
    model_labels <- as.character(model_labels)
  }

  MVBU_ModelList(
    models = models,
    model_labels = model_labels,
    metadata = as.list(metadata)
  )
}

#' Subset an MVBU_ModelList object
#'
#' @param x An [MVBU_ModelList] object.
#' @param i Index or logical vector.
#' @return A subsetted [MVBU_ModelList] object.
#' @export
`[.MVBU_ModelList` <- function(x, i) {
  sub_models <- x@models[i]
  sub_labels <- x@model_labels[i]
  as_model_list(sub_models, model_labels = sub_labels, metadata = x@metadata)
}

#' Extract a single model element from an MVBU_ModelList
#'
#' @param x An [MVBU_ModelList] object.
#' @param i Index or scalar element name.
#' @return The extracted model object.
#' @export
`[[.MVBU_ModelList` <- function(x, i) {
  x@models[[i]]
}

#' @rdname get_cue_labels
#' @export
S7::method(get_cue_labels, MVBU_ModelList) <- function(x, indices = NULL, ...) {
  get_cue_labels(x@models[[1]], indices = indices, ...)
}

#' @rdname get_category_labels
#' @export
S7::method(get_category_labels, MVBU_ModelList) <- function(x, indices = NULL, ...) {
  get_category_labels(x@models[[1]], indices = indices, ...)
}

#' @rdname get_metadata
#' @export
S7::method(get_metadata, MVBU_ModelList) <- function(x, ...) {
  x@metadata
}

#' S3 Print method for MVBU_ModelList
#' @param x An [MVBU_ModelList] object.
#' @param ... Additional arguments.
#' @export
print.MVBU_ModelList <- function(x, ...) {
  n <- length(x@models)
  cat(sprintf("MVBeliefUpdatr Model List containing %d model(s):\n", n))
  for (i in seq_len(min(n, 10L))) {
    m <- x@models[[i]]
    lbl <- if (i <= length(x@model_labels)) x@model_labels[i] else paste0("Model_", i)
    cat(sprintf("  [%d] %-15s : %s\n", i, lbl, get_model_family(m)))
  }
  if (n > 10L) {
    cat(sprintf("  ... and %d more models.\n", n - 10L))
  }
  invisible(x)
}

#' S3 Summary method for MVBU_ModelList
#' @param object An [MVBU_ModelList] object.
#' @param ... Additional arguments.
#' @export
summary.MVBU_ModelList <- function(object, ...) {
  print(object)
  n <- length(object@models)
  if (n > 0L) {
    first_m <- object@models[[1L]]
    has_moments <- tryCatch({
      !is.null(get_marginal_mu(first_m))
    }, error = function(e) FALSE)
    if (has_moments) {
      cat("\nModel Moments Summary across Steps:\n")
      rows <- list()
      for (i in seq_along(object@models)) {
        m <- object@models[[i]]
        lbl <- if (i <= length(object@model_labels)) object@model_labels[i] else paste0("Model_", i)
        cats <- tryCatch(get_category_labels(m), error = function(e) character(0))
        for (cat_name in cats) {
          exp_mu <- tryCatch(get_expected_mu(m, categories = cat_name), error = function(e) NA)
          marg_mu <- tryCatch(get_marginal_mu(m, categories = cat_name), error = function(e) NA)
          rows[[length(rows) + 1L]] <- data.frame(
            step = lbl,
            category = cat_name,
            expected_mu = paste(round(as.numeric(exp_mu), 3), collapse = ", "),
            marginal_mu = paste(round(as.numeric(marg_mu), 3), collapse = ", "),
            stringsAsFactors = FALSE
          )
        }
      }
      if (length(rows) > 0L) {
        mom_df <- do.call(rbind, rows)
        print(mom_df, row.names = FALSE)
      }
    }
  }
  invisible(object)
}

S7::method(summary, MVBU_ModelList) <- summary.MVBU_ModelList

# -----------------------------------------------------------------------------
# Model List: add_category_representation method
# -----------------------------------------------------------------------------

# Add a category representation to each cognitive model in the model list.
# Returns an updated MVBU_ModelList with the representation added to each model.
S7::method(add_category_representation, MVBU_ModelList) <- function(
  x,
  representation,
  name = NULL,
  category_prior = NULL,
  lapse_bias = NULL,
  ...
) {
  updated_models <- lapply(x@models, function(m) {
    add_category_representation(
      m,
      representation,
      name = name,
      category_prior = category_prior,
      lapse_bias = lapse_bias,
      ...
    )
  })
  MVBU_ModelList(
    models = updated_models,
    model_labels = x@model_labels,
    metadata = x@metadata
  )
}

