#' An S4 class for ideal_adaptor stanfit objects that use one of the ideal_adaptor Stan programs.
#'
#' @name ideal_adaptor_stanfit-class
#' @aliases NIW_ideal_adaptor_stanfit
#' @docType class
#'
#' @details
#' See \code{methods(class = "ideal_adaptor_stanfit")} for an overview of available methods.
#'
#' @slot data A \code{data.frame} containing the data used to fit the model.
#' @slot staninput A named \code{list} containing the data handed to rstan through
#'   \code{\link{make_staninput}}.
#' @slot stanvars A \code{\link{stanvars}} object.
#' @slot backend The name of the backend used to fit the model.
#' @slot stan_args Named list of additional control arguments that were passed
#'   to the Stan backend directly. NOT YET USED
#' @slot stanfit An object of class \code{\link[rstan]{stanfit}}
#'   among others containing the posterior draws.
#' @slot basis An object that contains a small subset of the Stan data
#'   created at fitting time, which is needed to process new data correctly. NOT YET USED
#' @slot transform_information list containing elements transform.parameters, transform.function, and
#'   untransform.function.
#' @slot criteria An empty \code{list} for adding model fit criteria
#'   after estimation of the model. NOT YET USED
#' @slot file Optional name of a file in which the model object was stored in
#'   or loaded from.
#' @slot version The versions of \pkg{MVBeliefUpdatr} and \pkg{rstan} with
#'   which the model was fitted.
#' @slot labels list
setClass(
  "ideal_adaptor_stanfit",
  slots = list(
    data = "data.frame",
    staninput = "list",
    stanvars = "ANY",
    backend = "character",
    save_pars = "ANY",
    stan_args = "list",
    stanfit = "ANY",
    basis = "ANY",
    transform_information = "ANY",
    criteria = "list",
    file = "character",
    version = "ANY"
  )
)

# ideal_adaptor_stanfit class constructor
ideal_adaptor_stanfit <- function(
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

  # Add check here that staninput is a valid staninput object. But don't confuse it with
  # is.ideal_adaptor_staninput, which is the currently confusingly named list of
  # staninput, data, and transform_information
  # assert_that(
  #   is.ideal_adaptor_staninput(staninput),
  #   msg = paste("staninput is not an acceptable input for ideal_adaptor_stanfit stan program."))

  version <- get_current_versions()

  new(
    "ideal_adaptor_stanfit",
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
    version = version
  )
}

methods(class = "ideal_adaptor_stanfit")

setMethod(
  "recover_types", "ideal_adaptor_stanfit",
  function(fit, staninput = NULL) {
    stanfit <- get_stanfit(fit)
    if (is.null(staninput)) staninput <- get_staninput(fit)

    # the levels information recovered below should probably should be stored in a more systematic way, either as attributes
    # to the model or as some part of a list
    stanfit %<>%
      recover_types(
        crossing(
          category = factor(colnames(staninput$z_test_counts), levels = colnames(staninput$z_test_counts)),
          group = factor(attr(staninput$y_test, "levels"), levels = attr(staninput$y_test, "levels")),
          cue = factor(dimnames(staninput$x_test)[[2]], levels = attr(staninput$x_test, "cues")),
          cue2 = cue))

    fit <- set_stanfit(fit, stanfit)
    return(fit)
  }
)




#' An S3 class for MVG objects
#'
#' @name MVG-class
#' @aliases MVG
#' @docType class
#'
#' @details
#' See \code{methods(class = "MVG")} for an overview of available methods.
#'
#' @slot category A list of category labels
#' @slot mu A list of mean vectors for each category
#' @slot Sigma A list of covariance matrices for each category
#' @slot dim A vector of names of the cue dimensions
#'
NULL

MVG <- function(category, mu, Sigma, dim) {
  x <-
    tibble(
      category = factor(category),
      mu = mu,
      Sigma = Sigma
    )

  assert_that(is.MVG(x))

  # For S4 class in the future
  # setClass(
  #   "MVG",
  #   contains = "tibble",
  #   package = "MVBeliefUpdatr")

  # To make S3 class, which however will require a lot of changes
  class(x) <- "MVB"

  x
}
