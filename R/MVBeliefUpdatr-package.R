#' @import rstantools
#' @import methods
#' @importFrom stats cov2cor density dmultinom dnorm plogis prcomp predict
#'   qlogis quantile rbinom rnorm runif sd var
#' @importFrom utils data globalVariables
#' @importFrom Rdpack reprompt
#' @importFrom magrittr %<>% %T>%
#' @importFrom Hmisc %nin%
#' @importFrom rlang !! !!! .data .env is_symbol sym syms expr as_name
#'   quo_is_null is_missing
#' @importFrom purrr map map2 pmap reduce
#' @importFrom dplyr %>% select filter mutate mutate_at summarise summarise_at
#'   left_join rename rename_at group_by ungroup between case_when pull
#' @importFrom tidyr complete crossing drop_na nest replace_na unnest
#' @importFrom tidyselect starts_with
#' @importFrom tibble tibble is_tibble
#' @importFrom rstan sampling
#' @importFrom LaplacesDemon is.positive.definite
#' @useDynLib MVBeliefUpdatr, .registration=TRUE
NULL

utils::globalVariables(".")

#' @section Overview:
#' `MVBeliefUpdatr` provides a unified, object-oriented framework for
#' Bayesian ideal observers, ideal adaptors, incremental belief updating,
#' and hierarchical Bayesian cognitive modeling.
#'
#' The package supports continuous distributional representations (univariate
#' Gaussian, multivariate Gaussian, multi-univariate Gaussian cue integration),
#' conjugate prior distributions over category parameters (Normal-Inverse-
#' Chisquare, Normal-Inverse-Wishart), and exemplar-based representations.
#'
#' Using analytical conjugate updating or hierarchical Bayesian inference
#' in Stan via \pkg{rstan}, `MVBeliefUpdatr` facilitates simulation, parameter
#' estimation, categorization prediction, and visualization of speech
#' perception and perceptual adaptation experiments.
#' @section Package Vignettes:
#' For detailed guides and worked examples, see the package vignettes:
#' \describe{
#'   \item{`vignette("s7-class-structure-and-workflows")`}{S7 Class Architecture,
#'     Constructors, and Workflows}
#'   \item{`vignette("visualizing-models-and-categories")`}{Visualizing Models and
#'     Categories}
#'   \item{`vignette("fitting-and-working-with-stanfit-models")`}{Fitting and
#'     Working with MVBeliefUpdatr Stanfit Models}
#' }
#'
#' @section S7 Core Class Architecture:
#' All model objects are organized under an explicit, compositional S7 class
#' hierarchy:
#' \describe{
#'   \item{`MVBU_CategoryRepresentation`}{Abstract base class for individual
#'     category representations, encapsulating category parameters, cue
#'     labels, and likelihood functions.}
#'   \item{`MVBU_CategoryRepresentationTemplate`}{Validated collection of
#'     category representations forming a perceptual inventory or template.}
#'   \item{`MVBU_CognitiveModel`}{Decision model pairing a category template
#'     with decision behavior (decision rule, category priors, lapse rate,
#'     lapse bias, and perceptual noise).}
#'   \item{`MVBU_IdealObserver`}{Subclass of `MVBU_CognitiveModel` for
#'     observers with fixed generative category distributions.}
#'   \item{`MVBU_IdealAdaptor`}{Subclass of `MVBU_CognitiveModel` for adaptors
#'     with conjugate prior beliefs that update incrementally from exposure.}
#'   \item{`MVBU_Stanfit`}{Wrapper for `rstan` `stanfit` objects inferring
#'     prior beliefs, lapse parameters, and noise from behavioral test data.}
#'   \item{`MVBU_StanfitInput`}{Container for exposure data, test data,
#'     priors, and transformation functions used to prepare Stan inputs.}
#'   \item{`MVBU_Staninput`}{Validated input data list passed to Stan models.}
#' }
#'
#' @section Supported Model Families:
#' `MVBeliefUpdatr` provides concrete S7 classes and constructors across seven
#' families:
#' \describe{
#'   \item{`UVG` / `NIX`}{Univariate Gaussian ideal observers
#'     ([UVG_IdealObserver]) and Normal-Inverse-Chisquare ideal adaptors
#'     ([NIX_IdealAdaptor]).}
#'   \item{`MUVG` / `MNIX`}{Multi-univariate Gaussian ideal observers
#'     ([MUVG_IdealObserver]) and independent Normal-Inverse-Chisquare cue
#'     integration adaptors ([MNIX_IdealAdaptor]).}
#'   \item{`MVG` / `NIW`}{Multivariate Gaussian ideal observers
#'     ([MVG_IdealObserver]) and Normal-Inverse-Wishart ideal adaptors
#'     ([NIW_IdealAdaptor]).}
#'   \item{`EXEMPLAR`}{Exemplar-based categorization models ([Exemplar_Model])
#'     with kernel density estimation over stored exemplars.}
#' }
#'
#' @section Acknowledgments:
#' Belief-updating formulations build upon conjugate Bayesian theory
#' \insertCite{murphy2012}{MVBeliefUpdatr}, phonetic sliding template models
#' \insertCite{nearey-assmann2007}{MVBeliefUpdatr}, Dave Kleinschmidt's
#' `BeliefUpdatr` \insertCite{kleinschmidt-jaeger2011}{MVBeliefUpdatr},
#' \insertCite{kleinschmidt-jaeger2012}{MVBeliefUpdatr},
#' \insertCite{kleinschmidt-jaeger2015}{MVBeliefUpdatr},
#' \insertCite{kleinschmidt-jaeger2016cogsci}{MVBeliefUpdatr},
#' and unsupervised adaptation modeling
#' \insertCite{yan:jaeger2018}{MVBeliefUpdatr}.
#'
#' Contributions and feedback from
#' Zach Burchill, Anna Persson, and Xin Xie are gratefully acknowledged.
#'
#' @references
#' \insertAllCited{}
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

get_current_versions <- function() {
  list(
    MVBeliefUpdatr = utils::packageVersion("MVBeliefUpdatr"),
    rstan = utils::packageVersion("rstan"),
    stanHeaders = utils::packageVersion("StanHeaders")
  )
}
