#' @import rstantools
#' @import methods
#' @import assertthat
#' @importFrom stats cov2cor density dmultinom dnorm plogis prcomp predict qlogis quantile rbinom rnorm runif sd var
#' @importFrom utils data globalVariables
#' @importFrom Rdpack reprompt
#' @importFrom magrittr %<>% %T>%
#' @importFrom Hmisc %nin%
#' @importFrom rlang !! !!! .data .env is_symbol sym syms expr as_name quo_is_null is_missing
#' @importFrom purrr map map2 pmap reduce
#' @importFrom dplyr %>% select filter mutate mutate_at summarise summarise_at left_join rename rename_at group_by ungroup between case_when
#' @importFrom tidyr crossing nest unnest
#' @importFrom tidyselect starts_with
#' @importFrom tibble tibble is_tibble
#' @importFrom rstan sampling
#' @importFrom LaplacesDemon is.positive.definite
#' @useDynLib MVBeliefUpdatr, .registration=TRUE

utils::globalVariables(".")

#' @section Overview:
#' This package provides convenience functions to model Bayesian ideal observers with multivariate Gaussian categories and
#' incremental Bayesian belief-updating for multivariate Gaussian categories. This includes conjugate belief-updating from
#' a Normal-Inverse-Wishart (NIW) prior based on exposure
#' data. Users can specify priors manually or based on existing data, prepare exposure data, update NIW beliefs under a
#' variety of assumptions (noise-free, noise added, etc.) both for labeled and unlabeled exposure data. Expected categories,
#' categorization functions, and categorizations under various decisions rules (e.g., criterion, proportional matching,
#' sampling) can be obtained and visualized.
#'
#' Additionally, the package provides a number of Stan programs that try to \emph{infer} NIW priors for multiple categories
#' from behavioral test responses. These functions use participants' categorization responses during test and, for example,
#' the sufficient statistics of the exposure data to infer a posterior distribution of the parameters of NIW priors for each
#' category. Users can either infer the strength of prior beliefs given a specified m and S parameter, or infer all four NIW
#' parameters together, although the latter requires test responses from multiple different exposure conditions.
#' Functions are provided to interact with stan through rstan, to (1) prepare data as input for the stan program
#' that implements the multivariate Bayesian belief-updating, and to (2) to summarize and visualize the prior and posterior
#' beliefs represented by the resulting fit.
#'
#' @section Basic class structure:
#' The package defines a number of new classes that are essentially tibbles with certain information.
#' \itemize{
#'   \item{\code{MVG}: }{one or more multivariate Gaussian categories, by default in long format with one category per row.
#'   Each row contains the mean mu and covariance matrix Sigma of the multivariate Gaussian.}
#'   \item{\code{MVG_ideal_observer}: }{an ideal observer with multivariate Gaussian categories, by default in long format
#'   with one category per row. In addition to the Gaussian categories the ideal observer contains the prior probability of
#'   each category and, optionally, a lapse rate, lapse bias, and/or perceptual noise matrix.}
#'   \item{\code{NIW_belief}: }{one or more Normal-Inverse-Wishart beliefs, by default in long format with one belief per row.
#'   A Normal-Inverse-Wishart belief specifies *uncertainty* about the about the location (i.e., mean mu) and shape (i.e.,
#'   covariance matrix Sigma) of a multivariate Gaussian. It does so in a specific way that makes assumptions about the
#'   way that the covariance of cues within a category relates to the covariance of the category means across contexts (e.g.,
#'   talkers). See the documentation for details.}
#'   \item{\code{NIW_ideal_adaptor}: }{an ideal adaptor with Normal-Inverse-Wishart beliefs, by default in long format with
#'   one row each for each belief. In addition to the Normal-Inverse-Wishart beliefs, the ideal adaptor contains the prior
#'   probability of each category (currently without uncertainty about those prior probabilities) and, optionally, a lapse
#'   rate, lapse bias, and/or perceptual noise matrix.}
#'   \item{\code{NIW_ideal_adaptor_MCMC}: }{a collection of MCMC samples, each of which constitutes an NIW ideal adaptor. In
#'   other words, this object describes uncertainty about the parameters of the NIW ideal adaptor (specifically, in the
#'   current implementation about the NIW beliefs and the lapse rate and lapse bias, but not yet about the category priors). This is
#'   used, for example, to represent the researchers uncertainty about the prior or posterior beliefs of an ideal adaptor.}
#'   \item{\code{NIW_ideal_adaptor_stanfit}: }{The stanfit resulting from inferring an \code{NIW_ideal_adaptor} from a collection
#'   of exposure and test data. This object contains an \code{NIW_ideal_adaptor_MCMC} object.}
#' }
#'
#' @section Acknowledgments:
#' The belief-updating formulas are taken from Murphy (2012). The package incorporates code
#' from Dave Kleinschmidt's BeliefUpdatr (Kleinschmidt and Jaeger, 2011, 2012, 2015, 2016) and Shaorong Yan's modeling of
#' unsupervised adaptation (Yan and Jaeger, 2018). Pull requests and suggestions from Zach Burchill, Anna Persson, and Xin Xie
#' are gratefully acknowledged.
#'
#' @section Package options: TBD
#'
#' @keywords internal
#' @references TBD
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

# This should be at the end of this file:
