



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
