#' @import assertthat
#' @importFrom mvtnorm rmvnorm
#' @importFrom purrr pmap
NULL

example_prior1 = function() {
  tibble(
    category = c("A", "B"),
    kappa = 10,
    nu = 30,
    M = list(c("cue1" = -2, "cue2" = -2), c("cue1" = 2, "cue2" = 2)),
    S = list(matrix(c(1, .3, .3, 1), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2"))),
             matrix(c(1, -.3, -.3, 1), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2")))),
    lapse = .05
  ) %>%
    mutate(category = factor(category))
}



#' Make multivariate Gaussian exposure data.
#'
#' Returns a tibble of observations drawn from multivariate Gaussians, with one observation per row. Each row
#' provides the category label and cue values. If \code{keep.parameters = T} then the parameters (\code{N, mean, sigma})
#' are also returned.
#'
#' The input is expected to be lists of parameters which the n-th element of each list specifying the n-th Gaussian and the
#' number of observations to be drawn from that Gaussian.
#'
#' @param Ns A list of integers, with each number specifying the number of observations to be drawn from the corresponding
#' Gaussian.
#' @param means List of mean vectors, each specifying the mean of a multivariate Gaussian.
#' @param sigmas List of covariance matrices, each specifying the covariance of a multivariate Gaussian.
#' @param category.labels List of category names, each specifying the category label of a multivariate Gaussian. If \code{NULL}
#' (default) then Gaussians will be numbered from 1:N.
#' @param cue.labels List of cue names. If \code{NULL} (default) then the cues will be numbered cue1, cue2, ...
#' @param keep.parameters Should the parameters handed to this function be included in the output? (default: FALSE)
#'
#' @return A tibble.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
#'
make_MV_exposure_data = function(Ns, means, Sigmas, category.labels = NULL, cue.labels = NULL, keep.parameters = F) {
  assert_that(!is.null(means), !is.null(Sigmas))
  assert_that(is.null(category.labels) | length(means) == length(category.labels),
              msg = "Number of category labels mismatch number of means.")
  assert_that(is.null(cue.labels) | length(means[[1]]) == length(cue.labels),
              msg = "Number of cue labels mismatches dimensionality of means.")

  if (is.null(category.labels)) category.labels = 1:length(means)
  if (is.null(cue.labels)) cue.labels = paste0("cue", 1:length(means[[1]]))

  x = tibble(n = unlist(Ns), mean = means, sigma = Sigmas) %>%
    mutate(data = pmap(.l = list(n, mean, sigma), .f = rmvnorm)) %>%
    mutate(category = unlist(category.labels)) %>%
    mutate(data = map(data, ~ .x %>% as.data.frame() %>% rename_all(~ cue.labels))) %>%
    unnest(data)

  if (keep.parameters) return(x) else return(x %>% select(category, everything(), -c(n, mean, sigma)))
}



#' Update prior beliefs about multivariate Gaussian category based on exposure data and the conjugate NIW prior.
#'
#' Returns updated/posterior beliefs about the Gaussian categories.
#'
#' The priors for the categories are specified through the \code{priors} argument. This is expected to be a tibble
#' of the same format as the posterior draws stored in an MV IBBU stanfit object. Each row of the tibble specifies
#' the prior for one category (specified in the \code{category} column). The four parameters of the NIW are the
#' pseudocounts indicating the strength of the prior beliefs into the mean (\code{kappa}) and covariance (\code{
#' nu}), as well as the prior mean of means (\code{M}, same as \code{M_0} in Murphy, 2012) and the prior scatter
#' matrix (\code{S}, same as \code{S_0} in Murphy, 2012).
#'
#' @param data \code{data.frame} or \code{tibble} with exposure data. Each row is assumed to contain one observation.
#' @param category Name of variable in \code{data} that contains the category information. (default: "category")
#' @param cues Name(s) of variables in \code{data} that contain the cue information.
#' @param priors A tibble with information about the prior. See Details for expected format of the \code{priors} argument.
#'
#' @return A tibble.
#'
#' @seealso TBD
#' @keywords TBD
#' @references Murphy, K. P. (2012). Machine learning: a probabilistic perspective. MIT press.
#' @examples
#' TBD
#' @export
#'
update_MV_beliefs <- function(
  data,
  priors,
  category = "category",
  cues = paste0("cue", 1:length(priors$M[[1]])),
  store.history = T
){
  assert_that(is.NIW_belief(priors),
              msg = "Priors must be NIW beliefs. Check is.NIW_belief().")

  # Number of dimensions/cues
  D = length(cues)
  if (any(priors$nu < D + 2))
    message(paste0("Prior for at least one category had nu smaller than allowed (", D, ").\n"))

  # Prepare data
  data %<>%
    mutate(cues = pmap(.l = list(!!! syms(cues)), .f = ~ c(...))) %>%
    select(category, cues)

  if (store.history)
    priors %<>%
      mutate(observation = 0)

  for (i in 1:nrow(data)) {
    posteriors = if (store.history) priors %>% filter(observation == i - 1) else priors

    current_category_index = which(posteriors$category == data[i,]$category)
    current_observation = unlist(data[i, "cues"])

    # Keep this order, see Murphy 2012, p. 134
    posteriors[current_category_index,]$S[[1]] =
      posteriors[current_category_index,]$S[[1]] +
      # Using centered versions, rather than uncentered sum of squares
      (posteriors[current_category_index,]$kappa / (posteriors[current_category_index,]$kappa + 1)) *
      matrix(current_observation - posteriors[current_category_index,]$M[[1]]) %*%
      t(matrix(current_observation - posteriors[current_category_index,]$M[[1]]))

    posteriors[current_category_index,]$M[[1]] =
      (posteriors[current_category_index,]$kappa / (posteriors[current_category_index,]$kappa + 1)) * posteriors[current_category_index,]$M[[1]] +
      (1 / (posteriors[current_category_index,]$kappa + 1)) * current_observation

    posteriors[current_category_index,]$kappa = posteriors[current_category_index,]$kappa + 1
    posteriors[current_category_index,]$nu = posteriors[current_category_index,]$nu + 1

    if (store.history) {
      posteriors %<>%
        mutate(observation = i)
      priors = rbind(priors, posteriors)
    } else priors = posteriors
  }

  return(priors)
}
