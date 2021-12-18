#' Add draws of posterior MCMC samples to a data.frame.
#'
#' DESCRIBE HERE
#' @param fit Stanfit object. (NOT CHECKED YET)
#' @param pars Names of parameters to be extracted or NULL if all parameters are to be returned. (default: NULL)
#' @param summarize Should quantiles be calculated (TRUE) or should all draws be returned (FALSE)?
#' @param ndraws Number of random draws to be return or NULL if all draws are to be returned. (default: NULL)
#' @param draws Vector with specific draw(s) to be returned, or NULL if no specific draws are to be returned. (default: NULL)
#' @param quantiles If summarize is TRUE, quantiles that are calculated; ignored if summarize if FALSE.
#'
#' @return tibble with posterior draws.
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
#'
#' # Draw random sample of parameters from posterior distribution of stanfit and
# return them as a data.frame. If summarize == T, then quantile summaries are
# provided. Otherwise sample.n random samples are returned, with one sample
# per column of the data.frame. Either way, each row will correspond to one
# parameter.
add_draws = function(fit,
                     pars = NULL,
                     summarize = if (is.null(draw) & is.null(draw.n)) TRUE else FALSE,
                     ndraws = NULL,
                     draws = NULL,
                     qi = if (summarize) c(.025, .5, .975) else NA
)
{
  # TO DO: make new object class
  # assert_that(is.mvbeliefs(fit))
  assert_that(is.null(pars) | is.character(pars),
              msg = "Argument pars must be be NULL or a vector of characters.")
  if (is.null(pars)) pars = get_params(fit)

  assert_that(is.null(ndraws) | is.count(ndraws),
              msg = "Argument ndraws must be NULL or a count.")
  assert_that(is.null(draws) |
                (is.numeric(draws) & all(between(draws, 1, get_number_of_draws(fit)))),
              msg = "Argument draws must be NULL or a vector of positive numbers.")
  assert_that(is.flag(summarize))
  assert_that(!all(!summarize, is.null(draws), is.null(ndraws)),
              msg = "If summarize is FALSE, draws or ndraws cannot both be NULL.")
  assert_that(!summarize | (summarize & is.numeric(qi)),
              msg = "If summarize is TRUE, at least one quantile must be specified.")

  if (summarize)
    dat = as.data.frame(rstan::summary(fit, pars = pars, probs = quantiles)$summary[, paste0(quantiles * 100, "%")])
  else {
    dat = as.data.frame(fit, pars = pars)

    if (!is.null(draws)) dat %<>% slice(draws)
    else dat %<>% sample_n(size = ndraws)

    dat %<>%
      t() %>%
      as.data.frame()
  }

  return(dat)
}
