% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-info-from-MVG-IO.R
\name{get_MVG_likelihood}
\alias{get_MVG_likelihood}
\alias{get_likelihood_from_MVG}
\title{Get likelihood}
\usage{
get_MVG_likelihood(
  x,
  mu,
  Sigma,
  Sigma_noise = NULL,
  noise_treatment = infer_default_noise_treatment(Sigma_noise),
  log = T
)

get_likelihood_from_MVG(
  x,
  model,
  noise_treatment = infer_default_noise_treatment(model$Sigma_noise),
  log = T,
  category = "category",
  category.label = NULL,
  wide = FALSE
)
}
\arguments{
\item{x}{Observations. Can be a vector with k elements for a single observation or a matrix with k
columns and n rows, in which case each row of the matrix is taken to be one observation. If x is a
tibble with k columns or a list of vectors of length k, it is reduced into a matrix with k columns.}

\item{mu}{The category mean mu. Should be a k x 1 or 1 x k
matrix, or vector of length k.}

\item{Sigma}{The category covariance matrix Sigma. Should be a square k x k matrix.}

\item{Sigma_noise}{Optionally, a covariance matrix describing the perceptual noise to be applied while
calculating the posterior predictive. (default: `NULL`)}

\item{noise_treatment}{Determines whether perceptual noise is considered during categorization, and how.
Can be "no_noise", "sample", or "marginalize". If "no_noise", no noise will be applied to the input,
and no noise will be assumed during categorization. If "marginalize", average noise (i.e., no noise)
will be added to the stimulus, and `Sigma_noise` is added to Sigma when calculating the likelihood.
This simulates the expected consequences for perceptual noise on categorization *in the limit*, i.e,
if the input was categorized infinitely many times. If "sample", then noise is sampled and applied to
the input, and `Sigma_noise` is added to Sigma when calculating the likelihood. This simulates the
consequence of perceptual noise *on a particular observation*. If "sample" or "marginalize" are chosen,
`Sigma_noise` must be a covariance matrix of appropriate dimensions. (default: "no_noise" if Sigma_noise
is NULL, "marginalize" otherwise).}

\item{log}{Should the log-transformed density be returned (`TRUE`)? (default: `TRUE`)}
}
\description{
Get likelihood of observation(s) x given the MVG parameters mu and Sigma. This is the density of
a multivariate normal distribution over k dimensions.
}
\seealso{
TBD
}
\keyword{TBD}
