% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make-objects.R
\name{sample_MVG_data}
\alias{sample_MVG_data}
\alias{sample_MVG_data_from_model}
\title{Sample multivariate Gaussian exposure data.}
\usage{
sample_MVG_data(
  Ns,
  mus,
  Sigmas,
  category.labels = NULL,
  cue.labels = NULL,
  randomize.order = F,
  keep.input_parameters = F
)

sample_MVG_data_from_model(
  Ns,
  model = NULL,
  randomize.order = F,
  keep.input_parameters = F
)
}
\arguments{
\item{Ns}{Integer vector, with each number specifying the number of observations to be drawn from the corresponding
Gaussian.}

\item{mus}{List of mean vectors, each specifying the mean of a multivariate Gaussian.}

\item{category.labels}{Character vector of category names, each specifying the category label of a multivariate Gaussian. If \code{NULL}
(default) then Gaussians will be numbered from 1:N.}

\item{cue.labels}{Character vector of cue names. If \code{NULL} (default) then the cues will be numbered cue1, cue2, ...}

\item{randomize.order}{Should the order of the data be randomized? (default: FALSE) This won't affect the final outcome of
NIW belief updating, but it will change the incremental updates (and thus, for example, visualizations of the update process).}

\item{keep.input_parameters}{Should the parameters handed to this function be included in the output? (default: FALSE)}

\item{model}{\code{\link{MVG}}, \code{\link{MVG_ideal_observer}}, \code{\link{NIW_belief}}, or \code{\link{NIW_ideal_adaptor}} object.}

\item{sigmas}{List of covariance matrices, each specifying the covariance of a multivariate Gaussian.}
}
\value{
A tibble.
}
\description{
Returns a tibble of observations drawn from multivariate Gaussians, with one observation per row. Each row
provides the category label and cue values. If \code{keep.input_parameters = T} then the parameters (\code{N, mean, sigma})
are also returned.
}
\details{
The input is expected to be lists/vectors of parameters with the n-th element of each list/vector specifying the
category label, number of observations, \code{mu}, and \code{Sigma} of the n-th Gaussian.
}
\seealso{
TBD
}
\keyword{TBD}
