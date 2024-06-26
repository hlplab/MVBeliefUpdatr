% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualize-basics.R
\name{add_exposure_data_to_1D_plot}
\alias{add_exposure_data_to_1D_plot}
\alias{add_test_data_to_1D_plot}
\alias{add_exposure_data_to_2D_plot}
\alias{add_test_data_to_2D_plot}
\alias{add_exposure_summary_to_1D_plot}
\alias{add_exposure_summary_to_2D_plot}
\title{Add exposure/test data to a plot as points or ellipse}
\usage{
add_exposure_data_to_1D_plot(
  data,
  cue.labels,
  category.ids,
  category.labels,
  category.colors
)

add_test_data_to_1D_plot(data, cue.labels)

add_exposure_data_to_2D_plot(
  data,
  cue.labels,
  category.ids,
  category.labels,
  category.colors
)

add_test_data_to_2D_plot(data, cue.labels)

add_exposure_summary_to_1D_plot(data)

add_exposure_summary_to_2D_plot(data, level = 0.95)
}
\arguments{
\item{data}{A `tibble` or `data.frame`. For plotting points, this table must contain columns
named in `cue.labels`. For ellipsis plots, the table must contain the columns `mean` and `cov`.}

\item{cue.labels}{Name of the columns in `data` that contain the cue values of the data to be
plotted.}

\item{level}{Level of probability mass for which ellipsis should be plotted (default: .95)}
}
\value{
A list of geoms.
}
\description{
Add exposure/test data as points or ellipse to a ggplot2. For the 1D functions, the x-axis is
assumed to be a phonetic cue. For the 2D functions, the x- and y-axes are assumed to be
phonetic cues. The difference between adding exposure and test data is whether the data are
colored by the category (exposure) or not (test).
}
\seealso{
TBD
}
\keyword{TBD}
