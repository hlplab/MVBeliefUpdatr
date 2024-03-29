% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/override.R
\name{rmultinom}
\alias{rmultinom}
\title{Overrides stats::rmultinom, reformatting its output. Specifically, the output of stats:rmultinom is transpased, so
that the returned value is a matrix in which each row (rather than column) corresponds to one observation and each
column (rather than row) corresponds to one of the categorical outcomes. Each cell represents the counts observed
for each outcomes on that observation.}
\usage{
rmultinom(n, size, prob)
}
\description{
Overrides stats::rmultinom, reformatting its output. Specifically, the output of stats:rmultinom is transpased, so
that the returned value is a matrix in which each row (rather than column) corresponds to one observation and each
column (rather than row) corresponds to one of the categorical outcomes. Each cell represents the counts observed
for each outcomes on that observation.
}
