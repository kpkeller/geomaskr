% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_sampling.R
\name{calc_radius}
\alias{calc_radius}
\alias{calc_k}
\title{calc_radius}
\usage{
calc_radius(k_anon, density, units_scale = 1)

calc_k(radius, r0 = 0, density)
}
\arguments{
\item{k_anon}{K-anonymity metric}

\item{density}{Population density}

\item{units_scale}{Scaling factor for distances; defaults to 1. See details.}

\item{radius}{Radius for calculating k_anonymity}

\item{r0}{Minimum radius}
}
\description{
Calculates radius for a specified k-anonymity and population density.
}
\details{
The scaling factor is applied to the radius at the end of the calculation. For example, if density is in persons/km2, then setting \code{units_scale=1000} in \code{calc_radius} would result in the necessary radius in meters.
}
\seealso{
\code{\link{st_sample_radius_bounded}}
}
