% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binning.R
\name{rand_rect_bin_1d}
\alias{rand_rect_bin_1d}
\title{Random 1d Rectangular Binning}
\usage{
rand_rect_bin_1d(xs, origin, width, output = "centers")
}
\arguments{
\item{xs}{Vector of numeric values to bin.}

\item{origin}{numeric value at which to begin standard rectangular binning. This should be selected at or below minimum value in \code{xs}.}

\item{width}{Positive numeric value at which to space standard rectangular bins.}

\item{output}{Parameter reserved for expanding output types in future package versions.}
}
\value{
bin centers as a vector corresponding to \code{xs} vector
}
\description{
Establish bounding bin centers for each observation
Then use Unif(0,1) draws compared to assignment probs to allocate
Reassignment for values below first center or above highest center
}
\examples{
rand_rect_bin_1d(1:10, origin=.5, width=2)
}

