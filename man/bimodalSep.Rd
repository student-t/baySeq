\name{bimodalSeparator}
\alias{bimodalSeparator}
\title{A function that, given a numeric vector, finds the value which
  splits the data into two sets of minimal total variance using Otsu's
  method.}
\description{This function takes a numeric vector and finds the value
  which splits the data into two sets of minimal total variance,
  weighted by the size of subsets (Otsu's method). It is
  principally intended to be a quick and easy way of separating
  bimodally distributed data.}
\usage{
bimodalSeparator(x, weights = NULL, minperc = 0.1, elbow = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A numeric vector containing the data to be split.}
  \item{weights}{Possible weightings on the values in x for calculating
    the variance.}
  \item{minperc}{The required minimum size of each of the two subsets,
    expressed as a percentage of the total size. See Details.}
  \item{elbow}{If set, finds the 'left' or 'right' elbow of variance, instead of the
    minimum; defaults to NULL. See Details.}
}

\details{This function is intended to give a quick and easy way of
  splitting bimodally distributed data. Where there are large outliers
  in the data, it may be that the value which minimises the variance
  does not split the bimodal data but isolates the outliers. The
  \code{'minperc'} parameter can be used to ensure that each subset of
  the split data will be of some minimum size, avoiding the outlier
  problem.

  If 'elbow = NULL' (the default) then the split occurs at the value
  that minimises the variance, x0. If 'elbow = left' then we attempt to find
  the elbow point to the left of the value that minimises the variance,
  if 'elbow = right' then we find the elbow point to the right of the
  value that minimises the variance. Elbow points are found by drawing a
  line from the first point (for the left elbow) or the last point (for
  the right elbow) to x0, and finding the location on the curve of
  summed variances which maximises the distance to that line.
}
\value{Numeric value which splits the data.}

\author{Thomas J. Hardcastle}

\examples{
bimodalSeparator(c(rnorm(200, mean = c(5,7), sd = 1)))
}

\keyword{models}
