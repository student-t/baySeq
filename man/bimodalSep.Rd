\name{bimodalSep}
\alias{bimodalSep}
\title{A function that, given a numeric vector, finds the value which
  splits the data into two sets of minimal total variance.}
\description{This function takes a numeric vector and finds the value
  which splits the data into two sets of minimal total variance. It is
  principally intended to be a quick and easy way of separating
  bimodally distributed data.}
\usage{
bimodalSep(z, weights = NULL, bQ = c(0,1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{A numeric vector containing the data to be split.}
  \item{weights}{Possible weightings on the values in z for calculating
    the variance.}
  \item{bQ}{Lower and upper limits on the quartile of z in which a
    separating value is sought. See Details.}
}

\details{This function is intended to give a quick and easy way of
  splitting bimodally distributed data. Where there are large outliers
  in the data, it may be that the value which minimises the variance
  does not split the bimodal data but isolates the outliers. The
  \code{'bQ'} parameter can be used to ensure that the split occurs within
  some range of quantiles of the data.
}
\value{Numeric.
}

\author{Thomas J. Hardcastle}

\examples{
bimodalSep(c(rnorm(200, mean = c(5,7), sd = 1)))
}

\keyword{models}
