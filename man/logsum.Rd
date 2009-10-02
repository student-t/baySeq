\name{logsum}
\alias{logsum}
\title{An internal function for the baySeq package that calculates the log of the sum of
  the exponential of two variables}
\description{
  If x = log(y), where y is a vector, this function returns log(sum(y))
  given x.
}
\usage{
logsum(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{Vector of log-values.}
}

\value{
  Numeric.
}
\author{Thomas J. Hardcastle}

\seealso{\code{\link{countDataPosterior}}}

\examples{

x <- logsum(log(c(3, 4, 5)))
exp(x) == sum(c(3,4,5))

}
\keyword{arith}
