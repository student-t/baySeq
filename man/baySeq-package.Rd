\name{baySeq-package}
\alias{baySeq-package}
\alias{baySeq}
\docType{package}
\title{Empirical Bayesian analysis of patterns of differential expression in count data.}
\description{
This package is intended to identify differential expression in
high-throughput 'count' data, such as that derived from next-generation
sequencing machines. We achieve this by empirical bayesian methods,
first bootstrapping to estimate prior parameters from the
data and then assessing posterior likelihoods of the models proposed.
}
\details{
\tabular{ll}{
Package: \tab baySeq\cr
Type: \tab Package\cr
Version: \tab 1.1.1\cr
Date: \tab 2009-16-05\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}
To use the package, construct a \code{\link{countData}} object and use
the functions documented in \link{getPriors} to empirically determine
priors on the data. Then use the functions documented in
\link{getLikelihoods} to establish posterior likelihoods for the models
proposed. A few convenience functions, \code{\link{getTPs}} and
\code{\link{topCounts}} are also included.

The package (optionally) makes use of the 'snow' package for
parallelisation of computationally intensive functions. This is highly
recommended for large data sets.

See the vignette for more details.
}
\author{
Thomas J. Hardcastle

Maintainer: Thomas J. Hardcastle <tjh48@cam.ac.uk>
}
\references{Hardcastle T.J., and Kelly, K (2009). Empirical Bayesian
  methods for differential expression in count data. In submission.}

\keyword{ package }
\examples{

# See vignette for more examples.

# load test data
data(simCount)
data(libsizes)

# define hypotheses on data
groups <- list(c(1,1,1,1,1,1,1,1,1,1), c(1,1,1,1,1,2,2,2,2,2))

# construct 'countData' object
CD <- new("countData", data = simCount, libsizes = libsizes, groups = groups)
CD[1:10,]

# estimate prior distributions on 'countData' object using Poisson
# method. Other methods are available - see getPriors
CDP.Poi <- getPriors.Pois(CD, samplesize = 20, takemean = TRUE, cl = NULL)

# estimate posterior likelihoods for each row of data belonging to each hypothesis
CDPost.Poi <- getLikelihoods.Pois(CDP.Poi, prs = c(0.5, 0.5), estimatePriors = TRUE, cl = NULL)

# display the rows of data showing greatest association with the second
# hypothesis (differential expression)
topCounts(CDPost.Poi, group = 2, number = 10)

# find true positive selection rate
getTPs(CDPost.Poi, group = 2, TPs = 1:100)[1:100]

}
