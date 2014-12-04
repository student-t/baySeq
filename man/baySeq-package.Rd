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
\references{Hardcastle T.J., and Kelly, K. baySeq: Empirical Bayesian
  Methods For Identifying Differential Expression In Sequence Count Data. BMC Bioinformatics (2010)}


\keyword{ package }
\examples{

# See vignette for more examples.

# load test data
data(simData)


# replicate structure of data
replicates <- c("simA", "simA", "simA", "simA", "simA", "simB", "simB", "simB", "simB", "simB")

# define hypotheses on data
groups <- list(NDE = c(1,1,1,1,1,1,1,1,1,1), DE = c(1,1,1,1,1,2,2,2,2,2))

# construct 'countData' object
CD <- new("countData", data = simData, replicates = replicates, groups =
groups)

#estimate library sizes for countData object
libsizes(CD) <- getLibsizes(CD)

# estimate prior distributions on 'countData' object using negative binomial
# method. Other methods are available - see getPriors
CDPriors <- getPriors.NB(CD, cl = NULL)

# estimate posterior likelihoods for each row of data belonging to each hypothesis
CDPost <- getLikelihoods(CDPriors, cl = NULL)

# display the rows of data showing greatest association with the second
# hypothesis (differential expression)
topCounts(CDPost, group = "DE", number = 10)

# find true positive selection rate
getTPs(CDPost, group = "DE", TPs = 1:100)[1:100]

}
