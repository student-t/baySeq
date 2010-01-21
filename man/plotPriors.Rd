\name{plotPriors}
\alias{plotPriors}
\title{Plots the density of the log values estimated for the mean rate
  in the prior data for the Negative Binomial approach to detecting
  differential expression}
\description{
  This function plots the density of the log values estimated for the
  mean rate in the data used to estimate a prior distribution for data
  under the assumption of a Negative Binomial distribution. This
  function is useful for looking for bimodality of the distributions,
  and thus determining whether we should try and identify data with no
  true expression.
}
\usage{
plotPriors(cD, group)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{cD}{\code{\link{countData}} object, for which priors have been
    estimated using the assumption of a Negative Binomial distribution
    (see \code{\link{getPriors.NB}}).}
  \item{group}{Which group should we plot the priors for? In general,
    should be the group that defines non-differentially expressed data.}
}
\details{
  If the plot of the data appears bimodal, then it may be sensible to
  try and look for data with no true expression by using the option
  \code{nullPosts = TRUE} in \code{\link{getLikelihoods.NBboot}}.
}
\value{
  Plotting function.
}
\author{Thomas J. Hardcastle}

\seealso{\code{\link{getPriors.NB}}, \code{\link{getLikelihoods.NBboot}}}

\examples{

# Create a {countData} object and estimate priors for the Poisson methods.
data(simSeg)
data(libsizes)
groups <- list(c(1,1,1,1,1,1,1,1,1,1), c(1,1,1,1,1,2,2,2,2,2))
CD <- new("countData", data = simSeg[,-1], seglens = simSeg[,1], libsizes = libsizes, groups = groups)
CDP.NB <- getPriors.NB(CD, samplesize = 1000, estimation = "ML", cl = NULL)

plotPriors(CDP.NB, group = 1)
}

\keyword{manip}
