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
    should be the group that defines non-differentially expressed
    data. Can be defined either as the number of the element in
    \code{'cD@groups'} or as a string which will be partially matched to
    the names of the \code{'cD@groups'} elements.} 
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

\seealso{\code{\link{getPriors.NB}}, \code{\link{getLikelihoods.NB}}}

\examples{

# We load in a `countData' object containing the estimated priors (see `getPriors').

data(CDPriors)

try(stopCluster(cl))

plotPriors(CDPriors, group = "NDE")
}

\keyword{hplot}
