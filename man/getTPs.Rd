\name{getTPs}
\alias{getTPs}
\title{Gets the number of true positives in the top n counts selected by
  ranked posterior likelihoods}
\description{
  If the true positives are known, this function will return a vector,
  the ith member of which gives the number of true positives identified
  if the top i counts, based on estimated posterior likelihoods, are chosen.
}
\usage{
getTPs(cDP, group, decreasing = TRUE, TPs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{cDP}{\code{\link{countDataPosterior}} object, containing
    posterior likelihoods for each group.}
  \item{group}{Which group should we give the counts for?}
  \item{decreasing}{Ordering on posterior likelihoods.}
  \item{TPs}{Known true positives.}
}
\details{
  In the rare (or simulated) cases where the true positives are known,
  this function will calculate the number of true positives selected at
  any cutoff.
}
\value{
  A vector, the ith member of which gives the number of true positives
  identified if the top i counts are chosen.
}
\author{Thomas J. Hardcastle}

\seealso{\code{\link{countDataPosterior}}}

\examples{

# Create a {countData} object and estimate priors for the Poisson methods.
data(simCount)
data(libsizes)
groups <- list(c(1,1,1,1,1,1,1,1,1,1), c(1,1,1,1,1,2,2,2,2,2))
CD <- new("countData", data = simCount, libsizes = libsizes, groups = groups)
CDP.Poi <- getPriors.Pois(CD, samplesize = 20, iterations = 1000,
takemean = TRUE)

# Get likelihoods for data with Poisson method
CDPost.Poi <- getLikelihoods.Pois(CDP.Poi, prs = c(0.5, 0.5),
estimatePriors = TRUE, cl = NULL)

# If the first hundred rows in the 'simCount' matrix are known to be
# truly differentially expressed (the second hypothesis defined in the
# 'groups' list) then we find the number of true positives for the top n
# genes selected as the nth member of

getTPs(CDPost.Poi, group = 2, decreasing = TRUE, TPs = 1:100)

}

\keyword{manip}
