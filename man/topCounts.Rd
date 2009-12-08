\name{topCounts}
\alias{topCounts}
\title{Get the top counts corresponding to some group from a
  'countData' object}
\description{
  Takes posterior likelihoods and returns the counts with highest (or
  lowest) likelihood of association with a given group.
}
\usage{
topCounts(cD, group, decreasing = TRUE, number = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cD}{\code{\link{countData}} object, containing
    posterior likelihoods for each group.}
  \item{group}{Which group should we give the counts for? See Details.}
  \item{decreasing}{Ordering on posterior likelihoods.}
  \item{number}{How many results should be returned?}
}

\value{
  A dataframe of the top counts associated with some model (group),
  described by annotation drawn from the '@annotation' slot of the
  'cD' object and the raw data from the '@data' slot, together with the
  posterior log-likelihoods.

  If group = NULL, then the function looks at the posterior likelihoods
  that the data have no true differential expression (if calculated).
}
\author{Thomas J. Hardcastle}
\seealso{\code{\link{countData}}}

\examples{

data(simCount)
data(libsizes)

# Make 'countData' object and calculate posterior likelihoods for each
# item belonging to each hypothesis.
groups <- list(c(1,1,1,1,1,1,1,1,1,1), c(1,1,1,1,1,2,2,2,2,2))
CD <- new("countData", data = simCount, libsizes = libsizes, groups = groups)
CDP.Poi <- getPriors.Pois(CD, samplesize = 20, cl = NULL)
CDPost.Poi <- getLikelihoods.Pois(CDP.Poi, prs = c(0.5, 0.5), estimatePriors = TRUE, cl = NULL)

# Report the top ten rows of data that have highest (log) likelihood of belonging to
# group 2 of the data (i.e., differentially expressed)

topCounts(CDPost.Poi, group = 2, number = 10)

}

\keyword{print}

