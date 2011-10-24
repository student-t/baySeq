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
getTPs(cD, group, decreasing = TRUE, TPs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{cD}{\code{\link{countData}} object, containing
    posterior likelihoods for each group.}
  \item{group}{Which group should we give the counts for? See Details.}
  \item{decreasing}{Ordering on posterior likelihoods.}
  \item{TPs}{Known true positives.}
}
\details{
  In the rare (or simulated) cases where the true positives are known,
  this function will calculate the number of true positives selected at
  any cutoff.

  The 'group' can be defined either as the number of the element in
  \code{'cD@groups'} or as a string which will be partially matched to
  the names of the \code{'cD@groups'} elements.If group = NULL, then
  the function looks at the posterior likelihoods that the data have
  no true differential expression (if calculated). 
}
\value{
  A vector, the ith member of which gives the number of true positives
  identified if the top i counts are chosen.
}
\author{Thomas J. Hardcastle}

\seealso{\code{\link{countData}}}

\examples{

# See vignette for more examples.

# If we do not wish to parallelise the functions we set the cluster
# object to NULL.

cl <- NULL

# Alternatively, if we have the 'snow' package installed we
# can parallelise the functions. This will usually (not always) offer
# significant performance gain.

try(library(snow))
try(cl <- makeCluster(4, "SOCK"))

# load test data
data(simData)

# Create a {countData} object from test data.

replicates <- c("typeA", "typeA", "typeA", "typeA", "typeA",
                "typeB", "typeB", "typeB", "typeB", "typeB")

groups <- list(NDE = c(1,1,1,1,1,1,1,1,1,1), DE = c(1,1,1,1,1,2,2,2,2,2))
CD <- new("countData", data = simData, replicates = replicates, groups = groups)

#estimate library sizes for countData object
CD@libsizes <- getLibsizes(CD)

# Get priors for negative binomial method
CDPriors <- getPriors.NB(CD, samplesize = 10^5, estimation = "QL", cl = cl)

# Get likelihoods for data with negative binomial method.

CDPost <- getLikelihoods.NBboot(CDPriors, prs = c(0.5, 0.5),
pET = "BIC", marginalise = FALSE, bootStraps = 1, cl = cl)

# If the first hundred rows in the 'simData' matrix are known to be
# truly differentially expressed (the second hypothesis defined in the
# 'groups' list) then we find the number of true positives for the top n
# genes selected as the nth member of

getTPs(CDPost, group = "DE", decreasing = TRUE, TPs = 1:100)

}

\keyword{manip}
