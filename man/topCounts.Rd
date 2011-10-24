\name{topCounts}
\alias{topCounts}
\title{Get the top counts corresponding to some group from a
  'countData' object}
\description{
  Takes posterior likelihoods and returns the counts with highest (or
  lowest) likelihood of association with a given group.
}
\usage{
topCounts(cD, group, decreasing = TRUE, number = 10, normaliseData = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cD}{\code{\link{countData}} object, containing
    posterior likelihoods for each group.}
  \item{group}{Which group should we give the counts for? See Details.}
  \item{decreasing}{Ordering on posterior likelihoods.}
  \item{number}{How many results should be returned?}
  \item{normaliseData}{Should the displayed counts be normalised by library
    size? Defaults to FALSE.}
}

\value{
  A dataframe of the top counts associated with some model (group),
  described by annotation drawn from the '@annotation' slot of the
  'cD' object and the raw data from the '@data' slot, together with the
  posterior likelihoods and false discovery rates.

  The argument 'group' can be specified either as a number, giving the
  index of an element in the \code{cD@groups} list, or as a character
  string identifying an element by name. Partial matching is allowed. If
  group = NULL, then the function looks at the posterior likelihoods
  that the data have no true differential expression (if calculated). 
}
\author{Thomas J. Hardcastle}
\seealso{\code{\link{countData}}}

\examples{

# Create a {countData} object.
data(simData)

replicates <- c("simA", "simA", "simA", "simA", "simA", "simB", "simB", "simB", "simB", "simB")
groups <- list(NDE = c(1,1,1,1,1,1,1,1,1,1), DE = c(1,1,1,1,1,2,2,2,2,2))
CD <- new("countData", data = simData, replicates = replicates, groups = groups)

#estimate library sizes for countData object
CD@libsizes <- getLibsizes(CD)

# If we do not wish to parallelise the functions we set the cluster
# object to NULL.

cl <- NULL

# Alternatively, if we have the 'snow' package installed we
# can parallelise the functions. This will usually (not always) offer
# significant performance gain.

try(library(snow))
try(cl <- makeCluster(4, "SOCK"))

CDPriors <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = cl)
CDPost <- getLikelihoods(CDPriors, cl = cl)

# Report the top ten rows of data that have highest likelihood of belonging to
# group 2 of the data (i.e., differentially expressed)

topCounts(CDPost, group = "DE", number = 10)

# equivalently...
topCounts(CDPost, group = 2, number = 10)

}

\keyword{print}

