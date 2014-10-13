\name{makeOrderings}
\alias{makeOrderings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Construct orderings for count data given a model structure and an
  ordering function.
}
\description{
  Given a model structure as defined in the `@groups' slot of a
  \code{\link{countData}} object containing more than one group, it is
  often possible to define an ordering on the groups for a given genomic
  event. To take a simple example, if the average expression of a gene
  is higher in sample set A then in sample set B, then we might impose
  an ordering A>B.
}
\usage{
makeOrderings(cD, orderingFunction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cD}{
    A \code{\link{countData}} object, or a descendant thereof.
  }
  \item{orderingFunction}{
    A function defining the orderings. If not given, will be taken from
    the `@densityFunction' slot of `cD'. See Details, and examples below.
}
}
\details{
  The orderingFunction takes 'dat' and 'observables' as arguments. 'dat'
  is equivalent to the `@data' slot of the `cD' object, and
  'observables' the combined data in the `@sampleObservables',
  `@rowObservables' and `@cellObservables' slots.
}

\value{
  A \code{\link{countData}} with populated `@orderings' slot.
}
\references{Hardcastle T.J., and Kelly, K. baySeq: Empirical Bayesian
  Methods For Identifying Differential Expression In Sequence Count
  Data. BMC Bioinformatics (2010)}
\author{Thomas J. Hardcastle}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# load test data
data(simData)

# Create a countData object from test data

replicates <- c("simA", "simA", "simA", "simA", "simA", "simB", "simB", "simB", "simB", "simB")
groups <- list(NDE = c(1,1,1,1,1,1,1,1,1,1), DE = c(1,1,1,1,1,2,2,2,2,2))
CD <- new("countData", data = simData, replicates = replicates, groups = groups)
libsizes(CD) <- getLibsizes(CD)

# order on expression normalised for library size (scaling factor) and gene length

CD <- makeOrderings(CD, orderingFunction = function(dat, observables) dat / observables$libsizes)

# orderings calculated for DE group
head(CD@orderings)

# load test (paired) data
data(pairData)  

# create a countData object from paired data
pairCD <- new("countData", data = list(pairData[,1:4], pairData[,5:8]),
                  replicates = c(1,1,2,2),
                  groups = list(NDE = c(1,1,1,1), DE = c(1,1,2,2)),
               densityFunction = bbDensity)
libsizes(pairCD) <- getLibsizes(pairCD)

# order on (log-)ratio of pairs, with fudge-factor on zeros.

pairCD <- makeOrderings(pairCD, orderingFunction = function(dat, observables) {
                                      data <- dat / observables$libsizes               
                                      adjmin <- min(data[data > 0]) / 10
                                      log(data[,,1] + adjmin) - log(data[,,2] + adjmin)
})

# orderings calculated for DE group
head(pairCD@orderings)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
