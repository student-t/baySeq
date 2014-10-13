\name{allModels}
\alias{allModels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to generate all possible models for a countData object based
  on the replicate data.
}
\description{
  This function populates the `@groups' slot of the supplied countData
  object with all possible models for equivalence/non-equivalence of
  expression between replicate groups.
}
\usage{
allModels(CD)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CD}{
    A countData object with a populated `@replicates' slot.
}
}
\details{
  Given a large number of different replicate groups, the total number
  of possible models listed in the `@groups' slot rises
  exponentially. This function will attempt to list them all. The use of
  consensus priors (see \code{\link{getPriors}}) is recommended if the
    number of models is high. 
}
\value{
  A \code{\link{countData}} with populated `@groups' slot.
}
\references{Hardcastle T.J., and Kelly, K. baySeq: Empirical Bayesian
  Methods For Identifying Differential Expression In Sequence Count
  Data. BMC Bioinformatics (2010)}
\author{Thomas J. Hardcastle}

\seealso{
  \code{\link{getPriors}}
}
\examples{
# load test data
data(simData)

# Create a {countData} object from test data, supposing that there are
# multiple experimental groups present.

replicates <- c("simA", "simA", "simB", "simC", "simC", "simD", "simE", "simE", "simF", "simG")
CD <- new("countData", data = simData, replicates = replicates)
CD <- allModels(CD)

# The total number of models generated is high.
length(CD@groups)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

