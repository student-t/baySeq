\name{marginalisePairwise}
\alias{marginalisePairwise}
\title{
Computes marginal likelihoods that two replicate groups differ, in a
paticular direction.
}
\description{
In cases where multiple models are simultaneously evaluated in the
'getLikelihoods' function, the posterior likelihoods for each model in
which one condition is greater than another  can be summed to give the
marginal likelihood of (directed) difference for all biomolecular events
(i.e., data rows). 
}
\usage{
marginalisePairwise(cD, greaterThan, lessThan)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cD}{
    A \code{\link{countData}} object with evaluated posterior
    likelihoods in the '@posteriors' slot.
}
  \item{greaterThan}{
    A defined group name (or vector of group names) to identify in the
    '@groups' slot of the \code{\link{countData}} object 'cD'; the
    function will identify all models in which these groups are
    equivalent and greater than that defined in the 'lessThan' variable.
}
  \item{lessThan}{
    A defined group name (or vector of group names) to identify in the
    '@groups' slot of the \code{\link{countData}} object 'cD'; the
    function will identify all models in which these groups are
    eqivalent and less than that defined in the 'greaterThan' variable.
}
}

\value{
  A vector of marginal posterior likelihoods defining the probability
  that the two group identifiers are (directionally) different for each row of the data.
}
\author{
Thomas J. Hardcastle
}
\seealso{
  \code{\link{allModels}}
  \code{\link{marginaliseEqual}}
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

# Priors and likelihoods acquired through standard means.

\dontrun{CD <- getPriors(CD, cl = cl)}
\dontrun{CD <- getLikelihoods(CD, cl = cl)}

# Marginal likelihood that 'simA' condition is greater than 'simD' group
# for each row of the data.

\dontrun{marginalisePairwise(CD, "simA", "simD")}

}

\keyword{manip}
