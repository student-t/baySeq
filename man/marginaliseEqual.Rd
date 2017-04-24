\name{marginaliseEqual}
\alias{marginaliseEqual}
\title{
Computes marginal likelihoods that two replicate groups are equal.
}
\description{
In cases where multiple models are simultaneously evaluated in the
'getLikelihoods' function, the posterior likelihoods for each model in
which two conditions are equivalent can be summed to give the marginal
likelihood of equivalence for all biomolecular events (i.e., data rows).
}
\usage{
marginaliseEqual(cD, r1, r2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cD}{
    A \code{\link{countData}} object with evaluated posterior
    likelihoods in the '@posteriors' slot.
}
  \item{r1}{
    A defined group name to identify in the '@groups' slot of the
    \code{\link{countData}} object 'cD'.
}
  \item{r2}{
    A defined group name to identify in the '@groups' slot of the
    \code{\link{countData}} object 'cD'.
}
}
\value{
  A vector of marginal posterior likelihoods defining the probability
  that the two group identifiers are equal for each row of the data.
}
\author{
Thomas J. Hardcastle
}
\seealso{
  \code{\link{allModels}}
  \code{\link{marginalisePairwise}}
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

# Marginal likelihood that 'simA' and 'simD' replicate groups are equal
# for each row of the data.

\dontrun{marginaliseEqual(CD, "simA", "simD")}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{mainip}

