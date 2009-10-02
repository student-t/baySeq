\name{PDgivenr}
\alias{PDgivenr}
\alias{PDgivenr.Dirichlet}
\alias{PDgivenr.Pois}
\alias{PDgivenr.PoisIndie}
\alias{PDgivenr.NBIndie}
\title{Functions for estimating the likelihood of the data given an
  equivalence relation for some assumption on the distribution.}
\description{
  These functions aim to calculate the likelihood of the data given an
  equivalence relation R for some assumed distribution. They are
  internal functions to \link{baySeq} and should not be called by
  the user.
}
\usage{
PDgivenr.Dirichlet(us, prior, group)
PDgivenr.Pois(us, ns, prior, group)
PDgivenr.PoisIndie(us, ns, prior, group)
PDgivenr.NBIndie(us, ns, prior, group, priorWeights)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{us}{Number of counts.}
  \item{ns}{Library size.}
  \item{prior}{Prior parameters on distribution.}
  \item{group}{Hypothesised grouping of the data based on equivalence
    relation.}
  \item{priorWeights}{Weightings for each prior to be used in numerical integration.}
}
\details{
  Not intended to be called by the user.
}
\value{
  Numeric value giving log-likelihood of data given distribution and grouping.
}
\references{Hardcastle T.J., and Kelly, K (2009). Empirical Bayesian
  methods for differential expression in count data. In submission.}
\author{Thomas J. Hardcastle}

\examples{

data(simCount)
data(libsizes)
groups <- list(c(1,1,1,1,1,1,1,1,1,1), c(1,1,1,1,1,2,2,2,2,2))
CD <- new("countData", data = simCount, libsizes = libsizes, groups =
groups)

# Estimate priors using Poisson method.
CDP.Poi <- getPriors.Pois(CD, samplesize = 20, iterations = 1000,
takemean = TRUE)

# Calculate likelihood of first row of data given the first hypothesis
# defined by the 'groups' slot of CDP.Poi by Poisson methods

PDgivenr.Pois(us = CDP.Poi@data[1,], ns = CDP.Poi@libsizes, prior =
CDP.Poi@priors$priors[[1]], group = CDP.Poi@groups[[1]])

# Calculate likelihood of each row of data given the first hypothesis
# defined by the 'groups' slot of CDP.Poi by Poisson methods

apply(CDP.Poi@data, 1, PDgivenr.Pois, ns = CDP.Poi@libsizes, prior =
CDP.Poi@priors$priors[[1]], group = CDP.Poi@groups[[1]])

}

\keyword{distribution}
\keyword{models}
