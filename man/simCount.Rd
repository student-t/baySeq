\name{simCount}
\docType{data}
\alias{simCount}

\title{Simulated data for testing the baySeq package methods; simulated
  counts from a pairwise differential expression analysis}

\description{This data set is a matrix of simulated counts from a simple
  pairwise expression analysis. It is simulated according to a negative
  binomial distribution with varying parameters for each row. The first
  hundred rows of the data are truly differentially expressed, the
  remainder have no differential expression. Library sizes for these
  data sets are given in \code{\link{libsizes}}.}

\usage{simCount}
\format{A matrix of which each of the ten columns represents a sample,
  and each row some discrete data (acquired by sequencing).}
\source{Simulation.}
\seealso{\code{\link{libsizes}}}
\references{Hardcastle T.J., and Kelly, K (2009). Empirical Bayesian
  methods for differential expression in count data. In submission.
}

\keyword{datasets}
    