\name{factCount}
\docType{data}
\alias{factCount}

\title{Simulated data for testing the baySeq package methods; simulated
  counts from a factorial design differential expression analysis}

\description{This data set is a matrix of simulated counts from a simple
  pairwise expression analysis. It is simulated according to a negative
  binomial distribution with varying parameters for each row. The first
  hundred rows of the data are truly differentially expressed between
  the first four samples and the second four samples. The second hundred
  rows of the data are truly differentially expressed between samples
  1,2,5,6 and samples 3,4,7,8.
}

\usage{factCount}
\format{A matrix of which each of the eight columns represents a sample,
  and each row some discrete data (acquired by sequencing).}
\source{Simulation.}
\seealso{\code{\link{factlibsizes}}}
\references{Hardcastle T.J., and Kelly, K (2010). Identifying Patterns
  of Differential Expression in Count Data. In submission.}


\keyword{datasets}
    