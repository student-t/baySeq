\name{simSeg}
\docType{data}
\alias{simSeg}

\title{Simulated data for testing the baySeq package methods; simulated
  counts and segment lengths from a pairwise differential expression analysis}

\description{This data set is a matrix of simulated counts from a simple
  pairwise expression analysis of genomic regions. The first column of
  the data is the length of the segment, simulated from a negative
  binomial distribution. The remaining columns of the data are simulated according to a negative
  binomial distribution with varying parameters for each row, such that
  the rate of expression is proportional to the length of the segment. The first
  hundred rows of the data are truly differentially expressed, the
  second hundred rows have no true expression of any kind, the
  remainder are expressed but not differentially expressed. Library sizes for these
  data are given in \code{\link{libsizes}}.}

\usage{simCount}
\format{A matrix of which the first column gives the length of a genomic
  region, and the following ten columns represents the discrete data
  (acquired by sequencing) observed at each region.}
\source{Simulation.}
\seealso{\code{\link{libsizes}}}
\references{Hardcastle T.J., and Kelly, K. baySeq: Empirical Bayesian
  Methods For Identifying Differential Expression In Sequence Count Data. BMC Bioinformatics (2010)}


\keyword{datasets}
    