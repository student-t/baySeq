\name{pairData}
\docType{data}
\alias{pairData}

\title{Simulated data for testing the baySeq package methods; simulated
  counts and library sizes from paired samples for differential expression
  analysis.}

\description{This data set is a matrix (`pairData') of simulated counts
  from paired samples for a pairwise expression analysis. It is
  simulated according to a negative binomial distribution with varying
  parameters for each row. The first hundred rows of the data are truly
  differentially expressed, the remainder have no differential expression.}

\usage{pairData}
\format{A matrix of which each of the eight columns represents a sample,
  and each row some discrete data (acquired by sequencing). The first
  four columns are paired with the second four columns,
  respectively.}
\source{Simulation.}

\keyword{datasets}
    