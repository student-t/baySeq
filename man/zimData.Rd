\name{zimData}
\docType{data}
\alias{zimData}

\title{Simulated data for testing the baySeq package methods}

\description{This data set is a matrix (`zimData') of zero-inflated
  simulated counts from a simple pairwise expression analysis. It is
  simulated according to a negative binomial distribution with varying
  parameters for each row, and with zero-inflation applied to each
  row. The first hundred rows of the data are truly differentially
  expressed, the remainder have no differential expression.}

\usage{zimData}
\format{A matrix of which each of the ten columns represents a sample,
  and each row some discrete data (acquired by sequencing).}
\source{Simulation.}

\keyword{datasets}
    