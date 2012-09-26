\name{pairData}
\docType{data}
\alias{pairData}

\title{Simulated data for testing the baySeq package methods for paired data}

\description{This data set is a matrix (`pairData') of simulated counts
  from a set of high-throughput sequencing data from a paired
  experimental design. The first four columns of data are to be paired
  with the second four columns of data respectively. The first two
  paired samples form one replicate group, the second two paired samples
  form another replicate group. The first hundred rows of the data are
  truly differentially expressed between replicate groups, the second
  hundred are differentially expressed between pairs, the remainder
  have no differential expression. 

  It is simulated according to a set of Poisson distributions whose
  parameters for each row are determined by a beta distribution on the
  relative proportions of data in each pairing.}

\usage{pairData}
\format{A matrix of which each of the eight columns represents a sample,
  and each row some discrete data (acquired by sequencing).}
\source{Simulation.}

\keyword{datasets}
    