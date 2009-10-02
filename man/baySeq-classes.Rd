\name{baySeq-classes}
\docType{class}
\alias{baySeq-classes}
\alias{baySeq-class}
\alias{countData}
\alias{countDataPosterior}
\alias{countData-class}
\alias{countDataPosterior-class}
\alias{dim,countData-method}
\alias{dim,countDataPosterior-method}
\alias{[,countData-method}
\alias{[,countDataPosterior-method}
\alias{show,countData-method}
\alias{show,countDataPosterior-method}

\title{baySeq - classes}

\description{
The \code{countData} class defines a class for summarising count data
and estiablishing prior parameters on distributions defined upon the
count data. The \code{countDataPosterior} extends this class by
including posterior likelihood estimates on the data. Only
\code{countData} objects should be established by the user
directly. \code{countDataPosterior} objects are established by the
functions described in \code{\link{getLikelihoods}}.}

\section{Slots}{
Objects of these class should contain the following list components:  
  \tabular{ll}{
  \code{data}:\tab Count data (matrix).\cr
  \code{libsizes}:\tab Vector of library size for each sample.\cr
  \code{groups}:\tab Group (model) structure to test on the data
  (list). \cr
  \code{annotation}:\tab Annotation data for each count (data.frame). \cr
  \code{priors}:\tab Prior parameter information. \cr
  \code{posteriors}:\tab Estimated posterior likelihoods for each group
  (matrix). 'countDataPosterior' class only. \cr
  \code{estProps}:\tab Estimated proportion of tags belonging to each
  group (numeric). 'countDataPosterior' class only. \cr
  }}

\section{Methods}{
Methods 'dim', '[' and 'show' have been defined for these classes.
}

\author{Thomas J. Hardcastle}

\examples{

library(baySeq)

data(simCount)
data(libsizes)

groups <- list(c(1,1,1,1,1,1,1,1,1,1), c(1,1,1,1,1,2,2,2,2,2))

#create new 'countData' object
CD <- new("countData", data = simCount, libsizes = libsizes, groups = groups)

CD[1:10,]
dim(CD)

}

\keyword{classes}
