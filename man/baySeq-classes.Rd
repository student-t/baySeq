\name{baySeq-classes}
\docType{class}
\alias{baySeq-classes}
\alias{baySeq-class}
\alias{countData}
\alias{countData-class}
\alias{dim,countData-method}
\alias{[,countData-method}
\alias{show,countData-method}

\title{baySeq - classes}

\description{
  The \code{countData} class is used to define summaries of count data
  and establishing prior and posterior parameters on distributions defined upon the
  count data.
}

\section{Slots}{
  Objects of the \code{'countData'} class should contain the following components:  
  \tabular{ll}{
    \code{data}:\tab Count data (matrix).\cr
    \code{libsizes}:\tab Vector of library size for each sample.\cr
    \code{groups}:\tab Group (model) structure to test on the data
    (list). \cr
    \code{annotation}:\tab Annotation data for each count (data.frame). \cr
    \code{priorType}:\tab Character string describing the type of prior
  information available in slot \code{'priors'}. \cr
    \code{priors}:\tab Prior parameter information. Calculated by the
    functions described in \code{\link{getPriors}}. \cr
    \code{posteriors}:\tab Estimated posterior likelihoods for each group
    (matrix). Calculated by the functions described in \code{\link{getLikelihoods}}. \cr
    \code{estProps}:\tab Estimated proportion of tags belonging to each
    group (numeric). Calculated by the functions described in
    \code{\link{getLikelihoods}}. \cr
    \code{nullPosts}:\tab If calculated, the posterior likelihoods for
  the data having no true expression of any kind. \cr
    \code{seglens}:\tab Lengths of segments containing the counts
    described in \code{data}. A matrix, but may be initialised with a
  vector, or ignored altogether. See Details. \cr
}}

\section{Details}{
The \code{seglens} slot describes, for each row of the \code{data}
object, the length of the 'segment' that contains the number of counts
described by that row. For example, if we are looking at the number of
hits matching genes, the \code{seglens} object would consist of
transcript lengths. Exceptionally, we may want to use different segment
lengths for different samples and so the slot takes the form of a
matrix. If the matrix has only one column, it is duplicated for all
samples. Otherwise, it should have the same number of columns as the
'@data' slot. If the slot is the empty matrix, then it is assumed that
all segments have the same length.
}

\section{Methods}{
  Methods 'new', 'dim', '[' and 'show' have been defined for this class.
}

\author{Thomas J. Hardcastle}

\examples{

library(baySeq)

data(simCount)
data(libsizes)

replicates <- c(1,1,1,1,1,2,2,2,2,2)
groups <- list(c(1,1,1,1,1,1,1,1,1,1), c(1,1,1,1,1,2,2,2,2,2))

#create new 'countData' object
CD <- new("countData", data = simCount, replicates = replicates, libsizes = libsizes, groups = groups)

CD[1:10,]
dim(CD)

}

\keyword{classes}
