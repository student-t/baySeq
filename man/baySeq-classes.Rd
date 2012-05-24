\name{baySeq-classes}
\docType{class}
\alias{baySeq-classes}
\alias{baySeq-class}
\alias{countData}
\alias{countData-class}
\alias{pairedData}
\alias{pairedData-class}
\alias{show,pairedData-method}
\alias{[,pairedData-method}
\alias{dim,countData-method}
\alias{[,countData-method}
\alias{show,countData-method}
\alias{rbind,countData-method}
\alias{replicates,countData-method}
\alias{replicates<-,countData-method}
\alias{replicates}
\alias{replicates<-}
\alias{groups,countData-method}
\alias{groups<-,countData-method}
\alias{groups}
\alias{groups<-}
\alias{rbind}

\title{baySeq - classes}

\description{
  The \code{countData} class is used to define summaries of count data
  and establishing prior and posterior parameters on distributions
  defined upon the count data.

  The \code{pairedData} class is based on the countData class and performs
  the same functions for paired count data.
}

\section{Slots}{
  Objects of these class contain the following components:  
  \tabular{ll}{
    \code{data}:\tab Count data (matrix).\cr
    \code{pairData}:\tab Count data (matrix) for paired samples; slot
  exists only in \code{'pairedData'} class.\cr
    \code{replicates}:\tab The replicate structure of the data. Stored
  as a factor, but can be given in any form.\cr
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
  The standard methods 'new', 'dim', '[', 'show' and 'rbind' have been
  defined for these classes. The methods 'groups', 'groups<-', 'replicates'
  and 'replicates<-' have also been defined in order to access and
  modify these slots, and their use is recommended.
}

\author{Thomas J. Hardcastle}

\examples{

#load test data
data(simData)

# Create a 'countData' object from test data.
replicates <- c("simA", "simA", "simA", "simA", "simA", "simB", "simB", "simB", "simB", "simB")
groups <- list(NDE = c(1,1,1,1,1,1,1,1,1,1), DE = c(1,1,1,1,1,2,2,2,2,2))
CD <- new("countData", data = simData, replicates = replicates, groups = groups)

#estimate library sizes for countData object
CD <- getLibsizes(CD)

CD[1:10,]
dim(CD)

}

\keyword{classes}
