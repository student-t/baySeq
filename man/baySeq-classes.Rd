\name{baySeq-classes}
\docType{class}
\alias{baySeq-classes}
\alias{baySeq-class}
\alias{countData}
\alias{countData-class}
\alias{dim,countData-method}
\alias{[,countData-method}
\alias{show,countData-method}
\alias{segData}
\alias{segData-class}
\alias{dim,segData-method}
\alias{[,segData-method}
\alias{show,segData-method}


\title{baySeq - classes}

\description{
  The \code{countData} class defines a class for summarising count data
  and estiablishing prior and posterior parameters on distributions defined upon the
  count data. The \code{segData} class extends the \code{countData} class
  by including segment lengths.}

\section{Slots}{
  Objects of these class should contain the following list components:  
  \tabular{ll}{
    \code{data}:\tab Count data (matrix).\cr
    \code{libsizes}:\tab Vector of library size for each sample.\cr
    \code{groups}:\tab Group (model) structure to test on the data
    (list). \cr
    \code{annotation}:\tab Annotation data for each count (data.frame). \cr
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
    described in \code{data}. Only in 'segData' class. See Details. \cr
}}

\section{Details}{
The \code{seglens} slot describes, for each row of the \code{data}
object, the length of the 'segment' that contains the number of counts
described by that row. For example, if we are looking at the number of
hits matching genes, the \code{seglens} object would consist of
transcript lengths.
}

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
