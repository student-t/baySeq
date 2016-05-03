\name{mobData}
\docType{data}
\alias{mobData}

\title{Data from a set of small RNA sequencing experiments carried out
  on grafts of Arabidopsis thaliana intended for differential expression
  analyses.}

\description{This data set is a matrix (`mobData') of counts acquired
  for three thousand small RNA loci from a set of Arabidopsis grafting
  experiments. Three different biological conditions exist within these
  data; one in which a Dicer 2,3,4 triple mutant shoot is grafted onto a
  Dicer 2,3,4 triple mutant root (SL236 and SL260), one
  in which a wild-type shoot is grafted onto a wild-type root
  (SL239 and SL240), and one in which a wild-type shoot
  is grafted onto a Dicer 2,3,4 triple mutant root (SL237 and
  SL238). Dicer 2,3,4 is required for the production of 22nt
  and 24nt small RNAs, as well as some 21nt ones. Consequently, if we
  detect differentially expressed  sRNA loci in the root stock of the
  grafts, we can make inferences about the mobility of small RNAs.

  The annotation of the loci from which these data derive is in data
  file `mobAnnotation'.
}

\usage{mobData}
\format{A matrix of which each of the six columns represents a sample,
  and each row an sRNA locus (acquired by sequencing).}
\source{Illumina sequencing.}
\seealso{\code{\link{mobAnnotation}}}

\references{Molnar A. and Melnyk C.W. et al. Small silencing RNAs in
  plants are mobile and direct epigenetic modification in recipient
  cells. Science (2010)}

\keyword{datasets}
    