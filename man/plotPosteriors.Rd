\name{plotPosteriors}
\alias{plotPosteriors}

\title{
Plots the posterior likelihoods estimated for a 'countData' object
against the log-ratios observed between two sets of sample data.
}
\description{
This function plots the posterior likelihoods estimated for a
'countData' object against the log-ratios observed between two sets of
sample data. For those data where the log-ratio is infinite (because in
one set of sample data all observed counts are zero), we plot instead
the log-values of the other group.
}
\usage{
plotPosteriors(cD, group, samplesA, samplesB, ...)
}
\arguments{
  \item{cD}{
    A \code{\link{countData}} object, for which posterior likelihoods have been
    estimated (see \code{\link{getPosteriors}}).
  }
  \item{group}{
    From which group (as defined in the \code{'cD@groups'} slot) should
    posterior likelihoods be shown? Can be defined either as the number
    of the element in \code{'cD@groups'} or as a string which will be
    partially matched to the names of the \code{'cD@groups'} elements.
  }
  \item{samplesA}{
    A numerical vector giving the columns of data in the
    \code{'countData'} object that forms sample set A.
  }
  \item{samplesB}{
    A numerical vector giving the columns of data in the
    \code{'countData'} object that forms sample set B.
  }
  \item{\dots}{
    Any other parameters to be passed to the \code{\link{plot}} function.
  }
}

\value{
  Plotting function.
}
\author{Thomas J. Hardcastle}
\seealso{\code{\link{getPosteriors}}}

\examples{

# We load in a `countData' object containing the estimated posterior
# likelihoods of expression (see `getLikelihoods').

data(CDPost)

plotPosteriors(CDPost, group = "DE", samplesA = 1:5, samplesB = 6:10)

# equivalent to plotPosteriors(CDPost, group = 2, samplesA = 1:5, samplesB = 6:10)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
