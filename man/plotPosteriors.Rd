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
plotPosteriors(cD, group = 1, samplesA, samplesB, ...)
}
\arguments{
  \item{cD}{
    A \code{\link{countData}} object, for which posterior likelihoods have been
    estimated (see \code{\link{getPosteriors}}).
  }
  \item{group}{
    From which group (as defined in the \code{'cD@groups'} slot) should
    posterior likelihoods be shown?
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
data(simCount)
data(libsizes)
replicates <- c(1,1,1,1,1,2,2,2,2,2)
groups <- list(c(1,1,1,1,1,1,1,1,1,1), c(1,1,1,1,1,2,2,2,2,2))
CD <- new("countData", data = simCount, replicates = replicates, libsizes = libsizes, groups = groups)

# Get priors for negative binomial method
CDP.NB <- getPriors.NB(CD, samplesize = 10^5, estimation = "QL", cl = NULL)

# Get likelihoods for data with negative binomial method

CDPost.NB <- getLikelihoods.NB(CDP.NB, prs = c(0.5, 0.5),
pET = "BIC", cl = NULL)

plotPosteriors(CDPost.NB, group = 1, samplesA = 1:5, samplesB = 6:10)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
