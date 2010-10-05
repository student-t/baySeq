\name{plotMA.CD}
\alias{plotMA.CD}

\title{
'MA'-plot for count data.
}
\description{
This function creates an MA-plot from two sets of samples. For those
data where the log-ratio is infinite (because in one set of sample data
all observed counts are zero), we plot instead the log-values of the
other group.
}
\usage{
plotMA.CD(cD, samplesA, samplesB, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cD}{
    A \code{\link{countData}} object.
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
\author{
  Thomas J. Hardcastle
}

\seealso{
\code{\link{countData}}
}
\examples{
data(simCount)
data(libsizes)
replicates <- c(1,1,1,1,1,2,2,2,2,2)
groups <- list(c(1,1,1,1,1,1,1,1,1,1), c(1,1,1,1,1,2,2,2,2,2))
CD <- new("countData", data = simCount, replicates = replicates, libsizes = libsizes, groups = groups)

plotMA.CD(CD, samplesA = 1:5, samplesB = 6:10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplots}
