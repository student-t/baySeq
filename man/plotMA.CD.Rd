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
plotMA.CD(cD, samplesA, samplesB, normaliseData = TRUE, scale = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cD}{
    A \code{\link{countData}} object.
  }
  \item{samplesA}{
    Either a character vector, identifying sample set A by either
    replicate name or sample name, or a numerical vector giving the columns of data in the
    \code{'countData'} object that forms sample set A. See Details.
}
\item{samplesB}{
  Either a character vector, identifying sample set B by either
  replicate name or sample name, or a numerical vector giving the columns of data in the
  \code{'countData'} object that forms sample set B. See Details.
  }
  \item{normaliseData}{Should the data be normalised by library size
    before computing log-ratios? Defaults to TRUE.}
  \item{scale}{If given, defines the scale on which the log-ratios will
    be plotted. Defaults to NULL, implying that the scale will be
    calculated by the function.}
  \item{\dots}{
    Any other parameters to be passed to the \code{\link{plot}} function.
  }
}
\details{
  The samples sets can be identified either by a numeric vector which
  specifies the columns of data from the \code{countData} object 'cD',
  or by a character vector. If a character vector is used, the members
  of the character vector will first be searched for in the
  \code{@replicates} slot of the 'cD' object. Any members of the vector not found
  in the replicates slot, will be searched for in the column names of the
  \code{@data} slot of the 'cD' object. Different classes of vector can
  be used for 'samplesA' and 'samplesB', as shown in the example below.
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
data(simData)

replicates <- c("simA", "simA", "simA", "simA", "simA", "simB", "simB", "simB", "simB", "simB")
groups <- list(NDE = c(1,1,1,1,1,1,1,1,1,1), DE = c(1,1,1,1,1,2,2,2,2,2))
CD <- new("countData", data = simData, replicates = replicates, groups = groups)

#estimate library sizes for countData object
CD@libsizes <- getLibsizes(CD)

#MA-plot comparing replicate groups
plotMA.CD(CD, samplesA = "simA", samplesB = 6:10)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplots}
