\name{getLibsizes}
\alias{getLibsizes}

\title{
Estimates library scaling factors (library sizes) for count data.
}
\description{
This function estimates the library scaling factors that should be used
for either a 'countData', or a matrix of counts and replicate information.
}
\usage{
getLibsizes(cD, data, replicates, subset = NULL, estimationType = c("quantile", "total", "edgeR"), quantile = 0.75, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cD}{
    A \code{\link{countData}} object.
  }
  \item{data}{
    A matrix of count values. Ignored if 'cD' is given.
  }
  \item{replicates}{
    A replicate structure for the data given in 'data'. Ignored if 'cD'
    is given.
  }
  \item{subset}{
    A numerical vector indicating the rows of the 'countData' object
    that should be used to estimate library scaling factors.
  }
  \item{estimationType}{
    One of 'quantile', 'total', or 'edgeR'. Partial matching is
    allowed. See Details.
  }
  \item{quantile}{A value between 0 and 1 indicating the level of
    trimming that should take place. See Details.}
  \item{\dots}{
    Additional parameters to be passed to the 'edgeR' calcNormFactors function.
  }
}
\details{This function estimates the library scaling factors (surrogates
  for library size) in one of several ways, depending on the
  'estimationType' argument. 'total' will give the library sizes by
  summing all counts in each sample. 'quantile' will give a library scaling
  factor by the method of Bullard et al (Bioinformatics 2010), summing
  all counts in each sample whose value below the qth quantile of
  non-zero counts for that sample. 'edgeR' uses the Trimmed Mean of
  M-vales (TMM) method of Robinson \& Oshlack (Genome Biology, 2010) via
  the 'edgeR' calcNormFactors function; other options are
  available through this function.

  If a \code{\link{countData}} object 'cD' is given, the library sizes
  will be inferred from this. Alternatively, a matrix of count values
  (columns are libraries) and a replicate structure (a vector defining
  which samples belong to which replicate group) can be given.
  }
\value{
  A numerical vector of library sizes (scaling factors) for each library
  in the data.
}
\author{
  Thomas J. Hardcastle
}

\seealso{
\code{\link{countData}}
}
\examples{
data(simData)
replicates <- c(1,1,1,1,1,2,2,2,2,2)
groups <- list(c(1,1,1,1,1,1,1,1,1,1), c(1,1,1,1,1,2,2,2,2,2))
CD <- new("countData", data = simData, replicates = replicates, groups = groups)

CD@libsizes <- getLibsizes(CD)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
