\name{plotNullPrior}
\alias{plotNullPrior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plots distribution of null function and shows the threshold separator.
}
\description{
  In sequencing expression of various genomic events, it is not uncommon
  to find a subset of genomic events that are qualitatively different
  from the remainder of the data. Thus, for some function of the
  estimated priors, we may observe bimodality or long tails which
  correlate to this subset.
}
\usage{
plotNullPrior(cD, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cD}{
    A \code{\link{countData}} object with a '@nullFunction' slot in its '@densityFunction'.
}
  \item{\dots}{
    Additional arguments to be passed to 'plot'
}
}
\value{
Invisibly, the numeric value of the threshold.
}

\author{
Thomas J. Hardcastle
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
