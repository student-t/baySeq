\name{methObservables}
\alias{methObservables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generation of intermediate values in likelihood estimation for methylation
data.
}
\description{
  Estimating prior and posterior values for methylation data that
  account for non-conversion rates is a time-consuming
  process. Significant increases in speed can be made by calculating in
  advance sets of data that will be re-used at several points of these
  analyses. This function populates the '@cellObservables' slot of a
  `countData' object that contains a `nonconversion' object in the
  `@sampleObservables' slot.
}
\usage{
methObservables(mD)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mD}{
    A \code{\link{countData}} object, or descendant, which contains a
    numeric vector 'nonconversion' in the `@sampleObservables' slot.
}
}

\value{
A \code{\link{countData}} object with the '@cellObservables' slot
populated with temporary values useful in the faster calculation of likelihoods.
}
\author{
  Thomas J. Hardcastle
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
