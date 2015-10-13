\name{summarisePosteriors}
\alias{summarisePosteriors}

\title{
Summarises expected number of genomic events given the calculated
posterior likelihoods of a countData object.
}
\description{
  Given posterior likelihoods for each model, we can calculate the
  expected number of genomic events corresponding to each model (and to
  each ordering within each model) by summing the posterior likelihoods.
}
\usage{
summarisePosteriors(cD, orderings = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cD}{
    A \code{\link{countData}} object.
}
  \item{orderings}{
    Indicates whether models should be split by orderings of the data
    under the model (defaults to TRUE).
}
}
\value{
  Numeric vector of expected number of genomic events belonging to each
  model (optionally, split by orderings).
}
\author{
Thomas J. Hardcastle
}
\seealso{\code{\link{topCounts}}, \code{\link{selectTop}}}

\examples{

# We load in a `countData' object containing the estimated posterior
# likelihoods of expression (see `getLikelihoods').

data(CDPost)

# summarise the expected number of genomic events in each category
summarisePosteriors(CDPost)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print}