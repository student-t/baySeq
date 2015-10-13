\name{topCounts}
\alias{topCounts}
\title{Get the top counts corresponding to some group from a
  'countData' object}
\description{
  Takes posterior likelihoods and returns the counts with highest (or
  lowest) likelihood of association with a given group.
}
\usage{
topCounts(cD, group, ordering, decreasing = TRUE, number = 10, likelihood, FDR,
FWER, normaliseData = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cD}{\code{\link{countData}} or \code{\link{pairedData}} object,
    containing posterior likelihoods for each group.}
  \item{group}{Which group should we give the counts for? See Details.}
  \item{ordering}{If specified, restricts the analysis to a particlar
    ordering on the group.}
  \item{decreasing}{Ordering on posterior likelihoods.}
  \item{number}{How many results should be returned?}
  \item{likelihood}{If given, ignores `number' and returns all results
    above a certain likelihood.}
  \item{FDR}{If given, ignores `number' and returns all results with an
    FDR lower than this threshold. Will be ignored if `likelihood' is
    given.}
  \item{FWER}{If given, ignores `number' and returns all results with an
    FWER lower than this threshold. Will be ignored if `likelihood' or
    FDR are given.}
  \item{normaliseData}{Should the displayed counts be normalised? See
    details. Defaults to FALSE.}
}

\details{  
  The argument 'group' can be specified either as a number, giving the
  index of an element in the \code{cD@groups} list, or as a character
  string identifying an element by name. Partial matching is allowed. If
  group = NULL, then the function looks at the posterior likelihoods
  that the data have no true differential expression (if calculated). 
  
  If a \code{\link{countData}} object is given, the returned dataframe
  will contain either the raw counts for that object, or (if
  'normaliseData = TRUE' the counts normalised by library size.

  If a \code{\link{pairedData}} object is given, the returned data frame
  will contain either the raw counts shown as a ratio, or (if
  'normaliseData = TRUE' the counts in the '@data' slot as a percentage
  of the total counts in the '@data' and '@pairData' slots.
}

\value{
  A dataframe of the top counts associated with some model (group),
  described by annotation drawn from the '@annotation' slot of the
  'cD' object and the raw data from the '@data' slot, together with the
  posterior likelihoods and false discovery rates.
}
\author{Thomas J. Hardcastle}
\seealso{\code{\link{countData}}}

\examples{

# We load in a `countData' object containing the estimated posterior
# likelihoods of expression (see `getLikelihoods').

data(CDPost)

# Report the top ten rows of data that have highest likelihood of belonging to
# group 2 of the data (i.e., differentially expressed)

topCounts(CDPost, group = "DE", number = 10)

# equivalently...
topCounts(CDPost, group = 2, number = 10)


# Report the top ten rows of data that have highest likelihood of belonging to
# group 2 of the data (i.e., differentially expressed), with group 1
# being overexpressed compared to group 2.

topCounts(CDPost, group = "DE", ordering = "1>2", number = 10)

}

\keyword{print}

