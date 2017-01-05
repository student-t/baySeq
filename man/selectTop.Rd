\name{selectTop}
\alias{selectTop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Selects the top genomic events, based on posterior likelihoods, from a
  `countData' object.
}
\description{
  This function subsets a countData object by selecting those events
  that best (or least) represent a model, based on the posterior
  likelihoods estimated for that model and some threshold. Selection can
  be done for a specific model (and ordering of the data under that
  model) or for all models (and all orderings).
}
\usage{
selectTop(cD, group, ordering, orderings = TRUE, decreasing = TRUE,
number = 10, likelihood, FDR, FWER, posteriors)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cD}{
    A \code{\link{countData}} object, with a populated `@posteriors'
    slot (unless `posteriors' is specified; see below).
}
  \item{group}{
    Optionally, the model of interest, as defined in the `@groups' slot
    of the \code{\link{countData}} object. If unspecified, subsets for
    all models will be returned as a list.
}
  \item{ordering}{
    If `group' is specified, a particular ordering of the data based on
    that group can also be specified.
}
  \item{orderings}{
    If no group is specified, should the selection of models also be
    split by the orderings of the data under the models? Defaults to TRUE.
}
  \item{decreasing}{
    If FALSE, considers the data with the lowest posterior likelihoods,
    rather than the greatest (i.e., selects those data least likely to
    conform to a particular model.
}
  \item{number}{
    If given, selects the top 'number' of genomic events for each model
    (and optionally, ordering). Ignored if another selection criteria is
    chosen, unless this criteria would return no values.
}
  \item{likelihood}{
    If given, selects all genomic events for a given model (and
    optionally, ordering) with posterior likelihood exceeding this
    value.
}
  \item{FDR}{
    If given, selects all genomic events for a given model (and
    optionally, ordering) with false discovery rate below this
    value. Ignored if likelihood is specified.
}
  \item{FWER}{
    If given, selects all genomic events for a given model (and
    optionally, ordering) with family-wise error rate below this
    value. Ignored if likelihood or FDR is specified.
  }
  \item{posteriors}{
    If given, a vector of log-posterior likelihoods to use instead of
    those in the `@posteriors' slot of the `cD' object.
    }
}

\value{
  Either a single \code{\link{countData}} object (if `group' is
  specified), or a named list of \code{\link{countData}} objects.
}

\author{
Thomas J. Hardcastle
}

\seealso{
  \code{\link{topCounts}}
}

\examples{
# We load in a `countData' object containing the estimated posterior
# likelihoods of expression (see `getLikelihoods').

data(CDPost)

# select from all models and orderings with FDR equal to or lower than 0.01.

selectTop(CDPost, FDR = 0.01)

}

\keyword{manip}
