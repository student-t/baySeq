\name{densityFunction-class}
\Rdversion{1.1}
\docType{class}
\alias{densityFunction-class}

\title{Class \code{"densityFunction"}}
\description{
  This function fills the '@densityFunction' slot of a `countData'
  object. It defines the distribution used to estimate posterior
  likelihoods, and associated values used in these calculations.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("densityFunction", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{description}:}{A description of the distribution defined.}
    \item{\code{density}:}{A \code{"function"}, defining the likelihood
      of a data array given observed data and hyperparameters.}
    \item{\code{initiatingValues}:}{A \code{"list"} of functions (may be
    supplied as numerics) that define initial values of numeric prior discovery.}
    \item{\code{equalOverReplicates}:}{A \code{"logical"}, describing
      which of the hyperparameters are equally marginally distributed
      over all groups, and which are not.}
    \item{\code{lower}:}{A \code{"function"}, required to define the
      lower limit of optimisation in the case where only one
      hyperparameter is *not* equally marginally distributed over all groups.}
    \item{\code{upper}:}{A \code{"function"}, required to define the
      upper limit of optimisation, as for `lower', above.}
    \item{\code{stratifyFunction}:}{An optional \code{"function"}, used
      to stratify the data for more accurate prior estimation.}
    \item{\code{stratifyBreaks}:}{An optional \code{"numeric"}, used to
      define the number of strata in a stratification.}
    \item{\code{nullFunction}:}{An optional \code{"function"} on the
      hyperparameters, used to generate a one-dimensional distribution
      which can be partitioned to identify `null' data.}
    \item{\code{orderingFunction}:}{An optional \code{"function"} for
      ordering the data between groups of a model.}
  }
}
\section{Methods}{
No methods defined with class "densityFunction" in the signature.
}
\author{
Thomas J. Hardcastle
}

\examples{
showClass("densityFunction")
}
\keyword{classes}
