\name{getPriors}
\alias{getPriors}
\alias{getPriors.BB}
\alias{getPriors.NB}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimates prior parameters for the underlying distributions of
  'count' data.}
\description{
  These functions estimate, via maximum or quasi-likelihood methods, the
  parameters of the underlying distributions for negative binomial
  distributions on count data, or for beta-binomial distributions on
  paired count data.
}
\usage{
getPriors.NB(cD, samplesize = 1e5, samplingSubset = NULL,
equalDispersions = TRUE, estimation = "QL", verbose = TRUE, zeroML =
FALSE, consensus = FALSE, cl, ...)
getPriors.BB(cD, samplesize = 1e5, samplingSubset = NULL, verbose = TRUE, cl, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cD}{A \code{\link{countData}} or \code{\link{pairedData}} object.}
  \item{samplesize}{How large a sample should be taken in estimating the
    priors?}
  \item{samplingSubset}{If given, the priors will be sampled only from
    the subset specified.}
  \item{equalDispersions}{Should we assume equal dispersions of data
    across all groups in the \code{'cD'} object? Defaults to TRUE; see Details.}
  \item{estimation}{Defaults to "QL", indicating quasi-likelihood
    estimation of priors. Currently, the only other possibilities are "ML",
    a maximum-likelihood method, and "edgeR", the moderated dispersion
    estimates produced by the 'edgeR' package. See Details.}
  \item{verbose}{Should status messages be displayed? Defaults to TRUE.}
  \item{zeroML}{Should parameters from zero data (rows that within a
    group are all zeros) be estimated using maximum likelihood
    methods (which will result in zeros in the parameters? See Details.}
  \item{consensus}{If TRUE, creates a consensus distribution rather than
    a separate distribution for each member of the groups structure in
    the `cD' object. See Details.}
  \item{cl}{A SNOW cluster object.}
  \item{...}{Additional parameters to be passed to the
    \code{\link[edgeR:edgeR-package]{estimateTagwiseDisp}} function if
    \code{'estimation = "edgeR"'}.}
}
\details{
  These functions empirically estimate prior parameters for the
  distributions used in estimating posterior likelihoods of each count
  belonging to a particular group. For unpaired count data, the
  distributions on the data are assumed to be negative binomial. For
  paired count data, the distributions on the data are assumed to be
  beta-binomial.

  For priors estimated for the negative binomial methods, three options
  are available. Differences in the options focus on the way in which
  the dispersion is estimated for the data. In simulation studies,
  quasi-likelihood methods (\code{'estimation = "QL"'}) performed best
  and so these are used by default. Alternatives are maximum-likelihood
  methods (\code{'estimation = "ML"'}), and the 'edgeR' packages
  moderated dispersion estimates (\code{'estimation = "edgeR"'}). 
  
  The priors estimated for the negative binomial methods
  (\code{'getPriors.NB'}) may assume that the dispersion of data for a
  given row is identical for all group structures defined in
  \code{'cD@groups'} (\code{'equalDispersions = TRUE'}). Alternatively,
  the dispersions may be estimated individually for each group structure
  (\code{'equalDispersions = FALSE'}). Unless there is a strong reason
  for believing that the data are differently dispersed between groups,
  \code{'equalDispersions = TRUE'} is recommended. If \code{'estimation
    = "edgeR"'} then this parameter is ignored and dispersion is assumed
  identical for all group structures.
  
  If all counts in a given row for a given group are zero, then maximum
  and quasi-likelihood estimation methods will result in a zero
  parameter for the mean. In analyses where segment length is a factor, this  makes
  it hard to differentiate between (for example) a region which contains
  no reads but is only ten bases long and one which likewise contains no
  reads but is ten megabases long. If \code{'zeroML'} is \code{FALSE},
  therefore, the dispersion is set to 1 and the mean estimated as the
  value that leaves the likelihood of zero data at fifty percent.

  If `consensus = TRUE', then a consensus distribution is created and
  used for each group in the 'cD' object. This allows faster computation
  of the priors and likelihoods, but with some degradation of
  accuracy.
  
  A 'cluster' object is recommended in order to estimate the priors for
  the negative binomial distribution. Passing NULL to this variable will
  cause the function to run in non-parallel mode.
}
\value{
  A \code{\link{countData}} object.
}
\references{Hardcastle T.J., and Kelly, K. baySeq: Empirical Bayesian
  Methods For Identifying Differential Expression In Sequence Count Data. BMC Bioinformatics (2010)}
\author{Thomas J. Hardcastle}

\seealso{\code{\link{countData}}, \code{\link{getLikelihoods}}}
\examples{
# See vignette for more examples.


# If we do not wish to parallelise the functions we set the cluster
# object to NULL.

cl <- NULL

# Alternatively, if we have the 'snow' package installed we
# can parallelise the functions. This will usually (not always) offer
# significant performance gain.

\dontrun{try(library(snow))}
\dontrun{try(cl <- makeCluster(4, "SOCK"))}

# load test data
data(simData)

# Create a {countData} object from test data.

replicates <- c("simA", "simA", "simA", "simA", "simA", "simB", "simB", "simB", "simB", "simB")
groups <- list(NDE = c(1,1,1,1,1,1,1,1,1,1), DE = c(1,1,1,1,1,2,2,2,2,2))
CD <- new("countData", data = simData, replicates = replicates, groups = groups)

#estimate library sizes for countData object
CD <- getLibsizes(CD)

# Get priors for negative binomial method
CDPriors <- getPriors.NB(CD, samplesize = 10^5, estimation = "QL", cl = cl)

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{models}
