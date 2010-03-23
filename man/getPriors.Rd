\name{getPriors}
\alias{getPriors}
\alias{getPriors.Dirichlet}
\alias{getPriors.Pois}
\alias{getPriors.NB}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimates prior parameters for the underlying distributions of
  'count' data.}
\description{
  These functions estimate, via maximum likelihood methods, the
  parameters of the underlying distributions for the different methods
  of describing the 'count' data.
}
\usage{
getPriors.Dirichlet(cD, samplesize = 10^5, perSE = 1e-1, maxit = 10^6,
verbose = TRUE)
getPriors.Pois(cD, samplesize = 10^5, perSE = 1e-1, takemean = TRUE,
maxit = 10^5, verbose = TRUE, cl)
getPriors.NB(cD, samplesize = 10^5, equalDispersions = TRUE, estimation
= "QL", verbose = TRUE, cl, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cD}{A \code{\link{countData}} object.}
  \item{samplesize}{How large a sample should be taken in estimating the
    priors?}
  \item{perSE}{What should the relative standard error of the estimated
    parameters fall below?}
  \item{maxit}{Over how many iterations (at most) should we take samples and
    re-estimate the priors in order to achieve convergence?}
  \item{takemean}{If TRUE (recommended), we take the mean of the
    estimated priors to define a gamma distribution. If FALSE, we use
    all estimated priors to define an empirical distribtion on the
    parameters of the gamma distribution.}
  \item{equalDispersions}{Should we assume equal dispersions of data
    across all groups in the \code{'cD'} object? Defaults to TRUE; see Details.}
  \item{estimation}{Defaults to "QL", indicating quasi-likelihood
    estimation of priors. Currently, the only other possibilities are "ML",
    a maximum-likelihood method, and "edgeR", the moderated dispersion
    estimates produced by the 'edgeR' package. See Details.}
  \item{verbose}{Should status messages be displayed? Defaults to TRUE.}
  \item{cl}{A SNOW cluster object.}
  \item{...}{Additional parameters to be passed to the
    \code{\link[edgeR:edgeR-package]{estimateTagwiseDisp}} function if
    \code{'estimation = "edgeR"'}.}
}
\details{
  These functions empirically estimate prior parameters for different
  distributions used in estimating posterior likelihoods of each count
  belonging to a particular group. The choice of which function to use
  for estimating the prior parameters will depend on the choice of which
  method is being used to estimate the posterior likelihoods (see
  \link{getLikelihoods}).

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
  
  A 'cluster' object is recommended in order to estimate the priors for
  the negative binomial distribution. Passing NULL to this variable will
  cause the function to run in non-parallel mode.

  getPriors.Dirichlet and getPriors.Pois will issue warnings if the
  estimation of any priors fails to achieve less than the relative
  standard error specified in the maximum number of iterations.
}
\value{
  A \code{\link{countData}} object.
}
\references{Hardcastle T.J., and Kelly, K (2010). Identifying Patterns
  of Differential Expression in Count Data. In submission.}
\author{Thomas J. Hardcastle}

\seealso{\code{\link{countData}}, \code{\link{getLikelihoods}}}
\examples{
# See vignette for more examples.

# Create a {countData} object.
data(simCount)
data(libsizes)
replicates <- c(1,1,1,1,1,2,2,2,2,2)
groups <- list(c(1,1,1,1,1,1,1,1,1,1), c(1,1,1,1,1,2,2,2,2,2))
CD <- new("countData", data = simCount, replicates = replicates, libsizes = libsizes, groups =
groups)

# If we have the 'snow' package installed we can parallelise the prior
# estimation. This will usually (depending on your parallelisation
# set-up) offer significant performance gains. 

cl <- NULL
try(library(snow))
try(cl <- makeCluster(4, "SOCK"))

# Estimate priors using Poisson method.
DP.Poi <- getPriors.Pois(CD, samplesize = 20, takemean = TRUE, cl = cl)

# Alternatively, get priors for negative binomial method.

CDP.NBML <- getPriors.NB(CD, samplesize = 10^5, estimation = "QL", cl = cl)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{models}
