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
getPriors.Dirichlet(cD, samplesize = 10^5, perSE = 1e-1, maxit = 10^6)
getPriors.Pois(cD, samplesize = 10^5, perSE = 1e-1, takemean = TRUE,
maxit = 10^5, cl)
getPriors.NB(cD, samplesize = 10^5, estimation = "ML", cl)
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
  \item{estimation}{Defaults to "ML", indicating maximum likelihood
    estimation of priors. Currently, the only other possibility is "QL",
  a quasi-likelihood method.}
  \item{cl}{A SNOW cluster object.}
}
\details{
  These functions empirically estimate prior parameters for different
  distributions used in estimating posterior likelihoods of each count
  belonging to a particular group. The choice of which function to use
  for estimating the prior parameters will depend on the choice of which
  method is being used to estimate the posterior likeihoods (see
  \link{getLikelihoods}).

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
\references{Hardcastle T.J., and Kelly, K (2009). Empirical Bayesian
  methods for differential expression in count data. In submission.}
\author{Thomas J. Hardcastle}

\seealso{\code{\link{countData}}, \code{\link{getLikelihoods}}}
\examples{
# See vignette for more examples.

# Create a {countData} object.
data(simCount)
data(libsizes)
groups <- list(c(1,1,1,1,1,1,1,1,1,1), c(1,1,1,1,1,2,2,2,2,2))
CD <- new("countData", data = simCount, libsizes = libsizes, groups =
groups)

# Estimate priors using Poisson method.
CDP.Poi <- getPriors.Pois(CD, samplesize = 20, takemean = TRUE, cl = NULL)

\dontrun{

# Alternatively, get priors for negative binomial method

CDP.NBML <- getPriors.NB(CD, samplesize = 10^5, estimation = "ML", cl = NULL)

# Alternatively, if we have the 'snow' package installed we
# can parallelise the prior estimation (for the negative binomial
# methods only). This will usually (not always) offer
# significant performance gain.

library(snow)
cl <- makeCluster(4, 'SOCK')

CDP.NBML <- getPriors.NB(CD, samplesize = 10^5, estimation = "ML", cl = cl)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{models}
