\name{getPosteriors}
\alias{getPosteriors}
\title{An internal function in the baySeq package for calculating
  posterior likelihoods given likelihoods of the data.}
\description{
  For likelihoods of the data given a set of models, this function
  calculates the posterior likelihoods of the models given the data. An
  internal function of baySeq, which should not in general be called by
  the user.
}
\usage{
getPosteriors(ps, prs, pET = "none", marginalise = FALSE, groups, priorSubset = NULL, maxit = 100, accuracy =
1e-5, cl = cl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ps}{A matrix containing likelihoods of the data for each count
    (rows) under each model (columns).}
  \item{prs}{(Initial) prior probabilities for each of the models.}
  \item{pET}{What type of prior re-estimation should be
    attempted? Defaults to "none"; "BIC" and "iteratively" are also
    available.}
  \item{marginalise}{Should an attempt be made to numerically
    marginalise over a prior distribution iteratively estimated from the
    posterior distribution? Defaults to FALSE, as in general offers
    little performance gain and increases computational cost considerably.}
  \item{groups}{Group structure from which likelihoods in \code{'ps'}
    were defined.}
  \item{priorSubset}{If \code{'estimatePriors = TRUE'}, what subset of
    the data should be used to re-estimate the priors? Defaults to NULL,
    implying all data will be used.}
  \item{maxit}{What is the maximum number of iterations that should be
    tried if we are bootstrapping prior probabilities from the data?}
  \item{accuracy}{How small should the difference in estimated priors be
    before we stop bootstrapping.}
  \item{cl}{A SNOW cluster object.}
}

\details{
  An internal function, that will not in general be called by the
  user. It takes the log-likelihoods of the data given the models being
  tested and returns the posterior likelihoods of the models.

  The function may attempt to estimate the prior likelihoods
  either by using the Bayesian Information Criterion (\code{'pET =
  "BIC"'}) to identify the proportion of the data best explained by each
  model and taking these proportions as prior. Alternatively, an
  iterative re-estimation of priors is possible (\code{'pET = "iteratively"'},
  in which an inital estimate for the prior likelihoods of the models is
  used to calculated the posteriors and then the priors are updated by
  taking the mean of the posterior likelihoods for each model across all
  data.
}
\value{
  A list containing
  posteriors: estimated posterior likelihoods of the model for each
  count (log-scale)
  priors: estimated (or given) prior probabilities of the model
}

\references{Hardcastle T.J., and Kelly, K (2010). Identifying Patterns
  of Differential Expression in Count Data. In submission.}
\author{Thomas J. Hardcastle}

\seealso{\code{\link{getLikelihoods}}}

\examples{

# Simulate some log-likeihoods of data given models (each model
# describes one column of the 'ps' object).
ps <- log(rbind(
                cbind(runif(10000, 0, 0.1), runif(10000, 0.3, 0.9)),
                cbind(runif(10000, 0.4, 0.9), runif(1000, 0, 0.2))))

# get posterior log-likelihoods of model, estimating prior likelihoods
# of each model from the data.

pps <- getPosteriors(ps, prs <- c(0.5, 0.5), pET = "none", cl =
NULL)

pps$priors

pps$posteriors[1:10,]

}

\keyword{models}
