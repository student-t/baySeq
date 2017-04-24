\name{getLikelihoods}
\alias{getLikelihoods}
\alias{getLikelihoods.NB}
\alias{getLikelihoods.BB}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Finds posterior likelihoods for each count or paired count as
  belonging to some model.}
\description{
  These functions calculate posterior probabilities for each of the
  rows in a `countData' object belonging to each
  of the models specified in the `groups' slot.
}
\usage{
getLikelihoods.NB(cD, prs, pET = "BIC", marginalise = FALSE, subset = NULL,
priorSubset = NULL, bootStraps = 1, conv = 1e-4, nullData = FALSE,
returnAll = FALSE, returnPD = FALSE, verbose = TRUE, discardSampling =
FALSE, cl, ...)
getLikelihoods.BB(cD, prs, pET = "BIC", marginalise = FALSE, subset =
NULL, priorSubset = NULL, bootStraps = 1, conv = 1e-04, nullData = FALSE,
returnAll = FALSE, returnPD = FALSE, verbose = TRUE, discardSampling =
FALSE, cl, ...)
getLikelihoods(cD, prs, pET = "BIC", marginalise = FALSE, subset = NULL,
priorSubset = NULL, bootStraps = 1, bsNullOnly = TRUE, conv = 1e-4, nullData = FALSE,
 weightByLocLikelihoods = TRUE, modelPriorSets = list(),
modelPriorValues = list(), returnAll = FALSE, returnPD = FALSE, verbose
= TRUE, discardSampling = FALSE, modelLikes = TRUE, cl = NULL, tempFile
= NULL, largeness = 1e+08)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cD}{An object of type \code{\link{countData}}, or descending
    from this class.}
  \item{prs}{(Initial) prior probabilities for each of the groups in the
    `cD' object. Should sum to 1, unless nullData is TRUE, in which
    case it should sum to less than 1.}
  \item{pET}{What type of prior re-estimation should be
    attempted? Defaults to "BIC"; "none" and "iteratively" are also
    available.}
  \item{marginalise}{Should an attempt be made to numerically
    marginalise over a prior distribution iteratively estimated from the
    posterior distribution? Defaults to FALSE, as in general offers
    little performance gain and increases computational cost considerably.}
  \item{subset}{Numeric vector giving the subset of counts for which
    posterior likelihoods should be estimated.}
  \item{priorSubset}{Numeric vector giving the subset of counts which may be
    used to estimate prior probabilities on each of the groups. See Details.}
  \item{bootStraps}{How many iterations of bootstrapping should be used
    in the (re)estimation of priors in the negative binomial method.}
  \item{bsNullOnly}{If TRUE (default, bootstrap hyper-parameters based
    on the likelihood of the null model and its inverse only; otherwise,
    on the likelihood of all models.}
  \item{conv}{If not null, bootstrapping iterations will cease if the
    mean squared difference between posterior likelihoods of consecutive
    bootstraps drops below this value.}
  \item{nullData}{If TRUE, looks for segments or counts with no true
    expression. See Details.}
  \item{weightByLocLikelihoods}{If a locLikelihoods slot is present in
    the `cD' object, and nullData = TRUE, then the initial weighting on
    nulls will be determined from the locLikelihoods slot. Defaults to TRUE.}
  \item{modelPriorSets}{If given, a list object, which defines subsets of
    the data for which different priors on the different models might be
    expected. See Details.}
  \item{modelPriorValues}{If given, a list object which defines priors
    on the different models. See Details.}
  \item{returnAll}{If TRUE, and bootStraps > 1, then instead of returning a
    single countData object, the function returns a list of countData
    objects; one for each bootstrap. Largely used for debugging
    purposes.}
  \item{returnPD}{If TRUE, then the function returns the (log)
    likelihoods of the data given the models, rather than the posterior
    (log) likelihoods of the models given the data. Not recommended for
    general use.}
  \item{verbose}{Should status messages be displayed? Defaults to TRUE.}
  \item{discardSampling}{If TRUE, discards information about which data
    rows are sampled to generate prior information. May slightly degrade
    the results but reduce computational time required. Defaults to
    FALSE.}
  \item{modelLikes}{If TRUE (default), returns likelihoods for each
    model. If FALSE, returns likelihoods for each hyper-parameter, from
    which the posterior joint distribution on hyper-parameters can be
    inferred.} 
  \item{cl}{A SNOW cluster object.}
  \item{tempFile}{Temporary file prefix for saving data
    likelihoods. Primarily for debugging purposes at this
    stage. Defaults to NULL, in which case no temporary data are saved.}
  \item{largeness}{The maximum size over which data likelihoods are
    calculated. Objects larger than this are split. This is most useful
    in combination with the saving of temporary files in the case of
    excessively large analyses.}
  \item{...}{Any additional information to be passed to the
    \code{'getLikelihoods'} function by the now deprecated functions.}
}

\details{
  These functions estimate, under the assumption of various
  distributions, the (log) posterior likelihoods that each count belongs
  to a group defined by the \code{@group} slot of the input
  object. The posterior likelihoods are stored on the natural log scale
  in the \code{@posteriors} slot of the \code{\link{countData}} object
  generated by this function. This is because the posterior likelihoods
  are calculated in this form, and ordering of the counts is better done
  on these log-likelihoods than on the likelihoods.  
  
  If \code{'pET = "none"'} then no attempt is made to re-estimate the
  prior likelihoods given in the \code{'prs'} variable. However, if
  \code{'pET = "BIC"'}, then the function will attempt to estimate the
  prior likelihoods by using the Bayesian Information Criterion to
  identify the proportion of the data best explained by each 
  model and taking these proportions as prior. Alternatively, an
  iterative re-estimation of priors is possible (\code{'pET = "iteratively"'}),
  in which an inital estimate for the prior likelihoods of the models is
  used to calculated the posteriors and then the priors are updated by
  taking the mean of the posterior likelihoods for each model across all
  data. This often works well, particularly if the 'BIC' method is used
  (see Hardcastle & Kelly 2010 for details). However, if the data are
  sufficiently non-independent, this approach may substantially
  mis-estimate the true priors. If it is possible to select a
  representative subset of the data by setting the variable
  \code{'subsetPriors'} that is sufficiently independent, then better
  estimates may be acquired. 

  In certain circumstances, it may be expected that certain subsets of the
  data are likely to behave differently to others; for example, if a set
  of genes are expected in advance to be differentially expressed, while
  the majority of the data are not. In this case, it may be advantageous
  (in terms of improving false discovery rates) to specify these different
  subsets in the modelPriorSets variable. However, care should be taken
  downstream to avoid confirmation bias.

  Filtering the data may be extremely advantageous in reducing run
  time. This can be done by passing a numeric vector to 'subset'
  defining a subset of the data for which posterior likelihoods are
  required.

  See Hardcastle & Kelly (2010) for a definition of the negative
  binomial methods.

  A 'cluster' object is strongly recommended in order to parallelise
  the estimation of posterior likelihoods, particularly for the
  negative binomial method. However, passing NULL to the \code{cl}
  variable will allow the functions to run in non-parallel mode.

  The `getLikelihoods.NB' and `getLikelihoods.BB' functions are now
  deprecated and will soon be removed.
}
\value{
  A \code{\link{countData}} object.
}
\references{Hardcastle T.J., and Kelly, K. baySeq: Empirical Bayesian
  Methods For Identifying Differential Expression In Sequence Count
  Data. BMC Bioinformatics (2010)}
\author{Thomas J. Hardcastle}


\seealso{\code{\link{countData}}, \code{\link{getPriors}},
  \code{\link{topCounts}}, \code{\link{getTPs}}}
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

# set negative binomial density function
densityFunction(CD) <- nbinomDensity

#estimate library sizes for countData object
libsizes(CD) <- getLibsizes(CD)

# Get priors for negative binomial method
\dontrun{CDPriors <- getPriors(CD, samplesize = 10^5, estimation = "QL", cl = cl)}

# To speed up the processing of this example, we have already created
# the `CDPriors' object.
data(CDPriors)

# Get likelihoods for data with negative binomial method.

CDPost <- getLikelihoods(CDPriors, pET = "BIC", cl = cl)

try(stopCluster(cl))

}

\keyword{distribution}
\keyword{models}
