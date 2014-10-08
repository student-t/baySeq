#
#.logRowSum <- function(z)
#  {
#    maxes <- do.call(pmax, c(as.list(data.frame(z)), list(na.rm = TRUE)))
#    pmax(maxes, maxes + log(rowSums(exp(z - maxes), na.rm = TRUE)), na.rm = TRUE)
#  }
#
`getPriors.BB` <-
function (cD, samplesize = 1e5, samplingSubset = NULL, verbose = TRUE,
          #zeroDispersion = FALSE, monoModal = FALSE,
          cl, ...)
{

  warning("This function is deprecated; using 'getPriors' instead.")
  densityFunction(cD) <- bbDensity
  return(getPriors(cD, samplesize = samplesize, samplingSubset = samplingSubset, verbose = TRUE, cl = cl))
  
  stop("You shouldn't ever see this message.")
}

#  
#  monoModal = FALSE
#  zeroDispersion = FALSE
#  dbetabinom <- function(x, n, prop, disp, monoModal, log = TRUE) {
#
#    ps <- rep(NA, (length(x)))
#    disp <- rep(disp, length(prop))    
#
#    smallDisp <- disp < 1e-15 & disp >= 0
#    largeDisp <- disp > 1e-15
#    
#    if(any(disp < 0 | disp >= 1 | prop > 1 | prop < 0)) ps[disp < 0 | disp >= 1 | prop > 1 | prop < 0] <- NA
#    if(any(smallDisp))
#      ps[smallDisp] <- dbinom(x[smallDisp], n[smallDisp], prob = prop[smallDisp], log = log)
#
#    if(any(largeDisp)) {
#      alpha = (1/disp - 1) * prop
#      beta = (1/disp - 1) * (1-prop)
#
#      if(monoModal)
#        if (any(alpha < 1) | any(beta < 1)) return(NA)
#
#      if(!log) {
#        ps[largeDisp] <- choose(n[largeDisp], x[largeDisp]) * beta((x + alpha)[largeDisp], (n - x + beta)[largeDisp]) / beta(alpha[largeDisp], beta[largeDisp])
#      } else ps[largeDisp] <- lchoose(n[largeDisp], x[largeDisp]) + lbeta((x + alpha)[largeDisp], (n - x + beta)[largeDisp]) - lbeta(alpha[largeDisp], beta[largeDisp])
#    }
#    return(ps)
#  }
#  
#  if(!inherits(cD, what = "countData") && !inherits(cD, what = "methData"))
#    stop("variable 'cD' must be of or descend from class 'countData' or 'methData'")
#
#  if(verbose) message("Finding priors...", appendLF = FALSE)
#
#  getDisps <- function(selRow, replicates, moderate, zeroDispersion, monoModal)
#    {
#      nzRep <- levels(replicates)[sapply(levels(replicates), function(rep) any(y[selRow, replicates == rep] != 0) & any(secondy[selRow, replicates == rep] != 0))]
#      nzreplicates <- replicates[replicates %in% nzRep]
#      if(length(nzreplicates) == 0) return(NA)
#      
#      propDispReplicates <- function(propDisp, z, secondz, subLibsizes, subPairLibsizes, monoModal)
#        {
#          props <- propDisp[match(nzreplicates, nzRep)]          
#          props <- props * subLibsizes / (props * subLibsizes + (1-props) * subPairLibsizes)
#
#          if(any(props < 0, na.rm = TRUE)| any(props > 1, na.rm = TRUE)) return(NA)          
#          
#          disp <- propDisp[length(propDisp)]
#
#          if(disp < 0 | disp > 1) return(NA)
#          
#          sum(dbetabinom(z, z + secondz, props, disp = disp, monoModal = monoModal))
#        }
#      
#      z <- y[selRow,]
#      secondz <- secondy[selRow,]
#
#      if(all(z == 0) | all(secondz == 0))
#        if(zeroDispersion) return(0) else return(NA)
#      
#      if(moderate)
#        if(all(secondz[z != 0] == 0) & all(z[secondz != 0] == 0) & any(z != 0) & any(secondz != 0))
#          {              
#            maxInZ <- max(z[secondz == 0] / libsizes[secondz == 0])
#            maxInSecondZ <- max(secondz[z == 0] / pairLibsizes[z == 0])
#            if(maxInZ >= maxInSecondZ) {
#              secondz[which(secondz == 0)[which.max(z[secondz == 0] / libsizes[secondz == 0])]] <- 1
#            } else {
#              z[which(z == 0)[which.max(secondz[z == 0] / pairLibsizes[z == 0])]] <- 1
#            }
#          }                
#      pars <- optim(par = c(rep(0.5, length(nzRep)), 0.01), propDispReplicates, z = z[replicates %in% nzRep], secondz = secondz[replicates %in% nzRep], 
#                    subLibsizes = libsizes[replicates %in% nzRep], subPairLibsizes = pairLibsizes[replicates %in% nzRep], monoModal = monoModal,
#                    control = list(fnscale = -1, trace = FALSE, reltol = 1e-50, maxit = 1000))$par
#        disp <- pars[length(pars)]
#      disp
#    }         
#
#  optimoverPriors <- function(selRow, selcts)
#    {
#      propOnly <- function(prop, disp, z, secondz, libz, pairLibz, nonconv)
#        {
#          props <- rep(prop, length(z))
#          if(all(is.na(nonconv))) props <- props * libz / (props * libz + (1-props) * pairLibz) else {
#            k = nonconv / (1 - nonconv) * secondz
#            props = (props * (z + secondz) + k) / (z + secondz)
#          }              
#          sum(dbetabinom(z, z + secondz, props, disp = disp, monoModal = FALSE))
#        }
#      
#      
#      pars <- c(optimize(propOnly, interval = c(0, 1),
#                         disp = disps[selRow],
#                         z = y[selRow,selcts], secondz = secondy[selRow,selcts],
#                         libz = libsizes[selcts], pairLibz = pairLibsizes[selcts], nonconv = nonconversion[selcts], tol = 1e-50, maximum = TRUE)$maximum, disps[selRow])
#
#    }
#
#  
#  if(!is.null(cl))
#    {
#      getPriorsEnv <- new.env(parent = .GlobalEnv)
#      environment(optimoverPriors) <- getPriorsEnv
#    }
#  
#  
#  if(is.null(samplingSubset))
#    samplingSubset <- 1:nrow(cD)
#
#  samplingSubset <- samplingSubset[rowSums(is.na(cD@data[samplingSubset,,drop = FALSE])) == 0]
#  samplingSubset <- samplingSubset[rowSums(cD@data[samplingSubset,,drop = FALSE]) > 0 | rowSums(cD@pairData[samplingSubset,,drop = FALSE]) > 0]
#
#  if(any(sapply(cD@groups, class) != "factor"))
#    {
#      cD@groups <- lapply(cD@groups, as.factor)
#      warning("Not all members of the '@groups' slot were factors; converting now.")
#    }
#  if(class(cD@replicates) != "factor")
#    {
#      cD@replicates <- as.factor(cD@replicates)
#      warning("The '@replicates' slot is not a factor; converting now.")
#    }      
#  
#  sD <- cD[samplingSubset,]
#
#  if(inherits(cD, what = "countData")) {
#    libsizes <- as.double(sD@libsizes)
#    pairLibsizes <- as.double(sD@pairLibsizes)
#    if(length(cD@libsizes) == 0) libsizes <- rep(1, ncol(sD))
#    if(length(cD@pairLibsizes) == 0) pairLibsizes <- libsizes
#    nonconversion <- rep(NA, ncol(sD))
#    y <- round(sD@data)
#    secondy <- round(sD@pairData)
#  } else if(inherits(cD, what = "countData")) {
#    pairLibsizes <- libsizes <- rep(NA, ncol(sD))
#    nonconversion <- as.double(sD@nonconversion)
#    if(length(nonconversion) == 0) nonconversion = rep(0, ncol(sD))
#    y <- round(sD@Cs)
#    secondy <- round(sD@Ts)
#  }
#      
#  groups <- sD@groups
#  replicates <- as.factor(sD@replicates)
#
#  if(nrow(sD@seglens) > 0) seglens <- sD@seglens[,,drop = TRUE] else seglens <- rep(1, nrow(sD@data))
#  if(is.vector(seglens)) lensameFlag <- TRUE else lensameFlag <- FALSE
#
#  #tupData <- log(sD@data) - log(sD@pairData)
#  
#  #carve out some stratified sampling here...
#  #also account for duplicates
#  
#  selrow <- which(rowSums(y != 0) > 0 | rowSums(secondy != 0) > 0)
#
#  if(length(selrow) > samplesize)
#    selrow <- sample(selrow, size = samplesize)
#  
#  disps <- c()
#
#  fixDispNA <- function(disps, selrow) {
#    disps[selrow[is.na(disps[selrow])]] <- sapply(selrow[is.na(disps[selrow])], function(dispna) {
#      meanNA <- mean(y[dispna,] / (y[dispna,] + secondy[dispna,]), na.rm = TRUE)
#      meanRat <- rowMeans(y / (y + secondy), na.rm = TRUE)
#      sample(disps[!is.na(disps)], 1, prob = dexp(abs(meanNA - meanRat)[!is.na(disps)], rate = 100))
#    })
#    disps
#  }
#  
#  if (!is.null(cl)) {
#    clustAssign <- function(object, name)
#      {
#        assign(name, object, envir = .GlobalEnv)
#        NULL
#      }
#    getPriorEnv <- new.env(parent = .GlobalEnv)
#    environment(clustAssign) <- getPriorEnv
#    environment(getDisps) <- getPriorEnv
#    environment(optimoverPriors) <- getPriorEnv
#    clusterCall(cl, clustAssign, y, "y")
#    clusterCall(cl, clustAssign, secondy, "secondy")
#    clusterCall(cl, clustAssign, libsizes, "libsizes")
#    clusterCall(cl, clustAssign, pairLibsizes, "pairLibsizes")
#    clusterCall(cl, clustAssign, nonconversion, "nonconversion")
#    clusterCall(cl, clustAssign, dbetabinom, "dbetabinom")
#
#    disps[selrow] <- parSapply(cl, selrow, getDisps, replicates = sD@replicates, moderate = moderate, zeroDispersion = zeroDispersion, monoModal = monoModal)
#    if(any(is.na(disps[selrow]))) ndisps <- fixDispNA(disps, selrow)
#    
#    clusterCall(cl, clustAssign, disps, "disps")
#    BBpar <- lapply(sD@groups, function(group)
#                       lapply(unique(levels(group)), function(uu){
#                         t(parSapply(cl, selrow, optimoverPriors, selcts = group == uu))
#                    }))    
#  } else {
#    disps[selrow] <-
#      sapply(selrow, getDisps, replicates = sD@replicates, moderate = moderate, zeroDispersion = zeroDispersion, monoModal = monoModal)
#    if(any(is.na(disps[selrow]))) disps <- fixDispNA(disps, selrow)
#    
#    BBpar <- lapply(sD@groups, function(group)
#                  lapply(unique(levels(group)), function(uu){
#                    t(sapply(selrow, optimoverPriors, selcts = group == uu))
#                  }))
#  }
#
#  sampled <- cbind(sampled = selrow, representative = 1:length(selrow))
#  weights <- rep(1, length(selrow))
#
#  if(verbose) message("done.")
#  
#  new(class(cD), cD, priorType = "BB", priors = list(sampled = sampled, weights = weights, priors = BBpar))
#}
#
#
#
`getLikelihoods.BB` <- function(cD, prs, pET = "BIC", marginalise = FALSE, subset = NULL, priorSubset = NULL, bootStraps = 1, conv = 1e-4, nullData = FALSE, returnAll = FALSE, returnPD = FALSE, verbose = TRUE, discardSampling = FALSE, cl, ...)
  {
    warning("This function is deprecated; using 'getLikelihoods' instead.")
    if(!missing(prs)) return(getLikelihoods(cD, prs = prs, pET = pET, marginalise = marginalise, subset = subset, priorSubset = priorSubset, bootStraps = bootStraps, conv = conv, nullData = nullData, returnAll = returnAll, returnPD = returnPD, verbose = verbose, discardSampling = discardSampling, cl = cl, ...)) else return(getLikelihoods(cD, pET = pET, marginalise = marginalise, subset = subset, priorSubset = priorSubset, bootStraps = bootStraps, conv = conv, nullData = nullData, returnAll = returnAll, returnPD = returnPD, verbose = verbose, discardSampling = discardSampling, cl = cl, ...))
    
    stop("You shouldn't ever see this message.")
  }    
#    constructWeights <- function(withinCluster = FALSE)
#      {
#        priorWeights <- lapply(1:length(BBpriors), function(ii)
#                               lapply(1:length(BBpriors[[ii]]), function(jj)
#                                      {                                      
#                                        samplings <- numintSamp[[ii]][[jj]]
#                                        samplings <- samplings[order(samplings[,2]),]
#                                        intervals <- findInterval(1:nrow(BBpriors[[ii]][[jj]]) - 0.5, samplings[,2]) + 1
#                                        intlens <- diff(c(intervals, nrow(samplings) + 1))                                                                       
#                                        weights <- c()
#                                        weights[intlens == 1] <- samplings[intervals[intlens == 1],3]
#                                        if(any(intlens > 1))
#                                          weights[which(intlens > 1)] <- sapply(which(intlens > 1), function(kk)
#                                                                                sum(samplings[intervals[kk] + 0:(intlens[kk] - 1),3]))
#                                        weights
#                                      }))
#                                                                     
#        if(withinCluster) {
#          assign("priorWeights", priorWeights, envir = .GlobalEnv)
#          return(invisible(NULL))
#        } else return(priorWeights)
#      }
#    
#    BBbootStrap <- function(row, groups)
#      {
#        `logsum` <-
#          function(x)
#            max(x, max(x, na.rm = TRUE) + log(sum(exp(x - max(x, na.rm = TRUE)), na.rm = TRUE)), na.rm = TRUE)
#        
#        PDgivenr.BB <- function (number, cts, secondcts, priors, group, wts, sampInfo)
#          {
#            dbetabinom <- function(x, n, prop, disp) {
#              
#              smallDisp <- disp < 1e-15 & disp >= 0
#              largeDisp <- disp > 1e-15
#
#              ps <- matrix(NA, ncol = ncol(prop), nrow = nrow(prop))
#              disp <- matrix(disp, ncol = ncol(prop), nrow = nrow(prop))
#              x <- matrix(x, ncol = ncol(prop), nrow = nrow(prop), byrow = TRUE)
#              n <- matrix(n, ncol = ncol(prop), nrow = nrow(prop), byrow = TRUE)
#              
#              if(any(largeDisp)) {
#
#                alpha = (1/disp - 1) * prop
#                beta = (1/disp - 1) * (1-prop)
#                
#                ps[largeDisp,] <- lchoose(n[largeDisp,,drop = FALSE], x[largeDisp,,drop = FALSE]) +
#                  lbeta((x + alpha)[largeDisp,,drop = FALSE], (n - x + beta)[largeDisp,, drop = FALSE]) -
#                    lbeta(alpha[largeDisp,, drop = FALSE], beta[largeDisp,, drop = FALSE]) 
#              }
#              if(any(smallDisp))
#                ps[smallDisp,] <- dbinom(x[smallDisp,,drop = FALSE], n[smallDisp,,drop = FALSE], prob = prop[smallDisp,,drop = FALSE], log = TRUE)
#              
#              return(ps)
#            }
#                                   
#            sum(
#                sapply(1:length(levels(group)), function(ll) {
#                  selcts <- group == levels(group)[ll] & !is.na(group)
#                  prior <- priors[[ll]]
#                  weightings <- wts[[ll]]
#                  nzWts <- weightings != 0 & !is.na(weightings)
#                  #wsInfo <- which(sampInfo[[ll]][,1] == number)
#                  #weightings[sampInfo[[ll]][wsInfo,2]] <- weightings[sampInfo[[ll]][wsInfo,2]] - sampInfo[[ll]][wsInfo,3]
#                  
#                  p = outer(prior[nzWts,1], libsizes[selcts]) / (outer(1-prior[nzWts,1], pairLibsizes[selcts]) + outer(prior[nzWts,1], libsizes[selcts]))
#                  disp = prior[nzWts,2]
#                  ps <- dbetabinom(cts[selcts], cts[selcts] + secondcts[selcts], prop = p, disp = disp)
#                  
#                  logsum(
#                         rowSums(ps) + log(weightings[nzWts])                                 
#                         )  - log(sum(weightings[nzWts]))
#                })
#                )
#          }
#        
#        sapply(1:length(BBpriors), function(gg)
#               PDgivenr.BB(number = row, cts = y[row,], secondcts = secondy[row,], priors = BBpriors[[gg]], group = groups[[gg]], wts = priorWeights[[gg]], sampInfo = numintSamp[[gg]])                           
#               )
#      }
#    
#    if(!is.null(cl))
#      {
#        clustAssign <- function(object, name)
#          {
#            assign(name, object, envir = .GlobalEnv)
#            NULL
#          }
#        
#        getLikelihoodsEnv <- new.env(parent = .GlobalEnv)
#        environment(clustAssign) <- getLikelihoodsEnv
#        environment(BBbootStrap) <- getLikelihoodsEnv
#      }
#
#    
#    listPosts <- list()
#    
#    if(!inherits(cD, what = "countData") && !inherits(cD, what = "methData"))
#      stop("variable 'cD' must be of or descend from class 'countData' or 'methData'")
#    
#    if(cD@priorType != "BB") stop("Incorrect prior type for this method of likelihood estimation")
#    
#    if(pET %in% c("none", "iteratively"))
#      {
#        if(length(prs) != length(cD@groups)) stop("'prs' must be of same length as the number of groups in the 'cD' object")
#        if(any(prs < 0))
#          stop("Negative values in the 'prs' vector are not permitted")
#        if(missing(nullProps) & sum(prs) != 1)
#          stop("If no values given for 'nullProps' then the 'prs' vector should sum to 1.")
#        
#        if(!missing(nullProps) & sum(prs) >= 1)
#          stop("If some value given for 'nullProps' then the 'prs' vector should sum to less than 1.")
#        
#      } else if(pET %in% c("BIC"))
#        prs <- rep(NA, length(cD@groups))
#    
#    if(!(class(subset) == "integer" | class(subset) == "numeric" | is.null(subset)))
#      stop("'subset' must be integer, numeric, or NULL")
#    
#    if(is.null(conv)) conv <- 0
#
#    groups <- cD@groups
#    BBpriors <- cD@priors$priors
#    numintSamp <- cD@priors$sampled
#    y <- round(cD@data)
#    secondy <- round(cD@pairData)
#    libsizes <- cD@libsizes
#    pairLibsizes <- cD@pairLibsizes
#
#    if(is.null(subset)) subset <- 1:nrow(cD)
#    subset <- subset[(rowSums(y[subset,,drop = FALSE]) != 0 | rowSums(secondy[subset,,drop = FALSE]) != 0)]
#    
#    if(is.null(priorSubset)) priorSubset <- subset
#
#
#    if(length(libsizes) == 0) libsizes <- rep(1, ncol(cD))
#    if(length(pairLibsizes) == 0) pairLibsizes <- libsizes
#    
#    if(!missing(nullProps)) {
#
#      ndelocGroup <- which(unlist(lapply(cD@groups, function(x) all(x[!is.na(x)] == x[!is.na(x)][1])))) 
#      if(length(ndelocGroup) == 0) stop("If 'nullProps' is given then there must exist some vector in groups whose members are all identical") else ndelocGroup <- ndelocGroup[1]      
#      
#      nulpPriors <- lapply(nullProps, function(x) list(cbind(x, BBpriors[[ndelocGroup]][[1]][,2])))
#      nulpGroups <- lapply(nullProps, function(x) as.factor(rep(1, ncol(y))))
#      names(nulpGroups) <- sapply(nullProps, function(x) paste("nde_", x, sep = ""))
#      groups <- c(nulpGroups, cD@groups)
#      BBpriors <- c(nulpPriors, BBpriors)      
#    }
#        
#    if(discardSampling) numintSamp[,1] <- NA
#
#    if(is.matrix(numintSamp))
#      numintSamp <- lapply(BBpriors, function(x) lapply(x, function(z) numintSamp))
#    if(is.null(numintSamp))
#      numintSamp <- lapply(BBpriors, function(x) lapply(x, function(z) cbind(sampled = rep(-1, nrow(z)), representative = 1:nrow(z))))
#
#    weights <- cD@priors$weights
#    if(is.null(weights))
#      weights <- lapply(numintSamp, function(x) lapply(x, function(z) weights = rep(1, nrow(z))))
#    if(is.numeric(weights))
#      weights <- lapply(numintSamp, function(x) lapply(x, function(z) weights = weights))
#
#    numintSamp <- lapply(1:length(numintSamp), function(ii) lapply(1:length(numintSamp[[ii]]), function(jj) cbind(numintSamp[[ii]][[jj]], weights = weights[[ii]][[jj]])))
#    priorWeights <- constructWeights()
#    
#    posteriors <- matrix(NA, ncol = length(groups), nrow = nrow(cD@data))
#    propest <- NULL
#    converged <- FALSE
#
#    if(!is.null(cl))
#      {
#        clusterCall(cl, clustAssign, BBpriors, "BBpriors")
#        clusterCall(cl, clustAssign, libsizes, "libsizes")
#        clusterCall(cl, clustAssign, pairLibsizes, "pairLibsizes")
#      }  
#
#    if(verbose) message("Finding posterior likelihoods...", appendLF = FALSE)
#
#    priorReps <- unique(unlist(sapply(numintSamp, function(x) as.integer(unique(sapply(x, function(z) z[,1]))))))
#    priorReps <- priorReps[priorReps > 0 & !is.na(priorReps)]
#
#    if(!all(priorReps %in% 1:nrow(cD@data)) & bootStraps > 1)
#      {
#        warning("Since the sampled values in the '@priors' slot are not available, bootstrapping is not possible.")
#        bootStraps <- 1
#      }
#    
#    postRows <- unique(c(priorReps, priorSubset, subset))
#
#    .fastUniques <- function(x){
#      if (nrow(x) > 1) {
#        return(c(TRUE, rowSums(x[-1L, , drop = FALSE] == x[-nrow(x),, drop = FALSE]) != ncol(x)))
#      } else return(TRUE)
#    }    
#
##    if(length(priorReps) == 0 || !any(priorReps %in% postRows))
##      {
##        orddat <- do.call("order", c(lapply(1:ncol(seglens), function(ii) seglens[postRows,ii]), lapply(1:ncol(cD@data), function(ii) cD@data[postRows,ii])))
##        whunq <- .fastUniques(cbind(seglens[postRows,], cD@data[postRows,])[orddat,,drop = FALSE])
##        postRows <- postRows
##      } else whunq <- TRUE
#
##    orddat <- 1:nrow(cD)
#
#    whunq <- TRUE
#
#    for(cc in 1:bootStraps)
#      {
#        if(cc > 1) numintSamp <- lapply(1:length(numintSamp), function(ii) lapply(1:length(numintSamp[[ii]]), function(jj) cbind(numintSamp[[ii]][[jj]][,1:2], weights = exp(posteriors[numintSamp[[ii]][[jj]][,1],ii]))))
#
#        if (is.null(cl)) {
#          if(cc > 1) priorWeights <- constructWeights()            
#          ps <- sapply(postRows[whunq], BBbootStrap, groups = groups)
#        } else {
#          environment(constructWeights) <- getLikelihoodsEnv
#          clusterCall(cl, clustAssign, numintSamp, "numintSamp")
#          clusterCall(cl, constructWeights, TRUE)
#          clusterCall(cl, clustAssign, y, "y")
#          clusterCall(cl, clustAssign, secondy, "secondy")
#
#          ps <- parSapply(cl, postRows[whunq], BBbootStrap, groups = groups)
#        }
#        
#        ps <- matrix(ps, ncol = length(groups), byrow = TRUE)
#        rps <- matrix(NA, ncol = length(groups), nrow = nrow(cD@data))
#        
#        rps[postRows,] <- ps
#
#        if(returnPD) {
#              if(verbose) message("done.")
#              return(rps)
#            }
#
#        if(pET != "none")
#          {
#            restprs <- getPosteriors(rps[priorSubset,, drop = FALSE], prs, pET = pET, marginalise = FALSE, groups = groups, priorSubset = NULL, cl = cl)$priors
#          } else restprs <- prs
#        
#        pps <- getPosteriors(rps[union(priorReps,subset),], prs = restprs, pET = "none", marginalise = marginalise, groups = groups, priorSubset = NULL, cl = cl)
#        
#        if(any(!is.na(posteriors)))
#          if(all(abs(exp(posteriors[union(priorReps,subset),]) - exp(pps$posteriors)) < conv)) converged <- TRUE
#        
#        posteriors[union(priorReps, subset),] <- pps$posteriors
#
#        prs <- pps$priors
#        propest <- rbind(propest, prs)
#
#        estProps <- pps$priors
#        names(estProps) <- names(cD@groups)
#
#        cat(".")        
#
#        if(returnAll | converged | cc == bootStraps)
#          {
#            retPosts <- posteriors
#            retPosts[priorReps[!(priorReps %in% subset)],] <- NA
#            
#            nullPosts <- matrix(nrow = 0, ncol = 0)
#            if(!missing(nullProps)) {
#              nullPosts <- retPosts[,1:length(nullProps), drop = FALSE]
#              colnames(nullPosts) <- sapply(nullProps, function(x) paste("nde_", x, sep = ""))
#              retPosts <- retPosts[,-(1:length(nullProps)), drop = FALSE]
#            }
#            estProps <- apply(exp(retPosts), 2, mean, na.rm = TRUE)
#            
#            colnames(retPosts) <- names(cD@groups)
#            listPosts[[cc]] <- (new(class(cD), cD, posteriors = retPosts, estProps = estProps, nullPosts = nullPosts))
#            listPosts[[cc]]@orderings <- makeOrderings(listPosts[[cc]])
#          }
#
#        if(converged)
#          break()
#      }
#
#    if(!is.null(cl))
#      clusterEvalQ(cl, rm(list = ls()))
#
#    if(verbose) message("done.")
#
#    if(!returnAll) return(listPosts[[cc]]) else {
#      if(length(listPosts) == 1) return(listPosts[[1]]) else return(listPosts)
#    }
#    
#  }
#
