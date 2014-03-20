.optimoverPriorsLN <- function(x, estimation, replicates, groups, libsizes, equalDispersions, lensameFlag, zeroML, consensus, disp)
    {
      x <- x[-1]
      likelnDisp <- function(vars, x, replicates, libsizes) {
        disp <- rep(vars[1], length(x)) 
        means <- rep(vars[-1], sapply(levels(replicates), function(rep) sum(replicates == rep)))
        if(any(disp < 0)) return(NA)
        sum(dnorm(x, mean = means * libsizes, sd = disp, log = TRUE))
      }
      likelnMeans <- function(mean, disp, x, libsizes) {
        if(any(disp < 0)) return(NA)
        sum(dnorm(x, mean = mean * libsizes, sd = disp, log = TRUE))
      }
      if(equalDispersions) {
        disp <- optim(par = c(0.1, 1, 1), fn = likelnDisp, x = x, replicates = replicates, control = list(fnscale = -1), libsizes = libsizes)$par[1]
        groupness <- lapply(groups, function(group)
                            lapply(levels(group), function(gg)                                   
                                   c(disp, optimise(likelnMeans, interval = c(min(x / libsizes) - 5, max(x / libsizes) + 5), x = x[group == gg], disp = disp, libsizes = libsizes[group == gg], maximum = TRUE)[[1]]))
                            )
      } else {
        groupness <- lapply(groups, function(group)
                            lapply(levels(group), function(gg) {
                              disp <- optim(par = c(0.1, rep(1, length(unique(replicates[group == gg])))), fn = likelnDisp, x = x[group == gg], replicates = replicates[group == gg], libsizes = libsizes[group == gg], control = list(fnscale = -1))$par[1]
                              c(disp, optimise(likelnMeans, interval = c(min(x / libsizes) - 1, max(x / libsizes) + 1), x = x[group == gg], disp = disp, maximum = TRUE)[[1]])
                            }))
      }
      groupness
    }


.optimoverPriorsGamma <- function(x, estimation, replicates, groups, libsizes, equalDispersions, lensameFlag, zeroML, consensus, disp)
    {
      x <- x[-1]
      print(x)
      likelnShape <- function(vars, x, replicates, libsizes) {
        shape <- rep(vars[1], length(x)) 
        scale <- rep(vars[-1], sapply(levels(replicates), function(rep) sum(replicates == rep)))
        if(any(vars < 0)) return(NA)
        sum(dgamma(x, shape = shape, scale = scale * libsizes, log = TRUE))
      }
      likelnMeans <- function(scale, shape, x, libsizes) {
        if(any(shape < 0) | any(scale < 0)) return(NA)
        sum(dgamma(x, shape = shape, scale = scale * libsizes, log = TRUE))
      }
      if(equalDispersions) {
        shape <- optim(par = c(1, 1, 1), fn = likelnShape, x = x, replicates = replicates, control = list(fnscale = -1), libsizes = libsizes)$par[1]
        groupness <- lapply(groups, function(group)
                            lapply(levels(group), function(gg)                                   
                                   c(shape, optimise(likelnMeans, interval = c(0, max(x) * 10 + 1), x = x[group == gg], shape = shape, libsizes = libsizes[group == gg], maximum = TRUE)[[1]]))
                            )
      } else {
        groupness <- lapply(groups, function(group)
                            lapply(levels(group), function(gg) {
                              shape <- optim(par = c(1, rep(1, length(unique(replicates[group == gg])))), fn = likelnShape, x = x[group == gg], replicates = replicates[group == gg], control = list(fnscale = -1))$par[1]
                              c(shape, optimise(likelnMeans, interval = c(0, max(x) * 10 + 1), x = x[group == gg], shape = shape, maximum = TRUE)[[1]])
                            }))
      }
      groupness
    }

.optimoverPriorsEXP <- function(x, estimation, replicates, groups, libsizes, equalDispersions, lensameFlag, zeroML, consensus, disp)
    {
      x <- x[-1]
      likelnExp <- function(vars, x, libsizes) {
        mean <- rep(vars[1], length(x)) 
        if(any(mean <= 0)) return(NA)
        z <- sum(dexp(x, rate = 1 / (mean * libsizes), log = TRUE))
        z
      }
      groupness <- lapply(groups, function(group)
                          lapply(levels(group), function(gg)                                   
                                 c(optimise(likelnExp, interval = c(0, 1 + max(x) * 5), x = x[group == gg], libsizes = libsizes[group == gg], maximum = TRUE)[[1]])))
                          
      groupness
    }



`getPriors` <-
function (cD, distribution = c("normal", "exponential", "negative-binomial", "gamma"), samplesize = 1e5, samplingSubset = NULL, equalDispersions = TRUE, estimation = "QL", verbose = TRUE, zeroML = FALSE, consensus = FALSE, cl, ...)
{
  distributions <- c("normal", "exponential", "negative-binomial", "gamma")
  distribution <- distributions[pmatch(distribution, distributions)]
  if(distribution == "negative-binomial")
    return(getPriors.NB(cD = cD, samplesize = samplesize, samplingSubset = samplingSubset, equalDispersions = equalDispersions, estimation = estimation, verbose = verbose, zeroML = zeroML, consensus = consensus, cl = cl))
  if(distribution == "normal") {
    optimFunction <- .optimoverPriorsLN
    estimation = "ML"
  }
  if(distribution == "exponential") {
    optimFunction <- .optimoverPriorsEXP
    estimation = "ML"
  }
  if(distribution == "gamma") {
    optimFunction <- .optimoverPriorsGamma
    estimation = "ML"
  }
  
  if(!inherits(cD, what = "countData"))
    stop("variable 'cD' must be of or descend from class 'countData'")
  
  if(length(cD@libsizes) == 0)
    {
      warning("'@libsizes' slot empty; inferring libsizes using default settings")
      cD@libsizes <- getLibsizes(cD)
    }
    
  if(verbose) message("Finding priors...", appendLF = FALSE)          
  
  if(!is.null(cl))
    {
      getPriorsEnv <- new.env(parent = .GlobalEnv)
      environment(optimFunction) <- getPriorsEnv
    }
  
  
  if(any(sapply(cD@groups, class) != "factor"))
    {
      cD@groups <- lapply(cD@groups, as.factor)
      warning("Not all members of the '@groups' slot were factors; converting now.")
    }
  if(class(cD@replicates) != "factor")
    {
      cD@replicates <- as.factor(cD@replicates)
      warning("The '@replicates' slot is not a factor; converting now.")
    }

  if(is.null(samplingSubset))
    samplingSubset <- 1:nrow(cD)
  
  samplingSubset <- samplingSubset[rowSums(do.call("cbind", lapply(cD@groups, function(x) do.call("cbind", lapply(levels(x), function(rep) rowSums(is.na(cD@data[samplingSubset,which(x == rep),drop = FALSE])) == length(which(x == rep))))))) == 0]

  sD <- cD[samplingSubset,]
  
  libsizes <- as.double(sD@libsizes)
  groups <- sD@groups
  replicates <- as.factor(sD@replicates)

  if(nrow(sD@seglens) > 0) seglens <- sD@seglens[,,drop = TRUE] else seglens <- rep(1, nrow(sD@data))
  if(is.vector(seglens)) lensameFlag <- TRUE else lensameFlag <- FALSE

  tupData <- cbind(sD@data, seglens)
  
  ordData <- do.call(order, as.data.frame(tupData))
  if(length(ordData) > 1) {
    dups <- c(1, which(rowSums(tupData[ordData[-1],, drop = FALSE] == (tupData)[ordData[-length(ordData)],, drop = FALSE], na.rm = TRUE) != rowSums(!is.na(tupData[ordData[-1],, drop=FALSE]))) + 1)
  } else dups <- 1  
  
  if(length(dups) <= samplesize) {
    y <- sD@data[ordData[dups],,drop = FALSE]
    copies <- diff(c(dups, nrow(sD@data) + 1))
    sampled <- (1:nrow(sD))[ordData]
    sy <- cbind(sampled = sampled, representative = rep(1:length(copies), copies))
    weights <- rep(1, nrow(sy))
    if(lensameFlag) seglensy <- seglens[ordData[dups]] else seglensy <- seglens[ordData[dups],]
  } else {
    sampData <- colSums(t(sD@data/seglens) / sD@libsizes, na.rm = TRUE) / ncol(sD)

    if(zeroML) nzSD <- sampData[sampData > 0] else nzSD <- sampData
      
    sqnum <- length(nzSD) / 1000
    squant <- quantile(nzSD, 1:sqnum / sqnum, na.rm = TRUE)
    sqdup <- c(1, which(diff(squant) > min(1 / libsizes) / 10))
    z <- cbind(as.numeric(squant[sqdup]), c(as.numeric(squant[sqdup[-1]]), max(nzSD, na.rm = TRUE)))
    if(zeroML) z[1,1] <- 0 else z[1,1] <- -Inf

    sy <- do.call("rbind",
                  lapply(1:nrow(z), function(ii) {
                    w <- z[ii,]
                    inbetweener <- which(sampData > w[1] & sampData <= w[2])
                    samplenum <- min(length(inbetweener), ceiling(samplesize / nrow(z)), na.rm = TRUE)
                    cbind(sample(inbetweener, size = samplenum, replace = FALSE), length(inbetweener) / samplenum)})
                )

    weights <- sy[,2]
    sy <- sy[,1]

    y <- sD@data[sy,,drop = FALSE]
    if(lensameFlag) seglensy <- seglens[sy] else seglensy <- seglens[sy,,drop = FALSE]

    ordData <- do.call(order, as.data.frame(cbind(y, seglensy)))
    y <- y[ordData,, drop = FALSE]
    weights <- weights[ordData]
    if(lensameFlag) seglensy <- seglensy[ordData] else seglensy <- seglensy[ordData,,drop = FALSE]
    
    dups <- c(1, which(rowSums((cbind(y,seglensy))[-1,] == (cbind(y, seglensy))[-nrow(y),], na.rm = TRUE) != rowSums(is.na(cbind(y,seglensy)[-1,]))) + 1)
    copies <- diff(c(dups, nrow(y) + 1))

    sy <- cbind(sampled = sy[ordData], representative = rep(1:length(copies), copies))

    y <- y[dups,, drop = FALSE]
    if(lensameFlag) seglensy <- seglensy[dups] else seglensy <- seglensy[dups,, drop = FALSE]
    
    if(any(sampData == 0) & zeroML)
      {
        sy <- rbind(sy, cbind(sampled = which(sampData == 0), representative = nrow(y) + 1))
        y <- rbind(y, rep(0, ncol(y)))
        if(lensameFlag) seglensy <- c(seglensy, 1) else seglensy <- rbind(seglensy, rep(1, ncol(seglensy)))
        weights <- c(weights, rep(1, sum(sampData == 0)))
      }
  }

  z <- cbind(seglensy, y)  

  LNpar <- list()

  if(is.null(cl)) {
    parEach <- apply(z, 1, optimFunction, estimation = estimation, replicates = replicates, groups = groups, libsizes = libsizes, equalDispersions = equalDispersions, lensameFlag = lensameFlag, zeroML = zeroML)
  } else {
    parEach <- parApply(cl, z, 1, optimFunction, estimation = estimation, replicates = replicates, groups = groups, libsizes = libsizes, equalDispersions = equalDispersions, lensameFlag = lensameFlag, zeroML)
  }
  
  LNpar <- lapply(1:length(groups), function(gg)
                  lapply(1:length(levels(groups[[gg]])), function(ii)
                         do.call("rbind", lapply(parEach, function(x) x[[gg]][[ii]]))))

  if(verbose) message("done.")

  sy[,1] <- samplingSubset[sy[,1]]
  names(LNpar) <- names(groups)
  LNpar <- list(sampled = sy, weights = weights, priors = LNpar)
  new(class(cD), cD, priorType = distribution, priors = LNpar)
}

`getLikelihoods` <-
function(cD, distribution = c("normal", "exponential", "negative-binomial", "gamma"), prs, pET = "BIC", marginalise = FALSE, subset = NULL, priorSubset = NULL, bootStraps = 1, conv = 1e-4, nullData = FALSE, returnAll = FALSE, returnPD = FALSE, verbose = TRUE, discardSampling = FALSE, cl)
  {
    distributions <- c("normal", "exponential", "negative-binomial", "gamma")
    distribution <- distributions[pmatch(distribution, distributions)]
    if(distribution == "negative-binomial")
      return(getLikelihoods.NB(cD = cD, prs, pET = "BIC", marginalise = marginalise, subset = subset, priorSubset = priorSubset, bootStraps = bootStraps, conv = conv, nullData = nullData, returnAll = returnAll, returnPD = returnPD, verbose = verbose, discardSampling = discardSampling, cl = cl))
    if(distribution == "normal") {
      distFunctionConsensus <- function(cts, priors, libsizes, seglens) {
        matrix(dnorm(rep(cts, each = nrow(priors)),
                     sd = priors[,1],
                     mean = rep(libsizes * seglen, each = nrow(priors)) * priors[,2]
                     , log = TRUE),
               ncol = length(cts))
      }
      distFunction <- function(cts, selcts, nzWts, prior, libsizes, seglen) {
        dnorm(rep(cts[selcts], each = sum(nzWts)),
              sd = prior[nzWts,1],
              mean = rep(libsizes[selcts] * seglen[selcts], each = sum(nzWts)) * prior[nzWts,2], log = TRUE)
      }
    }
    if(distribution == "exponential") {
      distFunctionConsensus <- function(cts, priors, libsizes, seglens) {
        matrix(dexp(rep(cts, each = nrow(priors)),
                    rate = 1 / (priors[,1] * libsizes * seglen) , log = TRUE),
               ncol = length(cts))
      }
      distFunction <- function(cts, selcts, nzWts, prior, libsizes, seglen) {
        dexp(rep(cts[selcts], each = sum(nzWts)),
             rate = 1 / (prior[nzWts,1] * libsizes * seglen), log = TRUE)
      }
    }
    if(distribution == "gamma") {
      distFunctionConsensus <- function(cts, priors, libsizes, seglens) {
        matrix(dgamma(rep(cts, each = nrow(priors)),
                      shape = priors[,1], scale = priors[,2] * libsizes * seglen , log = TRUE),
               ncol = length(cts))
      }
      distFunction <- function(cts, selcts, nzWts, prior, libsizes, seglen) {
        dgamma(rep(cts[selcts], each = sum(nzWts)),
               shape = prior[nzWts,1], scale = prior[nzWts,2] * libsizes * seglen, log = TRUE)
      }
    }

    
    constructWeights <- function(withinCluster = FALSE, consensus = FALSE)
      {
        priorWeights <- lapply(1:length(numintSamp), function(ii)
                               lapply(1:length(numintSamp[[ii]]), function(jj)
                                      {                                      
                                        samplings <- numintSamp[[ii]][[jj]]
                                        samplings <- samplings[order(samplings[,2]),]
                                        if(consensus) intervals <- findInterval(1:nrow(LNpriors) - 0.5, samplings[,2]) + 1 else intervals <- findInterval(1:nrow(LNpriors[[ii]][[jj]]) - 0.5, samplings[,2]) + 1
                                        intlens <- diff(c(intervals, nrow(samplings) + 1))                                                                       
                                        weights <- c()
                                        weights[intlens == 1] <- samplings[intervals[intlens == 1],3]
                                        if(any(intlens > 1))
                                          weights[which(intlens > 1)] <- sapply(which(intlens > 1), function(kk)
                                                                                sum(samplings[intervals[kk] + 0:(intlens[kk] - 1),3]))
                                        weights
                                      }))
                                                                     
        if(withinCluster) {
          assign("priorWeights", priorWeights, envir = .GlobalEnv)
          return(invisible(NULL))
        } else return(priorWeights)
      }
    
    LNbootStrap <- function(us, libsizes, groups, lensameFlag, consensus = FALSE, differentWeights = differentWeights) {
      `logsum` <-
        function(x)
          max(x, max(x, na.rm = TRUE) + log(sum(exp(x - max(x, na.rm = TRUE)), na.rm = TRUE)), na.rm = TRUE)


      PDgivenr.LNConsensus <- function (number, cts, seglen, libsizes, priors, groups, priorWeights, numintSamp, differentWeights)
        {
          dnorms <- distFunctionConsensus(cts, prior, libsizes, seglens)                           

          if(differentWeights)
            {
              ld <- sapply(1:length(groups), function(grpnum) {
                group <- groups[[grpnum]]
                wts <- priorWeights[[grpnum]]
                sampInfo <- numintSamp[[grpnum]]
                sum(sapply(1:length(levels(group)), function(gg) {
                  selcts <- group == levels(group)[gg] & !is.na(group)
                  weightings <- wts[[gg]]
                  nzWts <- weightings != 0
                  wsInfo <- which(sampInfo[[gg]][,1] == number)
                  weightings[sampInfo[[gg]][wsInfo,2]] <- weightings[sampInfo[[gg]][wsInfo,2]] - sampInfo[[gg]][wsInfo,3]
                  logsum(rowSums(dnorms[nzWts,selcts,drop = FALSE], na.rm = TRUE) + log(weightings[nzWts])) - log(sum(weightings[nzWts]))
                }))
              })
            } else {
              weightings <- priorWeights[[1]][[1]]
              wsInfo <- which(numintSamp[[1]][[1]][,1] == number)
              weightings[numintSamp[[1]][[1]][wsInfo,2]] <- weightings[numintSamp[[1]][[1]][wsInfo,2]] - numintSamp[[1]][[1]][wsInfo,3]
              nzWts <- weightings != 0            
              dnorms <- dnorms[nzWts,,drop = FALSE]
              lweight <- log(weightings[nzWts])
              lsumweight <- log(sum(weightings[nzWts]))
              
              ld <- sapply(1:length(groups), function(grpnum) {
                group <- groups[[grpnum]]
                sum(sapply(1:length(levels(group)), function(gg) {
                  selcts <- which(group == levels(group)[gg] & !is.na(group))
                  logsum(rowSums(dnorms[,selcts,drop = FALSE], na.rm = TRUE) + lweight) - lsumweight
                })) 
              })
              
            }
          ld
        }
      
      PDgivenr.LN <- function (number, cts, seglen, libsizes, priors, group, wts, sampInfo, differentWeights)
        {
          if(length(seglen) == 1 & length(cts) > 1)
            seglen <- rep(seglen, length(cts))
          
          sum(
              sapply(1:length(levels(group)), function(gg) {
                selcts <- group == levels(group)[gg] & !is.na(group)
                prior <- priors[[gg]]
                weightings <- wts[[c(1, gg)[differentWeights + 1]]]
                nzWts <- weightings != 0
                wsInfo <- which(sampInfo[[c(1, gg)[differentWeights + 1]]][,1] == number)
                weightings[sampInfo[[c(1, gg)[differentWeights + 1]]][wsInfo,2]] <- weightings[sampInfo[[c(1, gg)[differentWeights + 1]]][wsInfo,2]] - sampInfo[[c(1, gg)[differentWeights + 1]]][wsInfo,3]
                
                logsum(
                  rowSums(
                    matrix(distFunction(cts, selcts, nzWts, prior, libsizes, seglen), ncol = sum(selcts)), na.rm = TRUE) + log(weightings[nzWts])) - log(sum(weightings[nzWts]))                       
              })
              )
        }
      
      number <- us[1]
      us <- us[-1]
      if(lensameFlag)
        {
          cts <- us[-1]
          seglen <- us[1]
        } else {
          cts <- us[-(1:(length(us) / 2))]
          seglen <- us[1:(length(us) / 2)]
        }

      

      if(consensus) {
        PDgivenr.LNConsensus(number, cts, seglen = seglen, libsizes = libsizes, priors = LNpriors, groups = groups, priorWeights = priorWeights, numintSamp = numintSamp, differentWeights = differentWeights)
      } else {      
        sapply(1:length(LNpriors), function(gg)
               PDgivenr.LN(number = number, cts = cts, seglen = seglen, libsizes = libsizes, priors = LNpriors[[gg]], group = groups[[gg]],
                           wts = priorWeights[[c(1, gg)[differentWeights + 1]]],
                           sampInfo = numintSamp[[c(1, gg)[differentWeights + 1]]],
                           differentWeights = differentWeights)
               )
      }
    }

    if(!is.null(cl))
      {
        clustAssign <- function(object, name)
          {
            assign(name, object, envir = .GlobalEnv)
            NULL
          }
        
        getLikelihoodsEnv <- new.env(parent = .GlobalEnv)
        environment(clustAssign) <- getLikelihoodsEnv
        environment(LNbootStrap) <- getLikelihoodsEnv
        clusterCall(cl, clustAssign, distFunctionConsensus, "distFunctionConsensus")
        clusterCall(cl, clustAssign, distFunction, "distFunction")
      }

    
    listPosts <- list()
    
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")
    
#    if(cD@priorType != "LN") stop("Incorrect prior type for this method of likelihood estimation")
    
    if(pET %in% c("none", "iteratively"))
      {
        if(length(prs) != length(cD@groups)) stop("'prs' must be of same length as the number of groups in the 'cD' object")
        if(any(prs < 0))
          stop("Negative values in the 'prs' vector are not permitted")
        if(!nullData & sum(prs) != 1)
          stop("If 'nullData = FALSE' then the 'prs' vector should sum to 1.")
        
        if(nullData & sum(prs) >= 1)
          stop("If 'nullData = TRUE' then the 'prs' vector should sum to less than 1.")
        
      } else if(pET %in% c("BIC"))
        prs <- rep(NA, length(cD@groups))
    
    if(!(class(subset) == "integer" | class(subset) == "numeric" | is.null(subset)))
      stop("'subset' must be integer, numeric, or NULL")
    
    if(is.null(subset)) subset <- 1:nrow(cD)

    subset <- subset[rowSums(do.call("cbind", lapply(cD@groups, function(x) do.call("cbind", lapply(levels(x), function(rep) rowSums(is.na(cD@data[subset,which(x == rep),drop = FALSE])) == length(which(x == rep))))))) == 0]

    if(is.null(priorSubset)) priorSubset <- subset
    
    if(is.null(conv)) conv <- 0
    if(nrow(cD@seglens) > 0) seglens <- cD@seglens else seglens <- matrix(1, ncol = 1, nrow = nrow(cD@data))
    if(ncol(seglens) == 1) lensameFlag <- TRUE else lensameFlag <- FALSE    

    libsizes <- as.double(cD@libsizes)
    groups <- cD@groups
    LNpriors <- cD@priors$priors
    numintSamp <- cD@priors$sampled
    weights <- cD@priors$weights    

    if(is.matrix(LNpriors)) consensus <- TRUE else consensus <- FALSE
    
    if(is.numeric(weights) & is.matrix(numintSamp) & bootStraps == 1 & !nullData)
      {
        differentWeights <- FALSE
        numintSamp <- list(list(cbind(numintSamp, weights)))
        priorWeights <- constructWeights(consensus = consensus)
      } else {
        differentWeights <- TRUE
            
        if(discardSampling) numintSamp[,1] <- NA    
        if(is.matrix(numintSamp))
          numintSamp <- lapply(groups, function(x) lapply(1:length(levels(x)), function(z) numintSamp))
        if(is.null(numintSamp))
          numintSamp <- lapply(groups, function(x) lapply(1:length(levels(x)), function(z) cbind(sampled = rep(-1, nrow(z)), representative = 1:nrow(z))))
        
        if(is.null(weights))
          weights <- lapply(numintSamp, function(x) lapply(x, function(z) weights = rep(1, nrow(z))))
        if(is.numeric(weights))
          weights <- lapply(numintSamp, function(x) lapply(x, function(z) weights = weights))
        
        numintSamp <- lapply(1:length(numintSamp), function(ii) lapply(1:length(numintSamp[[ii]]), function(jj) cbind(numintSamp[[ii]][[jj]], weights = weights[[ii]][[jj]])))

        priorWeights <- constructWeights(consensus = consensus)
      }
    
    if(nullData)
      {
        ndelocGroup <- which(unlist(lapply(cD@groups, function(x) all(x[!is.na(x)] == x[!is.na(x)][1])))) 
        if(length(ndelocGroup) == 0) stop("If 'nullData = TRUE' then there must exist some vector in groups whose members are all identical") else ndelocGroup <- ndelocGroup[1]
        
        NZLs <- NZLpriors <- NULL
        
        ndePriors <- log(LNpriors[[ndelocGroup]][[1]][,1])
        weights <- priorWeights[[ndelocGroup]][[1]]
        sep <- bimodalSep(ndePriors[ndePriors > -Inf], weights[ndePriors > -Inf], bQ = c(0, 0.5))
        
        groups <- c(groups, groups[ndelocGroup])
        ndenulGroup <- length(groups)
        prs <- c(prs, 1 - sum(prs))

        LNpriors[[ndenulGroup]] <- LNpriors[[ndelocGroup]]
        priorWeights[[ndenulGroup]] <- priorWeights[[ndelocGroup]]
        priorWeights[[ndenulGroup]][[1]] <- priorWeights[[ndenulGroup]][[1]] * as.numeric(ndePriors <= sep)
        priorWeights[[ndelocGroup]][[1]] <- priorWeights[[ndelocGroup]][[1]] * as.numeric(ndePriors > sep)
        numintSamp[[ndenulGroup]] <- numintSamp[[ndelocGroup]]
        numintSamp[[ndelocGroup]][[1]][numintSamp[[ndelocGroup]][[1]][,2] %in% which(ndePriors <= sep),3] <- 0
        numintSamp[[ndenulGroup]][[1]][numintSamp[[ndenulGroup]][[1]][,2] %in% which(ndePriors > sep),3] <- 0
      }

    posteriors <- matrix(NA, ncol = length(groups), nrow = nrow(cD@data))
    propest <- NULL
    converged <- FALSE

    if(!is.null(cl)) {
      clusterCall(cl, clustAssign, LNpriors, "LNpriors")
    }

    if(verbose) message("Finding posterior likelihoods...", appendLF = FALSE)

    priorReps <- unique(unlist(sapply(numintSamp, function(x) as.integer(unique(sapply(x, function(z) z[,1]))))))
    priorReps <- priorReps[priorReps > 0 & !is.na(priorReps)]
    
    if(!all(priorReps %in% 1:nrow(cD@data)) & bootStraps > 1)
      {
        warning("Since the sampled values in the '@priors' slot are not available, bootstrapping is not possible.")
        bootStraps <- 1
      }

    postRows <- unique(c(priorReps, priorSubset, subset))   
    
    .fastUniques <- function(x){
      if (nrow(x) > 1) {
        return(c(TRUE, rowSums(x[-1L, , drop = FALSE] == x[-nrow(x),, drop = FALSE]) != ncol(x)))
      } else return(TRUE)
    }    

    if(length(priorReps) == 0 || !any(priorReps %in% postRows))
      {
        orddat <- do.call("order", c(lapply(1:ncol(seglens), function(ii) seglens[postRows,ii]), lapply(1:ncol(cD@data), function(ii) cD@data[postRows,ii])))
        whunq <- .fastUniques(cbind(seglens[postRows,], cD@data[postRows,])[orddat,,drop = FALSE])
        postRows <- postRows
      } else whunq <- TRUE
    
    for(cc in 1:bootStraps)
      {
        if(cc > 1) numintSamp <- lapply(1:length(numintSamp), function(ii) lapply(1:length(numintSamp[[ii]]), function(jj) cbind(numintSamp[[ii]][[jj]][,1:2], weights = exp(posteriors[numintSamp[[ii]][[jj]][,1],ii]))))          
        
        if (is.null(cl)) {
          if(cc > 1)
            priorWeights <- constructWeights()

          ps <- apply(cbind(1:nrow(cD@data), seglens, cD@data)[postRows[whunq],,drop = FALSE],
                      1, LNbootStrap, libsizes = libsizes, groups = groups, lensameFlag = lensameFlag, consensus = consensus, differentWeights = differentWeights)
        } else {
          environment(constructWeights) <- getLikelihoodsEnv
          clusterCall(cl, clustAssign, numintSamp, "numintSamp")
          #clusterCall(cl, clustAssign, priorWeights, "priorWeights")
          clusterCall(cl, constructWeights, withinCluster = TRUE, consensus = consensus)

          ps <- parRapply(cl[1:min(length(cl), length(postRows[whunq]))], cbind(1:nrow(cD@data), seglens, cD@data)[postRows[whunq],, drop = FALSE],
                          LNbootStrap, libsizes = libsizes, groups = groups, lensameFlag = lensameFlag, consensus = consensus, differentWeights = differentWeights)
        }

        ps <- matrix(ps, ncol = length(groups), byrow = TRUE)
        rps <- matrix(NA, ncol = length(groups), nrow = nrow(cD@data))

        rps[postRows[whunq],] <- ps

        if(length(priorReps) == 0 || !any(priorReps %in% postRows))
          {
            rps[postRows[orddat],] <- rps[postRows[orddat[rep(which(whunq), diff(c(which(whunq), length(whunq) + 1)))]],]
          }
        

        if(returnPD) {
              if(verbose) message("done.")
              return(rps)
            }
        
        if(pET != "none")
          {
            restprs <- getPosteriors(rps[priorSubset,, drop = FALSE], prs, pET = pET, marginalise = FALSE, groups = groups, priorSubset = NULL, cl = cl)$priors
          } else restprs <- prs
        
        pps <- getPosteriors(rps[union(priorReps,subset),,drop = FALSE], prs = restprs, pET = "none", marginalise = marginalise, groups = groups, priorSubset = NULL, cl = cl)

        if(any(!is.na(posteriors)))
          if(all(abs(exp(posteriors[union(priorReps,subset),,drop = FALSE]) - exp(pps$posteriors)) < conv)) converged <- TRUE
        
        posteriors[union(priorReps, subset),] <- pps$posteriors
        
        prs <- pps$priors
        propest <- rbind(propest, prs)

        estProps <- pps$priors
        names(estProps) <- names(cD@groups)

        cat(".")        

        if(returnAll | converged | cc == bootStraps)
          {
            retPosts <- posteriors
            retPosts[priorReps[!(priorReps %in% subset)],] <- NA
            
            nullPosts <- matrix(ncol = 0, nrow = 0)
            if(nullData) {
              nullPosts <- retPosts[,ndenulGroup,drop = FALSE]
              retPosts <- retPosts[,-ndenulGroup, drop = FALSE]
            }
            estProps <- apply(exp(retPosts), 2, mean, na.rm = TRUE)
            
            colnames(retPosts) <- names(cD@groups)
            listPosts[[cc]] <- (new(class(cD), cD, posteriors = retPosts, estProps = estProps, nullPosts = nullPosts))
            listPosts[[cc]]@orderings <- .makeOrderings(listPosts[[cc]])
          }

        if(converged)
          break()
      }

    if(!is.null(cl))
      clusterEvalQ(cl, rm(list = ls()))
    
    if(verbose) message("done.")

    if(!returnAll) return(listPosts[[cc]]) else {
      if(length(listPosts) == 1) return(listPosts[[1]]) else return(listPosts)
    }
    
  }
