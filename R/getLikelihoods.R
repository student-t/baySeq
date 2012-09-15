#'getLikelihoods' <- function(cD, prs, pET = "BIC", marginalise = FALSE, subset = NULL, priorSubset = NULL, verbose = TRUE, ..., cl)
#  {
#    type = cD@priorType
#    switch(type,
#           "Dir" = getLikelihoods.Dirichlet(cD = cD, prs = prs, pET = pET, marginalise = marginalise, subset = subset, priorSubset = priorSubset, verbose = verbose, cl = cl),
#           "Poi" = getLikelihoods.Pois(cD = cD, prs = prs, pET = pET, marginalise = marginalise, subset = subset, priorSubset = priorSubset, verbose = verbose, ..., cl = cl),
#           "NB" = getLikelihoods.NB(cD = cD, prs = prs, pET = pET, marginalise = marginalise, subset = subset, priorSubset = priorSubset, verbose = verbose, ..., cl = cl))
#  }


`getPosteriors` <-
function(ps, prs, pET = "none", marginalise = FALSE, groups, priorSubset = NULL, maxit = 100, accuracy = 1e-5, cl = cl)
  {            
    getPosts <- function(ps, prs)
      {
        posts <- t(t(ps) + log(prs))
        maxes <- do.call("pmax", data.frame(posts))
        posts <- posts - pmax(maxes, maxes + log(rowSums(exp(posts - maxes))))
      }
    
    if(is.null(priorSubset))
      priorSubset <- 1:nrow(ps)

    prs <- switch(pET,                  
                  iteratively = {
                    oldprs <- prs
                    for(ii in 1:maxit)
                      {
                        posteriors <- getPosts(ps[priorSubset,], prs)

                        prs <- colSums(exp(posteriors)) / nrow(posteriors)
                        if(all(abs(oldprs - prs) < accuracy)) break
                        oldprs <- prs
                      }
                    if(ii == maxit)
                      warning("Convergence not achieved to required accuracy.")
                    prs
                  },
                  BIC = {
                    sampleSize <- length(groups[[1]])
                    bicps <- t(-2 * t(ps[priorSubset, , drop = FALSE]) + (1 + (unlist(lapply(lapply(groups, levels), length)))) * log(sampleSize))
                    minbicps <- apply(bicps, 1, which.min)
                    prs <- sapply(1:length(groups), function(x) sum(minbicps == x, na.rm = TRUE))
                    if(any(prs == 0)) prs[prs == 0] <- 1
                    prs <- prs / sum(prs)
                    prs
                  },
                  none = prs
                  )

    posteriors <- getPosts(ps, prs)
    
    if(marginalise)
      {
        oldposts <- posteriors
        intmargins <- function(x)
          {
            z <- t(t(oldposts[setdiff(priorSubset, x[1]),]) + x[-1])
            maxes <- do.call("pmax", data.frame(z))
            z <- z - pmax(z, maxes + log(rowSums(exp(z - maxes))))
            .logRowSum(t(z)) - log(colSums(!is.na(z)))
          }
        
        if(!is.null(cl)) {
          clusterCall(cl, clustAssign, priorSubset, "priorSubset")
          environment(intmargins) <- getPostsEnv
        }
        for(ii in 1:100)
          {
            if(!is.null(cl)) {
              clusterCall(cl, clustAssign, oldposts, "oldposts")
              newposts <- matrix(parRapply(cl, cbind(1:nrow(ps), ps), intmargins), ncol = ncol(ps), byrow = TRUE)
            } else newposts <- t(apply(cbind(1:nrow(ps), ps), 1, intmargins))
            
            if(max(abs(exp(newposts) - exp(oldposts))) < 1e-3) break
            oldposts <- newposts
          }
        posteriors <- newposts
      }
            
    list(posteriors = posteriors, priors = prs)
  }

`getLikelihoods.NB` <-
function(cD, prs, pET = "BIC", marginalise = FALSE, subset = NULL, priorSubset = NULL, bootStraps = 1, conv = 1e-4, nullData = FALSE, returnAll = FALSE, returnPD = FALSE, verbose = TRUE, discardSampling = FALSE, cl)
  {
    constructWeights <- function(withinCluster = FALSE, consensus = FALSE)
      {
        priorWeights <- lapply(1:length(numintSamp), function(ii)
                               lapply(1:length(numintSamp[[ii]]), function(jj)
                                      {                                      
                                        samplings <- numintSamp[[ii]][[jj]]
                                        samplings <- samplings[order(samplings[,2]),]
                                        if(consensus) intervals <- findInterval(1:nrow(NBpriors) - 0.5, samplings[,2]) + 1 else intervals <- findInterval(1:nrow(NBpriors[[ii]][[jj]]) - 0.5, samplings[,2]) + 1
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
    
    NBbootStrap <- function(us, libsizes, groups, lensameFlag, consensus = FALSE, differentWeights = differentWeights) {
      `logsum` <-
        function(x)
          max(x, max(x, na.rm = TRUE) + log(sum(exp(x - max(x, na.rm = TRUE)), na.rm = TRUE)), na.rm = TRUE)


      PDgivenr.NBConsensus <- function (number, cts, seglen, libsizes, priors, groups, priorWeights, numintSamp, differentWeights)
        {
          dnbinoms <- matrix(
                             dnbinom(rep(cts, each = nrow(priors)),
                                     size = 1 / priors[,2],
                                     mu = rep(libsizes * seglen, each = nrow(priors)) * priors[,1]
                                     , log = TRUE),
                             ncol = length(cts))

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
                  logsum(rowSums(dnbinoms[nzWts,selcts,drop = FALSE], na.rm = TRUE) + log(weightings[nzWts])) - log(sum(weightings[nzWts]))
                }))
              })
            } else {
              weightings <- priorWeights[[1]][[1]]
              wsInfo <- which(numintSamp[[1]][[1]][,1] == number)
              weightings[numintSamp[[1]][[1]][wsInfo,2]] <- weightings[numintSamp[[1]][[1]][wsInfo,2]] - numintSamp[[1]][[1]][wsInfo,3]
              nzWts <- weightings != 0            
              dnbinoms <- dnbinoms[nzWts,,drop = FALSE]
              lweight <- log(weightings[nzWts])
              lsumweight <- log(sum(weightings[nzWts]))
              
              ld <- sapply(1:length(groups), function(grpnum) {
                group <- groups[[grpnum]]
                sum(sapply(1:length(levels(group)), function(gg) {
                  selcts <- which(group == levels(group)[gg] & !is.na(group))
                  logsum(rowSums(dnbinoms[,selcts,drop = FALSE], na.rm = TRUE) + lweight) - lsumweight
                })) 
              })
              
            }
          ld
        }
      
      PDgivenr.NB <- function (number, cts, seglen, libsizes, priors, group, wts, sampInfo, differentWeights)
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
                               matrix(
                                      dnbinom(rep(cts[selcts], each = sum(nzWts)),
                                              size = 1 / prior[nzWts,2],
                                              mu = rep(libsizes[selcts] * seglen[selcts], each = sum(nzWts)) * prior[nzWts,1]
                                              , log = TRUE),
                                      ncol = sum(selcts)), na.rm = TRUE) + log(weightings[nzWts])) - log(sum(weightings[nzWts]))                
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
        PDgivenr.NBConsensus(number, cts, seglen = seglen, libsizes = libsizes, priors = NBpriors, groups = groups, priorWeights = priorWeights, numintSamp = numintSamp, differentWeights = differentWeights)
      } else {      
        sapply(1:length(NBpriors), function(gg)
               PDgivenr.NB(number = number, cts = cts, seglen = seglen, libsizes = libsizes, priors = NBpriors[[gg]], group = groups[[gg]], wts = priorWeights[[c(1, gg)[differentWeights + 1]]], sampInfo = numintSamp[[c(1, gg)[differentWeights + 1]]], differentWeights = differentWeights)
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
        environment(NBbootStrap) <- getLikelihoodsEnv
      }

    
    listPosts <- list()
    
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")
    
    if(cD@priorType != "NB") stop("Incorrect prior type for this method of likelihood estimation")
    
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

    subset <- subset[rowSums(do.call("cbind", lapply(cD@groups, function(x) do.call("cbind", lapply(levels(x), function(rep) rowSums(is.na(cD@data[,x == rep,drop = FALSE])) == sum(x == rep)))))) == 0]
    
    if(is.null(priorSubset)) priorSubset <- subset
    
    if(is.null(conv)) conv <- 0
    if(nrow(cD@seglens) > 0) seglens <- cD@seglens else seglens <- matrix(1, ncol = 1, nrow = nrow(cD@data))
    if(ncol(seglens) == 1) lensameFlag <- TRUE else lensameFlag <- FALSE    

    libsizes <- as.double(cD@libsizes)
    groups <- cD@groups
    NBpriors <- cD@priors$priors
    numintSamp <- cD@priors$sampled
    weights <- cD@priors$weights    

    if(is.matrix(NBpriors)) consensus <- TRUE else consensus <- FALSE
    
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
        
        ndePriors <- log(NBpriors[[ndelocGroup]][[1]][,1])
        weights <- priorWeights[[ndelocGroup]][[1]]
        sep <- bimodalSep(ndePriors[ndePriors > -Inf], weights[ndePriors > -Inf], bQ = c(0, 0.5))
        
        groups <- c(groups, groups[ndelocGroup])
        ndenulGroup <- length(groups)
        prs <- c(prs, 1 - sum(prs))

        NBpriors[[ndenulGroup]] <- NBpriors[[ndelocGroup]]
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

    if(!is.null(cl))
      clusterCall(cl, clustAssign, NBpriors, "NBpriors")

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
        
        if (is.null(cl) | length(postRows[whunq]) < 2) {
          if(cc > 1)
            priorWeights <- constructWeights()

          ps <- apply(cbind(1:nrow(cD@data), seglens, cD@data)[postRows[whunq],,drop = FALSE],
                      1, NBbootStrap, libsizes = libsizes, groups = groups, lensameFlag = lensameFlag, consensus = consensus, differentWeights = differentWeights)
        } else {
          environment(constructWeights) <- getLikelihoodsEnv
          clusterCall(cl, clustAssign, numintSamp, "numintSamp")
          #clusterCall(cl, clustAssign, priorWeights, "priorWeights")
          clusterCall(cl, constructWeights, withinCluster = TRUE, consensus = consensus)
          ps <- parRapply(cl, cbind(1:nrow(cD@data), seglens, cD@data)[postRows[whunq],, drop = FALSE],
                          NBbootStrap, libsizes = libsizes, groups = groups, lensameFlag = lensameFlag, consensus = consensus, differentWeights = differentWeights)
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
          if(all(abs(exp(posteriors[union(priorReps,subset),]) - exp(pps$posteriors)) < conv)) converged <- TRUE
        
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
