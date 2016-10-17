.posteriorWeights <- function(cD, group, CDpriors, consensus) {
  numintSamp <- cD@priors$sampled
  weights <- cD@priors$weights
  numintSamp <- list(list(cbind(numintSamp, weights)))
  priorWeights <- .constructWeights(numintSamp = numintSamp, CDpriors = CDpriors, consensus = consensus) 
}

makeOrderings <- function(cD, orderingFunction)
  {
      if(missing(orderingFunction))
          orderingFunction <- lapply(cD@densityFunction, function(x) x@orderingFunction)
      if(!is.list(orderingFunction)) orderingFunction <- list(orderingFunction)
      if(length(orderingFunction) == 1) orderingFunction <- lapply(cD@groups, function(x) orderingFunction[[1]])

      if(any(sapply(orderingFunction, function(x) is.null(body(x))))) {
          warning("No valid ordering function available.")    
          cD@orderings <- do.call("data.frame", lapply(cD@groups, function(group) rep(NA, nrow(cD))))
      } else {
          observables <- .catObservables(cD)
      
          orderings <- do.call("data.frame", lapply(1:length(cD@groups), function(grp) {
              group <- cD@groups[[grp]]
              if(length(levels(group)) > 1) {
                  glev <- levels(group)
                  gorders <- sapply(glev, function(gg) {
                      orderon <- orderingFunction[[grp]](dat = .sliceArray(list(NULL, which(group == gg)), cD@data, drop = FALSE),
                                        observables = lapply(observables, function(obs) .sliceArray(list(NULL, which(group == gg)), obs, drop = FALSE)))
            if(is.matrix(orderon)) orderon <- rowMeans(orderon, na.rm = TRUE)
            orderon
                  })

          rownames(gorders) <- NULL

          matlev <- matrix(glev, ncol(gorders), nrow = nrow(gorders), byrow = TRUE)
          orderings <- rep("", nrow(gorders))
        
          for(ii in 1:ncol(gorders)) {
            maxes <- gorders == matrix(do.call("pmax", c(as.list(data.frame(gorders)), na.rm = TRUE)), ncol = ncol(gorders), nrow = nrow(gorders))
            maxes[is.na(maxes)] <- FALSE
            subord <- matrix(NA, ncol(gorders), nrow = nrow(gorders))
            subord[maxes] <- matlev[maxes]
            eqord <- do.call("paste", c(as.list(as.data.frame(subord)), sep = "="))
            eqord <- gsub("=*(NA)+=*", "=", eqord)
            eqord <- gsub("^=+", "", eqord)
            eqord <- gsub("=+$", "", eqord)
            eqord <- gsub("=+", "=", eqord)          
            orderings = paste(orderings, eqord, ">", sep = "")
            gorders[maxes] <- NA
          }
          orderings <- gsub(">+$", "", orderings)        
        } else orderings = rep("", nrow(cD))
        as.factor(orderings)
      }))                   
      
      colnames(orderings) <- names(cD@groups)
      cD@orderings <- orderings
    }
    cD
  }
    
`getPosteriors` <-
function(ps, prs, pET = "none", marginalise = FALSE, groups, priorSubset = NULL, maxit = 100, accuracy = 1e-5, eqOverRep = NULL, cl = cl)
  {
    if(nrow(ps) == 0) return(list(posteriors = ps, priors = prs))    
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
                    bicps <- do.call("cbind", lapply(1:ncol(ps), function(pp) -2 * ps[priorSubset,pp] + (sum(!eqOverRep[[pp]]) * nlevels(groups[[pp]]) + sum(eqOverRep[[pp]])) * log(sampleSize)))
                    minbicps <- apply(bicps, 1, which.min)
                    prs <- sapply(1:length(groups), function(x) sum(minbicps == x, na.rm = TRUE))
                    if(any(prs == 0)) prs[prs == 0] <- 1/length(priorSubset) / sum(prs == 0)
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
            clusterExport(cl, "priorSubset", envir = environment())
          
          getPostsEnv <- new.env(parent = .GlobalEnv)
          environment(intmargins) <- getPostsEnv
        }
        for(ii in 1:100)
          {
            if(!is.null(cl)) {
                clusterExport(cl, "oldposts", envir = environment())
                newposts <- matrix(parRapply(cl, cbind(1:nrow(ps), ps), intmargins), ncol = ncol(ps), byrow = TRUE)
            } else newposts <- t(apply(cbind(1:nrow(ps), ps), 1, intmargins))
            
            if(max(abs(exp(newposts) - exp(oldposts))) < 1e-3) break
            oldposts <- newposts
          }
        posteriors <- newposts
      }
            
    list(posteriors = posteriors, priors = prs)
  }

.constructWeights <- function(numintSamp, weights = weights, CDpriors, withinCluster = FALSE, consensus = FALSE)
  {
    if(is.numeric(weights)) weights <- list(list(weights))
    priorWeights <- lapply(1:length(numintSamp), function(ii)
                           lapply(1:length(numintSamp[[ii]]), function(jj)
                               {
                                    sampw <- cbind(numintSamp[[ii]][[jj]], exweight = weights[[ii]][[jj]])
                                    nS = numintSamp[[ii]][[jj]]; exweight = weights[[ii]][[jj]]
                                    save(sampw, nS, exweight, file = "~/sampw.RData")
                                    sampw <- sampw[order(sampw[,2]),]                                    
                                    sapply(split(sampw[,4] * sampw[,3], sampw[,2]), sum)
                                  }))
    
    if(withinCluster) {
      assign("priorWeights", priorWeights, envir = .GlobalEnv)
      return(invisible(NULL))
    } else return(priorWeights)
  }



`getLikelihoods` <-
function(cD, prs, pET = "BIC", marginalise = FALSE, subset = NULL, priorSubset = NULL, bootStraps = 1, bsNullOnly = TRUE, conv = 1e-4, nullData = FALSE, weightByLocLikelihoods = TRUE, modelPriorSets = list(), modelPriorValues = list(), returnAll = FALSE, returnPD = FALSE, verbose = TRUE, discardSampling = FALSE, modelLikes = TRUE, cl = NULL, tempFile = NULL, largeness = 1e8)#, restrictResults = FALSE)
  {
    .likeDataObs <- function(xdata, densityFunction, groups, consensus = FALSE, differentWeights = differentWeights, modelLikes = TRUE) { #, restrictResults) {
      `logsum` <-
        function(x) {
          max(x, max(x, na.rm = TRUE) + log(sum(exp(x - max(x, na.rm = TRUE)), na.rm = TRUE)), na.rm = TRUE)
        }
      
      PDgivenr.Consensus <- function (number, cts, xrobs, xcobs, sobs, priors, groups, priorWeights, numintSamp, differentWeights, densityFunction)
        {
          prior <- lapply(1:ncol(priors), function(jj) priors[,jj])
          repcts <- apply(cts, (1:length(dim(cts)))[-2], function(x) rep(x, nrow(priors)))

          xobs <- c(xrobs,
                    lapply(xcobs, function(obs)
                           array(obs, dim = c(ncol(cts) * nrow(priors), dim(obs)[-c(1:2)], 1))),
                    lapply(sobs, function(obs) {
                      slobs <- obs
                      if(is.vector(slobs) || length(dim(slobs)) == 1) {
                        return(rep(slobs, nrow(priors)))
                      } else apply(slobs, 2:length(dim(slobs)), function(x) rep(x, nrow(priors)))
                    })
                    )
          xobs <- c(xobs, list(dim = datdim))
          datalikes <- matrix(                         
                         densityFunction(repcts,
                                         observables = xobs,
                                         parameters = lapply(prior, function(priorpar) rep(priorpar, each = ncol(cts)))),
                         ncol = ncol(cts), byrow = TRUE)

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
#                  weightings[sampInfo[[gg]][wsInfo,2]] <- weightings[sampInfo[[gg]][wsInfo,2]] - sampInfo[[gg]][wsInfo,3]
                  logsum(rowSums(datalikes[nzWts,selcts,drop = FALSE], na.rm = TRUE) + log(weightings[nzWts])) - log(sum(weightings[nzWts]))
                }))
              })
            } else {
              weightings <- priorWeights[[1]][[1]]
              wsInfo <- which(numintSamp[[1]][[1]][,1] == number)
              weightings[numintSamp[[1]][[1]][wsInfo,2]] <- weightings[numintSamp[[1]][[1]][wsInfo,2]] - numintSamp[[1]][[1]][wsInfo,3]
              nzWts <- weightings != 0            
              datalikes <- datalikes[nzWts,,drop = FALSE]
              lweight <- log(weightings[nzWts])
              lsumweight <- log(sum(weightings[nzWts]))
              
              matld <- sapply(whmat, function(grp) logsum(rowSums(datalikes[,grp,drop = FALSE]) + lweight) - lsumweight)
              ld <- sapply(split(matld[matmat], zid), sum)
              
#              system.time(ld <- sapply(1:1000, function(grpnum) {                
#                group <- groups[[grpnum]]
#                sum(sapply(1:length(levels(group)), function(gg) {
#                  selcts <- which(group == levels(group)[gg] & !is.na(group))
#                  logsum(rowSums(datalikes[,selcts,drop = FALSE], na.rm = TRUE) + lweight) - lsumweight
#                }))
#              })
              
            }
          ld
        }
      
      PDgivenr <- function (number, cts, xrobs, xcobs, sobs, priors, group, wts, sampInfo, differentWeights, modelLikes = TRUE,densityFunction)
        {          
          glikes <- sapply(1:length(levels(group)), function(gg) {
            selcts <- group == levels(group)[gg] & !is.na(group)
            weightings <- wts[[c(1, gg)[differentWeights + 1]]]
            nzWts <- weightings != 0
            
            prior <- lapply(1:ncol(priors[[gg]]), function(jj) priors[[gg]][nzWts,jj])
            xobs <- c(xrobs,
                      lapply(xcobs, function(obs)
                             array(.sliceArray(list(NULL, which(selcts)), obs, drop = FALSE), dim = c(sum(nzWts) * sum(selcts), dim(obs)[-(1:2)], 1))),
                      lapply(sobs, function(obs) {
                        slobs <- .sliceArray(list(which(selcts)), obs, drop = FALSE)
                        if(is.vector(slobs) || length(dim(slobs)) == 1) {
                          return(rep(slobs, sum(nzWts)))
                        } else apply(slobs, 2:length(dim(slobs)), function(x) rep(x, sum(nzWts)))                            
                      })
                      )
            xobs <- c(xobs, list(dim = datdim))
            
            repcts <- apply(.sliceArray(list(NULL, which(selcts)), cts, drop = FALSE), (1:length(dim(cts)))[-2], function(x) rep(x, sum(nzWts)))
            
            wsInfo <- which(sampInfo[[c(1, gg)[differentWeights + 1]]][,1] == number)
#            weightings[sampInfo[[c(1, gg)[differentWeights + 1]]][wsInfo,2]] <- weightings[sampInfo[[c(1, gg)[differentWeights + 1]]][wsInfo,2]] - sampInfo[[c(1, gg)[differentWeights + 1]]][wsInfo,3]
            
            likeD <-
              rowSums(
                matrix(
                  densityFunction(repcts,
                                  observables = xobs,
                                  parameters = lapply(prior, function(priorpar) rep(priorpar, each = sum(selcts))))
                  ,ncol = sum(selcts),byrow = TRUE)
                , na.rm = TRUE) + log(weightings[nzWts]) - log(sum(weightings[nzWts]))            
            
            if(modelLikes) {
              mL <- logsum(likeD)
              return(mL)
            } else {
              lD <- rep(-Inf, nrow(priors[[gg]]))
              lD[nzWts] <- likeD
              return(lD)
            }
          })
          if(modelLikes) {
            return(sum(glikes))
          } else return(glikes)        
        }
      
      xid <- xdata$id
      cts <- xdata$data

      xrobs <- xdata$rowObs
      xcobs <- xdata$cellObs
      
      if(consensus) {
        PDlikes <- PDgivenr.Consensus(number = xid, cts = cts, xrobs = xrobs, xcobs = xcobs, sobs = sampleObservables, priors = CDpriors, groups = groups, priorWeights = priorWeights, numintSamp = numintSamp, differentWeights = differentWeights, densityFunction = densityFunction[[1]])
      } else {      
          PDlikes <- lapply(1:length(CDpriors), function(gg)
              PDgivenr(
                  number = xid, cts = cts, xrobs = xrobs, xcobs = xcobs, sobs = sampleObservables, priors = CDpriors[[gg]], group = groups[[gg]],
                  wts = priorWeights[[c(1, gg)[differentWeights + 1]]],
                  sampInfo = numintSamp[[c(1, gg)[differentWeights + 1]]],
                  differentWeights = differentWeights, modelLikes = modelLikes,
                  densityFunction <- densityFunction[[gg]])
                            )
        if(modelLikes) PDlikes <- unlist(PDlikes)
      }
#      if(restrictResults) {
#
#        restrictThresh <- 1 / (length(groups))
#        maxprior <- 1 / length(groups)
#        minpriors <- log(1 - maxprior) - log(length(groups))
#        valid <- (logsum(c(log(maxprior) + log(1 - restrictThresh), log(restrictThresh) + minpriors))) + PDlikes > log(restrictThresh) + logsum(PDlikes + minpriors)
#        
#        return(list(whval = which(valid), PDL = PDlikes[valid], residual = logsum(PDlikes[!valid])))
#      }
      return(PDlikes)
  }

    
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")   
    if(length(modelPriorValues) == 0 & !missing(prs)) modelPriorValues <- prs
        
    listPosts <- list()
    
    if(!(class(subset) == "integer" | class(subset) == "numeric" | is.null(subset)))
      stop("'subset' must be integer, numeric, or NULL")    
    if(is.null(subset)) subset <- 1:nrow(cD)    
    if(is.null(priorSubset)) priorSubset <- subset    
    if(length(modelPriorSets) == 0) modelPriorSets <- list(subset)
    if(any(duplicated(unlist(modelPriorSets)))) stop("Some element appears twice in the modelPriorSets list")
    if(!all(subset %in% unlist(modelPriorSets))) stop("The modelPriorSets list does not contain all the data specified in the `subset' parameter (or all data, if this parameter is not specified).")
    
    modelPriorSets <- lapply(modelPriorSets, function(x) x[x %in% subset])
    if(is.numeric(modelPriorValues))
      modelPriorValues <- lapply(modelPriorSets, function(x) modelPriorValues)

    if(is.list(modelPriorValues) && length(modelPriorValues) == 0) {
      modelPriorValues <- lapply(modelPriorSets, function(x) rep(NA, length(cD@groups)))
    }
      
    if(length(modelPriorValues) != length(modelPriorSets)) stop("The length of 'modelPriorValues' (if a list) must be identical to that of 'modelPriorSets' (or zero).")    
        
    if(pET %in% c("none", "iteratively"))
      lapply(modelPriorValues, function(prs) {
        if(length(prs) != length(cD@groups)) stop("All members of modelPriorValues must be of same length as the number of groups in the 'cD' object")
        if(any(prs < 0))
          stop("Negative values in all members of modelPriorValues are not permitted")
        if(!nullData & sum(prs) != 1)
          stop("If 'nullData = FALSE' then all members of modelPriorValues should sum to 1.")
        
        if(nullData & sum(prs) >= 1)
          stop("If 'nullData = TRUE' then all members of modelPriorValues should sum to less than 1.")
      })
    if(pET %in% c("iteratively", "BIC") && (any(sapply(modelPriorValues, length) < 100) & sapply(modelPriorSets, is.null))) warning("Some subsets contain fewer than one hundred members; estimation of model priors may be unstable")
    
    if(is.null(conv)) conv <- 0

    groups <- cD@groups
    CDpriors <- cD@priors$priors
    cdDF <- cD@densityFunction
    if(length(cdDF) == 1) cdDF <- lapply(1:length(groups), function(ii) cdDF[[1]])
    if(length(cD@densityFunction) == 1) {
        densityFunction <- lapply(1:length(groups), function(ii) cD@densityFunction[[1]]@density)
    } else densityFunction <- lapply(1:length(groups), function(ii) cD@densityFunction[[ii]]@density)
    numintSamp <- cD@priors$sampled
    weights <- cD@priors$weights; if(is.null(weights)) weights <- rep(1, nrow(cD@priors$sampled))

    nullWeights <- cD@priors$nullWeights;
    data <- cD@data
    datdim <- dim(cD)

    sampleObservables <- cD@sampleObservables
    cellObservables <- cD@cellObservables
    rowObservables <- cD@rowObservables
    
    if(!("seglens" %in% names(cellObservables) || "seglens" %in% names(rowObservables))) rowObservables <- c(rowObservables, list(seglens = rep(1, nrow(cD))))
    if(!("libsizes" %in% names(sampleObservables))) sampleObservables <- c(sampleObservables, list(seglens = rep(1, ncol(cD))))
    
    if(is.matrix(CDpriors)) consensus <- TRUE else consensus <- FALSE
    
    if(is.numeric(weights) & is.matrix(numintSamp) & bootStraps == 1 & !nullData)
      {
        differentWeights <- FALSE
        numintSamp <- list(list(numintSamp))
        priorWeights <- .constructWeights(numintSamp = numintSamp, weights = weights, CDpriors = CDpriors, consensus = consensus)
      } else {
        differentWeights <- TRUE
            
        if(discardSampling) numintSamp[,1] <- NA    
        if(is.matrix(numintSamp))
          numintSamp <- lapply(groups, function(x) lapply(1:length(levels(x)), function(z) numintSamp))        
        if(is.null(weights))
          weights <- lapply(numintSamp, function(x) lapply(x, function(z) weights = rep(1, nrow(z))))
        if(is.numeric(weights))
          weights <- lapply(numintSamp, function(x) lapply(x, function(z) weights = weights))
#        numintSamp <- lapply(1:length(numintSamp), function(ii) lapply(1:length(numintSamp[[ii]]), function(jj) cbind(numintSamp[[ii]][[jj]], weights = weights[[ii]][[jj]])))        
        priorWeights <- .constructWeights(numintSamp = numintSamp, weights = weights, CDpriors = CDpriors, consensus = consensus)
      }

    if(nullData)
      {

          ndelocGroup <- which(unlist(lapply(cD@groups, function(x) all(x[!is.na(x)] == x[!is.na(x)][1])))) 
          if(length(ndelocGroup) == 0) stop("If 'nullData = TRUE' then there must exist some vector in groups whose members are all identical") else ndelocGroup <- ndelocGroup[1]
          
          nullFunction <- cD@densityFunction[[ndelocGroup]]@nullFunction
          modifyNullPriors <- cD@densityFunction[[ndelocGroup]]@modifyNullPriors
          if(is.null(body(nullFunction)) & nullData) {
              warning("nullData cannot be TRUE if no nullFunction is specified within the supplied densityFunction object.")
              nullData <- FALSE
          } else {
            groups <- c(groups, null = groups[ndelocGroup])
            densityFunction <- c(densityFunction, densityFunction[ndelocGroup])
            cdDF <- c(cdDF, cdDF[ndelocGroup])
            ndenulGroup <- length(groups)
            numintSamp[[ndenulGroup]] <- numintSamp[[ndelocGroup]]
            CDpriors[[ndenulGroup]] <- modifyNullPriors(CDpriors[[ndelocGroup]], datdim)
            
            if(!consensus) ndePriors <- nullFunction(CDpriors[[ndelocGroup]][[1]]) else ndePriors <- nullFunction(CDpriors)        
            
            if(weightByLocLikelihoods && "locLikelihoods" %in% slotNames(cD) && nrow(cD@locLikelihoods) > 0){
                newts <- exp(rowSums(log(1 - exp(cD@locLikelihoods)), na.rm = TRUE))[cD@priors$sampled[,1]]
                
                weights <- lapply(groups, function(x) lapply(levels(x), function(jj) return(1 - newts)))
                weights[[ndenulGroup]][[1]] <- newts
                cD
            } else {                                      
                if(is.null(nullWeights)) {
                    nullweights <- priorWeights[[ndelocGroup]][[1]]
                    sep <- bimodalSeparator(ndePriors[ndePriors > -Inf], nullweights[ndePriors > -Inf])
                    modelPriorValues <- lapply(modelPriorValues, function(prs) c(prs, 1 - sum(prs)))          
                    priorWeights[[ndenulGroup]] <- priorWeights[[ndelocGroup]]
                    priorWeights[[ndenulGroup]][[1]] <- priorWeights[[ndenulGroup]][[1]] * as.numeric(ndePriors <= sep)
                    priorWeights[[ndelocGroup]][[1]] <- priorWeights[[ndelocGroup]][[1]] * as.numeric(ndePriors > sep)
                    weights[[ndenulGroup]] <- weights[[ndelocGroup]]
                    weights[[ndelocGroup]][[1]][numintSamp[[ndelocGroup]][[1]][,2] %in% which(ndePriors <= sep)] <- 0
                    weights[[ndenulGroup]][[1]][numintSamp[[ndenulGroup]][[1]][,2] %in% which(ndePriors > sep)] <- 0
                } else weights[[ndenulGroup]] <- nullWeights
            }
            priorWeights <- .constructWeights(numintSamp = numintSamp, weights = weights, CDpriors = CDpriors, consensus = consensus)
        }
    }
    
    propest <- NULL
    converged <- FALSE

    if(consensus) {
#      rg <- do.call("rbind", groups)
#      greps <- unique(lapply(1:ncol(rg), function(ii) which(colSums(rg == rg[,ii]) == nrow(rg))))
#      gr <- c(); gr[unlist(greps)] <- rep(1:length(greps), sapply(greps, length))
#      permat <- do.call("rbind", lapply(1:max(gr), function(kk)
#                                        chooseMatrix(max(gr), kk)))[,gr] == 1

      z <- lapply(groups, function(grp) split(1:ncol(cD), as.numeric(grp)))
      whmat <- unique(do.call("c",z))

                                        #      whmat <- apply(permat, 1, function(x) which(x))
      matmat <- match(do.call("c", z), whmat)
      zid <- rep(1:length(z), sapply(z, length))
    } else {
      matmat <- NULL; whmat <- NULL; zid <- NULL
    }
    
    if(!is.null(cl))
      clusterExport(cl, c("CDpriors", "datdim", "densityFunction", "sampleObservables", ".sliceArray", "matmat", "whmat", "zid"), envir = environment())
    
    if(verbose) message("Finding posterior likelihoods...", appendLF = FALSE)

    if(bootStraps > 1) {
      priorReps <- unique(unlist(sapply(numintSamp, function(x) as.integer(unique(sapply(x, function(z) z[,1]))))))
      priorReps <- priorReps[priorReps > 0 & !is.na(priorReps)]
      
      if(!all(priorReps %in% 1:nrow(cD@data)) & bootStraps > 1)
      {
        warning("Since the sampled values in the '@priors' slot are not available, bootstrapping is not possible.")
        bootStraps <- 1
      }
    } else priorReps <- c()
    
    postRows <- unique(c(priorReps, priorSubset, subset))   

    if(verbose) message("Preparing data...", appendLF = FALSE)
    sliceData <- list()
    sliceData[postRows] <- lapply(postRows, function(id) {
        if(verbose) if(sample(1:round(nrow(data) / 50), 1) == 1) message(".", appendLF = FALSE)
        list(id = id,
             data = asub(data, id, dims = 1, drop = FALSE),
             cellObs = lapply(cellObservables, function(cob) asub(cob, id, dims = 1, drop = FALSE)),
             rowObs = lapply(rowObservables, function(rob) asub(rob, id, dims = 1, drop = FALSE)))
    })
    if(verbose) message("done.")
                    


    
    .fastUniques <- function(x){
      if (nrow(x) > 1) {
        return(c(TRUE, rowSums(x[-1L, , drop = FALSE] == x[-nrow(x),, drop = FALSE]) != ncol(x)))
      } else return(TRUE)
    }    

    whunq <- TRUE

    ncuts <- ceiling(length(groups) * as.double(nrow(cD)) / largeness)
    if(ncuts == 1) splitRows <- list(postRows) else splitRows <- split(postRows, cut(postRows, ncuts, labels = FALSE))
    

    posteriors <- NULL
    
    for(cc in 1:bootStraps)
      {          
        if(cc > 1)
          {
            if(!bsNullOnly | !nullData) {
              weights <- lapply(1:length(numintSamp), function(ii) lapply(1:length(numintSamp[[ii]]), function(jj) exp(posteriors[numintSamp[[ii]][[jj]][,1],ii])))
            } else {
              weights <- lapply(1:length(numintSamp), function(ii)
                                lapply(1:length(numintSamp[[ii]]), function(jj)
                                       {
                                         if(ii == ndenulGroup) weights = exp(posteriors[numintSamp[[ii]][[jj]][,1],ii]) else weights = 1 - exp(posteriors[numintSamp[[ii]][[jj]][,1],ndenulGroup])
                                            weights
                                          }))
            }
          }

        ps <- NULL
        if (is.null(cl)) {
          if(cc > 1)
            priorWeights <- .constructWeights(numintSamp = numintSamp, weights = weights, CDpriors = CDpriors)

          for(ss in 1:length(splitRows)) {
              tps <- lapply(sliceData[splitRows[[ss]]],
                            .likeDataObs, densityFunction = densityFunction, groups = groups, consensus = consensus, differentWeights = differentWeights, modelLikes = modelLikes)
              if(!is.null(tempFile)) {
                  save(tps, file = paste(tempFile, "_", ss, ".RData", sep = ""))
                  if(verbose) message(".", appendLF = FALSE)
              }
              ps <- c(ps, tps)
          }
      } else {
          clusterExport(cl, "numintSamp", envir = environment())
          
          clusterCall(cl, .constructWeights, numintSamp = numintSamp, weights = weights, CDpriors = CDpriors, withinCluster = TRUE, consensus = consensus)
          getLikesEnv <- new.env(parent = .GlobalEnv)
          environment(.likeDataObs) <- getLikesEnv
          for(ss in 1:length(splitRows)) {
              tps <- parLapplyLB(cl[1:min(length(cl), length(postRows[whunq]))], sliceData[splitRows[[ss]]],
                                 .likeDataObs, densityFunction = densityFunction, groups = groups, consensus = consensus, differentWeights = differentWeights, modelLikes = modelLikes)#, restrictResults = restrictResults)
              if(!is.null(tempFile)) {
                  save(tps, file = paste(tempFile, "_", ss, ".RData", sep = ""))
                  if(verbose) message(".", appendLF = FALSE)
              }
              ps <- c(ps, tps)
          }
      }

        if(!is.null(tempFile)) return(paste(tempFile, "_", 1:length(splitRows), ".RData", sep = ""))
        
        if(!modelLikes) return(ps)

        rps <- matrix(NA, ncol = length(groups), nrow = nrow(cD@data))

        rps[postRows[whunq],] <- do.call("rbind", ps)

        if(returnPD) {
            message("done.")
            return(rps)
        }

        restprs <- lapply(1:length(modelPriorSets), function(pp) {
          pSub <- intersect(priorSubset, modelPriorSets[[pp]])
          prs <- modelPriorValues[[pp]]
          if(pET == "iterative" || (pET == "BIC" & all(is.na(modelPriorValues[[pp]]))))
            {
                if(length(pSub) == nrow(cD) && all(1:nrow(cD) == pSub)) pps <- rps else pps <- rps[pSub,,drop = FALSE]
                restprs <- getPosteriors(ps = pps, prs, pET = pET, marginalise = FALSE, groups = groups, priorSubset = NULL,
                                         eqOverRep = lapply(cdDF, function(x) x@equalOverReplicates(dim(cD))), cl = cl)$priors
            } else restprs <- prs
          restprs
        })
        restprs <- lapply(restprs, function(x) {
          names(x) <- names(groups)
          x
        })
        names(restprs) <- names(modelPriorSets)


        ppsPosts <- lapply(1:length(modelPriorSets), function(pp) {
          pSub <- intersect(union(priorReps, subset), modelPriorSets[[pp]])
          if(length(pSub) == nrow(cD) && all(1:nrow(cD) == pSub)) pps <- rps else pps <- rps[pSub,,drop = FALSE]
          pps <- getPosteriors(ps = pps, prs = restprs[[pp]], pET = "none", marginalise = marginalise, groups = groups, priorSubset = NULL, cl = cl)
          list(pps = pps, pSub = pSub)
        })

        if(is.null(posteriors)) posteriors <- matrix(NA, ncol = length(groups), nrow = nrow(cD@data))
        compPosts <- do.call("c", lapply(ppsPosts, function(x) x$pSub))
        
        newPosts <- posteriors

        for(ii in 1:length(ppsPosts)) {
          newPosts[ppsPosts[[ii]]$pSub,] <- ppsPosts[[ii]]$pps$posteriors
        }
        
        if(any(!is.na(posteriors)))
          if(all(abs(exp(posteriors[compPosts,,drop = FALSE]) - exp(newPosts[compPosts,,drop = FALSE])) < conv)) converged <- TRUE

        posteriors <- newPosts

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
            colnames(retPosts) <- names(cD@groups)            
            if(nullData) {
              cD@priors$weights <- weights[-ndenulGroup]
              cD@priors$nullWeights <- weights[[ndenulGroup]]
            } else cD@priors$weights <- weights
            listPosts[[cc]] <- (new(class(cD), cD, posteriors = retPosts,
                                    nullPosts = nullPosts, priorModels = restprs))
            #if(!is.null(body(cD@densityFunction@orderingFunction)))
            listPosts[[cc]] <- makeOrderings(listPosts[[cc]])
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

