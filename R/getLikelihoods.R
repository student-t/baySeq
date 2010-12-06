'getLikelihoods' <- function(cD, prs, pET = "BIC", marginalise = FALSE, subset = NULL, priorSubset = NULL, verbose = TRUE, ..., cl)
  {
    type = cD@priorType
    switch(type,
           "Dir" = getLikelihoods.Dirichlet(cD = cD, prs = prs, pET = pET, marginalise = marginalise, subset = subset, priorSubset = priorSubset, verbose = verbose, cl = cl),
           "Poi" = getLikelihoods.Pois(cD = cD, prs = prs, pET = pET, marginalise = marginalise, subset = subset, priorSubset = priorSubset, verbose = verbose, ..., cl = cl),
           "NB" = getLikelihoods.NB(cD = cD, prs = prs, pET = pET, marginalise = marginalise, subset = subset, priorSubset = priorSubset, verbose = verbose, ..., cl = cl))
  }
           
`getLikelihoods.Dirichlet` <-
function(cD, prs, pET = "BIC", marginalise = FALSE, subset = NULL, priorSubset = NULL, verbose = TRUE, cl)
  {
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")

    if(cD@priorType != "Dir") stop("Incorrect prior type for this method of likelihood estimation")
    if(class(cD) != "countData")
      stop("variable 'countData' in 'getGroupPriors' must be of class 'countData'")

    if(pET %in% c("none", "iteratively"))
      {
        if(length(prs) != length(cD@groups)) stop("'prs' must be of same length as the number of groups in the 'cD' object")
        if(any(prs < 0))
          stop("Negative values in the 'prs' vector are not permitted")
      } else if(pET %in% c("BIC"))
        prs <- rep(NA, length(cD@groups))

    
    if(!(class(subset) == "integer" | class(subset) == "numeric" | is.null(subset)))
      stop("'subset' must be integer, numeric, or NULL")

    if(is.null(subset)) subset <- 1:nrow(cD@data)

    if(!is.null(priorSubset))
      {
        priorSub <- rep(FALSE, nrow(cD@data))
        priorSub[priorSubset] <- TRUE
        priorSubset <- priorSub[subset]
      }
    
    `PrgivenD.Dirichlet` <-
      function(us, libsizes, priors, groups, priorgroups = c(0.99, 0.01))
        {
          `PDgivenr.Dirichlet` <-
            function(us, prior, group)
              {
                lbeta.over.beta <- function(us, alphas)
                  {
                    if(all(alphas > 0))
                      sum(lgamma(alphas + us)) + lgamma(sum(alphas)) - sum(lgamma(alphas)) - lgamma(sum(alphas) + sum(us)) else -Inf
                  }
                
                sum(sapply(1:nrow(us), function(ii) lfactorial(sum(us[ii,])) - sum(lfactorial(us[ii,])))) +
                  sum(sapply(unique(group), function(gg)
                             lbeta.over.beta(apply(matrix(us[group == unique(group)[gg],], ncol = ncol(us)), 2, sum),
                                             prior[[gg]])))
              }

          
          us <- cbind(us, libsizes - us)
          ps <- sapply(1:length(groups), function(ii)
                       PDgivenr.Dirichlet(us, priors[[ii]], groups[[ii]]))
          
          ps
        }

    if(!is.null(cl))
      {
        getLikelihoodsEnv <- new.env(parent = .GlobalEnv)
        environment(PrgivenD.Dirichlet) <- getLikelihoodsEnv
      }

    if(verbose) message("Finding posterior likelihoods...", appendLF = FALSE)
    
    if(is.null(cl)) {
      ps <- apply(cD@data[subset,, drop = FALSE], 1, PrgivenD.Dirichlet, cD@libsizes, cD@priors$priors, cD@groups)
    } else {
      ps <- parRapply(cl, cD@data[subset,, drop = FALSE], PrgivenD.Dirichlet, cD@libsizes, cD@priors$priors, cD@groups)
    }
    ps <- matrix(ps, ncol = length(cD@groups), byrow = TRUE)
    groups <- cD@groups

    if(!is.null(cl))
      clusterEvalQ(cl, rm(list = ls()))
    
    pps <- getPosteriors(ps = ps, prs = prs, pET = pET, marginalise = marginalise, groups = groups, priorSubset = priorSubset, cl = cl)
    
    posteriors <- matrix(NA, ncol = length(cD@groups), nrow(cD@data))
    posteriors[subset,] <- t(pps$posteriors)

    if(verbose) message("done.")
    
    colnames(posteriors) <- names(cD@groups)
    estProps <- apply(exp(posteriors), 2, mean)
    names(estProps) <- names(cD@groups)
    
    new(class(cD), cD, posteriors = posteriors, estProps = estProps)
  }


`getLikelihoods.Pois` <-
function(cD, prs, pET = "BIC", marginalise = FALSE, subset = NULL, priorSubset = NULL, distpriors = FALSE, verbose = TRUE, cl)
  {
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")

    if(cD@priorType != "Poi") stop("Incorrect prior type for this method of likelihood estimation")

    if(pET %in% c("none", "iteratively"))
      {
        if(length(prs) != length(cD@groups)) stop("'prs' must be of same length as the number of groups in the 'cD' object")
        if(any(prs < 0))
          stop("Negative values in the 'prs' vector are not permitted")
      } else if(pET %in% c("BIC"))
        prs <- rep(NA, length(cD@groups))

    
    if(!(class(subset) == "integer" | class(subset) == "numeric" | is.null(subset)))
      stop("'subset' must be integer, numeric, or NULL")

    if(is.null(subset)) subset <- 1:nrow(cD@data)

    if(!is.null(priorSubset))
      {
        priorSub <- rep(FALSE, nrow(cD@data))
        priorSub[priorSubset] <- TRUE
        priorSubset <- priorSub[subset]
      }
    
    if(nrow(cD@seglens) > 0) seglens <- cD@seglens else seglens <- matrix(1, ncol = 1, nrow = nrow(cD@data))
    
    if(ncol(seglens) == 1) lensameFlag <- TRUE else lensameFlag <- FALSE
    
    `PrgivenD.Pois` <-
      function(us, libsizes, priors, groups, distpriors = FALSE, lensameFlag)
        {
          `PDgivenr.Pois` <-
            function(cts, ns, seglen, prior, group)
              {
                if(length(seglen) == 1 & length(cts) > 1)
                  seglen <- rep(seglen, length(cts))
                
                  sum(sapply(unique(group), function(gg)
                         {
                           gprior <- matrix(unlist(prior[[gg]]), nrow = length(prior[[gg]]), byrow = TRUE)
                           selcts <- group == gg
                           lgamma(sum(cts[selcts]) + prior[[gg]][1]) - lgamma(prior[[gg]][1]) +
                             prior[[gg]][1] * log(prior[[gg]][2]) - (prior[[gg]][1] + sum(cts[selcts])) * log(prior[[gg]][2] + sum(ns[selcts] * seglen[selcts]))
                         }))
              }
          
          `PDgivenr.PoisIndie` <-
            function(cts, ns, seglen, prior, group)
              {
                if(length(seglen) == 1 & length(cts) > 1)
                  seglen <- rep(seglen, length(cts))
                `logsum` <-
                  function(x)
                    max(x, max(x, na.rm = TRUE) + log(sum(exp(x - max(x, na.rm = TRUE)), na.rm = TRUE)), na.rm = TRUE)
                
                sum(cts * log(ns * seglen)) - sum(lfactorial(cts)) +
                  sum(sapply(unique(group), function(gg)
                             {
                               selcts <- group == gg
                               logsum(lgamma(sum(cts[selcts]) + prior[[gg]][,1]) - lgamma(prior[[gg]][,1]) -
                                      sum(cts[selcts]) * log(sum(ns[selcts] * seglen[selcts]) + prior[[gg]][,2]) -
                                      prior[[gg]][,1] * log(1 + sum(ns[selcts] * seglen[selcts]) / prior[[gg]][,2])) - log(nrow(prior[[gg]]))
                             }))
              }


          if(lensameFlag)
            {
              cts <- us[-1]
              lens <- us[1]
            } else {
              cts <- us[-(1:(length(us) / 2))]
              lens <- us[1:(length(us) / 2)]
            }
          
          if(!distpriors) {
            ps <- sapply(1:length(groups), function(ii)
                         PDgivenr.Pois(cts = cts, ns = libsizes, seglen = lens, prior = priors[[ii]], group = groups[[ii]]))
          } else {
            ps <- sapply(1:length(groups), function(ii)
                         PDgivenr.PoisIndie(us = us, ns = libsizes, seglen = lens, prior = priors[[ii]], group = groups[[ii]]))
          }
        }

    if(!is.null(cl))
      {
        getLikelihoodsEnv <- new.env(parent = .GlobalEnv)
        environment(PrgivenD.Pois) <- getLikelihoodsEnv
      }

    if(verbose) message("Finding posterior likelihoods...", appendLF = FALSE)
    
    if(is.null(cl)) {
      ps <- apply(cbind(seglens, cD@data)[subset,, drop = FALSE], 1, PrgivenD.Pois, cD@libsizes, cD@priors$priors, cD@groups, distpriors, lensameFlag = lensameFlag)
    } else {
      ps <- parRapply(cl, cbind(seglens, cD@data)[subset,, drop = FALSE], PrgivenD.Pois, cD@libsizes, cD@priors$priors, cD@groups, distpriors, lensameFlag = lensameFlag)
      }                        
    ps <- matrix(ps, ncol = length(cD@groups), byrow = TRUE)
    groups <- cD@groups

    if(!is.null(cl))
      clusterEvalQ(cl, rm(list = ls()))
    
    pps <- getPosteriors(ps = ps, prs = prs, pET = pET, marginalise = marginalise, groups = groups, priorSubset = priorSubset, cl = cl)
    
    posteriors <- matrix(NA, ncol = length(cD@groups), nrow(cD@data))
    posteriors[subset, ] <- pps$posteriors

    if(verbose) message("done.")

    colnames(posteriors) <- names(cD@groups)
    estProps <- apply(exp(posteriors), 2, mean)
    names(estProps) <- names(cD@groups)
    
    new(class(cD), cD, posteriors = posteriors, estProps = estProps)
  }


getLikelihoods.NBboot <- function(...)
  {
    warning("This function has been deprecated and will soon be removed. Use 'getLikelihoods.NB' instead. Passing arguments to this function now.")
    getLikelihoods.NB(...)
  }


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
                    bicps <- t(-2 * t(ps[priorSubset, , drop = FALSE]) + (1 + (unlist(lapply(lapply(groups, unique), length)))) * log(sampleSize))
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
            `logsum` <- function(x)
              max(x, max(x, na.rm = TRUE) + log(sum(exp(x - max(x, na.rm = TRUE)), na.rm = TRUE)), na.rm = TRUE)
            z <- t(t(oldposts[setdiff(priorSubset, x[1]),]) + x[-1])
            maxes <- do.call("pmax", data.frame(z))
            z <- z - pmax(z, maxes + log(rowSums(exp(z - maxes))))
            apply(z, 2, logsum) - log(colSums(!is.na(z)))
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
function(cD, prs, pET = "BIC", marginalise = FALSE, subset = NULL, priorSubset = NULL, bootStraps = 1, conv = 1e-4, nullData = FALSE, returnAll = FALSE, returnPD = FALSE, verbose = TRUE, cl)
  {
    constructWeights <- function(withinCluster = FALSE)
      {
        priorWeights <- lapply(1:length(NBpriors), function(ii) lapply(1:length(NBpriors[[ii]]), function(jj)
                                                                       sapply(1:nrow(NBpriors[[ii]][[jj]]), function(z) sum(numintSamp[[ii]][[jj]][numintSamp[[ii]][[jj]][,2] == z,3]))))
        if(withinCluster) {
          assign("priorWeights", priorWeights, envir = .GlobalEnv)
          return(invisible(NULL))
        } else return(priorWeights)
      }
    
    NBbootStrap <- function(us, libsizes, groups, lensameFlag) {
      `logsum` <-
        function(x)
          max(x, max(x, na.rm = TRUE) + log(sum(exp(x - max(x, na.rm = TRUE)), na.rm = TRUE)), na.rm = TRUE)
      
      PDgivenr.NB <- function (number, cts, seglen, libsizes, priors, group, wts, sampInfo)
        {
          if(length(seglen) == 1 & length(cts) > 1)
            seglen <- rep(seglen, length(cts))
          
          sum(
              sapply(unique(group), function(gg) {
                selcts <- group == gg
                prior <- priors[[gg]]
                weightings <- wts[[gg]]
                wsInfo <- which(sampInfo[[gg]][,1] == number)
                weightings[sampInfo[[gg]][wsInfo,2]] <- weightings[sampInfo[[gg]][wsInfo,2]] - sampInfo[[gg]][wsInfo,3]
                
                logsum(
                       rowSums(
                               matrix(
                                      dnbinom(rep(cts[selcts], each = nrow(prior)),
                                              size = 1 / prior[,2],
                                              mu = rep(libsizes[selcts] * seglen[selcts], each = nrow(prior)) * prior[,1]
                                              , log = TRUE),
                                      ncol = sum(selcts))
                               ) + log(weightings)
                       ) - log(sum(weightings))
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
      
      sapply(1:length(NBpriors), function(gg)
             PDgivenr.NB(number = number, cts = cts, seglen = seglen, libsizes = libsizes, priors = NBpriors[[gg]], group = groups[[gg]], wts = priorWeights[[gg]], sampInfo = numintSamp[[gg]]))
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
    
    if(is.null(subset)) subset <- 1:nrow(cD@data)

    if(is.null(priorSubset)) priorSubset <- subset
    
#    if(!is.null(priorSubset))
#      {
#        priorSub <- rep(FALSE, nrow(cD@data))
#        priorSub[priorSubset] <- TRUE
#        priorSubset <- priorSub[subset]
#      }

    if(is.null(conv)) conv <- 0
    if(nrow(cD@seglens) > 0) seglens <- cD@seglens else seglens <- matrix(1, ncol = 1, nrow = nrow(cD@data))
    if(ncol(seglens) == 1) lensameFlag <- TRUE else lensameFlag <- FALSE    
    
    groups <- cD@groups
    NBpriors <- cD@priors$priors
    numintSamp <- cD@priors$sampled

    if(is.matrix(numintSamp))
      numintSamp <- lapply(NBpriors, function(x) lapply(x, function(z) numintSamp))
    if(is.null(numintSamp))
      numintsamp <- lapply(NBpriors, function(x) lapply(x, function(z) cbind(sampled = rep(-1, nrow(z)), representative = 1:nrow(z))))

    weights <- cD@priors$weights
    if(is.null(weights))
      weights <- lapply(numintSamp, function(x) lapply(x, function(z) weights = rep(1, nrow(z))))
    if(is.numeric(weights))
      weights <- lapply(numintSamp, function(x) lapply(x, function(z) weights = weights))

    numintSamp <- lapply(1:length(numintSamp), function(ii) lapply(1:length(numintSamp[[ii]]), function(jj) cbind(numintSamp[[ii]][[jj]], weights = weights[[ii]][[jj]])))
    priorWeights <- constructWeights()
    
    if(nullData)
      {
        ndelocGroup <- which(unlist(lapply(cD@groups, function(x) all(x == x[1])))) 
        if(length(ndelocGroup) == 0)
          stop("If 'nullData = TRUE' then there must exist some vector in groups whose members are all identical")
        
        NZLs <- NZLpriors <- NULL
        
        ndePriors <- log(NBpriors[[ndelocGroup]][[1]][,1])
        weights <- priorWeights[[ndelocGroup]][[1]]
        sep <- bimodalSep(ndePriors[ndePriors > -Inf], weights[ndePriors > -Inf], bQ = c(0, 1))
        
        groups <- c(groups, list(rep(1, ncol(cD@data))))
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
    
    for(cc in 1:bootStraps)
      {
        if (is.null(cl)) {
          if(cc > 1)
            {
              numintSamp <- lapply(1:length(numintSamp), function(ii) lapply(1:length(numintSamp[[ii]]), function(jj) cbind(numintSamp[[ii]][[jj]][,1:2], weights = exp(posteriors[numintSamp[[ii]][[jj]][,1],ii]))))
              priorWeights <- constructWeights()
            }
          ps <- apply(cbind(1:nrow(cD@data), seglens, cD@data)[postRows,,drop = FALSE],
                      1, NBbootStrap, libsizes = cD@libsizes, groups = groups, lensameFlag = lensameFlag)
        } else {
          environment(constructWeights) <- getLikelihoodsEnv
          clusterCall(cl, clustAssign, numintSamp, "numintSamp")
          clusterCall(cl, clustAssign, priorWeights, "priorWeights")
          clusterCall(cl, constructWeights, TRUE)
          ps <- parRapply(cl, cbind(1:nrow(cD@data), seglens, cD@data)[postRows,, drop = FALSE],
                          NBbootStrap, libsizes = cD@libsizes, groups = groups, lensameFlag = lensameFlag)
        }
        
        ps <- matrix(ps, ncol = length(groups), byrow = TRUE)
        rps <- matrix(NA, ncol = length(groups), nrow = nrow(cD@data))
        rps[postRows,] <- ps

        if(pET != "none")
          {
            restprs <- getPosteriors(rps[priorSubset,, drop = FALSE], prs, pET = pET, marginalise = FALSE, groups = groups, priorSubset = NULL, cl = cl)$priors
          } else restprs <- prs
        
        pps <- getPosteriors(rps[union(priorReps,subset),], prs = restprs, pET = "none", marginalise = marginalise, groups = groups, priorSubset = NULL, cl = cl)
        
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
            
            nullPosts <- numeric(0)
            if(nullData) {
              nullPosts <- retPosts[,ndenulGroup]
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

    if(returnPD) return(rps)
    
    if(!returnAll) return(listPosts[[cc]]) else {
      if(length(listPosts) == 1) return(listPosts[[1]]) else return(listPosts)
    }
    
  }
