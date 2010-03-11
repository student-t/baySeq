'getLikelihoods' <- function(cD, prs, pET = "BIC", subset = NULL, priorSubset = NULL, ..., cl)
  {
    type = cD@priorType
    switch(type,
           "Dir" = getLikelihoods.Dirichlet(cD = cD, prs = prs, pET = pET, subset = subset, priorSubset = priorSubset, cl = cl),
           "Poi" = getLikelihoods.Pois(cD = cD, prs = prs, pET = pET, subset = subset, priorSubset = priorSubset, ..., cl = cl),
           "NB" = getLikelihoods.NB(cD = cD, prs = prs, pET = pET, subset = subset, priorSubset = priorSubset, ..., cl = cl))
  }
           
`getLikelihoods.Dirichlet` <-
function(cD, prs, pET = "BIC", subset = NULL, priorSubset = NULL, cl)
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

    message("Finding posterior likelihoods...", appendLF = FALSE)
    
    if(is.null(cl)) {
      ps <- apply(cD@data[subset,], 1, PrgivenD.Dirichlet, cD@libsizes, cD@priors$priors, cD@groups)
    } else {
      ps <- parRapply(cl, cD@data[subset,], PrgivenD.Dirichlet, cD@libsizes, cD@priors$priors, cD@groups)
    }
    ps <- matrix(ps, ncol = length(cD@groups), byrow = TRUE)
    groups <- cD@groups

    pps <- getPosteriors(ps = ps, prs = prs, pET = pET, groups = groups, priorSubset = priorSubset, cl = cl)
    
    posteriors <- matrix(NA, ncol = length(cD@groups), nrow(cD@data))
    posteriors[subset,] <- t(pps$posteriors)

    message("done.")
    
    colnames(posteriors) <- names(cD@groups)
    estProps <- apply(exp(posteriors), 2, mean)
    names(estProps) <- names(cD@groups)
    
    new(class(cD), cD, posteriors = posteriors, estProps = estProps)
  }


`getLikelihoods.Pois` <-
function(cD, prs, pET = "BIC", subset = NULL, priorSubset = NULL, distpriors = FALSE, cl)
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
                
                sum(cts * log(ns * seglen)) - sum(lfactorial(cts)) + 
                  sum(sapply(unique(group), function(gg)
                         {
                           gprior <- matrix(unlist(prior[[gg]]), nrow = length(prior[[gg]]), byrow = TRUE)
                           selcts <- group == gg
                           lgamma(sum(cts[selcts]) + prior[[gg]][1]) - lgamma(prior[[gg]][1]) -
                             sum(cts[selcts]) * log(sum(ns[selcts] * seglen[selcts]) + prior[[gg]][2]) -
                               prior[[gg]][1] * log(1 + sum(ns[selcts] * seglen[selcts]) / prior[[gg]][2])
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

    message("Finding posterior likelihoods...", appendLF = FALSE)
    
    if(is.null(cl)) {
      ps <- apply(cbind(seglens, cD@data)[subset,], 1, PrgivenD.Pois, cD@libsizes, cD@priors$priors, cD@groups, distpriors, lensameFlag = lensameFlag)
    } else {
      ps <- parRapply(cl, cbind(seglens, cD@data)[subset,], PrgivenD.Pois, cD@libsizes, cD@priors$priors, cD@groups, distpriors, lensameFlag = lensameFlag)
      }                        
    ps <- matrix(ps, ncol = length(cD@groups), byrow = TRUE)
    groups <- cD@groups
    
    pps <- getPosteriors(ps = ps, prs = prs, pET = pET, groups = groups, priorSubset = priorSubset, cl = cl)
    
    posteriors <- matrix(NA, ncol = length(cD@groups), nrow(cD@data))
    posteriors[subset,] <- pps$posteriors

    message("done.")

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
function(ps, prs, pET = "none", groups, priorSubset = NULL, maxit = 100, accuracy = 1e-5, cl = cl)
  {
    getPosts <- function(x, prs)
      {
        `logsum` <-
          function(x)
            max(x, max(x, na.rm = TRUE) + log(sum(exp(x - max(x, na.rm = TRUE)), na.rm = TRUE)), na.rm = TRUE)

        posts <- x + log(prs)
        posts <- posts - logsum(posts)
        posts

      }
    
    if(is.null(priorSubset))
      priorSubset <- 1:nrow(ps)

    if(!is.null(cl))
      {
        getPostsEnv <- new.env(parent = .GlobalEnv)
        environment(getPosts) <- getPostsEnv
      }

    prs <- switch(pET,
                  iteratively = {
                    oldprs <- prs
                    for(ii in 1:maxit)
                      {
                        if(!is.null(cl)) {
                          posteriors <- matrix(parRapply(cl, ps[priorSubset,], getPosts, prs), ncol = ncol(ps), byrow = TRUE)
                        } else posteriors <- matrix(apply(ps[priorSubset,], 1, getPosts, prs), ncol = ncol(ps), byrow = TRUE)
                        prs <- colSums(exp(posteriors)) / nrow(posteriors)
                        if(all(abs(oldprs - prs) < accuracy)) break
                        oldprs <- prs
                      }
                    if(ii == maxit)
                      warning("Convergence not achieved to required accuracy.")
                    prs
                  },
                  BIC = {
                    sampleSize <- length(groups[[1]][[1]])
                    bicps <- t(-2 * t(ps[priorSubset,]) + (1 + (unlist(lapply(lapply(groups, unique), length)))) * log(sampleSize))
                    minbicps <- apply(bicps, 1, which.min)
                    prs <- sapply(1:length(groups), function(x) sum(minbicps == x))
                    if(any(prs == 0)) prs[prs == 0] <- 1
                    prs <- prs / sum(prs)
                    prs
                  },
                  none = prs
                  )
           
    list(posteriors = matrix(apply(ps, 1, getPosts, prs), ncol = ncol(ps), byrow = TRUE), priors = prs)
  }

`getLikelihoods.NB` <-
function(cD, prs, pET = "BIC", subset = NULL, priorSubset = NULL, bootStraps = 1, conv = 1e-4, nullData = FALSE, returnAll = FALSE, cl)
  {
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
    
    if(!is.null(priorSubset))
      {
        priorSub <- rep(FALSE, nrow(cD@data))
        priorSub[priorSubset] <- TRUE
        priorSubset <- priorSub[subset]
      }
    
    numintSamp <- cD@priors$sampled
    groups <- cD@groups

    NBpriors <- cD@priors$priors

    NBbootStrap <- function(us, libsizes, groups, lensameFlag) {
      `logsum` <-
        function(x)
          max(x, max(x, na.rm = TRUE) + log(sum(exp(x - max(x, na.rm = TRUE)), na.rm = TRUE)), na.rm = TRUE)
      
      PDgivenr.NB <- function (number, us, seglen, ns, prior, group, sampDist)
        {
          if(length(seglen) == 1 & length(us) > 1)
            seglen <- rep(seglen, length(us))
          sampDist <- sampDist[numintSamp != number]

          sum(sapply(unique(group), function(gg) {
            selus <- group == gg
            priors <- prior[[gg]]
            priors <- priors[numintSamp != number,,drop = FALSE]
            logsum(
                   rowSums(
                           matrix(
                                  dnbinom(rep(us[selus], each = nrow(priors)),
                                          size = 1 / priors[,2],
                                          mu = rep(ns[selus] * seglen[selus], each = nrow(priors)) * priors[,1]
                                          , log = TRUE),
                                  ncol = sum(selus))
                           ) +
                   log(sampDist) - log(1 / length(sampDist))
                   ) - log(length(sampDist))
          }))
        }
      
      number <- us[1]
      us <- us[-1]
      if(lensameFlag)
        {
          cts <- us[-1]
          lens <- us[1]
        } else {
          cts <- us[-(1:(length(us) / 2))]
          lens <- us[1:(length(us) / 2)]
        }
      
      sapply(1:length(NBpriors), function(group)
             PDgivenr.NB(number, cts, lens, libsizes, NBpriors[[group]], groups[[group]], sampPriors[[group]]))

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
    
    if(is.null(conv)) conv <- 0
    posteriors <- matrix(1 / length(cD@groups), ncol = length(cD@groups), nrow = nrow(cD@data))
    propest <- NULL

    if(nrow(cD@seglens) > 0) seglens <- cD@seglens else seglens <- matrix(1, ncol = 1, nrow = nrow(cD@data))
    
    if(ncol(seglens) == 1) lensameFlag <- TRUE else lensameFlag <- FALSE
    
    if(nullData)
      {
        ndelocGroup <- which(unlist(lapply(groups, function(x) all(x == x[1]))))
        if(length(ndelocGroup) == 0)
          stop("If 'nullData = TRUE' then there must exist some vector in groups whose members are all identical")
        groups <- c(groups, list(rep(1, ncol(cD@data))))
        ndenulGroup <- length(groups)
        prs <- c(prs, 1 - sum(prs))

        NBpriors[[ndenulGroup]] <- NBpriors[[ndelocGroup]]
        
        NZLs <- apply(cD@data[numintSamp,, drop = FALSE] != 0, 1, any)
        NZLpriors <- log(cD@priors$priors[[ndelocGroup]][[1]][NZLs,1])
        
        dsortp <- sort(NZLpriors, decreasing = TRUE)
        dcummeans <- cumsum(dsortp) / 1:length(dsortp)
        dcumsquare <- cumsum(dsortp ^ 2) / (1:length(dsortp) - 1)
        dcumvar <- dcumsquare - (dcummeans^2) * 1:length(dsortp) / (1:length(dsortp) - 1)
        dcumse <- sqrt(dcumvar / 1:length(dsortp))

        isortp <- sort(NZLpriors, decreasing = FALSE)
        icummeans <- cumsum(isortp) / 1:length(isortp)
        icumsquare <- cumsum(isortp ^ 2) / (1:length(isortp) - 1)
        icumvar <- icumsquare - (icummeans^2) * 1:length(isortp) / (1:length(isortp) - 1)
        icumse <- sqrt(icumvar / 1:length(isortp))

        meanvars <- (c(rev(dcumvar)[-1], NA) + icumvar) / 2
        
        varrise <- which(meanvars[-1] > meanvars[-length(meanvars)])
        sep <- min(varrise[varrise > length(NZLpriors) / 100])
        
        kstar <- c(mean(isortp[1:sep]), mean(isortp[(sep + 1):length(isortp)]))
        
        kmNZ <- kmeans(NZLpriors, centers = matrix(kstar, ncol = 1))
        km <- rep(1, length(numintSamp))
        km[NZLs] <- kmNZ$cluster

        sampPosteriors <- (sapply(1:length(groups), function(ii)
                                   c(1, 0)[as.numeric(km == c(1,2)[as.numeric(ii == ndenulGroup) + 1]) + 1]))
      } else sampPosteriors <- matrix(1, nrow = length(numintSamp), ncol = length(groups))


    posteriors <- matrix(NA, ncol = length(groups), nrow = nrow(cD@data))
    propest <- NULL
    converged <- FALSE

    if(!is.null(cl))
      {
        clusterCall(cl, clustAssign, NBpriors, "NBpriors")
        clusterCall(cl, clustAssign, numintSamp, "numintSamp")
      }

    message("Finding posterior likelihoods...", appendLF = FALSE)
    
    for(cc in 1:bootStraps)
      {
        sampPriors <-
          lapply(1:length(groups), function(ii)
                 sampPosteriors[,ii] / sum(sampPosteriors[,ii]))
                 

        if (is.null(cl)) {
          ps <- apply(cbind(1:nrow(cD@data), seglens, cD@data)[union(numintSamp, subset),,drop = FALSE],
                      1, NBbootStrap, libsizes = cD@libsizes, groups = groups, lensameFlag = lensameFlag)
        } else {
          clusterCall(cl, clustAssign, sampPriors, "sampPriors")
          ps <- parRapply(cl, cbind(1:nrow(cD@data), seglens, cD@data)[union(numintSamp, subset),, drop = FALSE],
                          NBbootStrap, libsizes = cD@libsizes, groups = groups, lensameFlag = lensameFlag)
        }
        
        ps <- matrix(ps, ncol = length(groups), byrow = TRUE)
        rps <- matrix(NA, ncol = length(groups), nrow = nrow(cD@data))
        rps[union(numintSamp, subset),] <- ps

        if(pET != "none")
          {
            restprs <- getPosteriors(rps[subset,], prs, pET = pET, groups = groups, priorSubset = priorSubset, cl = cl)$priors
          } else restprs <- prs
        
        pps <- getPosteriors(rps[union(numintSamp,subset),], prs = restprs, pET = "none", groups = groups, priorSubset = priorSubset, cl = cl)
        
        if(any(!is.na(posteriors)))
          if(all(abs(exp(posteriors[union(numintSamp,subset),]) - exp(pps$posteriors)) < conv)) converged <- TRUE
        
        posteriors[union(numintSamp, subset),] <- pps$posteriors

        sampPosteriors <- exp(posteriors[numintSamp,])
        
        prs <- pps$priors
        propest <- rbind(propest, prs)

        estProps <- pps$priors
        names(estProps) <- names(cD@groups)

        cat(".")        

        if(returnAll | converged | cc == bootStraps)
          {
            retPosts <- posteriors
            retPosts[numintSamp[!(numintSamp %in% subset)],] <- NA
            
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

    message("done.")
    
    if(!returnAll) return(listPosts[[cc]]) else {
      if(length(listPosts) == 1) return(listPosts[[1]]) else return(listPosts)
    }
    
  }
