`getLikelihoods.Dirichlet` <-
function(cD, prs, estimatePriors = TRUE, subset = NULL, cl)
  {
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")

    if(cD@priors$type != "Dir") stop("Incorrect prior type for this method of likelihood estimation")
    if(class(cD) != "countData")
      stop("variable 'countData' in 'getGroupPriors' must be of class 'countData'")
    if(length(prs) != length(cD@groups)) stop("'prs' must be of same length as the number of groups in the 'cD' object")
    if(!(class(subset) == "integer" | class(subset) == "numeric" | is.null(subset)))
      stop("'subset' must be integer, numeric, or NULL")

    if(is.null(subset)) subset <- 1:nrow(cD@data)
    
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
                
                pD <- 0
                for(ii in 1:nrow(us))
                  pD <- pD +
                    lfactorial(sum(us[ii,])) - sum(lfactorial(us[ii,]))
                for(gg in 1:length(unique(group)))
                  pD <- pD +
                    lbeta.over.beta(apply(matrix(us[group == unique(group)[gg],], ncol = ncol(us)), 2, sum),
                                    prior[gg,])
                pD
              }

          
          us <- cbind(us, libsizes - us)
          ps <- c()
          for(ii in 1:length(groups))
            ps[ii] <- PDgivenr.Dirichlet(us, priors[[ii]], groups[[ii]])
          
          ps
        }

    if(!is.null(cl))
      {
        getLikelihoodsEnv <- new.env(parent = .GlobalEnv)
        environment(PrgivenD.Dirichlet) <- getLikelihoodsEnv
      }
    
    if(is.null(cl)) {
      ps <- apply(cD@data[subset,], 1, PrgivenD.Dirichlet, cD@libsizes, cD@priors$priors, cD@groups)
    } else {
      ps <- parRapply(cl, cD@data[subset,], PrgivenD.Dirichlet, cD@libsizes, cD@priors$priors, cD@groups)
    }
    ps <- matrix(ps, ncol = length(cD@groups), byrow = TRUE)
    pps <- getPosteriors(ps, prs, estimatePriors = estimatePriors, cl = cl)

    posteriors <- matrix(NA, ncol = length(cD@groups), nrow(cD@data))
    posteriors[subset,] <- t(pps$posteriors)

    colnames(posteriors) <- names(cD@groups)
    estProps <- pps$priors
    names(estProps) <- names(cD@groups)
    
    new(class(cD), cD, posteriors = posteriors, estProps = estProps)
  }


`getLikelihoods.Pois` <-
function(cD, prs, estimatePriors = TRUE, subset = NULL, distpriors = FALSE, cl)
  {
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")

    if(cD@priors$type != "Poi") stop("Incorrect prior type for this method of likelihood estimation")
    if(length(prs) != length(cD@groups)) stop("'prs' must be of same length as the number of groups in the 'cD' object")
    
    if(!(class(subset) == "integer" | class(subset) == "numeric" | is.null(subset)))
      stop("'subset' must be integer, numeric, or NULL")

    if(is.null(subset)) subset <- 1:nrow(cD@data)

    if(class(cD) == "segData") seglens <- cD@seglens else seglens <- matrix(1, ncol = 1, nrow = nrow(cD@data))

    if(is.matrix(seglens))
      if(ncol(seglens) == 1) seglens <- matrix(seglens[,1], ncol = ncol(cD@data), nrow = nrow(cD@data))

    
    `PrgivenD.Pois` <-
      function(us, libsizes, priors, groups, distpriors = FALSE)
        {
          `PDgivenr.Pois` <-
            function(cts, ns, seglen, prior, group)
              {
                pD <- 0
                pD <- sum(cts * log(ns * seglen)) - sum(lfactorial(cts))
                for(gg in 1:length(unique(group)))
                  {
                    selcts <- group == gg
                    pD <- pD +
                      lgamma(sum(cts[selcts]) + prior[gg,1]) - lgamma(prior[gg,1]) -
                        sum(cts[selcts]) * log(sum(ns[selcts] * seglen[selcts]) + prior[gg,2]) -
                          prior[gg,1] * log(1 + sum(ns[selcts] * seglen[selcts]) / prior[gg,2])
                  }
                pD
              }
          
          `PDgivenr.PoisIndie` <-
            function(cts, ns, seglen, prior, group)
              {
                `logsum` <-
                  function(x)
                    max(x, max(x, na.rm = TRUE) + log(sum(exp(x - max(x, na.rm = TRUE)), na.rm = TRUE)), na.rm = TRUE)
                
                pD <- sum(cts * log(ns * seglen)) - sum(lfactorial(cts))
                for(gg in 1:length(unique(group)))
                  {
                    selcts <- group == gg
                    pD <- pD +
                      logsum(lgamma(sum(cts[selcts]) + prior[[gg]][,1]) - lgamma(prior[[gg]][,1]) -
                             sum(cts[selcts]) * log(sum(ns[selcts] * seglen[selcts]) + prior[[gg]][,2]) -
                             prior[[gg]][,1] * log(1 + sum(ns[selcts] * seglen[selcts]) / prior[[gg]][,2])) - log(nrow(prior[[gg]]))
                  }
                pD
              }

          
          cts <- us[-(1:(length(us) / 2))]
          lens <- us[1:(length(us) / 2)]
          
          ps <- c()
          for(ii in 1:length(groups))
            if(!distpriors) ps[ii] <-
              PDgivenr.Pois(cts = cts, ns = libsizes, seglen = lens, prior = priors[[ii]], group = groups[[ii]]) else ps[ii] <-
                PDgivenr.PoisIndie(us = us, ns = libsizes, seglen = lens, prior = priors[[ii]], group = groups[[ii]])
          ps
        }

    if(!is.null(cl))
      {
        getLikelihoodsEnv <- new.env(parent = .GlobalEnv)
        environment(PrgivenD.Pois) <- getLikelihoodsEnv
      }
    
    if(is.null(cl)) {
      ps <- apply(cbind(seglens, cD@data)[subset,], 1, PrgivenD.Pois, cD@libsizes, cD@priors$priors, cD@groups, distpriors)
    } else {
      ps <- parRapply(cl, cbind(seglens, cD@data)[subset,], PrgivenD.Pois, cD@libsizes, cD@priors$priors, cD@groups, distpriors)
      }                        
    ps <- matrix(ps, ncol = length(cD@groups), byrow = TRUE)
    pps <- getPosteriors(ps, prs, estimatePriors = estimatePriors, cl = cl)

    posteriors <- matrix(NA, ncol = length(cD@groups), nrow(cD@data))
    posteriors[subset,] <- pps$posteriors

    colnames(posteriors) <- names(cD@groups)
    estProps <- pps$priors
    names(estProps) <- names(cD@groups)
    
    new(class(cD), cD, posteriors = posteriors, estProps = estProps)
  }


`getLikelihoods.NBboot` <-
function(cD, prs, estimatePriors = TRUE, subset = NULL, bootStraps = 2, conv = 1e-4, nullData = FALSE, cl)
  {
    
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")

    if(cD@priors$type != "NB") stop("Incorrect prior type for this method of likelihood estimation")
    if(length(prs) != length(cD@groups)) stop("'prs' must be of same length as the number of groups in the 'cD' object")

    if(any(prs < 0))
      stop("Negative values in the 'prs' vector are not permitted")
    
    if(!nullData & sum(prs) != 1)
      stop("If 'nullData = FALSE' then the 'prs' vector should sum to 1.")

    if(nullData & sum(prs) >= 1)
            stop("If 'nullData = TRUE' then the 'prs' vector should sum to less than 1.")

    if(!(class(subset) == "integer" | class(subset) == "numeric" | is.null(subset)))
      stop("'subset' must be integer, numeric, or NULL")
    
    if(is.null(subset)) subset <- 1:nrow(cD@data)
    sy <- cD@priors$sampled
    groups <- cD@groups

    NBpriors <- cD@priors$priors

    NBbootStrap <- function(us, libsizes, groups) {
      PDgivenr.NB <- function (us, seglen, ns, prior, group, sampled)
          {
            `logsum` <-
              function(x)
                max(x, max(x, na.rm = TRUE) + log(sum(exp(x - max(x, na.rm = TRUE)), na.rm = TRUE)), na.rm = TRUE)

            pD <- 0
            for (gg in 1:length(unique(group))) {
              selus <- group == gg
              pD <- pD + logsum(
                                rowSums(matrix(
                                               dnbinom(rep(us[selus], each = nrow(prior[[gg]])),
                                                       size = 1 / prior[[gg]][,2],
                                                       mu = rep(ns[selus] * seglen[selus], each = nrow(prior[[gg]])) * prior[[gg]][,1], log = TRUE),
                                               ncol = sum(selus))
                                        ) + log(sampled))
            }
            pD
          }

      
      ps <- c()
      for (ii in 1:length(NBpriors))
        ps[ii] <- PDgivenr.NB(us[-(1:(length(us) / 2))], us[1:(length(us) / 2)], libsizes, NBpriors[[ii]], groups[[ii]], sampPriors[[ii]])
      ps
    }
    
    clustAssign <- function(object, name)
      {
        assign(name, object, envir = .GlobalEnv)
        NULL
      }

    if(!is.null(cl))
      {
        getLikelihoodsEnv <- new.env(parent = .GlobalEnv)
        environment(clustAssign) <- getLikelihoodsEnv
        environment(NBbootStrap) <- getLikelihoodsEnv
      }
    
    if(is.null(conv)) conv <- 0
    posteriors <- matrix(1 / length(cD@groups), ncol = length(cD@groups), nrow = nrow(cD@data))
    propest <- NULL

    if(class(cD) == "segData") seglens <- cD@seglens else seglens <- matrix(1, ncol = 1, nrow = nrow(cD@data))

    if(is.matrix(seglens))
      if(ncol(seglens) == 1) seglens <- matrix(seglens[,1], ncol = ncol(cD@data), nrow = nrow(cD@data))
    
    if(nullData)
      {
        ndelocGroup <- which(unlist(lapply(groups, function(x) all(x == x[1]))))
        if(length(ndelocGroup) == 0)
          stop("If 'nullData = TRUE' then there must exist some vector in groups whose members are all identical")
        groups <- c(groups, list(rep(1, ncol(cD@data))))
        ndenulGroup <- length(groups)
        prs <- c(prs, 1 - sum(prs))

        NBpriors[[ndenulGroup]] <- NBpriors[[ndelocGroup]]
        
        NZLs <- apply(cD@data[sy,, drop = FALSE] != 0, 1, any)
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
        km <- rep(1, length(sy))
        km[NZLs] <- kmNZ$cluster

        sampPosteriors <- matrix(NA, nrow = length(sy), ncol = length(groups))
        for(ii in 1:length(groups))
          sampPosteriors[,ii] <- c(1, 0)[as.numeric(km == c(1,2)[as.numeric(ii == ndenulGroup) + 1]) + 1]
      } else {
        sampPosteriors <- matrix(NA, nrow = length(sy), ncol = length(groups))
        for(ii in 1:length(groups))
          sampPosteriors[,ii] <- 1
      }

    posteriors <- matrix(NA, ncol = length(groups), nrow = nrow(cD@data))
    propest <- NULL
    converged <- FALSE

    if(!is.null(cl))
      clusterCall(cl, clustAssign, NBpriors, "NBpriors")
    
    for(cc in 1:bootStraps)
      {
        sampPriors <- list()
        for(ii in 1:length(groups))
          sampPriors[[ii]] <- sampPosteriors[,ii] / sum(sampPosteriors[,ii])

        if (is.null(cl)) {
          ps <- apply(cbind(seglens, cD@data)[union(sy, subset),,drop = FALSE],
                      1, NBbootStrap, libsizes = cD@libsizes, groups = groups)
        } else {
          clusterCall(cl, clustAssign, sampPriors, "sampPriors")
          ps <- parRapply(cl, cbind(seglens, cD@data)[union(sy, subset),, drop = FALSE],
                          NBbootStrap, libsizes = cD@libsizes, groups = groups)
        }
        
        ps <- matrix(ps, ncol = length(groups), byrow = TRUE)
        rps <- matrix(NA, ncol = length(groups), nrow = nrow(cD@data))
        rps[union(sy, subset),] <- ps
        
        pps <- getPosteriors(rps[subset,], prs, estimatePriors = TRUE, cl = cl)
        print(pps$priors)
        pps <- getPosteriors(rps[union(sy,subset),], pps$priors, estimatePriors = FALSE, cl = cl)

        if(any(!is.na(posteriors)))
          if(all(abs(exp(posteriors[union(sy,subset),]) - exp(pps$posteriors)) < conv)) converged <- TRUE
        
        posteriors[union(sy, subset),] <- pps$posteriors

        sampPosteriors <- exp(posteriors[sy,])
        
        prs <- pps$priors
        propest <- rbind(propest, prs)

        estProps <- pps$priors
        names(estProps) <- names(cD@groups)

        if(converged)
          break()
      }

    posteriors[sy[!(sy %in% subset)],] <- NA

    nullPosts <- numeric(0)
    if(nullData) {
      nullPosts <- posteriors[,ndenulGroup]
      estProps <- pps$priors[-ndenulGroup]
      posteriors <- posteriors[,-ndenulGroup, drop = FALSE]
    }
    colnames(posteriors) <- names(cD@groups)
    return(new(class(cD), cD, posteriors = posteriors, estProps = estProps, nullPosts = nullPosts))
  }


`getPosteriors` <-
function(ps, prs, estimatePriors = FALSE, maxit = 100, accuracy = 1e-5, cl = cl)
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

    if(!is.null(cl))
      {
        getPostsEnv <- new.env(parent = .GlobalEnv)
        environment(getPosts) <- getPostsEnv
      }
    
    if(estimatePriors)
      {
        oldprs <- prs
        for(ii in 1:maxit)
          {
            if(!is.null(cl)) {
              posteriors <- matrix(parRapply(cl, ps, getPosts, prs), ncol = ncol(ps), byrow = TRUE)
            } else posteriors <- matrix(apply(ps, 1, getPosts, prs), ncol = ncol(ps), byrow = TRUE)
            prs <- colSums(exp(posteriors)) / nrow(posteriors)
            if(all(abs(oldprs - prs) < accuracy)) break
            oldprs <- prs
          }
        if(ii == maxit)
          warning("Convergence not achieved to required accuracy.")
      }
list(posteriors = posteriors <- matrix(apply(ps, 1, getPosts, prs), ncol = ncol(ps), byrow = TRUE), priors = prs)
  }
