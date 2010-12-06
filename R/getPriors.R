`getPriors.Dirichlet` <-
function(cD, samplesize = 1e5, perSE = 1e-1, maxit = 1e6, verbose = TRUE)
  {
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")

    if(verbose) message("Finding priors...", appendLF = FALSE)
    
    `PgivenDir` <-
      function(alphas, us, ns)
        {
          if(any(alphas <= 0)) return(NA)
          sum(lgamma(alphas[1] + us) + lgamma(alphas[2] + (ns - us)) - lgamma(sum(alphas) + ns) -
              (sum(lgamma(alphas)) - lgamma(sum(alphas))))
        }
    
    priors <- list()
    libsizes <- as.double(cD@libsizes)
    y <- cD@data
    groups <- cD@groups
    priors <- lapply(1:length(groups), function(gg) {
      lapply(unique(groups[[gg]]), function(uu)
             {
               tempPriors <- NULL
               initial <- c(0.2, 1100)
               for(rr in 1:maxit)
                 {
                   if(!is.null(tempPriors))
                     initial <- apply(tempPriors, 2, mean)
                   
                   tempPriors <- rbind(tempPriors, optim(initial, PgivenDir,
                                                         control = list(fnscale = -1),
                                                         us = rowSums(data.frame(y[sample(1:nrow(y), samplesize, replace = FALSE),groups[[gg]] == unique(groups[[gg]])[uu]])),
                                                         ns = sum(libsizes[groups[[gg]] == unique(groups[[gg]])[uu]]))$par)
                   
                   if(nrow(tempPriors) > 1)
                     if(all((apply(tempPriors, 2, sd) / sqrt(nrow(tempPriors))) / apply(tempPriors, 2, mean) < perSE)) break()
                 }
               if(rr == maxit)
                 warning(paste("Convergence not achieved to required accuracy for model ", gg, ", group ", uu, sep = ""))
               if(verbose) message(".", appendLF = FALSE)
               apply(tempPriors, 2, mean)
             })
    })
    if(verbose) message("done.")
    names(priors) <- names(groups)
    new(class(cD), cD, priorType = "Dir", priors = list(priors = priors))
  }


`getPriors.Pois` <-
  function (cD, samplesize = 1e5, perSE = 1e-1
            , takemean = TRUE, maxit = 1e5, weights = NULL, verbose = TRUE, cl) 
{
  if(!inherits(cD, what = "countData"))
    stop("variable 'cD' must be of or descend from class 'countData'")

  if(verbose) message("Finding priors...", appendLF = FALSE)
  
  priorPars <- function(seluu, initial, lensameFlag)
    {
      PgivenPois <- function(priors, us, lens, ns, lensameFlag) {
        if (any(priors <= 0)) 
          return(NA)
        alpha <- priors[1]
        beta <- priors[2]
        sum(alpha * log(beta) - lgamma(alpha) + rowSums(log(outer(as.integer(lens), ns)) * us) - rowSums(lfactorial(us)) + lgamma(alpha + rowSums(us)) - (alpha + rowSums(us)) * log(beta + rowSums(outer(lens, ns))))
      }
      us <- 0
      while(all(us == 0))
        {
          sampcts <- sample((1:nrow(y))[rowSums(is.na(y)) == 0], size = samplesize, replace = FALSE, prob = weights)
          
          if(lensameFlag)
            {
              cts <- y[,-1,drop = FALSE]
              lens <- y[sampcts,1]
            } else {
              cts <- y[,-(1:(ncol(y) / 2)),drop = FALSE]
              len <- y[,1:(ncol(y) / 2), drop = FALSE]
              lens <- len[sampcts, seluu, drop = FALSE]
            }
          
          us <- cts[sampcts, seluu, drop = FALSE]
        }
      ns <- libsizes[seluu]
      
      if(is.null(initial))
        {
          m <- mean(apply(t(us / lens) / ns, 2, mean))
          v <- var(apply(t(us / lens) / ns, 2, mean))
          initial <- c(m^2 / v, m / v)
        }

      if(all(!is.na(initial) & initial > 0 & initial < Inf))
        {
          pars <- optim(initial, 
                        PgivenPois, control = list(fnscale = -1),
                        us = us, lens = lens,
                        ns = ns, lensameFlag = lensameFlag)$par
          return(pars)
        } else return(c(NA, NA))
    }

  if(nrow(cD@seglens) > 0) seglens <- cD@seglens else seglens <- matrix(1, ncol = 1, nrow = nrow(cD@data))
  
  if(ncol(seglens) == 1) lensameFlag <- TRUE else lensameFlag <- FALSE
  
  
  priors <- list()
  libsizes <- as.double(cD@libsizes)
  y <- cbind(seglens, cD@data)
  samplesize = min(samplesize, nrow(y))

  if(!is.null(cl))
    {
      clustAssignData <- function(y, libsizes, samplesize, weights)
        {
          assign("y", y, envir = .GlobalEnv)
          assign("libsizes", libsizes, envir = .GlobalEnv)
          assign("samplesize", samplesize, envir = .GlobalEnv)
          assign("weights", weights, envir = .GlobalEnv)
          NULL
        }
      
      getPriorEnv <- new.env(parent = .GlobalEnv)
      environment(clustAssignData) <- getPriorEnv
      environment(priorPars) <- getPriorEnv

      clusterCall(cl, clustAssignData, y, libsizes, samplesize, weights)
    }
    
    
  groups <- cD@groups
  priors <- lapply(1:length(groups), function(gg) {
    lapply(unique(groups[[gg]]), function(uu) {
      initial <- tempPriors <- NULL  
      seluu <- which(groups[[gg]] == unique(groups[[gg]])[uu])
      for(rr in 1:maxit) {
        if (!is.null(tempPriors)) 
          initial <- apply(tempPriors, 2, mean, na.rm = TRUE)
        if(!is.null(cl))
          {
            pars <- clusterCall(cl, priorPars, seluu = seluu, initial = initial, lensameFlag = lensameFlag)
            pars <- matrix(unlist(pars), ncol = 2, byrow = TRUE)
          } else {
            pars <- matrix(priorPars(seluu = seluu, initial = initial, lensameFlag = lensameFlag), ncol = 2, nrow = 1)
          }
        tempPriors <- rbind(tempPriors, pars[apply(!is.na(pars), 1, all),])
        if(nrow(tempPriors) > 1)
          if(all((apply(tempPriors, 2, sd) / sqrt(nrow(tempPriors))) / apply(tempPriors, 2, mean) < perSE)) break()
      }
      if(rr == maxit)
        warning(paste("Convergence not achieved to required accuracy for model ", gg, ", group ", uu, sep = ""))
      if(verbose) message(".", appendLF = FALSE)
      if(takemean) 
        return(apply(tempPriors, 2, mean))
      else return(tempPriors)
    })
  })
  
  if(verbose) message("done.")
  names(priors) <- names(groups)
  new(class(cD), cD, priorType = "Poi", priors = list(priors = priors))
}


`getPriors.NB` <-
function (cD, samplesize = 1e5, samplingSubset = NULL, equalDispersions = TRUE, estimation = "QL", verbose = TRUE, zeroML = FALSE, cl, ...)
{
  if(!inherits(cD, what = "countData"))
    stop("variable 'cD' must be of or descend from class 'countData'")

  if(verbose) message("Finding priors...", appendLF = FALSE)
  
  optimoverPriors <- function(x, estimation, replicates, groups, libsizes, equalDispersions, lensameFlag, zeroML)
    {
      dispalt <- function(dispersion, cts, mus, libsizes, len)
        abs(2 * (sum(cts[cts > 0] * log(cts[cts > 0]/(mus[cts > 0] * libsizes[cts > 0] * len[cts > 0]))) -
                 sum((cts + dispersion^-1) * log((cts + dispersion^-1) / (mus * libsizes * len + dispersion^-1)))) - (length(cts)-1))
      mualt <- function(mu, cts, dispersion, libsizes, len)
        sum(dnbinom(cts, size = 1/ dispersion, mu = mu * libsizes * len, log = TRUE))

      muZeros <- function(mu, cts, dispersion, libsizes, len)
        {
          if(mu == 0) return(NA) else abs(sum(pnbinom(cts, mu = mu * libsizes * len, size = 1 / dispersion, log = TRUE)) - log(0.5))
        }
      
      dispML <- function(alphas, replicates, y, libsizes, seglen)
        {
          if(any(alphas[-1] < 0) | alphas[1] <= 0) return(NA)
          
          dispersion <- alphas[1]
          mus <- alphas[-1]
          
          sum(sapply(unique(replicates), function(unqrep) sum(dnbinom(y[replicates == unqrep], size = 1 / dispersion,
                                                                      mu = mus[unique(replicates) == unqrep] * libsizes[unique(replicates) == unqrep] * seglen[unique(replicates) == unqrep], log = TRUE))))
        }


      findDisp.QL <- function(repdata)
        {
          replicates <- repdata$replicates
          cts <- repdata$cts
          len <- repdata$len
          libsizes <- as.double(repdata$libsizes)

          if(all(cts == 0))
            return(1)
          
          repmus <- sapply(unique(replicates), function(rep) mean(cts[replicates == rep] / (libsizes * len)[replicates == rep]))
          mus <- c()
          musmat <- matrix(unlist(lapply(unique(replicates), function(rep) rbind(which(replicates == rep), repmus[unique(replicates) == rep]))), ncol = 2, byrow = TRUE)
          mus[musmat[,1]] <- musmat[,2]
          
          newmus <- repmus
          newdisp <- 0
          for(ii in 1:1000) {
            disp <- optimise(dispalt, interval = c(0, 1e2), cts = cts, mus = mus, libsizes = libsizes, len = len, tol = 1e-30)$minimum
            newdisp <- c(newdisp, disp)
            repmus <- rep(0, length(unique(replicates)))
            nz <- sapply(unique(replicates), function(rep) any(cts[replicates == rep] != 0))
            repmus[nz] <- sapply(unique(replicates)[nz], function(rep) optimise(mualt, interval = c(0, 1000), cts = cts[replicates == rep], dispersion = disp, libsizes = libsizes[replicates == rep], len = len[replicates == rep], tol = 1e-50, maximum = TRUE)$maximum)

            musmat <- matrix(unlist(lapply(unique(replicates), function(rep) rbind(which(replicates == rep), repmus[unique(replicates) == rep]))), ncol = 2, byrow = TRUE)
            mus[musmat[,1]] <- musmat[,2]

            newmus <- rbind(newmus, repmus)
            
            if(any(abs(newdisp[-length(newdisp)] - disp) < 1e-50))
              {
                disp <- mean(newdisp[which(abs(newdisp[-length(newdisp)] - disp) < 1e-50):(length(newdisp) - 1)])
                break
              }
          }
          disp
        }

      findDisp.ML <- function(repdata)
        {
          replicates <- repdata$replicates
          cts <- repdata$cts
          len <- repdata$len
          libsizes <- as.double(repdata$libsizes)
          
          if(all(cts == 0))
            return(1)


          repmus <- sapply(unique(replicates), function(rep) mean(cts[replicates == rep] / (libsizes * len)[replicates == rep]))
          
          optim(par = c(0.01, repmus),
                fn = dispML, control = list(fnscale = -1), 
                y = cts, replicates = replicates, libsizes = libsizes, seglen = len)$par[1]
        }
  
      if(lensameFlag)
        {
          cts <- x[-1]
          len <- x[1]
        } else {
          cts <- x[-(1:(length(x) / 2))]
          len <- x[1:(length(x) / 2)]
        }
      if(length(len) == 1) len <- rep(len, length(cts))

      if(equalDispersions)
        {
          disp <- switch(estimation,
                         QL = findDisp.QL(list(replicates = replicates, cts = cts, len = len, libsizes = libsizes)),
                         ML = findDisp.ML(list(replicates = replicates, cts = cts, len = len, libsizes = libsizes))
                         )
          
          groupness <- lapply(groups, function(group) {
            dispersion <- rep(disp, length(unique(group)))
            mus <- sapply(unique(group), function(unqgrp)
              if(all(cts[group == unqgrp] == 0) & zeroML) {
                return(0)
              } else if(all(cts[group == unqgrp] == 0)) {
                return(
                       optimize(muZeros, interval = c(0, max(1 / libsizes[group == unqgrp] * len[group == unqgrp])),
                                cts = cts[group == unqgrp], dispersion = 1, libsizes = libsizes[group == unqgrp], len = len[group == unqgrp], tol = 1e-50, maximum = FALSE)$minimum
                       )
              } else return(optimise(mualt, interval = c(0, 1), cts = cts[group == unqgrp], dispersion = disp, libsizes = libsizes[group == unqgrp], len = len[group == unqgrp], tol = 1e-50, maximum = TRUE)$maximum))
            list(dispersion = dispersion, mus = mus)
          })
                        
          
        } else {
          repgroups <- lapply(groups, function(group) list(group = group, repData = lapply(unique(group), function(unqgrp) list(replicates = replicates[group == unqgrp], cts = cts[group == unqgrp], len = len[group == unqgrp], libsizes = libsizes[group == unqgrp]))))
          
          dispgroups <- switch(estimation,
                               QL = lapply(repgroups, function(repgroup) list(group = repgroup$group, dispersion = sapply(repgroup$repData, findDisp.QL))),
                               ML = lapply(repgroups, function(repgroup) list(group = repgroup$group, dispersion = sapply(repgroup$repData, findDisp.ML)))
                               )
          
          groupness <- lapply(dispgroups, function(dispgroup)
                              list(dispersion = dispgroup$dispersion,
                                   mus = sapply(unique(dispgroup$group),
                                     function(unqgrp)
                                     if(all(cts[dispgroup$group == unqgrp] == 0) & zeroML) {
                                       return(0)
                                     } else if(all(cts[dispgroup$group == unqgrp] == 0)) {
                                       return(
                                              optimize(muZeros, interval = c(0, max(1 / libsizes[dispgroup$group == unqgrp] * len[dispgroup$group == unqgrp])),
                                                       cts = cts[dispgroup$group == unqgrp], dispersion = 1, libsizes = libsizes[dispgroup$group == unqgrp], len = len[dispgroup$group == unqgrp], tol = 1e-50, maximum = FALSE)$minimum)
                                     } else return(optimise(mualt, interval = c(0, 1000), cts = cts[dispgroup$group == unqgrp], dispersion = dispgroup$dispersion[unique(dispgroup$group) == unqgrp], libsizes = libsizes[dispgroup$group == unqgrp], len = len[dispgroup$group == unqgrp], tol = 1e-50, maximum = TRUE)$maximum))))
        }
      
      groupness
      
    }
  
  if(!is.null(cl))
    {
      getPriorsEnv <- new.env(parent = .GlobalEnv)
      environment(optimoverPriors) <- getPriorsEnv
    }
  
  
  if(is.null(samplingSubset))
    samplingSubset <- 1:nrow(cD)
  
  sD <- cD[samplingSubset,]
  
  libsizes <- as.double(sD@libsizes)
  groups <- sD@groups
  replicates <- sD@replicates

  if(nrow(sD@seglens) > 0) seglens <- sD@seglens[,,drop = TRUE] else seglens <- rep(1, nrow(sD@data))
  if(is.vector(seglens)) lensameFlag <- TRUE else lensameFlag <- FALSE

  tupData <- cbind(sD@data, seglens)
  
  ordData <- do.call(order, as.data.frame(tupData))
  dups <- c(1, which(rowSums(tupData[ordData[-1],, drop = FALSE] == (tupData)[ordData[-length(ordData)],, drop = FALSE]) != ncol(tupData)) + 1)

  if(zeroML)
    {
      allzeros <- dups[which(rowSums(tupData[ordData[dups],-ncol(tupData), drop = FALSE]) == 0)]
      if(length(allzeros) > 0)
        dups <- c(1, dups[!(dups %in% allzeros)])
    }
  
  
  if(length(dups) <= samplesize) {
    y <- sD@data[ordData[dups],,drop = FALSE]
    copies <- diff(c(dups, nrow(sD@data) + 1))
    sampled <- (1:nrow(sD))[ordData]
    sy <- cbind(sampled = sampled, representative = rep(1:length(copies), copies))
    weights <- rep(1, nrow(sy))
    if(lensameFlag) seglensy <- seglens[ordData[dups]] else seglensy <- seglens[ordData[dups],]
  } else {
    sampData <- colSums(t(sD@data/seglens) / sD@libsizes) / ncol(sD)

    if(zeroML) nzSD <- sampData[sampData > 0] else nzSD <- sampData
      
    sqnum <- length(nzSD) / 1000
    squant <- quantile(nzSD, 1:sqnum / sqnum)
    sqdup <- c(1, which(diff(squant) > min(1 / libsizes) / 10))
    z <- cbind(as.numeric(squant[sqdup]), c(as.numeric(squant[sqdup[-1]]), max(nzSD)))
    if(zeroML) z[1,1] <- 0 else z[1,1] <- -Inf

    sy <- apply(z, 1, function(w) {
      inbetweener <- which(sampData > w[1] & sampData <= w[2])
      samplenum <- min(length(inbetweener), ceiling(samplesize / nrow(z)))
      rbind(sample(inbetweener, size = samplenum, replace = FALSE), length(inbetweener) / samplenum)
    })

    sy <- matrix(as.vector(sy), ncol = 2, byrow = TRUE)
    weights <- sy[,2]
    sy <- sy[,1]

    y <- sD@data[sy,,drop = FALSE]
    if(lensameFlag) seglensy <- seglens[sy] else seglensy <- seglens[sy,,drop = FALSE]
    
    ordData <- do.call(order, as.data.frame(cbind(y, seglensy)))
    y <- y[ordData,, drop = FALSE]
    weights <- weights[ordData]
    if(lensameFlag) seglensy <- seglensy[ordData] else seglensy <- seglensy[ordData,,drop = FALSE]
    
    dups <- c(1, which(rowSums((cbind(y,seglensy))[-1,] == (cbind(y, seglensy))[-nrow(y),]) != ncol(cbind(y, seglensy))) + 1)
    copies <- diff(c(dups, nrow(y) + 1))

    sy <- cbind(sampled = sy[ordData], representative = rep(1:length(copies), copies))

    y <- y[dups,, drop = FALSE]
    if(lensameFlag) seglensy <- seglens[dups] else seglensy <- seglens[dups,, drop = FALSE]
    
    if(any(sampData == 0) & zeroML)
      {
        sy <- rbind(sy, cbind(sampled = which(sampData == 0), representative = nrow(y) + 1))
        y <- rbind(y, rep(0, ncol(y)))
        if(lensameFlag) seglensy <- c(seglensy, 1) else seglensy <- rbind(seglensy, rep(1, ncol(seglensy)))
        weights <- c(weights, rep(1, sum(sampData == 0)))
      }
  }

  z <- cbind(seglensy, y)

  NBpar <- list()


  if(estimation == "edgeR")
    {
      dge <- new("DGEList")
      dge$counts = y
      dge$samples = data.frame(group = replicates, lib.size = libsizes)
      dge <- estimateCommonDisp(dge)
      dge <- estimateTagwiseDisp(dge, ...)
      disps <- dge$tagwise.dispersion

      parEach <- apply(cbind(disps, z), 1, function(x)
            {
              disp <- x[1]
              x <- x[-1]
              if(lensameFlag)
                {
                  cts <- x[-1]
                  len <- x[1]
                } else {
                  cts <- x[-(1:(length(x) / 2))]
                  len <- x[1:(length(x) / 2)]
                }
              if(length(len) == 1) len <- rep(len, length(cts))

              mualt <- function(mu, cts, dispersion, libsizes, len)
                sum(dnbinom(cts, size = 1/ dispersion, mu = mu * libsizes * len, log = TRUE))
              
              groupness <- lapply(groups, function(group) list(dispersion = rep(disp, length(unique(group))), mus = sapply(unique(group), function(unqgrp) optimise(mualt, interval = c(0, 1000), cts = cts[group == unqgrp], dispersion = disp, libsizes = libsizes[group == unqgrp], len = len[group == unqgrp], tol = 1e-50, maximum = TRUE)$maximum)))
            })
      
    } else {
      if(is.null(cl)) {
        parEach <- apply(z, 1, optimoverPriors, estimation = estimation, replicates = replicates, groups = groups, libsizes = libsizes, equalDispersions = equalDispersions, lensameFlag = lensameFlag, zeroML = zeroML)
      } else parEach <- parApply(cl, z, 1, optimoverPriors, estimation = estimation, replicates = replicates, groups = groups, libsizes = libsizes, equalDispersions = equalDispersions, lensameFlag = lensameFlag, zeroML)
    }

  if(!is.null(cl))
    clusterEvalQ(cl, rm(list = ls()))

  
  NBpar <- lapply(1:length(groups), function(gg)
                  lapply(unique(groups[[gg]]), function(ii) t(sapply(parEach, function(x) c(x[[gg]]$mus[ii], c(x[[gg]]$dispersion[ii], 1)[as.numeric(is.na(x[[gg]]$dispersion[ii])) + 1])))))
  
  if(verbose) message("done.")

  sy[,1] <- samplingSubset[sy[,1]]
  names(NBpar) <- names(groups)
  NBpar <- list(sampled = sy, weights = weights, priors = NBpar)
  new(class(cD), cD, priorType = "NB", priors = NBpar)
    }
