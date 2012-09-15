`getPriors.NB` <-
function (cD, samplesize = 1e5, samplingSubset = NULL, equalDispersions = TRUE, estimation = "QL", verbose = TRUE, zeroML = FALSE, consensus = FALSE, cl, ...)
{
  if(!inherits(cD, what = "countData"))
    stop("variable 'cD' must be of or descend from class 'countData'")

  if(length(cD@libsizes) == 0)
    {
      warning("'@libsizes' slot empty; inferring libsizes using default settings")
      cD@libsizes <- getLibsizes(cD)
    }
  
  
  if(verbose) message("Finding priors...", appendLF = FALSE)
  
  optimoverPriors <- function(x, estimation, replicates, groups, libsizes, equalDispersions, lensameFlag, zeroML, consensus, disp)
    {
      dispalt <- function(dispersion, cts, mus, libsizes, len)
        abs(2 * (sum(cts[cts > 0] * log(cts[cts > 0]/(mus[cts > 0] * libsizes[cts > 0] * len[cts > 0])), na.rm = TRUE) -
                 sum((cts + dispersion^-1) * log((cts + dispersion^-1) / (mus * libsizes * len + dispersion^-1)), na.rm = TRUE)) - (length(cts[!is.na(cts)])-1))
      mualt <- function(mu, cts, dispersion, libsizes, len)
        sum(dnbinom(cts, size = 1/ dispersion, mu = mu * libsizes * len, log = TRUE), na.rm = TRUE)

      muZeros <- function(mu, cts, dispersion, libsizes, len)
        {
          if(mu == 0) return(NA) else abs(sum(pnbinom(cts, mu = mu * libsizes * len, size = 1 / dispersion, log.p = TRUE), na.rm = TRUE) - log(0.5))
        }
      
      dispML <- function(alphas, replicates, y, libsizes, seglen)
        {
          if(any(alphas[-1] < 0) | alphas[1] <= 0) return(NA)
          
          dispersion <- alphas[1]
          
          mus <- alphas[-1]
          
          sum(sapply(levels(replicates), function(unqrep) sum(dnbinom(y[replicates == unqrep], size = 1 / dispersion,
                                                                      mu = mus[levels(replicates) == unqrep] * libsizes[levels(replicates) == unqrep] * seglen[levels(replicates) == unqrep], log = TRUE), na.rm = TRUE)))
        }

      findDisp.QL <- function(repdata)
        {
          replicates <- repdata$replicates
          cts <- repdata$cts
          len <- as.double(repdata$len)
          libsizes <- as.double(repdata$libsizes)
          searchReps <- levels(replicates)[levels(replicates) %in% replicates]
          
          if(all(cts == 0) | (sum(cts > 0) == 1 & sum(replicates %in% replicates[cts > 0]) == 1))
            return(1e-100)

          repmus <- sapply(searchReps, function(rep) return(mean(cts[replicates == rep] / (libsizes * len)[replicates == rep])))

          mus <- c()

          z <- unlist(lapply(searchReps, function(rep) rbind(which(replicates == rep), repmus[searchReps == rep])))
          
          musmat <- matrix(unlist(z), ncol = 2, byrow = TRUE)

          mus[musmat[,1]] <- musmat[,2]

          newmus <- repmus
          newdisp <- 0
          for(ii in 1:1000) {
            lowInt <- 0
            disp <- optimise(dispalt, interval = c(lowInt, 1e2), cts = cts, mus = mus, libsizes = libsizes, len = len, tol = 1e-50)$minimum

            newdisp <- c(newdisp, disp)
            repmus <- rep(0, length(searchReps))
            nz <- sapply(searchReps, function(rep) any(cts[replicates == rep] != 0))
            repmus[nz] <- sapply(searchReps[nz], function(rep) optimise(mualt, interval = c(0,  max(cts[replicates == rep] / libsizes[replicates == rep], na.rm = TRUE) * 2), cts = cts[replicates == rep], dispersion = disp, libsizes = libsizes[replicates == rep], len = len[replicates == rep], tol = 1e-50, maximum = TRUE)$maximum)

            musmat <- matrix(unlist(lapply(searchReps, function(rep) rbind(which(replicates == rep), repmus[searchReps == rep]))), ncol = 2, byrow = TRUE)
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
          searchReps <- levels(replicates)[levels(replicates) %in% replicates]
          
          if(all(cts == 0) | (sum(cts > 0) == 1 & sum(replicates == replicates[cts == 1]) == 1))
            return(1e-100)


          repmus <- sapply(searchReps, function(rep) mean(cts[replicates == rep] / (libsizes * len)[replicates == rep]))
          
          disp <- optim(par = c(1, repmus),
                fn = dispML, control = list(fnscale = -1), 
                y = cts, replicates = replicates, libsizes = libsizes, seglen = len)$par[1]

          disp
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

      getMu <- function(samples, disp) {
        if(all(cts[samples] == 0, na.rm = TRUE) & zeroML) {
          return(0)
        } else if(all(cts[samples] == 0, na.rm = TRUE)) {
          return(
                 optimize(muZeros, interval = c(0, max(1 / libsizes[samples] * len[samples], na.rm = TRUE)),
                          cts = cts[which(samples)], dispersion = 1, libsizes = libsizes[which(samples)], len = len[which(samples)], tol = 1e-50, maximum = FALSE)$minimum
                 )
        } else {
          return(optimise(mualt, interval = c(0, max(cts[samples] / libsizes[samples], na.rm = TRUE) * 2), cts = cts[which(samples)], dispersion = disp, libsizes = libsizes[which(samples)], len = len[which(samples)], tol = 1e-50, maximum = TRUE)$maximum)
        }
      }
      
      if(equalDispersions | consensus)
        {          
          disp <- switch(estimation,
                         QL = findDisp.QL(list(replicates = replicates, cts = cts, len = len, libsizes = libsizes)),
                         ML = findDisp.ML(list(replicates = replicates, cts = cts, len = len, libsizes = libsizes))
                         )
          if(consensus) {
            mu <- getMu(replicates == sample(levels(replicates), size = 1), disp = disp)
            groupness <- c(mu, disp)
          } else {
            groupness <- lapply(groups, function(group) {
              dispersion <- rep(disp, length(levels(group)))
              mus <- sapply(levels(group), function(unqgrp) getMu(group == unqgrp, disp = disp))
              list(dispersion = dispersion, mus = mus)
            })}
         } else {
           repgroups <- lapply(groups, function(group) list(group = group, repData = lapply(levels(group[!is.na(group)]), function(unqgrp) list(replicates = replicates[which(group == unqgrp)], cts = cts[which(group == unqgrp)], len = len[which(group == unqgrp)], libsizes = libsizes[which(group == unqgrp)]))))
           
          dispgroups <- switch(estimation,
                               QL = lapply(repgroups, function(repgroup) list(group = repgroup$group, dispersion = sapply(repgroup$repData, findDisp.QL))),
                               ML = lapply(repgroups, function(repgroup) list(group = repgroup$group, dispersion = sapply(repgroup$repData, findDisp.ML)))
                               )
          groupness <- lapply(dispgroups, function(dispgroup)
                              list(dispersion = dispgroup$dispersion,                                                                     
                                   mus = sapply(levels(dispgroup$group[!is.na(dispgroup$group)]), function(unqgrp) getMu(dispgroup$group == unqgrp, dispgroup$dispersion[levels(dispgroup$group) == unqgrp]))))
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

  samplingSubset <- samplingSubset[rowSums(do.call("cbind", lapply(cD@groups, function(x) do.call("cbind", lapply(levels(x), function(rep) rowSums(is.na(cD@data[,x == rep,drop = FALSE])) == sum(x == rep)))))) == 0]

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
  if(zeroML)
    {
      allzeros <- dups[which(rowSums(tupData[ordData[dups],-ncol(tupData), drop = FALSE], na.rm = TRUE) == 0)]
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

  NBpar <- list()

  if(estimation == "DSS")
    {
      if(length(levels(replicates)) != 2) stop("The DSS package only allows for pairwise comparisons at present. Use another method of parameter estimation")
      designs <- as.integer(replicates == levels(replicates)[2])
      colnames(y) <- NULL
      seqData <- newSeqCountSet(y, designs)
      seqData@normalizationFactor <- libsizes / min(libsizes)
      seqData <- estDispersion(seqData)
      dispersions <- seqData@dispersion

      parEach <- apply(cbind(dispersions, y), 1, function(z) {
        disp <- z[1]
        cts <- z[-1]

        mualt <- function(mu, cts, dispersion, libsizes, len = 1)
          sum(dnbinom(cts, size = 1/ dispersion, mu = mu * libsizes * len, log = TRUE), na.rm = TRUE)
        
        muZeros <- function(mu, cts, dispersion, libsizes, len = 1)
            if(mu == 0) return(NA) else abs(sum(pnbinom(cts, mu = mu * libsizes * len, size = 1 / dispersion, log.p = TRUE), na.rm = TRUE) - log(0.5))                  
        
        getMu <- function(samples, disp) {
          if(all(cts[samples] == 0, na.rm = TRUE) & zeroML) {
            return(0)
          } else if(all(cts[samples] == 0, na.rm = TRUE)) {
            return(
                   optimize(muZeros, interval = c(0, max(1 / libsizes[samples], na.rm = TRUE)),
                            cts = cts[which(samples)], dispersion = 1, libsizes = libsizes[which(samples)], tol = 1e-50, maximum = FALSE)$minimum
                   )
          } else {
            return(optimise(mualt, interval = c(0, max(cts[samples] / libsizes[samples], na.rm = TRUE) * 2), cts = cts[which(samples)], dispersion = disp, libsizes = libsizes[which(samples)], tol = 1e-50, maximum = TRUE)$maximum)
          }
        }
            
        groupness <- lapply(groups, function(group) {
          dispersion <- disp
          mus <- sapply(levels(group), function(unqgrp) getMu(group == unqgrp, disp = disp))
          list(dispersion = dispersion, mus = mus)
        })        
      })
      NBpar <- lapply(1:length(groups), function(gg)
                      lapply(1:length(levels(groups[[gg]])), function(ii) t(sapply(parEach, function(x) c(x[[gg]]$mus[ii], c(x[[gg]]$dispersion[ii], 1)[as.numeric(is.na(x[[gg]]$dispersion[ii])) + 1])))))    
    } else if(estimation == "edgeR")
    {
      if(!("edgeR" %in% loadedNamespaces()))
        library(edgeR)
      dge <- new("DGEList")
      dge$counts = y
      dge$samples = data.frame(group = as.numeric(replicates), lib.size = libsizes)
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
              
              groupness <- lapply(groups, function(group) list(dispersion = rep(disp, length(levels(group))), mus = sapply(levels(group), function(unqgrp) optimise(mualt, interval = c(0, 1000), cts = cts[which(group == unqgrp)], dispersion = disp, libsizes = libsizes[which(group == unqgrp)], len = len[which(group == unqgrp)], tol = 1e-50, maximum = TRUE)$maximum)))
            })
      
    } else {
      if(is.null(cl)) {
        parEach <- apply(z, 1, optimoverPriors, estimation = estimation, replicates = replicates, groups = groups, libsizes = libsizes, equalDispersions = equalDispersions, lensameFlag = lensameFlag, zeroML = zeroML, consensus = consensus)
      } else parEach <- parApply(cl, z, 1, optimoverPriors, estimation = estimation, replicates = replicates, groups = groups, libsizes = libsizes, equalDispersions = equalDispersions, lensameFlag = lensameFlag, zeroML, consensus = consensus)
    }

  if(!is.null(cl))
    clusterEvalQ(cl, rm(list = ls()))

  if(consensus) NBpar <- t(parEach) else {
    NBpar <- lapply(1:length(groups), function(gg)
                    lapply(1:length(levels(groups[[gg]])), function(ii) t(sapply(parEach, function(x) c(x[[gg]]$mus[ii], c(x[[gg]]$dispersion[ii], 1)[as.numeric(is.na(x[[gg]]$dispersion[ii])) + 1])))))
  }
  
  if(verbose) message("done.")

  sy[,1] <- samplingSubset[sy[,1]]
  names(NBpar) <- names(groups)
  NBpar <- list(sampled = sy, weights = weights, priors = NBpar)
  new(class(cD), cD, priorType = "NB", priors = NBpar)
}
