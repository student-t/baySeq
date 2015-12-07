
`getPriors` <- function (cD, samplesize = 1e5, samplingSubset = NULL, verbose = TRUE, consensus = FALSE, cl)
{

  if(is.null(body(cD@densityFunction@density))) {
    message("No density provided in cD object; using negative-binomial as default")
    densityFunction(cD) <- nbinomDensity
  }
  
  .getSinglePriors <- function(x, datdim, replicates, groups, optimFunction, eqovRep, initiatingValues, lower, upper, consensus)#, dataLL)
    {
      xid <- x$id
#      message(xid)
      dat <- x$data
#      cat(paste(xid, ",", sep = ""), file = "bugprior.txt", append = TRUE)

      xrobs <- x$rowObs
      xcobs <- x$cellObs
      
      repobs <- lapply(levels(replicates), function(rep) {
        repdim <- datdim
        repdim[2] <- sum(replicates == rep)
        c(xrobs,
          lapply(xcobs, function(obs) .sliceArray(list(1, which(replicates == rep)), obs)),
          lapply(sampleObservables, function(obs) .sliceArray(list(which(replicates == rep)), obs)),
          list(dim = repdim))
      })
      repdat <- lapply(levels(replicates), function(rep) .sliceArray(list(NULL, which(replicates == rep)), dat, drop = FALSE))
                             
      optimFixed <- function(vars, repdat, datdim, replicates, repobs)
        {    
            pars <- split(vars, splitOptimFixed)
            pars[!eqovRep] <- lapply(pars[!eqovRep], function(par) par[match(replicates, levels(replicates))])
            pars[eqovRep] <- lapply(pars[eqovRep], function(par) rep(par, datdim[2]))
            sum(
                sapply(levels(replicates), function(rep) {            
                    xobs <- repobs[[which(levels(replicates) == rep)]]
                    xobs$priorEstimator <- TRUE
                    sum(optimFunction(repdat[[which(levels(replicates) == rep)]], xobs, lapply(pars, function(x) x[replicates == rep])))
                })
            )
        }
                       
      optimGroup <- function(vars, fixed, gdat, xobs, ...)
        {
          pars <- list()
          pars[eqovRep] <- fixed
          pars[!eqovRep] <- split(vars, 1:length(vars))
          pars <- lapply(pars, function(par) rep(par, ncol(gdat)))
          xobs$priorEstimator <- TRUE
          sum(optimFunction(gdat, xobs, pars))
        }  

      parOptimFixed <- do.call("rbind", lapply(levels(replicates), function(rep) {
        repdat <- .sliceArray(list(NULL, which(replicates == rep)), dat, drop = FALSE)
        repobs <- c(xrobs,
                    lapply(xcobs, function(obs) .sliceArray(list(1, which(replicates == rep)), obs)),
                    lapply(sampleObservables, function(obs) .sliceArray(list(which(replicates == rep)), obs)),
                    list(dim = c(datdim[1], dim(repdat)[-1])))
        repInit <- initiatingValues(repdat, repobs)
      }))

      parOptimFixed[1,eqovRep] <- colMeans(parOptimFixed[,eqovRep,drop=FALSE])
      parOptimFixed[-1,eqovRep] <- NA
      parOptimFixed <- parOptimFixed[!is.na(parOptimFixed)]

      splitOptimFixed <- do.call("c", lapply(1:length(eqovRep), function(ii) if(eqovRep[ii]) ii else rep(ii, length(levels(replicates)))))

      if(any(eqovRep))
          fixed <- split(optim(par = parOptimFixed, fn = optimFixed, repdat = repdat, datdim = datdim, replicates = replicates, control = list(fnscale = -1, maxit = 5000, reltol = 1e-50), repobs = repobs)$par, splitOptimFixed)[eqovRep]

      groupValues <- function(gg, group) {
        gdat <- .sliceArray(list(NULL, which(group == gg)), dat, drop = FALSE)
        parGroup <- c()
        gpobs <- c(xrobs,
                   lapply(xcobs, function(obs) .sliceArray(list(1, which(group == gg)), obs)),
                   lapply(sampleObservables, function(obs) .sliceArray(list(which(group == gg)), obs)),
                   list(dim = c(datdim[1], dim(gdat)[-1])))
        parInitGroup <- initiatingValues(gdat, gpobs)[!eqovRep]
        if(sum(!eqovRep) > 1 || is.null(lower) || is.null(upper))
          parGroup[!eqovRep] <- optim(parInitGroup, fn = optimGroup, fixed = fixed, gdat = gdat,
                                      xobs = gpobs, control = list(fnscale = -1, maxit = 5000, reltol = 1e-50), method = "Nelder-Mead")$par
        if(sum(!eqovRep) == 1) {
          if(is.function(lower)) lowbound <- lower(gdat) else lowbound = lower
          if(is.function(upper)) uppbound <- upper(gdat) else uppbound = upper
          parGroup[!eqovRep] <- optimise(optimGroup, fixed = fixed, gdat = gdat,
                                                     xobs = gpobs, maximum = TRUE, lower = lowbound, upper = uppbound, tol = 1e-50)[[1]]
        }
        parGroup[eqovRep] <- unlist(fixed)
        parGroup
      }

      if(consensus) {
        groupness <- groupValues(gg = sample(levels(replicates), size = 1), group = replicates)
      } else {
        groupness <- lapply(groups, function(group) lapply(levels(group), groupValues, group = group))
      }

      groupness
    }
  
  
  if(!inherits(cD, what = "countData"))
    stop("variable 'cD' must be of or descend from class 'countData'")
  
  if(verbose) message("Finding priors...", appendLF = FALSE)          
  
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


#  if(nrow(cD@seglens) == 0) cD@seglens <- matrix(1, ncol = ncol(cD@data), nrow = nrow(cD@data))

  
  if(is.null(samplingSubset))
    samplingSubset <- 1:nrow(cD)
  
  sD <- cD[samplingSubset,]
  groups <- sD@groups
  replicates <- as.factor(sD@replicates)
  if(!all(levels(replicates) %in% replicates)) {
    warning("Not all replicate levels have corresponding samples; dropping these now")
    replicates <- factor(as.character(replicates), levels = levels(replicates)[levels(replicates) %in% replicates])
  }
  data <- sD@data  
  sampleObservables <- sD@sampleObservables
  cellObservables <- sD@cellObservables
  rowObservables <- sD@rowObservables
  densityFunction <- cD@densityFunction

  eqovRep <- densityFunction@equalOverReplicates(dim(cD))    

  if(!("seglens" %in% names(cellObservables) || "seglens" %in% names(rowObservables))) rowObservables <- c(rowObservables, list(seglens = rep(1, nrow(cD))))
  if(!("libsizes" %in% names(sampleObservables))) sampleObservables <- c(sampleObservables, list(seglens = rep(1, ncol(cD))))

  rowObservables <- lapply(rowObservables, function(x) if(length(dim(x)) == 1) as.vector(x) else x)
  sampleObservables <- lapply(sampleObservables, function(x) if(length(dim(x)) == 1) as.vector(x) else x)
  
  if(nrow(data) <= samplesize) {
    sampDat <- 1:nrow(data)
    weights <- rep(1, nrow(data))
  } else {  
    if(length(densityFunction@stratifyBreaks) > 0) breaks <- densityFunction@stratifyBreaks else breaks <- 1  
    if(!is.null(formals(densityFunction@stratifyFunction)) & breaks > 1) {
        if(length(data) == 1) stratDat <- densityFunction@stratifyFunction(data) else stratDat <- densityFunction@stratifyFunction(data)
        properCut <- function(x, breaks) if(breaks == 1) return(rep(1, length(x))) else return(cut(x, breaks, labels = FALSE))     
        cutdat <- properCut(stratDat, breaks)
        spldat <- split(1:nrow(data), cutdat)
#        spldat <- spldat[order(sapply(spldat, length), decreasing = FALSE)]
    } else spldat <- list(1:nrow(data))

    spllen <- sapply(spldat, length)

    if(any(spllen < ceiling(samplesize / spllen))) {
        sampType <- min(which(((samplesize - cumsum(spllen)) / (length(spldat):1 - 1))[-length(spllen)] <= spllen[-1]))
        sampLess <- ceiling((samplesize - sum(spllen[1:(sampType)])) / (length(spllen) - sampType))         
        sampDat <- c(spldat[1:sampType],
                     lapply(spldat[(sampType + 1):length(spldat)], function(x) sample(x, size = sampLess, replace = FALSE)))
        weights <- rep(sapply(spldat, length) / sapply(sampDat, length), sapply(sampDat, length))
        sampDat <- unlist(sampDat)
    } else {
      sampDat <- lapply(spldat, function(x) sample(x, size = min(length(x), ceiling(samplesize / length(spldat))), replace = FALSE))        
      weights <- rep(sapply(spldat, length) / sapply(sampDat, length), sapply(sampDat, length))
      sampDat <- unlist(sampDat)
    }
  }

  sliceData <- lapply(sampDat, function(id)
                      list(id = id,
                           data = asub(data, id, dims = 1, drop = FALSE),
                           cellObs = lapply(cellObservables, function(cob) asub(cob, id, dims = 1, drop = FALSE)),
                           rowObs = lapply(rowObservables, function(rob) asub(rob, id, dims = 1, drop = FALSE))))


#  dataLL = length(dim(data)) - 1
#  if(dataLL == 1) data <- data
#  if(dataLL > 1) data <- matrix(apply(z, 3, c), nrow = nrow(z))

  save(sliceData, file = "tempSlice.RData")
  
  if(is.null(cl)) {
      parEach <- lapply(sliceData, .getSinglePriors, replicates = replicates, datdim = dim(cD), groups = groups, optimFunction = densityFunction@density, eqovRep = eqovRep, initiatingValues = densityFunction@initiatingValues, lower = densityFunction@lower, upper = densityFunction@upper, consensus = consensus)
  } else {
      clusterExport(cl, c("sampleObservables", ".sliceArray"), envir = environment())

      getPriorsEnv <- new.env(parent = .GlobalEnv)
      environment(.getSinglePriors) <- getPriorsEnv
      parEach <- parLapplyLB(cl, sliceData, .getSinglePriors, replicates = replicates, datdim = dim(cD), groups = groups, optimFunction = densityFunction@density, eqovRep = eqovRep, initiatingValues = densityFunction@initiatingValues, lower = densityFunction@lower, upper = densityFunction@upper, consensus = consensus)
  }

  if(consensus) {
    LNpar <- do.call("rbind", parEach)
    names(LNpar) <- NULL  
  } else {
    LNpar <- lapply(1:length(groups), function(gg)
                    lapply(1:length(levels(groups[[gg]])), function(ii)
                           do.call("rbind", lapply(parEach, function(x) x[[gg]][[ii]]))))
  }

  if(verbose) message("done.")

  if(!is.null(cl))
    clusterEvalQ(cl, rm(list = ls()))

  sy <- cbind(sampled = samplingSubset[sampDat], representative = 1:length(sampDat))
  if(!consensus) names(LNpar) <- names(groups)
  LNpar <- list(sampled = cbind(sy, weights = weights / sum(weights)), priors = LNpar)
  cD@priorType <- "user-supplied function"; cD@priors = LNpar
  cD
}
