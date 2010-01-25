`getPriors.Dirichlet` <-
function(cD, samplesize = 10^5, perSE = 1e-1, maxit = 10^6)
  {
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")

    cat("Finding priors...")
    
    `PgivenDir` <-
      function(alphas, us, ns)
        {
          if(any(alphas <= 0)) return(NA)
          sum(lgamma(alphas[1] + us) + lgamma(alphas[2] + (ns - us)) - lgamma(sum(alphas) + ns) -
              (sum(lgamma(alphas)) - lgamma(sum(alphas))))
        }
    
    priors <- list()
    libsizes <- cD@libsizes
    y <- cD@data
    groups <- cD@groups
    for(gg in 1:length(groups))
      {
        priors[[gg]] <- list()
        for(uu in 1:length(unique(groups[[gg]])))
          {
            cat(".")
            tempPriors <- NULL
            for(rr in 1:maxit)
              {
                initial <- c(0.2, 1100)
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
            priors[[gg]][[uu]] <- apply(tempPriors, 2, mean)
          }
      }
    names(priors) <- names(groups)
    new(class(cD), cD, priors = new("priorData", type = "Dir", priors = priors))
  }


`getPriors.Pois` <-
  function (cD, samplesize = 10^5, perSE = 1e-1
            , takemean = TRUE, maxit = 10^5, cl) 
{
  if(!inherits(cD, what = "countData"))
    stop("variable 'cD' must be of or descend from class 'countData'")

  cat("Finding priors...")
  
  priorPars <- function(libsizes, samplesize, seluu, initial, lensameFlag)
    {
      PgivenPois <- function(priors, us, lens, ns, lensameFlag) {
        logProbDRepGivenR <- function(D, alpha, beta, nosamps, lensameFlag)
          {
            uses <- D[1:nosamps]
            ns <- D[1:nosamps + nosamps]
            if(lensameFlag) lambda <- D[2 * nosamps + 1] else lambda <- D[1:nosamps + 2 * nosamps]
            sum(uses * log(lambda * ns) - lfactorial(uses)) +
              alpha * log(beta) - lgamma(alpha) +
                lgamma(alpha + sum(uses)) - (alpha + sum(uses)) * log(sum(lambda * ns) + beta)
          }

        if (any(priors <= 0)) 
          return(NA)
        sum(apply(cbind(us, matrix(ns, nrow = nrow(us), ncol = ncol(us), byrow = TRUE), lens),
                  1, logProbDRepGivenR, alpha = priors[1], beta = priors[2], nosamps = ncol(us), lensameFlag = lensameFlag))
        
      }

      sampcts <- sample(1:nrow(y), size = samplesize, replace = FALSE)
      
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
      
      ns <- libsizes[seluu]
      
      if(is.null(initial))
        {
          m <- mean(apply(t(t(us / lens) / ns), 1, mean))
          v <- var(apply(t(t(us / lens) / ns), 1, mean))
          initial <- c(m^2 / v, m / v)
        }
      
      optim(initial, 
            PgivenPois, control = list(fnscale = -1),
            us = us, lens = lens,
            ns = ns, lensameFlag = lensameFlag)$par
    }

  if(nrow(cD@seglens) > 0) seglens <- cD@seglens else seglens <- matrix(1, ncol = 1, nrow = nrow(cD@data))
  
  if(ncol(seglens) == 1) lensameFlag <- TRUE else lensameFlag <- FALSE
  
  
  priors <- list()
  libsizes <- cD@libsizes
  y <- cbind(seglens, cD@data)

  if(!is.null(cl))
    {
      clustAssignData <- function(y)
        {
          assign("y", y, envir = .GlobalEnv)
          NULL
        }
      
      getPriorEnv <- new.env(parent = .GlobalEnv)
      environment(clustAssignData) <- getPriorEnv
      environment(priorPars) <- getPriorEnv

      clusterCall(cl, clustAssignData, y)
    }
    
    
  groups <- cD@groups
  for (gg in 1:length(groups)) {
    priors[[gg]] <- list()
    initial <- NULL
    for (uu in 1:length(unique(groups[[gg]]))) {
      tempPriors <- NULL
      seluu <- which(groups[[gg]] == unique(groups[[gg]])[uu])
      for(rr in 1:maxit) {
        if (!is.null(tempPriors)) 
          initial <- apply(tempPriors, 2, mean)
        if(!is.null(cl))
          {
            pars <- clusterCall(cl, priorPars, libsizes, samplesize, seluu, initial, lensameFlag = lensameFlag)
            pars <- matrix(unlist(pars), ncol = 2, byrow = TRUE)
          } else {
            pars <- matrix(priorPars(libsizes, samplesize, seluu, initial, lensameFlag = lensameFlag), ncol = 2, nrow = 1)
          }
        tempPriors <- rbind(tempPriors, pars[apply(!is.na(pars), 1, all),])
        if(nrow(tempPriors) > 1)
          if(all((apply(tempPriors, 2, sd) / sqrt(nrow(tempPriors))) / apply(tempPriors, 2, mean) < perSE)) break()
      }
      if(rr == maxit)
        warning(paste("Convergence not achieved to required accuracy for model ", gg, ", group ", uu, sep = ""))
      cat(".")
      if(takemean) 
        priors[[gg]][[uu]] <- apply(tempPriors, 2, mean)
      else priors[[gg]][[uu]] <- tempPriors
    }
  }
  names(priors) <- names(groups)
  new(class(cD), cD, priors = new("priorData", type = "Poi", priors = priors))
}


`getPriors.NB` <-
function (cD, samplesize = 10^5, estimation = "ML", cl)
{
  if(!inherits(cD, what = "countData"))
    stop("variable 'cD' must be of or descend from class 'countData'")

  cat("Finding priors...")
  
  ML.optimover <- function(x, libsizes, lensameFlag)
    {
      ML.NB <- function(alphas, y, libsizes, seglen)
        if (all(alphas > 0)) 
          sum(dnbinom(y, size = 1/alphas[2], mu = alphas[1] * libsizes * seglen, log = TRUE)) else return(NA)

      if(lensameFlag)
        {
          cts <- x[-1]
          len <- x[1]
        } else {
          cts <- x[-(1:(length(x) / 2))]
          len <- x[1:(length(x) / 2)]
        }
      
      optim(par = c(mean(cts/(libsizes * len)) + abs(rnorm(1, 0, 0.001)) * as.numeric(all(cts == 0)), 0.01),
            fn = ML.NB, control = list(fnscale = -1), 
            y = cts, libsizes = libsizes, seglen = len)$par
    }
  
  QL.optimover <- function(x, libsizes, lensameFlag)
    {
      dispalt <- function(dispersion, cts, mu, libsizes, len)
        abs(2 * (sum(cts[cts > 0] * log(cts[cts > 0]/(mu * libsizes[cts > 0] * len))) -
                 sum((cts + dispersion^-1) * log((cts + dispersion^-1) / (mu * libsizes * len + dispersion^-1)))) - (length(cts)-1))
      mualt <- function(mu, cts, dispersion, libsizes, len)
        sum(dnbinom(cts, size = 1/ dispersion, mu = mu * libsizes * len, log = TRUE))

      if(lensameFlag)
        {
          cts <- x[-1]
          len <- x[1]
        } else {
          cts <- x[-(1:(length(x) / 2))]
          len <- x[1:(length(x) / 2)]
        }

      mu <- mean(cts / (libsizes * len))
      disp <- 0
      newmu <- newdisp <- -1
      repeat {
        newdisp <- c(newdisp, optimise(dispalt, interval = c(0, 100), cts = cts, mu = mu, libsizes = libsizes, len = len, tol = 1e-30)$minimum)
        disp <- newdisp[length(newdisp)]
        newmu <- c(newmu, optimise(mualt, interval = c(0, 1000), cts = cts, dispersion = disp, libsizes = libsizes, len = len, maximum = TRUE)$maximum)
        if(any(abs(newmu[-length(newmu)] - newmu[length(newmu)]) < 1e-10 &
               abs(newdisp[-length(newdisp)] - newdisp[length(newdisp)]) < 1e-10))break
        mu <- newmu[length(newmu)]
      }
      c(newmu[length(newmu)], newdisp[length(newdisp)])
    }
  
  libsizes <- cD@libsizes
  y <- cD@data
  groups <- cD@groups

  if(nrow(y) < samplesize) samplesize <- nrow(y)
  sy <- sample(1:nrow(y), samplesize, replace = FALSE)
  sy <- sort(sy)

  if(nrow(cD@seglens) > 0) seglens <- cD@seglens else seglens <- matrix(1, ncol = 1, nrow = nrow(cD@data))
  
  if(ncol(seglens) == 1) lensameFlag <- TRUE else lensameFlag <- FALSE

  
  
  NBpar <- list()
  for (ii in 1:length(groups)) {
    NBpar[[ii]] <- list()
    for (jj in unique(groups[[ii]])) {
      NBpar.group <- c()
      if(lensameFlag)
        {
          z <- cbind(seglens[sy,], y[sy, groups[[ii]] == jj, drop = FALSE])
        } else {
          z <- cbind(seglens[sy,groups[[ii]] == jj], y[sy, groups[[ii]] == jj, drop = FALSE])
        }
      if (estimation == "ML") {
        if (is.null(cl)) 
          parEach <- t(apply(z, 1, ML.optimover, libsizes[groups[[ii]] == jj], lensameFlag = lensameFlag))
        else parEach <- t(parApply(cl, z, 1, ML.optimover, 
                                   libsizes[groups[[ii]] == jj], lensameFlag = lensameFlag))
        cat(".")
      } else if (estimation == "QL") {
        if(is.null(cl)) parEach <- t(apply(z, 1, QL.optimover, libsizes[groups[[ii]] == jj], lensameFlag = lensameFlag)) else parEach <- t(parApply(cl, z, 1, QL.optimover, libsizes[groups[[ii]] == jj], lensameFlag = lensameFlag))
      }
      rownames(parEach) <- sy
      NBpar[[ii]][[jj]] <- parEach
    }
  }
  names(NBpar) <- names(groups)
  NBpar <- new("priorData", type = "NB", sampled = sy, priors = NBpar)
  new(class(cD), cD, priors = NBpar)
}
