`getPriors.Dirichlet` <-
function(cD, samplesize = 10^5, iterations = 10^3)
  {
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")

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
        priors[[gg]] <- matrix(NA, ncol = 2, nrow = length(unique(groups[[gg]])))
        for(uu in 1:length(unique(groups[[gg]])))
          {
            cat(".")
            temp.priors <- NULL
            for(ii in 1:iterations)
              {
                initial <- c(0.2, 1100)
                if(!is.null(temp.priors))
                  initial <- apply(temp.priors, 2, mean)
                
                temp.priors <- rbind(temp.priors, optim(initial, PgivenDir,
                                                        control = list(fnscale = -1),
                                                        us = rowSums(data.frame(y[sample(1:nrow(y), samplesize, replace = FALSE),groups[[gg]] == unique(groups[[gg]])[uu]])),
                                                        ns = sum(libsizes[groups[[gg]] == unique(groups[[gg]])[uu]]))$par)
              }
            priors[[gg]][uu,] <- apply(temp.priors, 2, mean)
          }
      }
    names(priors) <- names(groups)
    new(class(cD), cD, priors = list(type = "Dir", priors = priors))
  }


`getPriors.Pois` <-
  function (cD, samplesize = 10^5, perSE = 1e-1
            , takemean = TRUE, cl) 
{
  if(!inherits(cD, what = "countData"))
    stop("variable 'cD' must be of or descend from class 'countData'")
  
  priorPars <- function(y, libsizes, samplesize, seluu, initial)
    {
      PgivenPois <- function(priors, us, lens = lens, ns) {
        logProbDRepGivenR <- function(D, alpha, beta, nosamps)
          {
            uses <- D[1:nosamps]
            ns <- D[1:nosamps + nosamps]
            lambda <- D[1:nosamps + 2 * nosamps]
            sum(uses * log(lambda * ns) - lfactorial(uses)) +
              alpha * log(beta) - lgamma(alpha) +
                lgamma(alpha + sum(uses)) - (alpha + sum(uses)) * log(sum(lambda * ns) + beta)
          }

        if (any(priors <= 0)) 
          return(NA)
        sum(apply(cbind(us, matrix(ns, nrow = nrow(us), ncol = ncol(us), byrow = TRUE), lens),
                  1, logProbDRepGivenR, priors[1], priors[2], ncol(us)))
        
      }

      cts <- y[,-(1:(ncol(y) / 2))]
      len <- y[,1:(ncol(y) / 2)]

      sampcts <- sample(1:nrow(y), size = samplesize, replace = FALSE)

      
      
      us <- cts[sampcts, seluu, drop = FALSE]
      lens <- len[sampcts, seluu, drop = FALSE]
      optim(initial, 
            PgivenPois, control = list(fnscale = -1),
            us = us, lens = lens,
            ns = libsizes[seluu])$par
    }
  
  clustAssignData <- function(y)
    {
      assign("y", y, envir = .GlobalEnv)
      NULL
    }
  
  getPriorEnv <- new.env(parent = .GlobalEnv)
  environment(clustAssignData) <- getPriorEnv
  environment(priorPars) <- getPriorEnv
  
  if(class(cD) == "segData") seglens <- cD@seglens else seglens <- matrix(1, ncol = 1, nrow = nrow(cD@data))
  
  if(is.matrix(seglens))
    if(ncol(seglens) == 1) seglens <- matrix(seglens[,1], ncol = ncol(cD@data), nrow = nrow(cD@data))
  
  
  priors <- list()
  libsizes <- cD@libsizes
  y <- cbind(seglens, cD@data)

  if(!is.null(cl))
    clusterCall(cl, clustAssignData, y)
    
  groups <- cD@groups
  for (gg in 1:length(groups)) {
    if (takemean) {
      priors[[gg]] <- matrix(NA, ncol = 2, nrow = length(unique(groups[[gg]])))
    } else priors[[gg]] <- list()
    initial <- c(0.2, 1100)
    for (uu in 1:length(unique(groups[[gg]]))) {
      tempPriors <- NULL
      seluu <- which(groups[[gg]] == unique(groups[[gg]])[uu])
      repeat {
        if (!is.null(tempPriors)) 
          initial <- apply(tempPriors, 2, mean)
        if(!is.null(cl))
          {
            pars <- clusterCall(cl, priorPars, y, libsizes, samplesize, seluu, initial)
            pars <- matrix(unlist(pars), ncol = 2, byrow = TRUE)
          } else {
            pars <- matrix(priorPars(y, libsizes, samplesize, seluu, initial), ncol = 2, nrow = 1)
          }
        tempPriors <- rbind(tempPriors, pars[apply(!is.na(pars), 1, all),])
        if(nrow(tempPriors) > 1)
          if(all((apply(tempPriors, 2, sd) / sqrt(nrow(tempPriors))) / apply(tempPriors, 2, mean) < perSE)) break()
      }
      cat(".")
      if(takemean) 
        priors[[gg]][uu, ] <- apply(tempPriors, 2, mean)
      else priors[[gg]][[uu]] <- tempPriors
    }
  }
  names(priors) <- names(groups)
  new(class(cD), cD, priors = list(type = "Poi", priors = priors))
}


`getPriors.NB` <-
function (cD, samplesize = 10^5, estimation = "ML", cl)
{
  if(!inherits(cD, what = "countData"))
    stop("variable 'cD' must be of or descend from class 'countData'")


  ML.optimover <- function(x, libsizes)
    {
      ML.NB <- function(alphas, y, libsizes, seglen)
        if (all(alphas > 0)) 
          sum(dnbinom(y, size = 1/alphas[2], mu = alphas[1] * libsizes * seglen, log = TRUE)) else return(NA)
      
      cts <- x[-(1:(length(x) / 2))]
      len <- x[1:(length(x) / 2)]
      
      optim(par = c(mean(cts/(libsizes * len)) + abs(rnorm(1, 0, 0.001)) * as.numeric(all(cts == 0)), 0.01),
            fn = ML.NB, control = list(fnscale = -1), 
            y = cts, libsizes = libsizes, seglen = len)$par
    }
  
  QL.optimover <- function(x, libsizes)
    {
      dispalt <- function(dispersion, cts, mu, libsizes, len)
        abs(2 * (sum(cts[cts > 0] * log(cts[cts > 0]/(mu * libsizes[cts > 0] * len))) -
                 sum((cts + dispersion^-1) * log((cts + dispersion^-1) / (mu * libsizes * len + dispersion^-1)))) - (length(cts)-1))
      mualt <- function(mu, cts, dispersion, libsizes, len)
        sum(dnbinom(cts, size = 1/ dispersion, mu = mu * libsizes * len, log = TRUE))

      cts <- x[-(1:(length(x) / 2))]
      len <- x[1:(length(x) / 2)]

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

  if(class(cD) == "segData") seglens <- cD@seglens else seglens <- matrix(1, ncol = 1, nrow = nrow(cD@data))
  
  if(is.matrix(seglens))
    if(ncol(seglens) == 1) seglens <- matrix(seglens[,1], ncol = ncol(cD@data), nrow = nrow(cD@data))

  
  NBpar <- list()
  for (ii in 1:length(groups)) {
    NBpar[[ii]] <- list()
    for (jj in unique(groups[[ii]])) {
      NBpar.group <- c()
      z <- cbind(seglens[sy,groups[[ii]] == jj], y[sy, groups[[ii]] == jj, drop = FALSE])
      if (estimation == "ML") {
        if (is.null(cl)) 
          parEach <- t(apply(z, 1, ML.optimover, libsizes[groups[[ii]] == jj]))
        else parEach <- t(parApply(cl, z, 1, ML.optimover, 
                                   libsizes[groups[[ii]] == jj]))
        cat(".")
      } else if (estimation == "QL") {
        if(is.null(cl)) parEach <- t(apply(z, 1, QL.optimover, libsizes[groups[[ii]] == jj])) else parEach <- t(parApply(cl, z, 1, QL.optimover, libsizes[groups[[ii]] == jj]))
      }
      rownames(parEach) <- sy
      NBpar[[ii]][[jj]] <- parEach
    }
  }
  names(NBpar) <- names(groups)
  NBpar <- list(type = "NB", sampled = sy, priors = NBpar)
  new(class(cD), cD, priors = NBpar)
}
