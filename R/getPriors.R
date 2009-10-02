`getPriors.Dirichlet` <-
function(cD, samplesize = 10^5, iterations = 10^3)
  {
    if(class(cD) != "countData")
      stop("variable 'cD' must be of class 'countData'")

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
    new("countData", cD, priors = list(type = "Dir", priors = priors))
  }


`getPriors.Pois` <-
  function (cD, samplesize = 10^5, iterations = 10^3, takemean = TRUE) 
{
    if (class(cD) != "countData") 
        stop("variable 'cD' must be of class 'countData'")
    PgivenPois <- function(alphas, us, ns) {
        if (any(alphas <= 0)) 
            return(NA)
        alpha <- alphas[1]
        beta <- alphas[2]
        Uq <- rowSums(us)
        Nq <- sum(ns)
        sum(rowSums(t(t(us) * log(ns)) - lfactorial(us)) + alpha * 
            log(beta) + lgamma(Uq + alpha) - ((Uq + alpha) * 
            log(Nq + beta) + lgamma(alpha)))
    }
    priors <- list()
    libsizes <- cD@libsizes
    y <- cD@data
    groups <- cD@groups
    for (gg in 1:length(groups)) {
        if (takemean) 
            priors[[gg]] <- matrix(NA, ncol = 2, nrow = length(unique(groups[[gg]])))
        else priors[[gg]] <- list()
        for (uu in 1:length(unique(groups[[gg]]))) {
            cat(".")
            temp.priors <- NULL
            for (ii in 1:iterations) {
                seluu <- which(groups[[gg]] == unique(groups[[gg]])[uu])
                initial <- c(0.2, 1100)
                if (!is.null(temp.priors)) 
                  initial <- apply(temp.priors, 2, mean)
                temp.priors <- rbind(temp.priors, optim(initial, 
                  PgivenPois, control = list(fnscale = -1), us = matrix(y[sample(1:nrow(y), 
                    samplesize, replace = FALSE), seluu], ncol = length(seluu)), 
                  ns = libsizes[seluu])$par)
            }
            if (takemean) 
                priors[[gg]][uu, ] <- apply(temp.priors, 2, mean)
            else priors[[gg]][[uu]] <- temp.priors
        }
    }
    names(priors) <- names(groups)
    new("countData", cD, priors = list(type = "Poi", priors = priors))
}


`getPriors.NB` <-
function (cD, samplesize = 10^5, estimation = "ML", cl)
{
  if(class(cD) != "countData")
    stop("variable 'cD' must be of class 'countData'")

  ML.optimover <- function(x, libsizes)
    {
      ML.NB <- function(alphas, y, libsizes)
        if (all(alphas > 0)) sum(dnbinom(y, size = 1/alphas[2], mu = alphas[1] * libsizes, log = TRUE)) else return(NA)
      optim(par = c(mean(x/libsizes) + abs(rnorm(1, 0, 0.001)) * as.numeric(all(x == 0)), 0.01), fn = ML.NB, control = list(fnscale = -1), y = x, libsizes = libsizes)$par
    }
  
  QL.optimover <- function(x, libsizes)
    {
      dispalt <- function(dispersion, x, mu, libsizes)
        abs(2 * (sum(x[x > 0] * log(x[x > 0]/(mu * libsizes[x > 0]))) -
                 sum((x + dispersion^-1) * log((x + dispersion^-1) / (mu * libsizes + dispersion^-1)))) - (length(x)-1))
      mualt <- function(mu, x, dispersion, libsizes)
        sum(dnbinom(x, size = 1/ dispersion, mu = mu * libsizes, log = TRUE))
      
      mu <- mean(x / libsizes)
      disp <- 0
      newmu <- newdisp <- -1
      repeat {
        newdisp <- c(newdisp, optimise(dispalt, interval = c(0, 100), x = x, mu = mu, libsizes = libsizes, tol = 1e-30)$minimum)
        disp <- newdisp[length(newdisp)]
        newmu <- c(newmu, optimise(mualt, interval = c(0, 1000), x = x, dispersion = disp, libsizes = libsizes, maximum = TRUE)$maximum)
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
  
  NBpar <- list()
  for (ii in 1:length(groups)) {
    NBpar[[ii]] <- list()
    for (jj in unique(groups[[ii]])) {
      NBpar.group <- c()
      z <- y[sy, groups[[ii]] == jj, drop = FALSE]
      if(estimation == "ML") {
        if(is.null(cl)) parEach <- t(apply(z, 1, ML.optimover, libsizes[groups[[ii]] == jj])) else parEach <- t(parApply(cl, z, 1, ML.optimover, libsizes[groups[[ii]] == jj]))
      } else if (estimation == "QL") {
        if(is.null(cl)) parEach <- t(apply(z, 1, QL.optimover, libsizes[groups[[ii]] == jj])) else parEach <- t(parApply(cl, z, 1, QL.optimover, libsizes[groups[[ii]] == jj]))
      }
      rownames(parEach) <- sy
      NBpar[[ii]][[jj]] <- parEach
    }
  }
  names(NBpar) <- names(groups)
  NBpar <- list(type = "NB", sampled = sy, priors = NBpar)
  new("countData", cD, priors = NBpar)
}
