
.betaBinomialFunction <- function(dat, observables, parameters) {  
dbetabinom <- function(x, n, prop, disp, log = TRUE) {  
  ps <- rep(NA, (length(x)))
  #disp <- rep(disp, length(prop))    
  
  smallDisp <- disp < 1e-15 & disp >= 0
  largeDisp <- disp > 1e-15
  
  if(any(disp < 0 | disp >= 1 | prop > 1 | prop < 0)) ps[disp < 0 | disp >= 1 | prop > 1 | prop < 0] <- NA
  if(any(smallDisp))
    ps[smallDisp] <- dbinom(x[smallDisp], n[smallDisp], prob = prop[smallDisp], log = log)
  
  if(any(largeDisp)) {
    alpha = (1/disp - 1) * prop
    beta = (1/disp - 1) * (1-prop)
    
    if(!log) {
      ps[largeDisp] <- choose(n[largeDisp], x[largeDisp]) * beta((x + alpha)[largeDisp], (n - x + beta)[largeDisp]) / beta(alpha[largeDisp], beta[largeDisp])
    } else ps[largeDisp] <- lchoose(n[largeDisp], x[largeDisp]) + lbeta((x + alpha)[largeDisp], (n - x + beta)[largeDisp]) - lbeta(alpha[largeDisp], beta[largeDisp])
  }  
  return(ps)
}
  
  if(any(sapply(parameters, function(par) any(par < 0) | any(par > 1)))) return(NA)
  props <- parameters[[1]]
  props <- props * observables$libsizes[,1] / (props * observables$libsizes[,1] + (1-props) * observables$libsizes[,2])
  (dbetabinom(dat[,,1], dat[,,1] + dat[,,2], prop = props, disp = parameters[[2]]))
}

.multiDirichletFunction <- function(dat, observables, parameters) {
  if(any(sapply(parameters, function(par) any(par < 0) | any(par > 1)))) return(NA)

  ddirmult <- function(x, prop, disp) {  
    
    if(nrow(x) != nrow(prop)) stop("Dimension of x does not match dimension of proportions")  
    
    ps <- rep(NA, (nrow(x)))
                                        #  disp <- rep(disp, nrow(prop))    
    
    alphas = (1/disp - 1) * prop
    
    ps <- lfactorial(rowSums(x)) - rowSums(lfactorial(x)) + lgamma(rowSums(alphas)) - lgamma(rowSums(alphas) + rowSums(x)) + rowSums(lgamma(x + alphas) - lgamma(alphas))
    return(ps)
  }

  disp <- parameters[[1]]
  props <- do.call("cbind", parameters[-1])
  props <- cbind(do.call("cbind", parameters[-1]), 1 - rowSums(props))
#  props <- cbind(props, 1 - props)

  if(any(props >= 1) | any(props <= 0) | any(disp <= 2e-15)) return(NA)
  
  #props <- props * observables$libsizes[,1] / (props * observables$libsizes[,1] + (1-props) * observables$libsizes[,2])

  propmod <- props * observables$libsizes
  propmod <- propmod / rowSums(propmod)

  x <- dat[,,,drop = TRUE]
  disp <- array(disp, dim = dim(x))
  
  ddirmult(x = x, prop = propmod, disp = disp)
}



.multiDirichletFunction2 <- function(dat, observables, parameters) {
  if(any(sapply(parameters, function(par) any(par < 0) | any(par > 1)))) return(NA)

  ddirmult <- function(x, prop, disp) {  
    
    if(nrow(x) != nrow(prop)) stop("Dimension of x does not match dimension of proportions")  
    
    ps <- rep(NA, (nrow(x)))
                                        #  disp <- rep(disp, nrow(prop))    
    
    alphas = (1/disp - 1) * prop
    
    ps <- lfactorial(rowSums(x)) - rowSums(lfactorial(x)) + lgamma(rowSums(alphas)) - lgamma(rowSums(alphas) + rowSums(x)) + rowSums(lgamma(x + alphas) - lgamma(alphas))
    return(ps)
  }

  x <- dat[,,,drop = TRUE]
  maxp <- order(colSums(x / observables$libsizes), decreasing = TRUE)[1:2]
  
  disp <- parameters[[1]]
  props <- matrix(NA, nrow = nrow(x), ncol = observables$dim[3])
  props[,maxp] <- do.call("cbind", parameters[-1])

  if(!all(parameters[[2]] >= parameters[[3]])) return(NA)
  
  props[is.na(props)] <- (1 - rowSums(props, na.rm = TRUE)) / (observables$dim[3] - 2)
    
  if(any(props >= 1) | any(props[,maxp[1]] < 1/observables$dim[3]) | any(props <= 0) | any(disp <= 2e-15)) return(NA)
  
  #props <- props * observables$libsizes[,1] / (props * observables$libsizes[,1] + (1-props) * observables$libsizes[,2])

  propmod <- props * observables$libsizes
  propmod <- propmod / rowSums(propmod)

  disp <- array(disp, dim = dim(x))  
  ddirmult(x = x, prop = propmod, disp = disp)
}

.multiDirichletFunction3 <- function(dat, observables, parameters) {
  if(any(sapply(parameters, function(par) any(par < 0) | any(par > 1)))) return(NA)

  ddirmult <- function(x, prop, disp) {  
    
    if(nrow(x) != nrow(prop)) stop("Dimension of x does not match dimension of proportions")  
    
    ps <- rep(NA, (nrow(x)))
                                        #  disp <- rep(disp, nrow(prop))    
    
    alphas = (1/disp - 1) * prop
    
    ps <- lfactorial(rowSums(x)) - rowSums(lfactorial(x)) + lgamma(rowSums(alphas)) - lgamma(rowSums(alphas) + rowSums(x)) + rowSums(lgamma(x + alphas) - lgamma(alphas))
    return(ps)
  }

  x <- dat[,,,drop = TRUE]
  maxp <- order(colSums(x / observables$libsizes), decreasing = TRUE)[1:3]
  
  disp <- parameters[[1]]
  props <- matrix(NA, nrow = nrow(x), ncol = observables$dim[3])
  props[,maxp] <- do.call("cbind", parameters[-1])

  if(!all(parameters[[2]] >= parameters[[3]])) return(NA)
  
  props[is.na(props)] <- (1 - rowSums(props, na.rm = TRUE)) / (observables$dim[3] - 3)
    
  if(any(props >= 1) | any(props[,maxp[1]] < 1/observables$dim[3]) | any(props <= 0) | any(disp <= 2e-15)) return(NA)
  
  #props <- props * observables$libsizes[,1] / (props * observables$libsizes[,1] + (1-props) * observables$libsizes[,2])

  propmod <- props * observables$libsizes
  propmod <- propmod / rowSums(propmod)

  disp <- array(disp, dim = dim(x))  
  ddirmult(x = x, prop = propmod, disp = disp)
}



.multiSymDirichletFunction <- function(dat, observables, parameters) {
  if(any(sapply(parameters, function(par) any(par < 0) | any(par > 1)))) return(NA)

  ddirmult <- function(x, conc) {  
    
    #if(nrow(x) != nrow(prop)) stop("Dimension of x does not match dimension of proportions")  
    
    ps <- rep(NA, (nrow(x)))
                                        #  disp <- rep(disp, nrow(prop))    
    
    alphas = 1/conc
    
    ps <- lfactorial(rowSums(x)) - rowSums(lfactorial(x)) + lgamma(rowSums(alphas)) - lgamma(rowSums(alphas) + rowSums(x)) + rowSums(lgamma(x + alphas) - lgamma(alphas))
    return(ps)
  }

  conc <- parameters[[1]]
#  props <- do.call("cbind", parameters[-1])
#  props <- cbind(do.call("cbind", parameters[-1]), 1 - rowSums(props))
#  props <- cbind(props, 1 - props)

  if(any(conc <= 2e-15)) return(NA)
  
  #props <- props * observables$libsizes[,1] / (props * observables$libsizes[,1] + (1-props) * observables$libsizes[,2])

#  propmod <- props * observables$libsizes
#  propmod <- propmod / rowSums(propmod)

  x <- dat[,,,drop = TRUE]
  conc <- array(conc, dim = dim(x))
  
  ddirmult(x = x, conc)
}



.nbinomDens <- function(dat, observables, parameters) {
  if(any(sapply(parameters, function(par) any(par < 0)))) return(NA)
#  if(any(zinf < 1-0.5^(1/observables$dim[2]))) return(NA)
#  if("priorEstimator" %in% names(observables) & sum(pnbinom(rep(0, length(dat)), mu = parameters[[1]] * observables$libsizes * observables$seglens, size = 1 / parameters[[2]], log.p = TRUE), na.rm = TRUE) < log(0.5)) return(NA)
#  else {    
    dnbinom(dat, mu = parameters[[1]] * observables$libsizes * observables$seglens, size = 1 / parameters[[2]], log = TRUE)
#  }
}


.normDensityFunction <- function(dat, observables, parameters) {
  if(any(parameters[[2]] < 0)) return(NA)
  dnorm(dat, parameters[[1]] * observables$libsizes, parameters[[2]], log = TRUE)
}

.dZINB <- function(dat, observables, parameters) {
  if(any(sapply(parameters, function(par) any(par < 0)))) return(NA)  
  zinf <- parameters[[3]]
  if(any(zinf >= 1)) return(NA)
  if(any(zinf < 1-0.5^(1/observables$dim[2]))) return(NA)
  
  nzero <- dnbinom(dat, mu = parameters[[1]] * observables$libsizes * observables$seglens, size = 1 / parameters[[2]], log = TRUE) + log(1 - zinf)
  zero <- c(-Inf, log(zinf))[as.integer(dat == 0) + 1]
  pmz <- pmax(nzero, zero)
  log(exp(zero - pmz) + exp(nzero - pmz)) + pmz
}

.betaBinomialNCFunction <- function(dat, observables, parameters) {

  .logsum <- function(x)
    max(x, max(x, na.rm = TRUE) + log(sum(exp(x - max(x, na.rm = TRUE)), na.rm = TRUE)), na.rm = TRUE)
  .fastUniques <- function(x){
    if (nrow(x) > 1) {
      return(c(TRUE, rowSums(x[-1L, , drop = FALSE] == x[-nrow(x),, drop = FALSE]) != ncol(x)))
    } else return(TRUE)
  }    


  dbetabinom <- function(x, n, prop, disp, lchoo, alpha, beta, lbaba, log = TRUE) {

    ps <- c()
                                        #disp <- rep(disp, length(prop))    
    
    smallDisp <- which(disp < 1e-15 & disp >= 0)
    largeDisp <- which(disp > 1e-15)
    
    #if(any(disp < 0 | disp >= 1 | prop > 1 | prop < 0)) ps[disp < 0 | disp >= 1 | prop > 1 | prop < 0] <- NA
    
    if(length(smallDisp) > 0)
      ps[smallDisp] <- dbinom(x[smallDisp], n[smallDisp], prob = prop[smallDisp], log = log)
    
    if(length(largeDisp) > 0) {
        if(missing(alpha))
            alpha <- (1/disp - 1) * prop
        if(missing(beta))
            beta <- (1/disp - 1) * (1-prop)
        if(missing(lbaba))
            lbaba <- lbeta(alpha[largeDisp], beta[largeDisp])
      
      if(!log) {
        ps[largeDisp] <- choose(n[largeDisp], x[largeDisp]) * beta((x + alpha)[largeDisp], (n - x + beta)[largeDisp]) / beta(alpha[largeDisp], beta[largeDisp])
    } else ps[largeDisp] <- lchoo[largeDisp] + lbeta((x + alpha)[largeDisp], (n - x + beta)[largeDisp]) - lbaba
    }  
    return(ps)
  }
    
#  if(any(sapply(parameters[1:2], function(par) any(par <= 0) | any(par >= 1)))) return(NA)

  if(any(parameters[[1]] <= 0) || any(parameters[[1]] >= 1)) return(NA)
  if(any(parameters[[2]] < 1e-15) || any(parameters[[2]] >= 1)) return(NA)
  
  repness <- observables$ncrange
  largeness <- 1e7
  if(round(sum(repness) / largeness) > 1) {
      repsplit <- split(1:nrow(repness), cut(cumsum(rowSums(repness)), breaks = round(sum(rowSums(repness)) / largeness), labels = FALSE))
  } else repsplit <- list(1:nrow(repness))
                                        #  names(parameters)[1:2] <- c("prop", "disp")
  ll <- do.call("c", lapply(repsplit, function(repid) {      
      datid <- dat[repid,,,drop = FALSE]
      obsid <- list()
      obsid$ncseq <- observables$ncseq[repid,,drop = FALSE]
      obsid$lchoose <- observables$lchoose[repid,,drop = FALSE]
      obsid$ncll <- observables$ncll[repid,,drop = FALSE]
      repnessid <- repness[repid,,drop = FALSE]
      tmat <- matrix(FALSE, nrow = nrow(dat), ncol = ncol(dat)); tmat[repid,] <- TRUE
      parid <- lapply(parameters, function(par) par[which(as.vector(tmat))])

      
      if(length(parid) == 2) {
          likes <- dbetabinom(
              x = rep(datid[,,1], repnessid) - unlist(obsid$ncseq),
              n = rep(datid[,,1] + datid[,,2], repnessid),
              lchoo = unlist(obsid$lchoose),
              prop = rep(parid[[1]], repnessid),
              disp = rep(parid[[2]], repnessid)) + unlist(obsid$ncll)
      } else {
          likes <- dbetabinom(
              x = rep(datid[,,1], repnessid) - unlist(obsid$ncseq),
              n = rep(datid[,,1] + datid[,,2], repnessid),
              lchoo = unlist(obsid$lchoose),
              prop = rep(parid[[1]], repnessid),
              disp = rep(parid[[2]], repnessid),
              alpha = rep(parid[[3]], repnessid),
              beta = rep(parid[[4]], repnessid),
              lbaba = rep(parid[[5]], repnessid))+ unlist(obsid$ncll)
      }

      ll <- sapply(
          split(likes, rep(1:length(repnessid), repnessid)),                
          .logsum)
  }))
                           

ll
}
