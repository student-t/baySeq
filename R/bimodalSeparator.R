
.bimodalKernel <- function(z, weights = NULL)
  {
                                        # forces bimodality by adjusting bandwidth until it is achieved    
    adjdens <- c(1)
    modes <- Inf   
    
    while(TRUE) {
      curadj <- adjdens[length(adjdens)]
      denz <- density(z, weights = weights, adjust = curadj)
      modes <- sum(denz$y > c(denz$y[-1], 0) & denz$y > c(0, denz$y[-length(denz$y)]))
      if(modes > 2)
        if(curadj == max(adjdens)) adjdens <- c(adjdens, max(adjdens) * 2) else adjdens <- c(adjdens, (curadj + min(adjdens[adjdens > curadj])) / 2)
      if(modes <= 2)
        if(curadj == min(adjdens)) adjdens <- c(adjdens, min(adjdens) / 2) else adjdens <- c(adjdens, (curadj + max(adjdens[adjdens < curadj])) / 2)
      if(modes == 2) {
        if(length(adjdens) > 1000) warning("Maximum number of iterations exceeded")
        if(abs(curadj - adjdens[length(adjdens)]) / curadj < 1e-5) break
      }
    }
    return(denz)
  }


bimodalSeparator <- function(x, weights = NULL, minperc = 0.1) {
  
  if(is.null(weights)) weights <- rep(1, length(x))

  xnif <- (x < Inf & x > -Inf)
  z <- x[xnif]
  weights <- weights[xnif]

#  weights <- round(weights/min(weights))
#  weights <- round(weights / max(weights) * 1000)
#  z <- sort(rep(x, weights))

  rordz <- order(z)
  z <- z[rordz]
  weights <- weights[rordz]

  weighted.var <- function(x, w)
    sum(w * (x - weighted.mean(x, w))^2) / (sum(w) - sum(w^2) / sum(w))
  
  zsel <- rep(TRUE, length(z))


  while(TRUE) {
    varest <- sapply(1:(length(z[zsel]) - 1), function(index)
                     index * weighted.var(z[zsel][1:index], weights[zsel][1:index]) + (length(z[zsel]) - index) * weighted.var(z[zsel][(index + 1):length(z[zsel])], weights[zsel][(index + 1):length(z[zsel])]))

    threshold <- mean(z[zsel][which.min(varest) + 0:1])
    if(sum(z > max(threshold)) < (minperc / 100) * length(z)) zsel[union(max(which(zsel)), which(z > max(threshold)))] <- FALSE
    if(sum(z < min(threshold)) < (minperc / 100) * length(z)) zsel[union(min(which(zsel)), which(z < min(threshold)))] <- FALSE
    if(sum(z > max(threshold)) > (minperc / 100) * length(z) & sum(z < max(threshold)) > (minperc / 100) * length(z)) break
  }
  return(threshold)
                                        #return(list(denz = denz, threshold = threshold))
}
  
