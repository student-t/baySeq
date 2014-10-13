
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


bimodalSeparator <- function(x, weights = NULL, minperc = 10) {
  
  if(is.null(weights)) weights <- rep(1, length(x))

  xnif <- (x < Inf & x > -Inf)
  x <- x[xnif]
  weights <- weights[xnif]

  weights <- round(weights/min(weights))  
  z <- sort(rep(x, weights))
  
  
  zsel <- rep(TRUE, length(z))
#  if(require("diptest")) dip <- dip.test(z)$p
#  if(dip > 0.05) warning(paste("Bimodality not (p = ", round(dip, digits = 3), ") present; null data estimation may be inaccurate.", sep = ""))
  
#  while(TRUE) {    
#    denz <- .bimodalKernel(z[zsel])
    
    # discard outlier cases
    
#    thresholds <- denz$x[which(denz$y > c(denz$y[-1], 0) & denz$y > c(0, denz$y[-length(denz$y)]))]
    
#    if(sum(z > max(thresholds)) < 0.05 * length(z)) zsel[union(max(which(zsel)), which(z > max(thresholds)))] <- FALSE
#    if(sum(z < min(thresholds)) < 0.05 * length(z)) zsel[union(min(which(zsel)), which(z < min(thresholds)))] <- FALSE
#    if(sum(z > max(thresholds)) > 0.05 * length(z) & sum(z < max(thresholds)) > 0.05 * length(z)) break
#  }
  
  # selects within that region using Otsu's graylevel conversion method


  while(TRUE) {
    varest <- sapply(1:(length(z[zsel]) - 1), function(index)
                     index * var(z[zsel][1:index]) + (length(z[zsel]) - index) * var(z[zsel][(index + 1):length(z[zsel])]))
  
    threshold <- mean(z[zsel][which.min(varest) + 0:1])
    if(sum(z > max(threshold)) < (minperc / 100) * length(z)) zsel[union(max(which(zsel)), which(z > max(threshold)))] <- FALSE
    if(sum(z < min(threshold)) < (minperc / 100) * length(z)) zsel[union(min(which(zsel)), which(z < min(threshold)))] <- FALSE
    if(sum(z > max(threshold)) > (minperc / 100) * length(z) & sum(z < max(threshold)) > (minperc / 100) * length(z)) break
  }
  return(threshold)
                                        #return(list(denz = denz, threshold = threshold))
}
  
