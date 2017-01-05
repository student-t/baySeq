
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


bimodalSeparator <- function(x, weights = NULL, minperc = 0.1, elbow = NULL) {
  
  if(is.null(weights)) weights <- rep(1, length(x))

  nas <- is.na(x) | is.na(weights)
  x <- x[!nas]; weights <- weights[!nas]
  
  xnif <- (x < Inf & x > -Inf)
  z <- x[xnif]
  weights <- weights[xnif]

#  weights <- round(weights/min(weights))
#  weights <- round(weights / max(weights) * 1000)
#  z <- sort(rep(x, weights))

  rordz <- order(z)
  z <- z[rordz]
  weights <- weights[rordz]
  weights <- weights / sum(weights)
  
  weighted.var <- function(x, w)
    sum(w * (x - weighted.mean(x, w))^2) / (sum(w) - sum(w^2) / sum(w))
  
  zsel <- rep(TRUE, length(z))

  cumwMean <- function(x, weights)
      cumsum(x * weights) / cumsum(weights)
  cumwVar <- function(x, weights) {
      wm <- cumwMean(x, weights)
      (cumsum(weights * x^2) + wm^2 * cumsum(weights) - 2 * wm * cumsum(weights * x)) / (cumsum(weights) - cumsum(weights^2) / cumsum(weights))
  }
  
  while(TRUE) {
      varest <- (1:length(z[zsel]) * cumwVar(z[zsel], weights[zsel]))[-sum(zsel)] +
          c(rev(1:length(z[zsel]) * cumwVar(rev(z[zsel]), rev(weights[zsel])))[c(-1, -sum(zsel))], NaN)

      minvar <- which.min(varest)
      threshold <- mean(z[zsel][which.min(varest) + 0:1])

      if(!is.null(elbow) && elbow == "left") {
          x0 <- z[zsel][2]; y0 <- varest[2]; x1 <- threshold; y1 <- min(varest, na.rm = TRUE)
          A = y1 - y0; B = -(x1 - x0); C = + (x1 - x0) * y0 - (y1 - y0) * x0
          threshold <- z[zsel][which.max((abs(A * z[zsel][-sum(zsel)] + B * varest + C) / sqrt(A^2 + B^2))[1:minvar])]
      } else if(!is.null(elbow) && elbow == "right") {
          x0 <- z[zsel][length(varest) - 1]; y0 <- varest[length(varest) - 1]; x1 <- threshold; y1 <- min(varest, na.rm = TRUE)
          A = y1 - y0; B = -(x1 - x0); C = + (x1 - x0) * y0 - (y1 - y0) * x0
          threshold <- z[zsel][
                           which.max((abs(A * z[zsel][-sum(zsel)] + B * varest + C) / sqrt(A^2 + B^2))[minvar:length(varest)]) + minvar - 1]
      }
      
    if(sum(z > max(threshold)) < (minperc / 100) * length(z)) zsel[union(max(which(zsel)), which(z > max(threshold)))] <- FALSE
    if(sum(z < min(threshold)) < (minperc / 100) * length(z)) zsel[union(min(which(zsel)), which(z < min(threshold)))] <- FALSE
    if(sum(z > max(threshold)) > (minperc / 100) * length(z) & sum(z < max(threshold)) > (minperc / 100) * length(z)) break
  }
  return(threshold)
                                        #return(list(denz = denz, threshold = threshold))
}
  
