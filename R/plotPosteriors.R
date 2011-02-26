plotPosteriors <- function(cD, group = 1, samplesA, samplesB, ...)
{
  if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")
  
  if(nrow(cD@posteriors) > 0)
    {
      Adata <- colSums(t(cD@data[,samplesA]) / cD@libsizes[samplesA]) / length(samplesA)
      Bdata <- colSums(t(cD@data[,samplesB]) / cD@libsizes[samplesB]) / length(samplesB)

      Azeros <- which(Adata == 0)
      Bzeros <- which(Bdata == 0)
      
      bexp <- log2(Bdata[Azeros] * mean(cD@libsizes[c(samplesA, samplesB)]))
      aexp <- log2(Adata[Bzeros] * mean(cD@libsizes[c(samplesA, samplesB)]))      

      minZeros <- floor(min(aexp[aexp > -Inf], bexp[bexp > -Inf]))
      maxZeros <- ceiling(max(aexp, bexp))
      
      logData <- log2(Adata / Bdata)
      infRatios <- which(abs(logData) == Inf | is.na(logData))

      nonInfMax <- ceiling(max(abs(logData[-infRatios]))) + 4

      logData[Bzeros] <- aexp + nonInfMax - minZeros
      logData[Azeros] <- -bexp - nonInfMax + minZeros
      
      plot(x = logData,
             y = exp(cD@posteriors[,group]),
             ylim = c(-0.2,1),
           xlim = c(-1,1) * nonInfMax + c(-1, 1) * (maxZeros - minZeros), ylab = "Posterior likelihood", xlab = "", axes = FALSE, ...)
      abline(v = c(-nonInfMax + 3, nonInfMax - 3), col = "orange", lty = 4)
      axis(side = 2, at = c(0:5 / 5))
      axis(side = 1, at = c(-maxZeros - nonInfMax, -nonInfMax, -nonInfMax + 3, -round(nonInfMax / 5), 0, round(nonInfMax / 5), nonInfMax - 3, nonInfMax, nonInfMax + maxZeros),
           labels = c(maxZeros, minZeros, -Inf, -round(nonInfMax / 5), 0, round(nonInfMax / 5), Inf, minZeros, maxZeros))
      text(x = c(-maxZeros - nonInfMax, 0, nonInfMax + maxZeros),
           y = -0.1, labels = c("log B", "log ratio", "log A"))
           
    }
  else stop("No posterior data found in 'cD' object!")
}
