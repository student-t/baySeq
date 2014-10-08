plotPosteriors <- function(cD, group, samplesA, samplesB, ...)
{
  if(length(dim(cD)) > 2) stop("This function is currently only applicable to 2-dimensional countData objects.")
  
  if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")

  if(is.character(group))
    group <- pmatch(group, names(cD@groups))
  
  if(missing(samplesA))
    samplesA <- which(cD@groups[[group]] == 1)
  
  if(missing(samplesB))
    samplesB <- which(cD@groups[[group]] == 2)

  if(length(samplesA) == 0 | length(samplesB) == 0)
    stop("Sample information is missing, and cannot be inferred.")

  if("libsizes" %in% names(cD@sampleObservables)) libsizes <- as.vector(cD@sampleObservables$libsizes) else libsizes <- rep(1, ncol(cD))
      
  if(nrow(cD@posteriors) > 0)
    {
      Adata <- colSums(t(cD@data[,samplesA]) / libsizes[samplesA]) / length(samplesA)
      Bdata <- colSums(t(cD@data[,samplesB]) / libsizes[samplesB]) / length(samplesB)

      Azeros <- which(Adata == 0)
      Bzeros <- which(Bdata == 0)
      
      bexp <- log2(Bdata[Azeros] * mean(libsizes[c(samplesA, samplesB)]))
      aexp <- log2(Adata[Bzeros] * mean(libsizes[c(samplesA, samplesB)]))      

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
