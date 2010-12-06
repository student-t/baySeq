plotMA.CD <- function(cD, samplesA, samplesB, normaliseData = TRUE, ...)
{
  if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")
  
  Adata <- colSums(t(cD@data[,samplesA])) / length(samplesA)
  Bdata <- colSums(t(cD@data[,samplesB])) / length(samplesB)

  if(normaliseData) {
    Adata <- Adata / cD@libsizes[samplesA] * mean(cD@libsizes[c(samplesA, samplesB)])
    Bdata <- Bdata / cD@libsizes[samplesB] * mean(cD@libsizes[c(samplesA, samplesB)])
  }
                   
  Azeros <- which(Adata == 0)
  Bzeros <- which(Bdata == 0)
  infRatio <- ceiling(max(abs((log2(Adata) - log2(Bdata))[-union(Azeros, Bzeros)])))

  M <- log2(Adata) - log2(Bdata)
  M[Azeros] <- -infRatio - 2
  M[Bzeros] <- infRatio + 2

  A <- (log2(Adata) + log2(Bdata)) / 2
  A[Azeros] <- log2(Bdata[Azeros])
  A[Bzeros] <- log2(Adata[Bzeros])
  
  plot(y = M, x = A, ylim = c(-infRatio - 3, infRatio + 3), axes = FALSE, ylab = "M", xlab = "A", ...)
  axis(side = 1)

  maxis <- axTicks(side = 2)
  maxis <- maxis[maxis <= infRatio + 1 - max(diff(maxis)) & maxis >= -infRatio - 1 + max(diff(maxis))]
  
  axis(side = 2, at = c(-infRatio - 1, maxis, infRatio + 1), labels = c(-Inf, maxis, Inf))
       
  abline(h = c(-1,1) * (1 + infRatio), col = "orange", lty = 3)

}
