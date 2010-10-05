plotMA.CD <- function(cD, samplesA, samplesB, ...)
{
  if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")
  
  Adata <- colSums(t(cD@data[,samplesA]) / cD@libsizes[samplesA]) / length(samplesA)
  Bdata <- colSums(t(cD@data[,samplesB]) / cD@libsizes[samplesA]) / length(samplesB)

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
  axis(side = 2, at = c(-infRatio - 1, axTicks(side = 2), infRatio + 1), labels = c(-Inf, axTicks(side = 2), Inf))
       
  abline(h = c(-1,1) * (1 + infRatio), col = "orange", lty = 3)

}
