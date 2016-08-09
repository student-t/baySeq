.ciEstimator <- function(posteriors, n = 10000, ci = 0.05) {
  zz <- matrix(rbinom(length(posteriors) * n, 1, posteriors), ncol = n)
  quantile(colSums(zz), probs = c(ci / 2, 1 - ci / 2))
}

summarisePosteriors <- function(cD, orderings = TRUE)
  {
    if(orderings & nrow(cD@orderings) == 0) {
      warning("No orderings contained in countData object.")
      orderings = FALSE
    }
    if(length(cD@nullPosts) > 0) nullPosts <- c(null = sum(exp(cD@nullPosts))) else nullPosts <- NULL
    if(!orderings) {
      return(c(nullPosts, colSums(exp(cD@posteriors), na.rm= TRUE)))
    } else {
      sumOrds <- do.call("c", lapply(1:ncol(cD@orderings),function(ii) {
        sumord <- sapply(levels(cD@orderings[,ii]), function(ord) sum(exp(cD@posteriors[cD@orderings[,ii] == ord,ii]), na.rm = TRUE))
        names(sumord) <- paste(colnames(cD@orderings)[ii], ":", levels(cD@orderings[,ii]), sep = "")
        names(sumord) <- gsub(":$", "", names(sumord))
        sumord
      }))
      return(c(nullPosts, sumOrds))
    }
  }

.controlFDR <- function(likes, FDR) {
  selnum <- max(which(cumsum((1 - sort(likes, decreasing = TRUE)) / 1:length(likes)) < FDR))
  if(selnum > 0) sellikes <- sort(order(likes, decreasing = TRUE)[1:selnum]) else sellikes <- integer()
  sellikes
}

.controlFWER <- function(likes, FWER) {
  llsum <- likes
  selnum <- max(which(1 - cumprod(sort(llsum, decreasing = TRUE)) < FWER))
  if(selnum > 0) sellikes <- sort(order(llsum, decreasing = TRUE)[1:selnum]) else sellikes <- integer()
  sellikes
}

.selectPosteriors <- function(cD, likelihood, FDR, FWER, perReplicate = TRUE) {  
                       #, returnBool = FALSE) {

  returnBool <- FALSE
  if(!missing(likelihood)) {
    selLoc <- cD@posteriors > log(likelihood)
    if(returnBool) return(selLoc) else selLoc <- which(rowSums(selLoc) > 0)
  } else {
    if(!missing(FDR)) {
      controlFunction <- .controlFDR
      controlCrit <- FDR
    } else if(!missing(FWER)) {
      controlFunction <- .controlFWER
      controlCrit <- FWER
    } else stop ("No criterion for locus selection given.")    
    if(perReplicate) {
      selRep <- lapply(1:ncol(cD@posteriors), function(jj) controlFunction(exp(cD@posteriors[,jj]), controlCrit))
      if(returnBool) {
        bool <- do.call("cbind", lapply(1:length(selRep), function(ii) {          
          selBool <- rep(FALSE, nrow(cD))
          if(length(selRep[[ii]]) > 0) selBool[selRep[[ii]]] <- TRUE
          selBool
        }))
        colnames(bool) <- colnames(cD@posteriors)
        return(bool)
      }
      selLoc <- sort(unique(unlist(selRep)))
    } else {
      selLoc <- controlFunction(1 - exp(rowSums(log(1 - exp(cD@posteriors)))), controlCrit)
      if(returnBool) {
        bool <- rep(FALSE, nrow(cD))
        bool[selLoc] <- TRUE
        return(bool)
      }
    }
  }
  if(length(selLoc) == 0) stop("No loci found for the given selection criterion.")
  cD[selLoc,]
}
