.catObservables <- function(cD) {
  cellObservables <- cD@cellObservables
  sampleObservables <- cD@sampleObservables
  rowObservables <- cD@rowObservables
  if(!("seglens" %in% names(cellObservables) || "seglens" %in% names(rowObservables))) rowObservables <- c(rowObservables, list(seglens = rep(1, nrow(cD))))
  if(!("libsizes" %in% names(sampleObservables))) sampleObservables <- c(sampleObservables, list(seglens = rep(1, ncol(cD))))
  
  observables <- c(cellObservables,
                   lapply(sampleObservables, function(obs)
                          array(rep(obs, each = nrow(cD)), dim = dim(cD))),
                   lapply(rowObservables, function(obs)
                          array(rep(obs, ncol(cD)), dim = dim(cD))))
  observables
}

.fastUniques <- function(x){
  if (nrow(x) > 1) {
    return(c(TRUE, rowSums(x[-1L, , drop = FALSE] == x[-nrow(x),, drop = FALSE]) != ncol(x)))
  } else return(TRUE)
}    


.showData <- function(data)
  {
    if(is.vector(data) || length(dim(data)) <= 2) return(data)
    dimsep <- c(":", "|")
    dimlen <- length(dim(data))    
    if(length(dim(data)) > 4)
      dimsep <- c(dimsep, sapply(2:(dimlen - 2), function(x) paste(rep("|", x), collapse = "")))
    dimsep <- c("", dimsep)
    dimsep <- dimsep[1:(dimlen - 1)]

    dimsep <- rev(dimsep)
    dimdat <- data    

    pasteDat <- function(x, dimnum) {
      if(length(dim(x)) > 2) {
        padat <- t(apply(x, 1, function(xx) paste(pasteDat(xx, dimnum = dimnum + 1), collapse = dimsep[dimnum])))
      } else {
        padat <- (apply(x, 1, function(z) paste(z, collapse = ":")))
      }
      return(padat)
    }
    pastemat <- t(apply(data, 1, pasteDat, dimnum = 1))
    pastemat
  }
      

.sliceArray <- function(slices, array, drop = FALSE) {
  if((is.vector(array) & sum(!sapply(slices, is.null)) > 1) || (is.array(array) & length(slices) > length(dim(array)))) warning("dimensions of slice exceed dimensions of array")

  sarray <- asub(array, slices, dims = 1:length(slices), drop = drop)
  
  sarray
}


.logsum <- function(x)
  max(x, max(x, na.rm = TRUE) + log(sum(exp(x - max(x, na.rm = TRUE)), na.rm = TRUE)), na.rm = TRUE)

.logRowSum <- function(z)
  {
    maxes <- do.call(pmax, c(as.list(data.frame(z)), list(na.rm = TRUE)))
    pmax(maxes, maxes + log(rowSums(exp(z - maxes), na.rm = TRUE)), na.rm = TRUE)
  }
