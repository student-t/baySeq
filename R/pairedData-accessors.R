setMethod("[", "pairedData", function(x, i, j, ..., drop = FALSE) {
  x <- callNextMethod(x, i, j, ..., drop = FALSE)
  if(missing(j))
    j <- 1:ncol(x@data)
  if(missing(i))
    i <- 1:nrow(x@data)

  x@pairData <- x@pairData[i,j, drop = FALSE]
  if(length(x@pairLibsizes) > 0) x@pairLibsizes <- x@pairLibsizes[j]
  x
})

setMethod("show", "pairedData", function(object) {
  cat(paste('An object of class "', class(object), '"\n', sep = ""))
  cat(paste(nrow(object), 'rows and', ncol(object), 'columns\n'))
  
  cat('\nSlot "replicates"\n')
  print(object@replicates)

  if(length(object@libsizes) > 0) {
    cat('\nSlot "libsizes"\n')
    print(object@libsizes)
  }

  if(length(object@groups) > 0) {
    cat('\nSlot "groups":\n')
    print(object@groups)
  }

  cat('\nSlot "data":\n')
  if(nrow(object) > 5)
    {
      print(object@data[1:5,])
      cat(paste(nrow(object) - 5), "more rows...\n")
    } else print(object@data)

  cat('\nSlot "pairData":\n')
  if(nrow(object) > 5)
    {
      print(object@pairData[1:5,])
      cat(paste(nrow(object) - 5), "more rows...\n")
    } else print(object@pairData)

  if(nrow(object@annotation) > 0) {
    cat('\nSlot "annotation":\n')
    if(nrow(object@annotation) > 5 & ncol(object@annotation) > 0)
      {
        print(object@annotation[1:5,])
        cat(paste(nrow(object) - 5), "more rows...\n")
      } else print(object@annotation)
  }
    
  if(nrow(object@posteriors) > 0)
    {
      cat('Slot "posteriors":\n')
      if(nrow(object@posteriors) > 5)
        {
          print(exp(object@posteriors[1:5,]))
          cat(paste(nrow(object) - 5), "more rows...\n")
        } else print(exp(object@posteriors))
    }
  if(length(object@estProps) > 0)
    {
      cat('\nSlot "estProps":\n')
      print(object@estProps)
    }
  if(length(object@priorType) > 1)
    {
      cat('Slot "priors":\n')
      cat(paste('Priors are of type:', object@priorType), '\n')
    }
})


setMethod("initialize", "pairedData", function(.Object, ...) {
  .Object <- callNextMethod(.Object, ...)
  if(any(dim(.Object@pairData) != dim(.Object@data)))
    stop("Dimensions of '@pairData' slot must be the same as that of the '@data' slot.")
  .Object
})


setMethod("libsizes<-", signature = "pairedData", function(x, value) {
  if(!is.list(value) | length(value) != 2) stop("Library sizes should be given as a list of length 2 for a pairedData object.")
  if(any(sapply(value, function(lib) !is.numeric(lib)))) stop("All members of 'libsizes' for a pairedData object must be numeric.")
  if(any(sapply(value, length) != ncol(x))) stop("Length of libsizes must be identical to the number of columns of the countData object.")
  if(any(sapply(value, function(lib) any(lib <= 0)))) stop("Library sizes less than or equal to zero make no sense to me!")
  x@libsizes <- value[[1]]
  x@pairLibsizes <- value[[2]]
  x
})

setMethod("libsizes", signature = "pairedData", function(x) {
  list(libsizes = x@libsizes, pairLibsizes = x@pairLibsizes)  
})
