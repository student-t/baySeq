setMethod("[", "pairedData", function(x, i, j, ..., drop = FALSE) {
  x <- callNextMethod(x, i, j, ..., drop = FALSE)
  if(missing(j))
    j <- 1:ncol(x@data)
  if(missing(i))
    i <- 1:nrow(x@data)

  x@pairData <- x@pairData[i,j, drop = FALSE]
  x
})

setMethod("show", "pairedData", function(object) {
  cat(paste('An object of class "', class(object), '"\n', sep = ""))
  cat(paste(nrow(object), 'rows and', ncol(object), 'columns\n'))
  
  cat('\nSlot "replicates"\n')
  print(object@replicates)
  
  cat('\nSlot "libsizes"\n')
  print(object@libsizes)

  cat('\nSlot "groups":\n')
  print(object@groups)

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
    } else print(object@data)

  
  cat('\nSlot "annotation":\n')
  if(nrow(object@annotation) > 5 & ncol(object@annotation) > 0)
    {
      print(object@annotation[1:5,])
      cat(paste(nrow(object) - 5), "more rows...\n")
    } else print(object@annotation)

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
