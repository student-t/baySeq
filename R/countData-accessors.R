setMethod("initialize", "countData", function(.Object, ...) {
  .Object <- callNextMethod()
  if(length(.Object@libsizes) != ncol(.Object@data))
    stop("Length of '@libsizes' slot must equal number of columns of '@data' slot.")
  if(nrow(.Object@annotation) > 0 & nrow(.Object@annotation) != nrow(.Object@data))
    warning("Number of rows of '@annotation' slot not same as '@data' slot.")

  if(any(lapply(.Object@groups, length) != ncol(.Object@data)))
    stop("All vectors in '@groups' slot must equal number of columns of '@data' slot.")

  if(ncol(.Object@posteriors) != length(.Object@groups) & ncol(.Object@posteriors) != 0)
    stop("Number of columns in '@posteriors' slot must equal length of '@groups' slot.")
  
  if(length(.Object@nullPosts) != nrow(.Object@data) & length((.Object@nullPosts) != 0))
    stop("Number of rows in '@data' slot must equal length of '@nullPosts' slot.")
  
  if(length(.Object@estProps) != length(.Object@groups) & length(.Object@estProps) != 0)
    stop("Length of '@estProps' slot must equal length of '@groups' slot.")
  .Object
})

setMethod("initialize", "segData", function(.Object, ..., seglens) {
  .Object <- callNextMethod(.Object, ...)
  if(!missing(seglens))
    {
      if(is.vector(seglens))
        seglens <- matrix(seglens, ncol = 1)
      .Object@seglens <- seglens
    }
  if(nrow(.Object@seglens) != nrow(.Object@data))
    stop("Number of rows (or length if submitting as vector) of '@seglens' slot must equal number of rows of '@data' slot.")
  .Object
})

setMethod("[", "countData", function(x, i, j, ..., drop = FALSE) {
  if(missing(j))
    j <- 1:ncol(x@data)
  if(!missing(j))
    if(!all(1:ncol(x@data) %in% j))
      {
        newgroups <- list()
        for(gg in 1:length(x@groups))
            newgroups[[gg]] <- x@groups[[gg]][j]
        warning("Selection of samples (columns) may adversely affect the argument in slot 'groups' and invalidate the values calculated in slot 'posteriors'.")
        x@groups <- newgroups
      }
  if(missing(i))
    i <- 1:nrow(x@data)
  
  x@data <- x@data[i,j, drop = FALSE]
  x@libsizes <- x@libsizes[j]
  x@annotation <- x@annotation[i,, drop = FALSE]
  if(nrow(x@posteriors) > 0)
    x@posteriors <- x@posteriors[i,]
  x
})

setMethod("[", "segData", function(x, i, j, ..., drop = FALSE) {
  if(missing(i))
    i <- 1:nrow(x@data)
  if(missing(j))
    j <- 1:ncol(x@data)
  x <- callNextMethod(x, i, j, ..., drop = FALSE)

  if(ncol(x@seglens) == 1) {
    x@seglens <- x@seglens[i,, drop = FALSE]
  } else x@seglens <- x@seglens[i, j, drop = FALSE]

  x
})


setMethod("dim", "countData", function(x) {
  dim(x@data)
})

setMethod("show", "countData", function(object) {
  cat(paste('An object of class "', class(object), '"\n', sep = ""))
  cat(paste(nrow(object), 'rows and', ncol(object), 'columns\n'))
  cat('\nSlot "data":\n')
  if(nrow(object) > 10)
    {
      print(object@data[1:10,])
      cat(paste(nrow(object) - 10), "more rows...\n")
    } else print(object@data)
  cat('\nSlot "libsizes":\n')
  print(object@libsizes)
  cat('\nSlot "groups":\n')
  print(object@groups)
  cat('\nSlot "annotation":\n')
  if(nrow(object@annotation) > 10 & ncol(object@annotation) > 0)
    {
      print(object@annotation[1:10,])
      cat(paste(nrow(object) - 10), "more rows...\n")
    } else print(object@annotation)

  if(nrow(object@posteriors) > 0)
    {
      cat('Slot "posteriors":\n')
      if(nrow(object@posteriors) > 10)
        {
          print(object@posteriors[1:10,])
          cat(paste(nrow(object) - 10), "more rows...\n")
        } else print(object@posteriors)
    }
  if(length(object@estProps) > 0)
    {
      cat('\nSlot "estProps":\n')
      print(object@estProps)
    }
  if(length(object@priors@type) > 1)
    {
      cat('Slot "priors":\n')
      cat(paste('Priors are of type:', object@priors@type), '\n')
    }
})
