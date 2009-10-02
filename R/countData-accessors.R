
setMethod("[", "countData", function(x, i, j, ..., drop = FALSE) {
  if(missing(j))
    j <- 1:ncol(x@data)
  if(!missing(j))
    if(!all(1:ncol(x@data) %in% j))
      {
        newgroups <- list()
        for(gg in 1:length(x@groups))
            newgroups[[gg]] <- x@groups[[gg]][j]
        warning("Selection of samples (columns) may adversely affect the argument in slot 'groups'.")
        x@groups <- newgroups
      }
  if(missing(i))
    i <- 1:nrow(x@data)
  x@data <- x@data[i,j, drop = drop]; x@libsizes <- x@libsizes[j]; x@annotation <- x@annotation[i,, drop = FALSE]
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
  cat('\nSlot "priors":\n')
  print(object@priors)
})
