setGeneric("rbind", function(..., deparse.level=1) standardGeneric("rbind"), signature = "...")

setMethod("rbind", "countData", function(x, ..., deparse.level = 1) {
  print(nargs())
  if(nargs() < 4) rbind2(x, ...) else rbind2(x, Recall(...))
})


setGeneric("groups<-", function(x, value) standardGeneric("groups<-"))
setMethod("groups<-", signature = "countData", function(x, value) {
  x@groups <- lapply(value, as.factor)
  x
})

setGeneric("groups", function(x) standardGeneric("groups"))
setMethod("groups", signature = "countData", function(x) {
  x@groups
})

setGeneric("libsizes<-", function(x, value) standardGeneric("libsizes<-"))
setMethod("libsizes<-", signature = "countData", function(x, value) {
  if(!is.numeric(value)) stop("All members of libsizes for a countData object must be numeric.")
  if(length(value) != ncol(x)) stop("Length of libsizes must be identical to the number of columns of the countData object.")
  if(any(value <= 0)) stop("Library sizes less than or equal to zero make no sense to me!")
  x@libsizes <- value
  x
})

setGeneric("libsizes", function(x) standardGeneric("libsizes"))
setMethod("libsizes", signature = "countData", function(x) {
  x@libsizes
})

setGeneric("replicates<-", function(x, value) standardGeneric("replicates<-"))
setMethod("replicates<-", signature = "countData", function(x, value) {
  x@replicates <- as.factor(value)
  x
})

setGeneric("replicates", function(x) standardGeneric("replicates"))
setMethod("replicates", signature = "countData", function(x) {
  x@replicates
})

setMethod("rbind2", "countData", function(x, y) {
  if(ncol(x) != ncol(y)) stop("Column numbers are not identical across the objects")

  data <- rbind(x@data, y@data)
  if(!all(x@libsizes == y@libsizes)) warning("'@libsizes' are not identical across objects")
  if(!all(x@replicates == y@replicates)) warning("'@replicates' are not identical across objects")

  xann <- x@annotation
  yann <- y@annotation

  incols <- intersect(colnames(xann), colnames(yann))

  zann <- rbind(subset(xann, select = incols), subset(yann, select = incols))

  unxann <- setdiff(colnames(xann), incols)
  if(length(unxann) > 0)
    zann <- cbind(zann, rbind(subset(xann, select = unxann), matrix(NA, ncol = length(unxann), nrow = nrow(yann), dimnames = list(NULL, unxann))))
  
  unyann <- setdiff(colnames(yann), incols)
  if(length(unyann) > 0)
    zann <- cbind(zann, rbind(matrix(NA, ncol = length(unyann), nrow = nrow(xann), dimnames = list(NULL, unyann)), subset(yann, select = unyann)))                

  if(nrow(x@posteriors) > 0 & nrow(y@posteriors) > 0)
    {
      if(length(x@groups) == length(y@groups)) {
        if(all(sapply(1:length(x@groups), function(ii) x@groups[[ii]] == y@groups[[ii]]))) {
          posteriors <- rbind(x@posteriors, y@posteriors)
        } else warning("'@groups' slots are not identical; posterior likelihoods will be discarded.")
      } else warning("'@groups' slots are not identical; posterior likelihoods will be discarded.")
    }
      



  z <- new(class(x), data = data, annotation = zann, posteriors = posteriors, libsizes = x@libsizes, replicates = x@replicates)
  
  if("groups" %in% slotNames(y) & "groups" %in% slotNames(x)) {
    z@groups <- c(x@groups, y@groups)
    if(length(z@groups) > 0)
      z@groups <- z@groups[!duplicated(groups)]
  }
  if("groups" %in% slotNames(x)) z@groups <- x@groups
      
  z
})
  


setMethod("initialize", "countData", function(.Object, ..., replicates, seglens) {
  .Object <- callNextMethod(.Object, ...)
  if(length(.Object@libsizes) != ncol(.Object@data) & length(.Object@libsizes) != 0)
    stop("Length of '@libsizes' slot, if provided, must equal number of columns of '@data' slot.")
  if(nrow(.Object@annotation) > 0 & nrow(.Object@annotation) != nrow(.Object@data))
    warning("Number of rows of '@annotation' slot not same as '@data' slot.")
  
  if(any(lapply(.Object@groups, length) != ncol(.Object@data)))
    stop("All vectors in '@groups' slot must equal number of columns of '@data' slot.")

  if(ncol(.Object@posteriors) != length(.Object@groups) & ncol(.Object@posteriors) != 0)
    stop("Number of columns in '@posteriors' slot must equal length of '@groups' slot.")

  if(length(.Object@nullPosts) != 0) {
    if(nrow(.Object@nullPosts) != nrow(.Object@data) & nrow((.Object@nullPosts) != 0))
      stop("Number of rows in '@data' slot must equal number of rows of '@nullPosts' slot.")
  } else nullPosts <- matrix(ncol = 0, nrow = nrow(.Object@data))
  
  if(length(.Object@estProps) != length(.Object@groups) & length(.Object@estProps) != 0)
    stop("Length of '@estProps' slot must equal length of '@groups' slot.")

  .Object@groups <- lapply(.Object@groups, as.factor)

  if(!missing(seglens))
    {
      if(is.vector(seglens))
        seglens <- matrix(seglens, ncol = 1)
      .Object@seglens <- seglens
      
  if(nrow(.Object@seglens) != nrow(.Object@data) & nrow(.Object@seglens) != 0)
    stop("If 'seglens' specified, the number of rows (or length if submitting as vector) of '@seglens' slot must equal number of rows of '@data' slot.")
  if(ncol(.Object@seglens) != 1 & ncol(.Object@seglens) != ncol(.Object@data))
    stop("If 'seglens' specified, it must either be a vector, a matrix with one column, or a matrix with the same number of columns as 'data'.")
    }
  if(!missing(replicates)) {
    if(length(replicates) != ncol(.Object@data))
      stop("The length of the '@replicates' slot must equal number of columns of '@data' slot.")
    .Object@replicates <- as.factor(replicates)
  } #else if(length(.Object@replicates) == 0) stop("Replicate data must be provided")

  if(length(colnames(.Object@data)) == 0) colnames(.Object@data) <- make.unique(c(as.character(unique(.Object@replicates)), as.character(.Object@replicates)))[-(1:(length(unique(.Object@replicates))))]
  if(length(.Object@libsizes) != 0)
    names(.Object@libsizes) <- colnames(.Object@data)
  .Object
})

setMethod("[", "countData", function(x, i, j, ..., drop = FALSE) {
  if(missing(j))
    j <- 1:ncol(x@data)
  if(!missing(j))
    if(!all(1:ncol(x@data) %in% j))
      {
        replicates(x) <- as.character(x@replicates[j])

        if(length(x@groups) > 0)
          {
            newgroups <- list()
            newgroups <- lapply(x@groups, function(x) {
              x[j]
              rep(1:length(unique(x[j])), sapply(unique(x[j]), function(z) sum(x[j] == z)))[unlist(sapply(unique(x[j]), function(z) which(x[j] == z)))]
            })
            x@groups <- newgroups[!duplicated(newgroups) | duplicated(x@groups)]
          }
            
        if(length(x@posteriors) > 0)
          {
            warning("Selection of samples (columns) will invalidate the values calculated in slot 'posteriors', and so these will be discarded.")
            x@posteriors <- matrix(nrow = 0, ncol = 0)
          }
      }
  if(missing(i))
    i <- 1:nrow(x@data)

  x@data <- x@data[i,j, drop = FALSE]
  if(length(x@libsizes) > 0) x@libsizes <- x@libsizes[j]
  x@annotation <- x@annotation[i,, drop = FALSE]
  if(nrow(x@posteriors) > 0)
    x@posteriors <- x@posteriors[i,, drop = FALSE]

  if(nrow(x@nullPosts) > 0)
    x@nullPosts <- x@nullPosts[i,,drop = FALSE]

  if(nrow(x@seglens) > 0)
    {
      if(ncol(x@seglens) == 1) {
        x@seglens <- x@seglens[i,, drop = FALSE]
      } else x@seglens <- x@seglens[i, j, drop = FALSE]
    }
  
  x
})

setMethod("dim", "countData", function(x) {
  dim(x@data)
})

setMethod("show", "countData", function(object) {

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

