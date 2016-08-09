#setGeneric("rbind", function(..., deparse.level=1) standardGeneric("rbind"), signature = "...")

#setMethod("rbind", "countData", function(x, ..., deparse.level = 1) {
#  print(nargs())
#  if(nargs() < 4) rbind2(x, ...) else rbind2(x, Recall(...))
#})

setMethod("summary", "countData", function(object, ...) {
  desc <- cat(paste('An object of class "', class(object), '"\n', sep = ""),
              paste(nrow(object), 'rows and', ncol(object), 'columns\n'))

  if(nrow(object@posteriors) > 0) {
    .getDEC <- function(cd, type = "expectations") {
      metric <- do.call("c", lapply(1:length(cd@groups), function(gg) {
        ords <- levels(cd@orderings[,gg])
        if(type == "expectations") sums <- sapply(ords, function(ord) sum(exp(cd@posteriors[,gg]) * (cd@orderings[,gg] ==  ord)))
        if(type == "FDR") {
            sums <- sapply(ords, function(ord) {
                if(sum(cd@orderings[,gg] == ord) > 0) {
                    tc <- suppressWarnings(topCounts(cd[cd@orderings[,gg] == ord,], gg, FDR = 0.05))
                    suppressWarnings(sum(tc$FDR < 0.05))
                } else return(0)
          })
        }
        names(sums) <- paste(names(cd@groups)[gg], ords, sep = ":")
        sums
      }))
      metric
    }  
    
    expecs <- round(.getDEC(object), 2)
    fdrsum <- .getDEC(object, type = "FDR")
    summarySum <- cbind("Expected" = expecs, "FDR (< 0.05)" = fdrsum)
    rownames(summarySum) <- names(expecs)
    print(summarySum)
  } else summarySum <- NULL
  invisible(list(desc = desc, expFDR = summarySum))
})
  
setGeneric("groups<-", function(x, value) standardGeneric("groups<-"))
setMethod("groups<-", signature = "countData", function(x, value) {
    if(any(sapply(value, length) != ncol(x))) stop(paste(sum(sapply(value, length) != ncol(x)), "vector(s) in the groups structure are the wrong length."))
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

  if(is.vector(value)) {
    if(length(value) != ncol(x)) stop("Length of libsizes must be identical to the number of columns of the countData object.")
    value <- array(value, length(value))
  } else if(is.array(value))
    if(any(dim(x@data)[-1] != dim(value))) stop("Dimension of libsizes must be identical to the dimension of the countData object (after dropping the first dimension).")
                                        
  if(any(value <= 0)) stop("Library sizes less than or equal to zero make no sense to me!")
  x@sampleObservables$libsizes <- value
  x
})

setGeneric("seglens<-", function(x, value) standardGeneric("seglens<-"))
setMethod("seglens<-", signature = "countData", function(x, value) {
  if(!is.numeric(value)) stop("All members of seglens for a countData object must be numeric.")

  if(inherits(value, "numeric")) {
    if(length(value) != ncol(x)) stop("Length of seglens must be identical to the number of columns of the countData object.")
    value <- matrix(value, ncol = 1)
  } else if(is.array(value))
    if(any(dim(x@data)[-1] != dim(value))) stop("Dimension of seglens must be identical to the dimension of the countData object (after dropping the first dimension).")
                                        
  if(any(value <= 0)) stop("Library sizes less than or equal to zero make no sense to me!")
  x@rowObservables$seglens <- value
  x
})


setGeneric("densityFunction", function(x) standardGeneric("densityFunction"))
setMethod("densityFunction", signature = "countData", function(x) {
  x@densityFunction
})

setGeneric("densityFunction<-", function(x, value) standardGeneric("densityFunction<-"))
setMethod("densityFunction<-", signature = "countData", function(x, value) {
    if(length(value) != 1 & length(value) != length(x@groups)) stop("The given value must be of length 1 or equal to the number of groups of the object.")
    if(is.list(value)) {
        if(any(!sapply(value, function(x) inherits(x, "densityFunction")))) stop("All members of the list must be of (or inherit from) the 'densityFunction' class") else x@densityFunction <- value
    } else if(!inherits(value, "densityFunction")) stop("The given value must be of (or inherit from) the 'densityFunction' class, or be a list object containing only elements of this class") else x@densityFunction <- list(value)        
  x
})

setGeneric("libsizes", function(x) standardGeneric("libsizes"))
setMethod("libsizes", signature = "countData", function(x) {
  x@sampleObservables$libsizes
})

setGeneric("seglens", function(x) standardGeneric("seglens"))
setMethod("seglens", signature = "countData", function(x) {
  if("seglens" %in% names(x@rowObservables)) return(x@rowObservables$seglens)
  if("seglens" %in% names(x@cellObservables)) return(x@cellObservables$seglens)
  return(matrix(rep(1, nrow(x)), ncol = 1))
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

#setMethod("rbind2", "countData", function(x, y) {
#  if(ncol(x) != ncol(y)) stop("Column numbers are not identical across the objects")

#  data <- rbind(x@data, y@data)
#  if(!all(x@libsizes == y@libsizes)) warning("'@libsizes' are not identical across objects")
#  if(!all(x@replicates == y@replicates)) warning("'@replicates' are not identical across objects")

#  xann <- x@annotation
#  yann <- y@annotation

#  incols <- intersect(colnames(xann), colnames(yann))

#  zann <- rbind(subset(xann, select = incols), subset(yann, select = incols))

#  unxann <- setdiff(colnames(xann), incols)
#  if(length(unxann) > 0)
#    zann <- cbind(zann, rbind(subset(xann, select = unxann), matrix(NA, ncol = length(unxann), nrow = nrow(yann), dimnames = list(NULL, unxann))))
  
#  unyann <- setdiff(colnames(yann), incols)
#  if(length(unyann) > 0)
#    zann <- cbind(zann, rbind(matrix(NA, ncol = length(unyann), nrow = nrow(xann), dimnames = list(NULL, unyann)), subset(yann, select = unyann)))                

#  if(nrow(x@posteriors) > 0 & nrow(y@posteriors) > 0)
#    {
#      if(length(x@groups) == length(y@groups)) {
#        if(all(sapply(1:length(x@groups), function(ii) x@groups[[ii]] == y@groups[[ii]]))) {
#          posteriors <- rbind(x@posteriors, y@posteriors)
#        } else warning("'@groups' slots are not identical; posterior likelihoods will be discarded.")
#      } else warning("'@groups' slots are not identical; posterior likelihoods will be discarded.")
#    }
      



#  z <- new(class(x), data = data, annotation = zann, posteriors = posteriors, libsizes = x@libsizes, replicates = x@replicates)
  
#  if("groups" %in% slotNames(y) & "groups" %in% slotNames(x)) {
#    z@groups <- c(x@groups, y@groups)
#    if(length(z@groups) > 0)
#      z@groups <- z@groups[!duplicated(groups)]
#  }
#  if("groups" %in% slotNames(x)) z@groups <- x@groups
      
#  z
#})
  


setMethod("initialize", "countData", function(.Object, ..., data, replicates, libsizes, seglens, densityFunction) {
  
  .Object <- callNextMethod(.Object, ...)
  
  if(!missing(data) && is.array(data)) .Object@data <- data
  if(!missing(data) && is.list(data)) .Object@data <- array(do.call("c", data), c(dim(data[[1]]), length(data)))  
  if(missing(replicates)) replicates <- .Object@replicates
  .Object@replicates <- as.factor(replicates)
  if(!missing(densityFunction) && inherits(densityFunction, "densityFunction")) densityFunction(.Object) <- list(densityFunction)

  if(length(dim(.Object@data)) == 1) .Object@data <- array(.Object@data, dim = c(dim(.Object@data), max(c(0, length(replicates), length(.Object@replicates)))))
  
  if(length(colnames(.Object@data)) == 0) colnames(.Object@data) <- make.unique(c(as.character(unique(.Object@replicates)), as.character(.Object@replicates)))[-(1:(length(unique(.Object@replicates))))]
  
  if(nrow(.Object@annotation) > 0 & nrow(.Object@annotation) != nrow(.Object@data))
    warning("Number of rows of '@annotation' slot not same as '@data' slot.")
  
  if(any(lapply(.Object@groups, length) != ncol(.Object@data)))
    stop("All vectors in '@groups' slot must equal number of columns of '@data' slot.")

  if(ncol(.Object@posteriors) != length(.Object@groups) & ncol(.Object@posteriors) != 0)
    stop("Number of columns in '@posteriors' slot must equal length of '@groups' slot.")

  if(length(.Object@densityFunction) > 1 & length(.Object@densityFunction) != length(.Object@groups))
      stop("Length of list of densityFunctions in '@densityFunction' slot must be 1 or equal to the length of the '@groups' slot.")
  
  if(length(.Object@nullPosts) != 0) {
    if(nrow(.Object@nullPosts) != nrow(.Object@data) & nrow((.Object@nullPosts) != 0))
      stop("Number of rows in '@data' slot must equal number of rows of '@nullPosts' slot.")
  } else nullPosts <- matrix(ncol = 0, nrow = nrow(.Object@data))
  
  if(length(.Object@estProps) != length(.Object@groups) & length(.Object@estProps) != 0)
    stop("Length of '@estProps' slot must equal length of '@groups' slot.")

  .Object@groups <- lapply(.Object@groups, as.factor)

  if(!missing(libsizes)) {
    if(is.array(libsizes) && (any(dim(libsizes) != dim(.Object@data)[-1])) || (is.vector(libsizes) & length(libsizes) != ncol(.Object@data)))       
      stop("If provided, the 'libsizes' variable must be a vector of equal length to the columns of the `@data' array or an array of equal dimension to a row of the `@data' array")
    if(is.array(libsizes) && is.null(colnames(libsizes))) colnames(libsizes) <- colnames(.Object@data)
    if(is.vector(libsizes) && is.null(names(libsizes))) names(libsizes) <- colnames(.Object@data)
    .Object@sampleObservables$libsizes <- libsizes
  }
  
  if(!missing(seglens))
    {
      if(is.vector(seglens)) {
        if(length(seglens) != nrow(.Object@data)) stop("If 'seglens' specified, and is a vector, the length of this variable must equal the number of rows of '@data' slot.")          
        .Object@rowObservables$seglens <- seglens
      }
      if(is.array(seglens)) {        
        if(length(dim(.Object@data)) != length(dim(seglens)) || (any(dim(.Object@data) != dim(seglens)))) stop("If 'seglens' specified, and is an array, the dimensions of this variable must equal the dimensions of the '@data' slot.")
        .Object@cellObservables$seglens <- seglens
      }
    }  

  if(length(.Object@rowObservables) > 0) {
    notRow <- sapply(.Object@rowObservables, length) != nrow(.Object@data)
    if(any(notRow)) stop(paste("The following '@rowObservables' elements have an incorrect length:", paste(names(notRow)[notRow], collapse = ",")))
  }
  if(length(.Object@sampleObservables) > 0) {
    notSample <- sapply(.Object@sampleObservables, function(x)
                        (is.vector(x) && length(x) != ncol(.Object@data)) | (is.array(x) && ((length(dim(x)) != length(dim(.Object@data)) - 1) | any(dim(x) != dim(.Object@data)[-1]))))
    
    if(any(notSample)) stop(paste("The following '@sampleObservables' elements have an incorrect length:", paste(names(notSample)[notSample], collapse = ",")))
  }
  if(length(.Object@cellObservables) > 0) {
    notCell <- sapply(.Object@cellObservables, function(oco) any(dim(oco)[1:2] != dim(.Object@data)[1:2]))
    if(any(notCell)) stop(paste("The following '@cellObservables' elements have incorrect dimensions:", paste(names(notCell)[notCell], collapse = ",")))
  }
  
    if(length(replicates) != 0 && length(replicates) != ncol(.Object@data))
      stop("The length of the '@replicates' slot must equal number of columns of '@data' slot.")

  .Object
})

setMethod("[", "countData", function(x, i, j, ..., drop = FALSE) {
  if(missing(j)) {
    j <- 1:ncol(x@data)
  } else {
    if(is.logical(j)) j <- which(j)
    
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
        if(length(x@orderings) > 0)
          {
            warning("Selection of samples (columns) will invalidate the values calculated in slot 'orderings', and so these will be discarded.")
            x@orderings <- data.frame()
          }

      }
  }

  
  if(missing(i))
    i <- 1:nrow(x@data)
  if(is.logical(i)) i <- which(i)

  x@data <- .sliceArray(list(i, j), x@data)
  
  x@annotation <- x@annotation[i,, drop = FALSE]
  if(nrow(x@posteriors) > 0)
    x@posteriors <- x@posteriors[i,, drop = FALSE]
  if(nrow(x@orderings) > 0)
    x@orderings <- x@orderings[i,, drop = FALSE]  
  if(length(x@nullPosts) > 0)
    x@nullPosts <- x@nullPosts[i,,drop = FALSE]

  x@rowObservables <- lapply(x@rowObservables, function(z) .sliceArray(list(i),z, drop = FALSE))
  x@sampleObservables <- lapply(x@sampleObservables, function(z) .sliceArray(list(j), z, drop = FALSE))
  x@cellObservables <- lapply(x@cellObservables, function(z) .sliceArray(list(i,j), z, drop = FALSE))

#  if(length(x@libsizes) > 0) x@libsizes <- .sliceArray(list(j), x@libsizes, drop = FALSE)
#  if(nrow(x@seglens) > 0)
#    {
#      if(ncol(x@seglens) == 1) {
#        x@seglens <- x@seglens[i,, drop = FALSE]
#      } else x@seglens <- x@seglens[i, j, drop = FALSE]
#    }
  
  x
})

setMethod("dim", "countData", function(x) {
  dim(x@data)
})

.pasteMatrixColon <- function(x, y) matrix(paste(x, y, sep = ":"), ncol = ncol(x), nrow = nrow(x))

setMethod("show", "countData", function(object) {

  cat(paste('An object of class "', class(object), '"\n', sep = ""))
  cat(paste(nrow(object), 'rows and', ncol(object), 'columns\n'))
  
  cat('\nSlot "replicates"\n')
  cat(as.character(object@replicates))
  
#  cat('\nSlot "libsizes"\n')
#  print(object@libsizes)

  cat('\nSlot "groups":\n')
  print(object@groups)

  cat('\nSlot "data":\n')

  if(nrow(object) > 5)
    {
      print(.showData(.sliceArray(list(1:5), object@data)), quote = FALSE)        
      cat(paste(nrow(object) - 5), "more rows...\n")
    } else print(.showData(object@data))
  
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

