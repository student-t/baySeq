'getLibsizes' <- function(cD, data, replicates, subset = NULL, estimationType = c("quantile", "total", "edgeR"), quantile = 0.75, ...)
  {
    if(!missing(cD)) {
      if(inherits(cD, what = "pairedData"))
        {
          data <- cbind(cD@data, cD@pairData)
          replicates <- paste(as.character(rep(cD@replicates, 2)), rep(c("a", "b"), each = ncol(cD)), sep = "")
        } else if(inherits(cD, what = "alignmentData")) {
          libSet <- !duplicated(as.character(values(cD@alignments)$tag))
          data <- sapply(1:ncol(cD), function(ii) as.integer(cD@data[libSet,ii]))
          replicates <- cD@replicates
        } else {
          data <- cD@data
          replicates <- cD@replicates
        }
      
    }
    
    if(missing(subset)) subset <- NULL
    if(is.null(subset)) subset <- 1:nrow(data)      
    estimationType = match.arg(estimationType)
    if(is.na(estimationType)) stop("'estimationType' not known")

    estLibs <- function(data, replicates)
      {
        libsizes <- switch(estimationType,
                           total = colSums(data[subset,,drop = FALSE], na.rm = TRUE),
                           quantile = apply(data[subset,, drop = FALSE], 2, function(z) {
                             x <- z[z > 0]
                             sum(x[x <= quantile(x, quantile, na.rm = TRUE)], na.rm = TRUE)
                               }),
                           
                           edgeR = {                             
                             if(!("edgeR" %in% loadedNamespaces()))
                                 requireNamespace("edgeR", quietly = TRUE)
                             d <- DGEList(counts = data[subset,, drop = FALSE], group = replicates, lib.size = colSums(data, na.rm = TRUE))
                             d <- calcNormFactors(d, ...)
                             d$samples$norm.factors * d$samples$lib.size
                           })
        names(libsizes) <- colnames(data)
        libsizes
      }

    if(length(dim(data)) == 2) estLibsizes <- estLibs(data, replicates)
    if(length(dim(data)) == 3) {
      combData <- do.call("cbind", lapply(1:dim(data)[3], function(kk) data[,,kk]))
      combReps <- paste(as.character(rep(replicates, dim(data)[3])), rep(c("a", "b"), each = ncol(data)), sep = "")
      estLibsizes <- estLibs(combData, combReps)
      estLibsizes <- do.call("cbind",
                             split(estLibsizes, cut(1:length(estLibsizes), breaks = dim(data)[3], labels = FALSE)))
    }
                            
    if(!missing(cD))
      if(inherits(cD, what = "pairedData")) return(list(estLibsizes[1:ncol(cD)], estLibsizes[1:ncol(cD) + ncol(cD)]))

    if(length(dim(data)) > 2) estLibsizes <- array(estLibsizes, dim = dim(cD@data)[-1])
    
    return(estLibsizes)
  }
