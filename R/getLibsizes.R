'getLibsizes' <- function(cD, data, replicates, subset = NULL, estimationType = c("quantile", "total", "edgeR"), quantile = 0.75, ...)
  {
    if(!missing(cD)) {
      data <- cD@data
      replicates <- cD@replicates
      returnCD = TRUE
    } else returnCD = FALSE
    
    if(missing(subset)) subset <- NULL
    if(is.null(subset)) subset <- 1:nrow(data)
      
    
    estimationType = match.arg(estimationType)
    if(is.na(estimationType)) stop("'estimationType' not known")
    estLibs <- function(data)
      {
        libsizes <- switch(estimationType,
                           total = colSums(data[subset,,drop = FALSE], na.rm = TRUE),
                           quantile = apply(data[subset,, drop = FALSE], 2, function(z) {
                             x <- z[z > 0]
                             sum(x[x <= quantile(x, quantile, na.rm = TRUE)], na.rm = TRUE) }),
                           edgeR = {
                             if(!("edgeR" %in% loadedNamespaces()))
                               library(edgeR)
                             d <- DGEList(counts = data[subset,, drop = FALSE], group = replicates)
                             d <- calcNormFactors(d, ...)
                             d$samples$norm.factors
                           })
        names(libsizes) <- colnames(data)
        libsizes
      }
    libsizes <- estLibs(data)
    
    if(returnCD)
      {
        cD@libsizes <- libsizes
        if(inherits(cD, what = "pairedData")) cD@pairLibsizes <- estLibs(cD@pairData)
        return(cD)
      } else return(libsizes)
  }

