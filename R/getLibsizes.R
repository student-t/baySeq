'getLibsizes' <- function(cD, data, replicates, subset = NULL, estimationType = c("quantile", "total", "edgeR"), quantile = 0.75, ...)
  {
    if(!missing(cD)) {
      data <- cD@data
      replicates <- cD@replicates
    }
    
    if(missing(subset)) subset <- NULL
    if(is.null(subset)) subset <- 1:nrow(data)
      
    
    estimationType = match.arg(estimationType)
    if(is.na(estimationType)) stop("'estimationType' not known")
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
    libsizes <- libsizes
    names(libsizes) <- colnames(data)
    libsizes
  }

