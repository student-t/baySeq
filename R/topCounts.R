`topCounts` <-
function(cD, group, decreasing = TRUE, number = 10, likelihood, FDR, normaliseData = FALSE)
  {
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")
    if(nrow(cD@posteriors) == 0)
      stop("The '@posteriors' slot of cD is empty!")

    if(is.character(group))
      group <- pmatch(group, names(cD@groups))
    if(!is.null(group) && is.na(group)) stop("Can't match this group name.")
    
    if(is.null(group)) {
      if(length(cD@nullPosts) == 0)
        stop("The '@nullPosts' slot of cD is empty - you can't use 'group = NULL'.")
      likes <- cD@nullPosts        
    } else likes <- cD@posteriors[,group,drop = FALSE]    

    if(!missing(likelihood)) cutNumber <- sum(likes > log(likelihood), na.rm = TRUE)
    if(missing(likelihood) & !missing(FDR)) cutNumber <- sum(cumsum(1 - exp(sort(likes[,1], decreasing = decreasing))) / 1:sum(!is.na(likes[,1])) < FDR, na.rm = TRUE)
        if(!missing(likelihood) | !missing(FDR))
      if(cutNumber == 0) warning("No features were found using the cutoffs for likelihood or FDR specified; using the 'number' argument instead") else number <- cutNumber

    number <- min(number, nrow(likes))
    
    selTags <- order(likes[,1], decreasing = decreasing)[1:number]

    if(length(cD@libsizes) == 0) {
      if(normaliseData) warning("No library sizes provided in '@libsizes' slot but normalise option selected.")
      cD@libsizes <- rep(1, ncol(cD)) 
    }
    if(inherits(cD, what = "pairedData") && length(cD@pairLibsizes) == 0) {
      if(normaliseData) warning("No library sizes provided in '@pairLibsizes' slot but normalise option selected.")
      cD@pairLibsizes <- rep(1, ncol(cD))
    }
    
    if(inherits(cD, what = "pairedData"))
      {        
        if(normaliseData) {          
          data <- round(t(t(cD@data) / cD@libsizes) * exp(mean(log(c(cD@libsizes, cD@pairLibsizes)))))[selTags,]
          pairData <- round(t(t(cD@pairData) / cD@pairLibsizes) * exp(mean(log(c(cD@libsizes, cD@pairLibsizes)))))[selTags,]
        } else {
          data <- cD@data[selTags,,drop = FALSE]
          pairData <- cD@pairData[selTags,,drop = FALSE]
        }
        data <- matrix(paste(round(data), round(pairData), sep = ":"), ncol = ncol(data), nrow = nrow(data))
      } else {
        if(normaliseData) data <- round(t(t(cD@data) / cD@libsizes) * exp(mean(log(cD@libsizes))))[selTags,,drop = FALSE] else data <- cD@data[selTags,,drop = FALSE]
      }
    ordering <- cD@orderings[selTags, group, drop = FALSE]
    
    colnames(data) <- colnames(cD@data)
    
    if(nrow(cD@annotation) == 0) annotation <- data.frame(rowID = selTags) else annotation <- cD@annotation[selTags,]    
    
    if(inherits(cD, what = "lociData") | inherits(cD, what = "methData"))
      annotation <- cbind(annotation, GenomicRanges::as.data.frame(cD@coordinates[selTags])) else annotation <- annotation

    topTags <- data.frame(annotation, data, Likelihood = exp(likes[selTags,]), ordering = ordering, FDR = cumsum(1 - exp(likes[selTags,1])) / 1:number)
    names(topTags)[names(topTags) == "FDR"] <- paste("FDR", names(cD@groups)[group[1]], sep = ".")
    rownames(topTags) <- rownames(cD@data)[selTags]
    topTags
  }
