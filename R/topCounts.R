selectTop <- function(cD, group, ordering, orderings = TRUE, decreasing = TRUE, number = 10, likelihood, FDR, FWER) {
  if(missing(likelihood)) likelihood <- NULL
  if(missing(FDR)) FDR <- NULL
  if(missing(FWER)) FWER <- NULL
  if(missing(ordering)) ordering <- NULL

  if(!missing(group)) {
    st <- .selectTags(cD, group = group, ordering = ordering, decreasing = decreasing, number = number, likelihood = likelihood, FDR = FDR, FWER = FWER)
    if(!is.null(st)) return(cD[st,]) else return(NULL)
  }

  if(length(cD@nullPosts) > 0) {
    st <- .selectTags(cD, NULL, ordering = NULL, decreasing = decreasing, number = number, likelihood, FDR, FWER)
    if(!is.null(st)) nullCD <- cD[st,] else nullCD <- NULL
  } else nullCD <- NULL
  
  
  
  if(!orderings) {
    selOrd <- lapply(1:length(cD@groups), function(ii) cD[.selectTags(cD, ii, ordering = NULL, decreasing = decreasing, number = number, likelihood, FDR, FWER),])                     
    names(selOrd) <- names(cD@groups)
    } else {
      selOrd <- do.call("c", lapply(1:ncol(cD@orderings),function(ii) {
        selord <- lapply(levels(cD@orderings[,ii]), function(ord) {
          st <- .selectTags(cD, ii, ordering = ord, decreasing = decreasing, number = number, likelihood, FDR, FWER)
          if(!is.null(st)) cD[st,] else NULL
        })
        names(selord) <- paste(colnames(cD@orderings)[ii], ":", levels(cD@orderings[,ii]), sep = "")
        names(selord) <- gsub(":$", "", names(selord))
        selord
      }))
    }
  if(length(cD@nullPosts) > 0) return(c(list(null = nullCD), selOrd)) else return(selOrd)
}
  
#.selectTop <- function(cD, group, ordering, decreasing = TRUE, number = 10, likelihood, FDR, FWER) #{
#  if(missing(likelihood)) likelihood <- NULL
#  if(missing(FDR)) FDR <- NULL
#  if(missing(FWER)) FWER <- NULL
#  if(missing(ordering)) ordering <- NULL  
  
#  selTags <- .selectTags(cD, group, ordering, decreasing = decreasing, number = number, likelihood, FDR, FWER)
#  cD[selTags,]
#}

.selectTags <- function(cD, group, ordering, decreasing = TRUE, number = 10, likelihood, FDR, FWER) {
  if(!inherits(cD, what = "countData"))
    stop("variable 'cD' must be of or descend from class 'countData'")
  if(nrow(cD@posteriors) == 0)
    stop("The '@posteriors' slot of cD is empty!")
  
  if(is.character(group))
    group <- pmatch(group, names(cD@groups))
  if(!is.null(group) && is.na(group)) stop("Can't match this group name.")
  
  if(!is.null(ordering)) {
    ordCD <- which(cD@orderings[,group] == ordering)
    cD <- cD[ordCD,]
  }
  
  if(is.null(group)) {
    if(length(cD@nullPosts) == 0)
      stop("The '@nullPosts' slot of cD is empty - you can't use 'group = NULL'.")
    likes <- cD@nullPosts
    neglikes <- cD@posteriors
  } else {
    likes <- cD@posteriors[,group,drop = FALSE]    
    if(length(cD@nullPosts) > 0) neglikes <- cbind(cD@posteriors[,-group,drop=FALSE], cD@nullPosts) else neglikes <- cD@posteriors[,-group, drop = FALSE]
  }
  
  ordgroups <- order(.logRowSum(neglikes), decreasing = !decreasing)

  cutNumber <- c()
  if(!is.null(likelihood))
      cutNumber <- c(cutNumber, sum(likes > log(likelihood), na.rm = TRUE))
  if (!is.null(FDR))
      cutNumber <- c(cutNumber, sum(cumsum(1 - exp(likes[ordgroups,1]))/ 1:sum(!is.na(likes[,1])) < FDR, na.rm = TRUE))
  if (!is.null(FWER))
      cutNumber <- c(cutNumber, sum(1 - cumprod(exp(likes[ordgroups,1])) < FWER, na.rm = TRUE))
  
  if(!is.null(likelihood) | !is.null(FDR) | !is.null(FWER)) {
      number <- min(cutNumber)
      if(cutNumber == 0) warning("No features were found using the cutoffs for likelihood, FDR or FWER specified")
   }
  number <- min(number, nrow(likes))

  if(number > 0) {
    selTags <- ordgroups[1:number]
    if(!is.null(ordering)) selTags <- ordCD[selTags]
    return(selTags)
  } else return(NULL)
}
  
`topCounts` <-
function(cD, group, ordering, decreasing = TRUE, number = 10, likelihood, FDR, FWER, normaliseData = FALSE)
  {
    if(missing(likelihood)) likelihood <- NULL
    if(missing(FDR)) FDR <- NULL
    if(missing(FWER)) FWER <- NULL
    if(missing(ordering)) ordering <- NULL

    if(is.character(group))
      group <- pmatch(group, names(cD@groups))
    if(!is.null(group) && is.na(group)) stop("Can't match this group name.")
    
    selTags <- .selectTags(cD, group, ordering, decreasing = decreasing, number = number, likelihood, FDR, FWER)
    

    if(!is.null(group)) likes <- cD@posteriors[selTags,group] else likes <- cD@nullPosts[selTags,1]

    if(all(is.null(selTags))) return(NULL)
    
    selData <- .sliceArray(list(selTags), cD@data)
    if(normaliseData) {
        observables <- .catObservables(cD[selTags,])
        
        meanLibs <- array(apply(
            matrix(apply(observables$libsizes, setdiff(1:length(dim(cD@data)), 2), function(x) exp(mean(log(x)))), nrow = nrow(observables$libsizes))
          , 2, function(x) matrix(x, nrow = length(x), ncol = dim(cD@data)[2])), dim(observables$libsizes))
        selData <- round(selData / observables$libsizes * meanLibs / observables$seglens * exp(mean(log(observables$seglens))))
      }
    showData <- .showData(selData)
    colnames(showData) <- colnames(cD@data)
    
    if(length(cD@orderings) > 0 && !is.null(group)) ordering <- cD@orderings[selTags, group, drop = TRUE] else ordering <- rep("", length(selTags))
    if(all(ordering == "")) noorder <- TRUE else noorder <- FALSE
    
    
    if(nrow(cD@annotation) == 0) annotation <- data.frame(rowID = selTags) else annotation <- cD@annotation[selTags,]    
    
    if(inherits(cD, what = "lociData") | inherits(cD, what = "methData"))
      annotation <- cbind(annotation, GenomicRanges::as.data.frame(cD@coordinates[selTags])) else annotation <- annotation    

    
    
    topTags <- data.frame(annotation, showData, Likelihood = exp(likes),
                          ordering = ordering,
                          FDR = cumsum(1 - exp(likes)) / 1:length(selTags),
                          FWER = 1 - cumprod(exp(likes)))
    names(topTags)[names(topTags) == "FDR"] <- paste("FDR", names(cD@groups)[group[1]], sep = ".")
    names(topTags)[names(topTags) == "FWER"] <- paste("FWER", names(cD@groups)[group[1]], sep = ".")
    if(noorder) topTags <- topTags[,-which(colnames(topTags) == "ordering")]
    rownames(topTags) <- rownames(cD@data)[selTags]
    topTags
  }
