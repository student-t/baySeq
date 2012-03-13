`topCounts` <-
function(cD, group, decreasing = TRUE, number = 10, likelihood, FDR, normaliseData = FALSE)
  {
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")
    if(nrow(cD@posteriors) == 0)
      stop("The '@posteriors' slot of cD is empty!")

    if(normaliseData)
      data <- round(t(t(cD@data) / cD@libsizes) * (prod(cD@libsizes))^(1/length(cD@libsizes))) else data <- cD@data

    if(is.character(group))
      group <- pmatch(group, names(cD@groups))

    if(nrow(cD@annotation) == 0) annotation <- data.frame(rowID = paste("row" , 1:nrow(cD), sep = "_")) else annotation <- cD@annotation
    
    if(class(cD) == "lociData") annotation <- cbind(data.frame(chr = as.character(seqnames(cD@coordinates)), start = as.numeric(start(cD@coordinates)), end = as.numeric(end(cD@coordinates))), annotation) else annotation <- annotation

    if(is.null(group)) {
      if(length(cD@nullPosts) == 0)
        stop("The '@nullPosts' slot of cD is empty - you can't use 'group = NULL'.")
      likes <- cD@nullPosts        
    } else likes <- cD@posteriors[,group]
    
    if(!missing(likelihood)) cutNumber <- sum(likes > log(likelihood))      
    if(missing(likelihood) & !missing(FDR)) cutNumber <- sum(cumsum(1 - exp(sort(likes, decreasing = decreasing))) / 1:length(likes) < FDR)

    if(!missing(likelihood) | !missing(FDR))
      if(cutNumber == 0) warning("No features were found using the cutoffs for likelihood or FDR specified; using the 'number' argument instead") else number <- cutNumber
      
    
    selTags <- order(likes, decreasing = decreasing)[1:number]
    topTags <- data.frame(annotation[selTags,, drop = FALSE], data[selTags,,drop = FALSE], Likelihood = exp(likes[selTags]), FDR = cumsum(1 - exp(likes[selTags])) / 1:number)

    rownames(topTags) <- rownames(cD@data)[selTags]
    topTags
  }

