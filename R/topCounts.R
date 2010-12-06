`topCounts` <-
function(cD, group, decreasing = TRUE, number = 10, normaliseData = FALSE)
  {
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")
    if(nrow(cD@posteriors) == 0)
      stop("The '@posteriors' slot of cD is empty!")

    if(normaliseData)
      data <- round(t(t(cD@data) / cD@libsizes) * (prod(cD@libsizes))^(1/length(cD@libsizes))) else data <- cD@data
    
    if(is.null(group))
      {
        if(length(cD@nullPosts) == 0)
          stop("The '@nullPosts' slot of cD is empty - you can't use 'group = NULL'.")
        selTags <- order(cD@nullPosts, decreasing = decreasing)[1:number]
        topTags <- data.frame(cD@annotation[selTags,, drop = FALSE], data[selTags,,drop = FALSE], Likelihood = exp(cD@nullPosts[selTags]), FDR = cumsum(1 - exp(cD@nullPosts[selTags])) / 1:number)
      } else
    {
      selTags <- order(cD@posteriors[,group], decreasing = decreasing)[1:number]
      topTags <- data.frame(cD@annotation[selTags,, drop = FALSE], data[selTags,,drop = FALSE], Likelihood = exp(cD@posteriors[selTags, group]), FDR = cumsum(1 - exp(cD@posteriors[selTags, group])) / 1:number)
    }
    rownames(topTags) <- rownames(cD@data)[selTags]
    topTags
  }

