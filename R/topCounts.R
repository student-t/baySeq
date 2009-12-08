`topCounts` <-
function(cD, group, decreasing = TRUE, number = 10)
  {
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")
    if(nrow(cD@posteriors) == 0)
      stop("The '@posteriors' slot of cD is empty!")
    if(is.null(group))
      {
        if(length(cD@nullPosts) == 0)
          stop("The '@nullPosts' slot of cD is empty - you can't use 'group = NULL'.")
        selTags <- order(cD@nullPosts, decreasing = decreasing)[1:number]
        topTags <- data.frame(cD@annotation[selTags,, drop = FALSE], cD@data[selTags,], logP = cD@nullPosts[selTags])
      } else
    {
      selTags <- order(cD@posteriors[,group], decreasing = decreasing)[1:number]
      topTags <- data.frame(cD@annotation[selTags,, drop = FALSE], cD@data[selTags,], logP = cD@posteriors[selTags, group])
    }
    rownames(topTags) <- rownames(cD@data)[selTags]
    topTags
  }

