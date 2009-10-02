`topCounts` <-
function(cDP, group, decreasing = TRUE, number = 10)
  {
    if(class(cDP) != "countDataPosterior")
      stop("variable 'cDP' must be of class 'countDataPosterior'")
    selTags <- order(cDP@posteriors[,group], decreasing = decreasing)[1:number]
    topTags <- data.frame(cDP@annotation[selTags,, drop = FALSE], cDP@data[selTags,], logP = cDP@posteriors[selTags, group])
    rownames(topTags) <- rownames(cDP@data)[selTags]
    #rownames(topTags) <- selTags
    topTags
  }

