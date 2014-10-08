getTPs <- function(cD, group, decreasing = TRUE, TPs)
  {
    if(is.character(group))
      group <- pmatch(group, names(cD@groups))
    
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")
    if(nrow(cD@posteriors) == 0)
      stop("The '@posteriors' slot of cD is empty!")

    nulls <- cD@nullPosts
    
    if(is.null(group))
      {
        if(length(cD@nullPosts) == 0)
          stop("The '@nullPosts' slot of cD is empty - you can't use 'group = NULL'.")
        posteriors <- cD@nullPosts
        pord <- order(posteriors, -apply(cD@posteriors, 1, .logsum), decreasing = decreasing)
      } else {
        if(length(nulls) == 0 || nrow(nulls) == 0) nulls <- matrix(ncol = 0, nrow = nrow(cD))
        posteriors <- cD@posteriors[,group]
        pord <- order(posteriors, -.logRowSum(cbind(cD@posteriors[,-group, drop = FALSE], nulls)), decreasing = decreasing)
      }
    cumsum(pord[1:length(posteriors)] %in% TPs)
  }



