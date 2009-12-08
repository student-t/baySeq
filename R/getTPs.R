getTPs <- function(cD, group, decreasing = TRUE, TPs)
  {
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")
    if(nrow(cD@posteriors) == 0)
      stop("The '@posteriors' slot of cD is empty!")
    countTPs <- function(selnum, ordering, TPs)
      sum(ordering[1:selnum] %in% TPs)
    if(is.null(group))
      {
        if(length(cD@nullPosts) == 0)
          stop("The '@nullPosts' slot of cD is empty - you can't use 'group = NULL'.")
        posteriors <- cD@nullPosts
      } else posteriors <- cD@posteriors[,group]
    pord <- order(posteriors, decreasing = decreasing)
    sapply(1:length(posteriors), countTPs, pord, TPs)
  }



