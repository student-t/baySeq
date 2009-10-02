getTPs <- function(cDP, group, decreasing = TRUE, TPs)
  {
    countTPs <- function(selnum, ordering, TPs)
      sum(ordering[1:selnum] %in% TPs)
    posteriors <- cDP@posteriors[,group]
    pord <- order(posteriors, decreasing = decreasing)
    sapply(1:length(posteriors), countTPs, pord, TPs)
  }



