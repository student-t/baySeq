`plotPriors` <-
function(cD, group)
  {
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")

    if(cD@priorType != "NB") stop("Incorrect prior type for this method of likelihood estimation")

    ungroup <- unique(cD@groups[[group]])
    
    par(mfrow = c(1, length(ungroup)))
    for(ii in ungroup)
      plot(density(rep(log(cD@priors$priors[[group]][[1]]), each = cD@priors$copies)), main = paste("Log prior means for group:", group))
  }

