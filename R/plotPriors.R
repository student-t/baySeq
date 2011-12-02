`plotPriors` <-
function(cD, group)
  {
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")

    if(is.character(group))
        group <- pmatch(group, names(cD@groups))

    if(cD@priorType != "NB") stop("Incorrect prior type for this method of likelihood estimation")

    ungroup <- 1:length(levels(cD@groups[[group]]))
    
    par(mfrow = c(1, length(ungroup)))
    for(ii in ungroup)
      plot(density(log(cD@priors$priors[[group]][[ii]][cD@priors$sampled[,2],1]), weights = cD@priors$weights / sum(cD@priors$weights)), main = paste("Log prior means for group:", group))
  }

