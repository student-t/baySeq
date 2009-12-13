`plotPriors` <-
function(cD, group)
  {
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")

    if(cD@priors$type != "NB") stop("Incorrect prior type for this method of likelihood estimation")

    ungroup <- unique(cD@groups[[group]])
    
    par(mfrow = c(1, length(ungroup)))
    for(ii in ungroup)
      plot(density(log(cD@priors$priors[[group]][[ii]][apply(cD@data[,cD@groups[[group]] == ii, drop = FALSE] != 0, 1, any),1])), main = paste("Log prior means for group:", group))
  }

