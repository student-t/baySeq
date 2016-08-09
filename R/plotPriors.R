`plotPriors` <-
function(cD, group, par = 1)
  {        
    if(!inherits(cD, what = "countData"))
      stop("variable 'cD' must be of or descend from class 'countData'")

    if(is.character(group))
        group <- pmatch(group, names(cD@groups))

#    if(cD@priorType != "NB" && cD@priorType != "NB-QL") stop("Incorrect prior type for this method of likelihood estimation")

    ungroup <- 1:length(levels(cD@groups[[group]]))
    
    par(mfrow = c(1, length(ungroup)))
    for(ii in ungroup)
      plot(density(log(cD@priors$priors[[group]][[ii]][cD@priors$sampled[,2],par]), weights = cD@priors$sampled[,3] / sum(cD@priors$sampled[,3])), main = paste("Log prior means for group:", group))
  }

plotNullPrior <- function(cD, ...)
  {
    nF <- densityFunction(cD)[[1]]@nullFunction
    if(is.null(body(nF))) stop("nullFunction of cD object not specified.")
    
    aleq <- sapply(groups(cD), function(x) all(x == x[1]))
    if(!any(aleq)) stop("No non-differentially expressed group exists.")
    
    nZ <- nF(cD@priors$priors[[which(aleq)[1]]][[1]][cD@priors$sampled[,2],,drop = FALSE])

    denz <- density(nZ, weights = cD@priors$sampled[,3] / sum(cD@priors$sampled[,3]))
    bS <- bimodalSeparator(nZ, weights = cD@priors$sampled[,"weights"])
    
#    if(require("diptest")) submain = paste("unimodality p:", dip.test(nZ[nZ > -Inf & nZ < -Inf])$p)

    main = "Density of null function"
    if("main" %in% names(list(...))) main = list(...)$main
    
    args <- modifyList(list(x = denz, main = main), list(...))    
    
    do.call("plot", args)
    
#    plot(bS$denz, sub = submain, main = main, ...)
    abline(v = bS, col = "red")
    invisible(bS)
  }
