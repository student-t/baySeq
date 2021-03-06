setClass("densityFunction", representation(description = "character", density = "function", initiatingValues = "function", equalOverReplicates = "function", lower = "function", upper = "function", stratifyFunction = "function", stratifyBreaks = "numeric", nullFunction = "function", orderingFunction = "function", modifyNullPriors = "function"))

setClass("countData", representation(data = "array", replicates = "factor", groups = "list", rowObservables = "list", sampleObservables = "list", cellObservables = "list", annotation = "data.frame", priorModels = "list", priorType = "character", densityFunction = "list", priors = "list", posteriors = "matrix", nullPosts = "matrix", estProps = "numeric", orderings = "data.frame"))

#seglens = "array", 
#libsizes = "array"

#setClass("pairedData", representation(pairData = "matrix", pairLibsizes = "numeric"), contains = "countData")
