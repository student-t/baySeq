setClass("countData", representation(data = "matrix", replicates = "factor", libsizes = "numeric", groups = "list", annotation = "data.frame", priorType = "character", priors = "list", posteriors = "matrix", nullPosts = "numeric", estProps = "numeric", seglens = "matrix"))

setClass("pairedData", representation(pairData = "matrix", pairLibsizes = "numeric"), contains = "countData")
