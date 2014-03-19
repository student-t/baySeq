setClass("countData", representation(data = "matrix", replicates = "factor", libsizes = "numeric", groups = "list", annotation = "data.frame", priorType = "character", priors = "list", posteriors = "matrix", nullPosts = "matrix", estProps = "numeric", seglens = "matrix", orderings = "data.frame"))

setClass("pairedData", representation(pairData = "matrix", pairLibsizes = "numeric"), contains = "countData")
