setClass("priorData", representation(type = "character", priors = "list", sampled = "numeric"))

setClass("countData", representation(data = "matrix", libsizes = "numeric", groups = "list", annotation = "data.frame", priors = "priorData", posteriors = "matrix", nullPosts = "numeric", estProps = "numeric"))

setClass("segData", contains = "countData", representation(seglens = "matrix"))
