setClass("countData", representation(data = "matrix", libsizes = "numeric", groups = "list", annotation = "data.frame", priors = "list", posteriors = "matrix", nullPosts = "numeric", estProps = "numeric"))

setClass("segData", contains = "countData", representation(seglens = "matrix"))
