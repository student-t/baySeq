setClass("countData", representation(data = "matrix", libsizes = "numeric", groups = "list", annotation = "data.frame", priors = "list"))
setClass("countDataPosterior", contains = "countData", representation(posteriors = "matrix", estProps = "numeric"))
