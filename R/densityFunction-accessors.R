setMethod("initialize", "densityFunction", function(.Object, ...) {
  dotlist <- list(...)
  if("initiatingValues" %in% names(dotlist)) dotlist$initiatingValues <- as.list(dotlist$initiatingValues)
  
  .Object <- do.call("callNextMethod", c(list(.Object), dotlist))
  .Object@initiatingValues <-
    lapply(.Object@initiatingValues, function(x) {
      if(is.function(x)) return(x)
      if(is.numeric(x)) return(eval(parse(text = paste("function(dat, observables)", x))))
      stop("All initiatingValues must be either numeric or a function generating a numeric.")
    })
  .Object
})
