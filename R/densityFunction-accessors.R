setMethod("initialize", "densityFunction", function(.Object, ...) {
  dotlist <- list(...)
  if("initiatingValues" %in% names(dotlist))
    if(is.numeric(dotlist$initiatingValues))
      dotlist$initiatingValues <- (eval(parse(text = paste("function(dat, observables)", "c(", paste(dotlist$initiatingValues, collapse = ", "), ")"))))
  if("equalOverReplicates" %in% names(dotlist))
    if(is.logical(dotlist$equalOverReplicates))
      dotlist$equalOverReplicates <- (eval(parse(text = paste("function(dat, observables)", "c(", paste(dotlist$equalOverReplicates, collapse = ", "), ")"))))
  
  .Object <- do.call("callNextMethod", c(list(.Object), dotlist))
  .Object
})
