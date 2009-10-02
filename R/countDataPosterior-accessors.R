
setMethod("[", "countDataPosterior", function(x, i, j, ..., drop = FALSE) {
  x <- callNextMethod(x, i = i, j = j, ..., drop = drop)
  x@posteriors <- x@posteriors[i,]
  x
})

setMethod("dim", "countDataPosterior", function(x) {
  x <- callNextMethod(x)
  x
})

setMethod("show", "countDataPosterior", function(object) {
  x <- callNextMethod(object)
  cat('Slot "posteriors":\n')
  if(nrow(object) > 10)
    {
      print(object@posteriors[1:10,])
      cat(paste(nrow(object) - 10), "more rows...\n")
    } else print(object@posteriors)
  cat('\nSlot "estProps":\n')
  print(object@estProps)
})
