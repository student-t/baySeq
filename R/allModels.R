
allModels <- function(CD) {
  replevs <- 1:length(levels(CD@replicates))
  
  x <- matrix(1, nrow = 1, ncol = 1)
  while(max(x) != max(replevs))
    {
      maxx <- apply(x, 1, max) + 1
      x <- cbind(x[rep(1:nrow(x), maxx),,drop = FALSE], sequence(maxx))
    }

  lenreps <- sapply(levels(CD@replicates), function(rep) sum(CD@replicates == rep))
  repmat <- matrix(rep(as.vector(t(x)), rep(lenreps, nrow(x))), nrow = nrow(x), byrow = TRUE)
  repmat[,unlist(lapply(levels(CD@replicates), function(rep) which(CD@replicates == rep)))] <- repmat
  
  groups(CD) <- split(repmat, 1:nrow(repmat))

  names(CD@groups) <-
    sapply(CD@groups, function(grp) {
    groupDesc <- lapply(levels(grp), function(gg) {
      grprep <- unique(CD@replicates[grp == gg])
      paste(grprep, collapse = ",")
    })
    groupDesc <- lapply(groupDesc, function(x) paste("{", x, "}", sep = ""))
    paste(groupDesc, collapse = ",")
  })
  
  CD
}
