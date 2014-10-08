densityFunctions <- function() {
  denfuncs <- Filter( function(x) 'densityFunction' %in% class( get(x) ), ls("package:baySeq") )
  cbind(densityFunction = denfuncs, description = sapply(denfuncs, function(df) eval(parse(text = paste(df, "@description", sep = "")))))
}

mdDensity <- new("densityFunction",
                 description = "A density function based on the multinomial-Dirichlet distribution, designed to assess genomic events that generate multiple counts (e.g. RNA-Seq from multiple tissue types taken from the same organism).",
                 density = .multiDirichletFunction,
                 initiatingValues = c(0.01, NA),
                 equalOverReplicates = c(TRUE, FALSE),
                 lower = function(data) 0,
                 upper = function(data) 1,
                 stratifyFunction = function(data) rowMeans(data[,,1] / rowSums(data), na.rm = TRUE),
                 stratifyBreaks = 10,
                 nullFunction = function(pars) apply(abs(pars[,-1] - 1 / dim(pars)), 1, max))



bbDensity <- new("densityFunction",
                 description = "A density function based on the betaBinomial distribution, designed to assess genomic events that generate paired data",
                 density = .betaBinomialFunction,
                 initiatingValues = c(function(dat, observables) {
                   if(sum(dat[,,1] + dat[,,2]) == 0) return(0.5)
                   if(sum(dat[,,1]) == 0) return(0.001)
                   if(sum(dat[,,2]) == 0) return(1 - 0.001)
                   mean(dat[,,1] / (dat[,,1] + dat[,,2]), na.rm = TRUE)
                 }, 0.01),
                 equalOverReplicates = c(FALSE, TRUE),
                 lower = function(data) 0,
                 upper = function(data) 1,
                 stratifyFunction = function(data) rowMeans(data[,,1] / (data[,,1] + data[,,2]), na.rm = TRUE),
                 stratifyBreaks = 10,
                 nullFunction = function(pars) abs(0.5 - pars[,1]),
#                 nullQuantiles = c(0,1),
                 orderingFunction = function(dat, observables) {
                   data <- dat / observables$libsizes               
                   adjmin <- min(data[data > 0]) / 10
                   log(.sliceArray(list(NULL, NULL, 1), data + adjmin, drop = TRUE)) -
                     log(.sliceArray(list(NULL, NULL, 2), data + adjmin, drop = TRUE))                   
                 }
                 )
                 
normDensity <- new("densityFunction",
                   description = "A density function based on the normal distribution.",
                   density = .normDensityFunction,
                   initiatingValues = c(function(dat, observables) mean(dat) / observables$libsizes, 1),
                   equalOverReplicates = c(FALSE, TRUE),
                   lower = function(dat) 0, upper = function(dat) 1 + max(dat) * 2,
                   stratifyFunction = rowMeans, stratifyBreaks = 10,
                   nullFunction = function(pars) pars[,1])


ZINBDensity <- new("densityFunction",
                   description = "A density function based on the zero-inflated negative-binomial distribution.",
                   density = .dZINB,
                   initiatingValues = c(function(dat, observables)
                     max((0.1 / observables$libsizes / observables$seglens), mean(dat / observables$libsizes / observables$seglens)),                         
                     0.1,
                     function(dat, observables) min(0.9, max(0.1, sum(dat == 0)/length(dat)))),
                   equalOverReplicates = c(FALSE, TRUE, TRUE),
                   lower = function(dat) 0, upper = function(dat) 1 + max(dat) * 2,
                   stratifyFunction = rowMeans, stratifyBreaks = 10)


nbinomDensity <- new("densityFunction",
                     description = "A density function based on the negative-binomial distribution.",
                     density = .nbinomDens,
                     initiatingValues = c(function(dat, observables) mean(pmax(dat,1) / observables$libsizes / observables$seglens), 0.01),
                     equalOverReplicates = c(FALSE, TRUE),
                     lower = function(dat) 0, upper = function(dat) 1 + max(dat) * 2,
                     stratifyFunction = rowMeans, stratifyBreaks = 10,
                     nullFunction = function(pars) log(pars[,1]),
                     orderingFunction = function(dat, observables) dat / observables$libsizes / observables$seglens)




bbNCDist <- new("densityFunction",
                description = "A density function based on the beta-binomial distribution, intended for analysis of methylation data with an observed non-conversion rate in each sample.",
                density = .betaBinomialNCFunction,
                initiatingValues = c(function(dat, observables) {
                  if(sum(dat[,,1] + dat[,,2]) == 0) return(0.5)
                  if(sum(dat[,,1]) == 0) return(0.001)
                  if(sum(dat[,,2]) == 0) return(1 - 0.001)
                  mean(dat[,,1] / (dat[,,1] + dat[,,2]), na.rm = TRUE)
                }, 0.01),
                equalOverReplicates = c(FALSE, TRUE),
                lower = function(data) 0,
                upper = function(data) 1,
                stratifyFunction = function(data) rowMeans(data[,,1] / (data[,,1] + data[,,2]), na.rm = TRUE),
                stratifyBreaks = 10,
                nullFunction = function(pars) abs(0.5 - pars[,1]))



methObservables <- function(mD) {
  cdat <- mD@data[,,1]
  upplim <- matrix(qbinom(0.005, as.vector(cdat), rep(mD@sampleObservables$nonconversion, each = nrow(mD)), lower.tail = FALSE), nrow = nrow(mD), ncol = ncol(mD))
  lowlim <- matrix(qbinom(0.005, as.vector(cdat), rep(mD@sampleObservables$nonconversion, each = nrow(mD)), lower.tail = TRUE), nrow = nrow(mD), ncol = ncol(mD))
  
  widlim <- upplim - lowlim + 1
  z <- split(dbinom(sequence(as.vector(widlim)) - 1 + rep(lowlim, widlim),
                    size = rep(as.vector(mD@data[,,1]), widlim),
                    prob = rep(rep(mD@sampleObservables$nonconversion, each = nrow(mD)), widlim), log = TRUE),
             rep(1:length(upplim), widlim))
  ncll <- array(z, dim = dim(cdat))
  ncseq <- array(split(sequence(as.vector(widlim)) - 1, rep(1:length(upplim), widlim)), dim = dim(cdat))
  ncmin <- array(lowlim, dim = dim(cdat))
  ncrange <- array(widlim, dim = dim(cdat))
  
  lchoo <- array(split(
                   lchoose(n = rep(as.vector(cdat) + as.vector(mD@data[,,2]), widlim), k = rep(as.vector(cdat), widlim) - (sequence(as.vector(widlim)) - 1 + rep(lowlim, widlim)))
                   , rep(1:length(upplim), widlim)), dim = dim(cdat))
  
  modls <- array(sapply(split(unlist(ncll) + unlist(lchoo), rep(1:length(upplim), widlim)), max), dim = dim(cdat))
  
  mD@cellObservables$ncll <- ncll
  mD@cellObservables$ncseq <- ncseq
  mD@cellObservables$ncrange <- ncrange
  mD@cellObservables$lchoose <- lchoo
  mD@cellObservables$ncmin <- ncmin
  mD@cellObservables$modls <- modls
  mD
}
