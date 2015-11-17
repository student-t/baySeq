densityFunctions <- function() {
  denfuncs <- Filter( function(x) 'densityFunction' %in% class( get(x) ), ls("package:baySeq") )
  cbind(densityFunction = denfuncs, description = sapply(denfuncs, function(df) eval(parse(text = paste(df, "@description", sep = "")))))
}

mdDensity <- new("densityFunction",
                 description = "A density function based on the multinomial-Dirichlet distribution, designed to assess genomic events that generate multiple counts (e.g. RNA-Seq from multiple tissue types taken from the same organism).",
                 density = .multiDirichletFunction,
                 initiatingValues = function(dat, observables) {
                   dat <- dat / array(observables$libsizes, dim = c(1, dim(observables$libsizes)))
                   fudgeFactor <- 0.1/observables$dim[3] / observables$dim[2] / max(observables$libsizes)
                   props <- rep(1/observables$dim[3], observables$dim[3])
                   if(sum(dat) > 0) {
                     sumdat <- apply(dat, 3, sum)
                     scaledat <- apply(dat, 2, sum)
                     props[sumdat == 0] <- min(fudgeFactor / scaledat)
                     props[sumdat > 0] <- apply(dat[,,sumdat > 0,drop = FALSE], 3, function(x) mean(x / scaledat, na.rm = TRUE)) - min(fudgeFactor / scaledat) * sum(sumdat == 0) / sum(sumdat != 0)
                   }
                     c(0.01, props[-observables$dim[3]])
                 },
                 equalOverReplicates = function(datdim) c(TRUE, rep(FALSE, datdim[3] - 1)),
                 lower = function(data) 0,
                 upper = function(data) 1,
                 stratifyFunction = function(data) rowMeans(data[,,1] / rowSums(data), na.rm = TRUE),
                 stratifyBreaks = 10,
                 nullFunction = function(pars) rowSums(abs(cbind(pars[,-1], 1 - rowSums(pars[,-1,drop = FALSE])) - 1/ncol(pars))))


md2Density <- new("densityFunction",
                 description = "A density function based on an R[2]-multinomial-Dirichlet distribution, in which the two highest proportion values in the multinomial fit are evaluated separately, and the remaining proportions are assumed to be equal. This density function is designed to assess genomic events that generate multiple counts (e.g. RNA-Seq from multiple tissue types taken from the same organism).",
                  density = .multiDirichletFunction2,
                 initiatingValues = function(dat, observables) {
                   dat <- dat / array(observables$libsizes, dim = c(1, dim(observables$libsizes)))
                   fudgeFactor <- 0.1/observables$dim[3] / observables$dim[2] / max(observables$libsizes)
                   props <- rep(1/observables$dim[3], observables$dim[3])
                   if(sum(dat) > 0) {
                     sumdat <- apply(dat, 3, sum)
                     scaledat <- apply(dat, 2, sum)
                     props[sumdat == 0] <- min(fudgeFactor / scaledat)
                     props[sumdat > 0] <- apply(dat[,,sumdat > 0,drop = FALSE], 3, function(x) mean(x / scaledat, na.rm = TRUE)) - min(fudgeFactor / scaledat) * sum(sumdat == 0) / sum(sumdat != 0)
                   }
                   c(0.01, sort(props, decreasing = TRUE)[1:2])
                 },
                  equalOverReplicates = function(datdim) c(TRUE, FALSE, FALSE),
                  lower = function(data) 0,
                  upper = function(data) 1,
                  stratifyFunction = function(data) rowMeans(data[,,1] / rowSums(data), na.rm = TRUE),
                  stratifyBreaks = 10,
                  nullFunction = function(pars) pars[,2],
                  modifyNullPriors = function(x, datdim) {x[[1]][,2:3] <- 1/datdim[2]; x})


md3Density <- new("densityFunction",
                 description = "A density function based on an R[3]-multinomial-Dirichlet distribution, in which the two highest proportion values in the multinomial fit are evaluated separately, and the remaining proportions are assumed to be equal. This density function is designed to assess genomic events that generate multiple counts (e.g. RNA-Seq from multiple tissue types taken from the same organism).",
                  density = .multiDirichletFunction3,
                 initiatingValues = function(dat, observables) {
                   dat <- dat / array(observables$libsizes, dim = c(1, dim(observables$libsizes)))
                   fudgeFactor <- 0.1/observables$dim[3] / observables$dim[2] / max(observables$libsizes)
                   props <- rep(1/observables$dim[3], observables$dim[3])
                   if(sum(dat) > 0) {
                     sumdat <- apply(dat, 3, sum)
                     scaledat <- apply(dat, 2, sum)
                     props[sumdat == 0] <- min(fudgeFactor / scaledat)
                     props[sumdat > 0] <- apply(dat[,,sumdat > 0,drop = FALSE], 3, function(x) mean(x / scaledat, na.rm = TRUE)) - min(fudgeFactor / scaledat) * sum(sumdat == 0) / sum(sumdat != 0)
                   }
                   c(0.01, sort(props, decreasing = TRUE)[1:3])
                 },
                  equalOverReplicates = function(datdim) c(TRUE, FALSE, FALSE, FALSE),
                  lower = function(data) 0,
                  upper = function(data) 1,
                  stratifyFunction = function(data) rowMeans(data[,,1] / rowSums(data), na.rm = TRUE),
                  stratifyBreaks = 10,
                  nullFunction = function(pars) pars[,2],
                  modifyNullPriors = function(x, datdim) {x[[1]][,2:4] <- 1/datdim[2]; x})


#msdDensity <- new("densityFunction",
#                  description = "A density function based on the multinomial-*symmetric*-Dirichlet distribution, designed to assess gen#omic events that generate multiple counts (e.g. RNA-Seq from multiple tissue types taken from the same organism).",
#                  density = .multiSymDirichletFunction,
#                  initiatingValues = 1,
#                  equalOverReplicates = TRUE,
#                  lower = function(data) 0,
#                  upper = function(data) 100,
#                  stratifyFunction = function(data) rowMeans(data[,,1] / rowSums(data), na.rm = TRUE),
#                  stratifyBreaks = 10,
#                  nullFunction = function(pars) pars[,1])



bbDensity <- new("densityFunction",
                 description = "A density function based on the betaBinomial distribution, designed to assess genomic events that generate paired data",
                 density = .betaBinomialFunction,
                 initiatingValues = function(dat, observables) {                   
                   if(sum(dat[,,1] + dat[,,2]) == 0) {
                     prop <- 0.5
                   } else if(sum(dat[,,1]) == 0) {
                     prop <- 0.001
                   } else if(sum(dat[,,2]) == 0) {
                     prop <- (1 - 0.001)
                   } else prop <- mean(dat[,,1] / (dat[,,1] + dat[,,2]), na.rm = TRUE)
                   c(prop, 0.1)
                 },
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
                 },
                 modifyNullPriors = function(x, datdim) {x[[1]][,1] <- 0.5; x})
                 
normDensity <- new("densityFunction",
                   description = "A density function based on the normal distribution.",
                   density = .normDensityFunction,
                   initiatingValues = function(dat, observables)
                       c(mean(dat / observables$libsizes), max(sd(dat), 1e-4)),
                   equalOverReplicates = c(FALSE, TRUE),
                   lower = function(dat) min(dat) - diff(range(dat)) - 1,
                   upper = function(dat) max(dat) + diff(range(dat)) + 1,
                   stratifyFunction = rowMeans, stratifyBreaks = 10,
                   nullFunction = function(pars) pars[,1])


ZINBDensity <- new("densityFunction",
                   description = "A density function based on the zero-inflated negative-binomial distribution.",
                   density = .dZINB,
                   initiatingValues = function(dat, observables)
                   c(mean(pmax(dat,0.1) / observables$libsizes / observables$seglens),
                     0.1,
                     min(0.9, max(1 - 0.5^(1/observables$dim[2]), sum(dat == 0)/length(dat)))),
                   equalOverReplicates = c(FALSE, TRUE, TRUE),
                   lower = function(dat) 0, upper = function(dat) 1 + max(dat) * 2,
                   stratifyFunction = rowMeans, stratifyBreaks = 10)


nbinomDensity <- new("densityFunction",
                     description = "A density function based on the negative-binomial distribution.",
                     density = .nbinomDens,
                     initiatingValues = function(dat, observables) c(mean(pmax(dat,0.1) / observables$libsizes / observables$seglens), 0.01),
                     equalOverReplicates = c(FALSE, TRUE),
                     lower = function(dat) 0, upper = function(dat) 1 + max(dat) * 2,
                     stratifyFunction = rowMeans, stratifyBreaks = 10,
                     nullFunction = function(pars) log(pars[,1]),
                     orderingFunction = function(dat, observables) dat / observables$libsizes / observables$seglens)




bbNCDist <- new("densityFunction",
                description = "A density function based on the beta-binomial distribution, intended for analysis of methylation data with an observed non-conversion rate in each sample.",
                density = .betaBinomialNCFunction,
                initiatingValues = function(dat, observables) {                   
                  if(sum(dat[,,1] + dat[,,2]) == 0) {
                    prop <- 0.5
                  } else if(sum(dat[,,1]) == 0) {
                    prop <- 0.001
                  } else if(sum(dat[,,2]) == 0) {
                    prop <- (1 - 0.001)
                  } else prop <- mean(dat[,,1] / (dat[,,1] + dat[,,2]), na.rm = TRUE)
                  c(prop, 0.1)
                },
                equalOverReplicates = c(FALSE, TRUE),
                lower = function(data) 0,
                upper = function(data) 1,
                                        #                stratifyFunction = function(data) rowMeans(data[,,1] / (data[,,1] + data[,,2]), na.rm = TRUE),
                stratifyFunction = function(data) rowSums(data[,,1]) / rowSums(data),                
                stratifyBreaks = 10,
                nullFunction = function(pars) abs(0.5 - pars[,1]),
                orderingFunction = function(dat, observables) {
                  Cs <- dat[,,1]
                  Ts <- dat[,,2]
                  ks <- Ts * (observables$nonconversion[,,1]) / (1 - (observables$nonconversion[,,1]))
                  (Cs - ks) / (Cs + Ts)
                })


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
