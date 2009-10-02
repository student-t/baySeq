`getLikelihoods.Dirichlet` <-
function(cDP, prs, estimatePriors = TRUE, subset = NULL, cl)
  {
    if(class(cDP) != "countData")
      stop("variable 'cDP' must be of class 'countData'")
    if(cDP@priors$type != "Dir") stop("Incorrect prior type for this method of likelihood estimation")
    if(class(cDP) != "countData")
      stop("variable 'countData' in 'getGroupPriors' must be of class 'countData'")
    if(length(prs) != length(cDP@groups)) stop("'prs' must be of same length as the number of groups in the 'cDP' object")
    if(!(class(subset) == "integer" | class(subset) == "numeric" | is.null(subset)))
      stop("'subset' must be integer, numeric, or NULL")

    if(is.null(subset)) subset <- 1:nrow(cDP@data)
    
    `PrgivenD.Dirichlet` <-
      function(us, libsizes, priors, groups, priorgroups = c(0.99, 0.01))
        {
          us <- cbind(us, libsizes - us)
          ps <- c()
          for(ii in 1:length(groups))
            ps[ii] <- PDgivenr.Dirichlet(us, priors[[ii]], groups[[ii]])
          
          ps
        }
    
    if(is.null(cl)) {
      ps <- apply(cDP@data[subset,], 1, PrgivenD.Dirichlet, cDP@libsizes, cDP@priors$priors, cDP@groups)
    } else {
      ps <- parRapply(cl, cDP@data[subset,], PrgivenD.Dirichlet, cDP@libsizes, cDP@priors$priors, cDP@groups)
    }
    ps <- matrix(ps, ncol = length(cDP@groups), byrow = TRUE)
    pps <- getPosteriors(ps, prs, estimatePriors = estimatePriors)

    posteriors <- matrix(NA, ncol = length(cDP@groups), nrow(cDP@data))
    posteriors[subset,] <- t(pps$posteriors)

    colnames(posteriors) <- names(cDP@groups)
    estProps <- pps$priors
    names(estProps) <- names(cDP@groups)
    
    new("countDataPosterior", cDP, posteriors = posteriors, estProps = estProps)
  }


`getLikelihoods.Pois` <-
function(cDP, prs, estimatePriors = TRUE, subset = NULL, distpriors = FALSE, cl)
  {
    if(class(cDP) != "countData")
      stop("variable 'cDP' must be of class 'countData'")
    if(cDP@priors$type != "Poi") stop("Incorrect prior type for this method of likelihood estimation")
    if(length(prs) != length(cDP@groups)) stop("'prs' must be of same length as the number of groups in the 'cDP' object")

    if(!(class(subset) == "integer" | class(subset) == "numeric" | is.null(subset)))
      stop("'subset' must be integer, numeric, or NULL")

    if(is.null(subset)) subset <- 1:nrow(cDP@data)

    `PrgivenD.Pois` <-
      function(us, libsizes, priors, groups, distpriors = FALSE)
        {
          ps <- c()
          for(ii in 1:length(groups))
            if(!distpriors) ps[ii] <-
              PDgivenr.Pois(us, libsizes, priors[[ii]], groups[[ii]]) else ps[ii] <-
                PDgivenr.PoisIndie(us, libsizes, priors[[ii]], groups[[ii]])
          ps
        }

    if(is.null(cl)) {
      ps <- apply(cDP@data[subset,], 1, PrgivenD.Pois, cDP@libsizes, cDP@priors$priors, cDP@groups, distpriors)
    } else {
      ps <- parRapply(cl, cDP@data[subset,], PrgivenD.Pois, cDP@libsizes, cDP@priors$priors, cDP@groups, distpriors)
      }                        
    ps <- matrix(ps, ncol = length(cDP@groups), byrow = TRUE)
    pps <- getPosteriors(ps, prs, estimatePriors = estimatePriors)

    posteriors <- matrix(NA, ncol = length(cDP@groups), nrow(cDP@data))
    posteriors[subset,] <- t(pps$posteriors)

    colnames(posteriors) <- names(cDP@groups)
    estProps <- pps$priors
    names(estProps) <- names(cDP@groups)
    
    new("countDataPosterior", cDP, posteriors = posteriors, estProps = estProps)
  }


`getPosteriors` <-
function(ps, prs, estimatePriors = FALSE, maxit = 100, accuracy = 1e-5)
  {
    getPosts <- function(x, prs)
      {
        posts <- x + log(prs)
        posts <- posts - logsum(posts)
        posts
      }
    if(estimatePriors)
      {
        oldprs <- prs
        for(ii in 1:maxit)
          {
            posteriors <- apply(ps, 1, getPosts, prs)
            prs <- rowSums(exp(posteriors)) / ncol(posteriors)
            if(all(abs(oldprs - prs) < accuracy)) break
            oldprs <- prs
          }
        if(ii == maxit)
          warning("Convergence not achieved to required accuracy.")
      }
    list(posteriors = posteriors <- apply(ps, 1, getPosts, prs), priors = prs)
  }


`getLikelihoods.NBboot` <-
function(cDP, prs, estimatePriors = TRUE, subset = NULL, bootStraps = 2, conv = 1e-4, cl)
  {
    
    if(class(cDP) != "countData")
      stop("variable 'cDP' must be of class 'countData'")
    if(cDP@priors$type != "NB") stop("Incorrect prior type for this method of likelihood estimation")
    if(length(prs) != length(cDP@groups)) stop("'prs' must be of same length as the number of groups in the 'cDP' object")

    if(!(class(subset) == "integer" | class(subset) == "numeric" | is.null(subset)))
      stop("'subset' must be integer, numeric, or NULL")
    
    if(is.null(subset)) subset <- 1:nrow(cDP@data)
    sy <- cDP@priors$sampled
    
    NBbootStrap <- function (us, groups, libsizes) 
      {
        ps <- c()
        for (ii in 1:length(groups))
          ps[ii] <- PDgivenr.NBIndie(us, libsizes, NBpriors[[ii]], groups[[ii]], sampPriors[[ii]])
        ps
      }    
    clustAssignPrior <- function(priors, sampPriors)
      {
        assign("NBpriors", priors, envir = .GlobalEnv)
        assign("sampPriors", sampPriors, envir = .GlobalEnv)
        NULL
      }
    if(is.null(conv)) conv <- 0
    posteriors <- matrix(1 / length(cDP@groups), ncol = length(cDP@groups), nrow = nrow(cDP@data))
    propest <- NULL
   
    for(cc in 1:bootStraps)
      {
        sampPriors <- list()
        for(ii in 1:length(cDP@groups))
            sampPriors[[ii]] <- posteriors[sy,ii] / sum(posteriors[sy,ii])

        NBpriors <- cDP@priors$priors
        if(is.null(cl)) {
          ps <- 
            apply(cDP@data[union(sy, subset),,drop = FALSE], 1, NBbootStrap,
                  groups = cDP@groups, libsizes = cDP@libsizes)
        } else {
          clusterCall(cl, clustAssignPrior, NBpriors, sampPriors)
          ps <- 
            parRapply(cl, cDP@data[union(sy, subset),, drop = FALSE], NBbootStrap,
                      groups = cDP@groups, libsizes = cDP@libsizes)
        }
        ps <- matrix(ps, ncol = length(cDP@groups), byrow = TRUE)
        pps <- getPosteriors(ps, prs, estimatePriors = TRUE)

        posteriors <- matrix(NA, ncol = length(cDP@groups), nrow = nrow(cDP@data))
        posteriors[union(sy, subset),] <- t(pps$posteriors)
        
        posteriors <- exp(posteriors)
        prs <- pps$priors
        propest <- rbind(propest, prs)
        print(prs)
    
    
        retposteriors <- matrix(NA, ncol = length(cDP@groups), nrow = nrow(cDP@data))
        retposteriors[union(sy, subset),] <- t(pps$posteriors)
        if(any(!(sy %in% subset)))
          retposteriors[sy[!(sy %in% subset)],] <- NA
        
        colnames(retposteriors) <- names(cDP@groups)
        estProps <- pps$priors
        names(estProps) <- names(cDP@groups)

        if(all(abs((posteriors[union(sy,subset),] - exp(t(pps$posteriors)))) < conv) | cc > bootStraps) break()
      }

    return(new("countDataPosterior", cDP, posteriors = retposteriors, estProps = estProps))
  }


`getPosteriors` <-
function(ps, prs, estimatePriors = FALSE, maxit = 100, accuracy = 1e-5)
  {
    getPosts <- function(x, prs)
      {
        posts <- x + log(prs)
        posts <- posts - logsum(posts)
        posts
      }
    if(estimatePriors)
      {
        oldprs <- prs
        for(ii in 1:maxit)
          {
            posteriors <- apply(ps, 1, getPosts, prs)
            prs <- rowSums(exp(posteriors)) / ncol(posteriors)
            if(all(abs(oldprs - prs) < accuracy)) break
            oldprs <- prs
          }
        if(ii == maxit)
          warning("Convergence not achieved to required accuracy.")
      }
    list(posteriors = posteriors <- apply(ps, 1, getPosts, prs), priors = prs)
  }

