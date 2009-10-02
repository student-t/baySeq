`PDgivenr.Dirichlet` <-
function(us, prior, group)
  {
    lbeta.over.beta <- function(us, alphas)
        {
          if(all(alphas > 0))
            sum(lgamma(alphas + us)) + lgamma(sum(alphas)) - sum(lgamma(alphas)) - lgamma(sum(alphas) + sum(us)) else -Inf
        }
    
    pD <- 0
    for(ii in 1:nrow(us))
      pD <- pD +
        lfactorial(sum(us[ii,])) - sum(lfactorial(us[ii,]))
    for(gg in 1:length(unique(group)))
      pD <- pD +
        lbeta.over.beta(apply(matrix(us[group == unique(group)[gg],], ncol = ncol(us)), 2, sum),
                        prior[gg,])
    pD
  }

`PDgivenr.Pois` <-
function(us, ns, prior, group)
  {
    pD <- 0
    pD <- sum(us * log(ns)) - sum(lfactorial(us))
    for(gg in 1:length(unique(group)))
      {
        selus <- group == gg
        pD <- pD +
          lgamma(sum(us[selus]) + prior[gg,1]) - lgamma(prior[gg,1]) -
            sum(us[selus]) * log(sum(ns[selus]) + prior[gg,2]) -
              prior[gg,1] * log(1 + sum(ns[selus]) / prior[gg,2])
      }
    pD
  }

`PDgivenr.PoisIndie` <-
  function(us, ns, prior, group)
  {
    pD <- sum(us * log(ns)) - sum(lfactorial(us))
    for(gg in 1:length(unique(group)))
      {
        selus <- group == gg
        pD <- pD +
          logsum(lgamma(sum(us[selus]) + prior[[gg]][,1]) - lgamma(prior[[gg]][,1]) -
                 sum(us[selus]) * log(sum(ns[selus]) + prior[[gg]][,2]) -
                 prior[[gg]][,1] * log(1 + sum(ns[selus]) / prior[[gg]][,2])) - log(nrow(prior[[gg]]))
      }
    pD
  }


`PDgivenr.NBIndie` <-
function(us, ns, prior, group, priorWeights)
  {
    pD <- 0
    for(gg in 1:length(unique(group)))
      {
        selus <- group == gg
        pD <- pD +
          logsum(
                 rowSums(matrix(
                                dnbinom(rep(us[selus], each = nrow(prior[[gg]])),
                                        size = 1 / prior[[gg]][,2],
                                        mu = rep(ns[selus], each = nrow(prior[[gg]])) * prior[[gg]][,1], log = TRUE),
                                ncol = sum(selus))
                 ) + log(priorWeights)
        )
      }
    pD
  }

