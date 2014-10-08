.bimodalSep <- function(z, weights = NULL, bQ = c(0,1))
  {
    if(is.null(weights))
      weights <- rep(1, length(z))

    decord <- order(z, decreasing = TRUE)
    dsortp <- z[decord]
    dwts <- weights[decord]
    dcummeans <- cumsum(dsortp * dwts) / cumsum(dwts)
    dcumsquare <- cumsum(dsortp ^ 2 * dwts) / (cumsum(dwts) - 1)
    dcumvar <- dcumsquare - (dcummeans^2) * cumsum(dwts) / (cumsum(dwts) - 1)

    incord <- order(z, decreasing = FALSE)
    isortp <- z[incord]
    iwts <- weights[incord]
    icummeans <- cumsum(isortp * iwts) / cumsum(iwts)
    icumsquare <- cumsum(isortp ^ 2 * iwts) / (cumsum(iwts) - 1)
    icumvar <- icumsquare - (icummeans^2) * cumsum(iwts) / (cumsum(iwts) - 1)
    
    meanvars <- (c(rev(dcumvar)[-1], NA) + icumvar) / 2

    if(round(bQ[1] * length(meanvars)) > 0)
      meanvars[1:round(bQ[1] * length(meanvars))] <- NA
    if(round(bQ[2] * length(meanvars)) < length(meanvars))
      meanvars[length(meanvars):round(bQ[2] * length(meanvars))] <- NA

    
    mean(isortp[which.min(meanvars) + 0:1], na.rm = TRUE)
  }
