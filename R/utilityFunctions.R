.logsum <- function(x)
  max(x, max(x, na.rm = TRUE) + log(sum(exp(x - max(x, na.rm = TRUE)), na.rm = TRUE)), na.rm = TRUE)

.logRowSum <- function(z)
  {
    maxes <- do.call(pmax, c(as.list(data.frame(z)), list(na.rm = TRUE)))
    pmax(maxes, maxes + log(rowSums(exp(z - maxes), na.rm = TRUE)), na.rm = TRUE)
  }
