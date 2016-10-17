marginaliseEqual <- function(cD, r1, r2) {
    sameGroups <- sapply(strsplit(names(cD@groups), "\\},\\{"), function(x) {
        z <- gsub("[\\{\\}]", "", x)
        any(
            sapply(strsplit(z, ","), function(zz)
                all(c(r1, r2) %in% zz)))
    })
    rowSums(exp(cD@posteriors[,sameGroups, drop = FALSE]))
}    



marginalisePairwise <- function (cD, greaterThan, lessThan) 
{
    diffGroups <- sapply(strsplit(names(cD@groups), "\\},\\{"), 
                         function(x) {
                             z <- gsub("[\\{\\}]", "", x)
                             sum(sapply(strsplit(z, ","), function(zz) xor(all(lessThan %in% zz), all(greaterThan %in% zz)))) == 2
                         })
    diffPosts <- sapply(which(diffGroups), function(gg) {
        cD@orderings[, gg]
        gnames <- strsplit(names(groups(cD)[gg]), "\\},\\{")[[1]]
        gord <- grep(paste(grep(greaterThan[1], gnames), ".*>.*", 
                           grep(lessThan[1], gnames), sep = ""), cD@orderings[, 
                                                                           gg])
        posts <- rep(NA, nrow(cD))
        posts[gord] <- cD@posteriors[gord, gg]
        posts
    })
    rowSums(exp(diffPosts), na.rm = TRUE)
}
