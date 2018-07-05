toCPXMatrix <- function(Amat){

    ## if (is.null(Amat)) {
    ##     return(NULL)
    ## } else if (is(Amat, "sparseMatrix")) {
        Amat <- as(Amat, "dgCMatrix")
        matbeg <- Amat@p
        matcnt <- diff(c(Amat@p, length(Amat@x)))
        matind <- Amat@i
        matval <- Amat@x
    ## } else if (is(Amat, "simple_triplet_matrix")) {
    ##     matbeg <- c(0L, cumsum(tabulate(Amat$j, Amat$ncol)))
    ##     matcnt <- tabulate(Amat$j, Amat$ncol)
    ##     matind <- Amat$i - 1L
    ##     matval <- Amat$v
    ## } else {
    ##     matbeg <- (0L:(ncol(Amat) - 1L)) * nrow(Amat)
    ##     matcnt <- rep(nrow(Amat), ncol(Amat))
    ##     matind <- rep(0L:(nrow(Amat) - 1L), ncol(Amat))
    ##     matval <- as.vector(Amat)
    ## }
    
    return(list(matbeg = as.integer(matbeg), matcnt = as.integer(matcnt),
        matind = as.integer(matind), matval = as.double(matval)))
}
