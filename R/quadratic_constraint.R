####################################################################
## Code for handling quadratic constraints

## Constructor for building quadratic constraints
quadratic_constraint <- function( QC, dir, b ) {
    ## from refcallablelibrary.pdf: only 'L' and 'G' are allowed in
    ## Q constraints
    stopifnot( all(dir %in% c("L", "G")) )
    len_Q <- length(QC$Q)
    len_L <- length(QC$L)
    if(!len_L){
        len_L <- len_Q
        QC$L <- vector("list", len_L)
    }
    stopifnot( all(c(len_Q, len_L, length(dir)) == length(b)) )
    QC$Q <- lapply( QC$Q, as.simple_triplet_matrix )
    QC$Q <- lapply( QC$Q, function(x) {
        x$i <- x$i - 1L
        x$j <- x$j - 1L
        x 
    })
    ## FIXME: Currently only dense vectors as input supported
    QC$L <- lapply( QC$L, as.sparse_vector)
    QC$linnzcount <- sapply( QC$L, function(x) {
                                         as.integer(length(x$v))} )
    if( !length(QC$linnzcount) ){
        QC$linnzcount <- rep(0L, len_Q)
    }

    ##  if( any( QC$linnzcount > 0 ) )
    ##     stop( "linear terms in quadratic constraints not supported in Rcplex." )

    QC$quadnzcount <- as.integer(lapply(QC$Q, function(x) length(x$v)))

    out <- list( QC, dir = as.character(dir), b = as.double(b) )

    class(out) <- "quadratic_constraint"
    out
}

as.sparse_vector <- function(x){
    UseMethod("as.sparse_vector")
}

as.sparse_vector.NULL <- function(x) {
    list(i = as.integer(NULL), v = as.double(NULL))
}

as.sparse_vector.numeric <- function(x){
    ind <- which(x != 0)
    val <- x[ind]
    list(i = as.integer(ind), v = as.double(val))
}

as.sparse_vector.simple_triplet_matrix<- function(x){
    stopifnot( dim(x)[1] == 1 )
    list(i = as.integer(x$j), v = as.double(x$v))
}

as.quadratic_constraint <- function(x, ...){
    UseMethod("as.quadratic_constraint")
}

as.quadratic_constraint.list <- function(x){
    quadratic_constraint(x$QC, x$dir, x$b)
}

as.quadratic_constraint.quadratic_constraint <- function(x){
    x
}

is.quadratic_constraint <- function(x){
    inherits(x, "quadratic_constraint")
}
