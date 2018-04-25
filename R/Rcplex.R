###########################################################################
## The following function solves a linear or quadratic program (LP|QP).

#' @name Rcplex2
#' @rdname internal
#' @title Rcplex2
#' @description
#'
#' Modification of Rcplex which takes in mipcontrol parameters.
#'
#' 
#' @export
Rcplex2 <- function(cvec, Amat, bvec, Qmat = NULL, lb = 0, ub = Inf,
                   control = list(), objsense = c("min", "max"), sense = "L",
                   vtype = NULL, n = 1)
{
    ## check constraints matrix A
    stopifnot((is(Amat, "matrix") && is.double(Amat)) ||
              is(Amat, "dsparseMatrix") ||
              is(Amat, "simple_triplet_matrix"))

    numrows <- nrow(Amat)
    numcols <- ncol(Amat)

    ## check Q matrix if given
    if (!is.null(Qmat)) {
        stopifnot((is(Qmat, "matrix") && is.double(Qmat)) ||
                 is(Qmat, "dsparseMatrix") ||
                 is(Qmat, "simple_triplet_matrix"),
                nrow(Qmat) == numcols, ncol(Qmat) == numcols)
    }

    ## check data dimensions
    stopifnot(length(cvec) == numcols, is.double(cvec), length(bvec) == numrows, is.double(bvec))


    ## check bounds
    if (length(lb) == 1L){
        lb <- rep(lb,numcols)
    }

    if (length(ub) == 1L){
        ub <- rep(ub,numcols)
    }

    stopifnot(length(lb) == numcols, is.double(lb), length(ub) == numcols, is.double(ub))

    ## check and set objective sense
    if(missing(objsense)){
        objsense <- "min"
    }

    stopifnot(objsense %in% c("min","max"))

    if (objsense == "min") {
        objsensei <- 1L
    } else {
        objsensei <- -1L
    }

    ## check constraints sense
    if (length(sense) == 1L){
        sense <- rep(sense, numrows)
    }

    stopifnot(length(sense) == numrows, is.character(sense), all(sense %in% c('L', 'G', 'E')))

    ## check variable type if needed
    if (!is.null(vtype)) {
        if (length(vtype) == 1L){
            vtype <- rep(vtype, numcols)
        }

        stopifnot(length(vtype) == numcols, is.character(vtype),
                all(vtype %in% c('C', 'I', 'B')))
        isMIP <- ifelse(any(vtype != "C"), 1L, 0L)
    } else{
        isMIP <- 0L
    }

    ## check number of solutions
    n <- as.integer(n)
    ## if NA then find all solutions (at most max. integer solutions)
    ## as Rcplex has troubles with .Machine$integer.max
    ## we take the value from the User's manual
    if(is.na(n)){
        n <- 2.1e+9L
    }

    ## coerce Amat and Qmat to CPX matrices
    Acpx <- toCPXMatrix(Amat)
    Qcpx <- toCPXMatrix(Qmat)

    isQP <- ifelse(is.null(Qcpx), 0L, 1L)

    control <- check.Rcplex.control(control, isQP)
    control <- split.control.list  (control)

    on.exit(.C("Rcplex_free"))

    ## Call the solver
    res <- .Call("Rcplex2",
                as.integer(numcols),
                as.integer(numrows),
                as.integer(objsensei),
                as.double(cvec),
                as.double(bvec),
                Acpx,
                Qcpx,
                as.double(lb),
                as.double(ub),
                as.character(sense),
                as.character(vtype),
                as.integer(isQP),
                as.integer(isMIP),
                as.integer(n),
                control$C,
                as.integer(control$R$maxcalls))

    if (isMIP) {
        intvars <- which(vtype != 'C')
        .canonicalize <- function(x){
            names(x) <- c("xopt", "obj", "status", "extra")
            names(x$extra) <- c("nodecnt", "slack")
            if(control$R$round){
                x$xopt[intvars] <- round(x$xopt[intvars])
            }
            x
        }
        res <- if(n > 1L){
            lapply(res, .canonicalize)
        } else{
            .canonicalize(res)
        }
    } else {
        names(res) <- c("xopt", "obj", "status", "extra")
        names(res$extra) <- c("lambda", "slack")
    }
    return(res)
}



