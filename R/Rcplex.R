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

###########################################################################
## The following function solves a given linear program (LP).

Rcplex_solve_LP <- function( cvec, Amat, bvec, lb = 0, ub = Inf, sense = "L",
                             objsense = c("min", "max"), vtype = NULL, n = 1,
                             control = list() ) {

  ## check and set objective sense
  objsense <- match.arg( objsense )

  ## solve the mathematical program
  Rcplex_solve_QP(cvec = cvec, Amat = Amat, bvec = bvec, Qmat = NULL,
                  lb = lb, ub = ub, sense = sense,objsense = objsense,
                  vtype = vtype, n = n, control = control )
}

###########################################################################
## The following function solves a given quadratic program (QP).
## It is similar to Rcplex() but forces the user to provide Qmat

Rcplex_solve_QP <- function( cvec, Amat, bvec, Qmat, lb = 0, ub = Inf,
                             sense = "L", objsense = c("min", "max"),
                             vtype = NULL, n = 1, control = list() ) {

  ## check and set objective sense
  objsense <- match.arg( objsense )

  ## solve the mathematical program
  Rcplex_solve_QCP( cvec = cvec, Amat = Amat, bvec = bvec, Qmat = Qmat,
                    QC = NULL, lb = lb, ub = ub, sense = sense,
                    objsense = objsense, vtype = vtype, n = n,
                    control = control )
}

###########################################################################
## The following function solves a quadratically constrained program (QCP)

Rcplex_solve_QCP <- function( cvec, Amat, bvec, Qmat = NULL, QC, lb = 0,
                              ub = Inf, sense = "L",
                              objsense = c("min", "max"), vtype = NULL,
                              n = 1, control = list() ) {

  ## check and set objective sense
  objsense <- match.arg(objsense)

  ## number objective variables, number constraints
  numrows <- nrow(Amat)
  numcols <- ncol(Amat)

  ## bounds
  if (length(lb) == 1L){
    lb <- rep(lb, numcols)
  }
  if (length(ub) == 1L){
    ub <- rep(ub, numcols)
  }

  ## check constraints sense
  if (length(sense) == 1L){
    sense <- rep(sense, numrows)
  }

    ## vtypes
    if( !is.null(vtype) ){
        if( length(vtype) == 1L ){
            vtype <- rep( vtype, numcols )
        }
    }

    ## Input validation
    Rcplex_check_input_for_sanity( cvec = cvec, Amat = Amat, bvec = bvec,
                                 Qmat = Qmat, QC = QC, lb = lb, ub = ub,
                                 sense = sense, objsense = objsense,
                                 vtype = vtype, nvars = numcols,
                                 nconstr = numrows )

    .Rcplex_solve( cvec = cvec, Amat = Amat, bvec = bvec,
                 Qmat = Qmat, QC = QC, lb = lb, ub = ub,
                 sense = sense, objsense = objsense,
                 vtype = vtype, n = n, control = control,
                 nvars = numcols, nconstr = numrows )
}

###########################################################################
## Calls the solver

.Rcplex_solve <- function( cvec, Amat, bvec, Qmat, QC, lb, ub,
                           sense, objsense, vtype, n, control, nvars, nconstr) {

  ## prepare date for solver

  ## coerce Amat and Qmat to CPX matrices
  Acpx <- toCPXMatrix(Amat)
  Qcpx <- toCPXMatrix(Qmat)

  if (objsense == "min") {
      objsensei <- 1L
  } else {
      objsensei <- -1L
  }

  ## do we have a mixed integer program (MIP)?
  if (!is.null(vtype)) {
      isMIP <- ifelse(any(vtype != "C"), 1L, 0L)
  } else{
      sMIP <- 0L
  }

  ## do we have a QP?
  isQP <- ifelse(is.null(Qcpx), 0L, 1L)

  ## do we have a QCP?
  isQCP <- ifelse(is.null(QC), 0L, 1L)

  ## check number of solutions
  n <- as.integer(n)
  ## if NA then find all solutions (at most max. integer solutions)
  ## as Rcplex has troubles with .Machine$integer.max
  ## we take the value from the User's manual
  if(is.na(n)){
      n <- 2.1e+9L
  }

  ## check control list
  control <- check.Rcplex.control(control, isQP)
  control <- split.control.list  (control)

  on.exit(.C("Rcplex_free"))

  ## Quadratic constraints
  if(isQCP){
      QC <- as.quadratic_constraint(QC)
  }
  nQC <- length(QC$dir)

  ## Call cplex interface
  res <- .Call( "Rcplex_QCP",
                as.integer(nvars),
                as.integer(nconstr),
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
                as.integer(control$R$maxcalls),
                as.integer(isQCP),
                QC,
                as.integer(nQC) )

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
    res
}

####################################################################
## Sanity checks

Rcplex_check_input_for_sanity <- function(cvec, 
    Amat, 
    bvec,
    Qmat,
    QC, 
    lb, 
    ub, 
    sense,
    objsense, 
    vtype, 
    nvars,
    nconstr) 
{

    ## check constraints matrix A
    stopifnot( Rcplex_check_Amat_for_sanity(Amat) )

    ## check Q matrix if given
    if( !is.null(Qmat) ) {
        stopifnot( Rcplex_check_Qmat_for_sanity(Qmat, nvars) )
    }

    ## check data dimensions
    stopifnot( Rcplex_check_cvec_for_sanity(cvec, nvars), Rcplex_check_bvec_for_sanity(bvec, nconstr) )

    stopifnot( Rcplex_check_bounds_for_sanity(lb, nvars), Rcplex_check_bounds_for_sanity(ub, nvars) )

    stopifnot( Rcplex_check_row_sense_for_sanity(sense, nconstr) )

    stopifnot( Rcplex_check_vtype_for_sanity(vtype, nvars) )

    invisible(TRUE)
}

Rcplex_check_cvec_for_sanity <- function(x, nvars){
    all( c(length(x) == nvars, is.double(x)) )
}

Rcplex_check_bvec_for_sanity <- function(x, nconstr){
    all( c(length(x) == nconstr, is.double(x)) )
}

Rcplex_check_matrix_for_sanity <- function(x){
    ( is( x, "matrix" ) && is.double(x) ) ||
    is( x, "dsparseMatrix" )            ||
    is( x, "simple_triplet_matrix" )
}

Rcplex_check_Amat_for_sanity <- function(x){
    Rcplex_check_matrix_for_sanity(x)
}

Rcplex_check_Qmat_for_sanity <- function(x, nvars){
    all(  c(Rcplex_check_matrix_for_sanity(x), nrow(x) == nvars, ncol(x) == nvars) )
}

Rcplex_check_row_sense_for_sanity <- function(x, nconstr){
    all( c(length(x) == nconstr, is.character(x), x %in% c('L', 'G', 'E')) )
}

Rcplex_check_bounds_for_sanity <- function(x, nvars){
    all( c(length(x) == nvars, is.double(x)) )
}

Rcplex_check_vtype_for_sanity <- function(x, nvars) {
    all( c(length(x) == nvars, is.character(x), x %in% c('C', 'I', 'B')) ) || is.null(x)
}

Rcplex.close <- function() {
    invisible(.C("Rcplex_close"))
}
