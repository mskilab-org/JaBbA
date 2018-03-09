// The actual solving procedure is called using the function Rcplex 
#include "Rcplex2.h"

SEXP Rcplex2(SEXP numcols_p,
	    SEXP numrows_p,
	    SEXP objsen_p,
	    SEXP cvec,
	    SEXP bvec,
	    SEXP Amat,
	    SEXP Qmat,
	    SEXP lb_p,
	    SEXP ub_p,
	    SEXP Rsense,
	    SEXP Rvtype,
	    SEXP isQP_p,
	    SEXP isMIP_p,
	    SEXP num_poplim,
	    SEXP control,
	    SEXP maxcalls)
{
  return(numcols_p);
}

