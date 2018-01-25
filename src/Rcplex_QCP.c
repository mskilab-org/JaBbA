// The actual solving procedure is called using the function Rcplex 
#include "Rcplex.h"

SEXP Rcplex_QCP(SEXP numcols_p,
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
		SEXP maxcalls,
		SEXP isQCP_p,
		SEXP QC,
		SEXP nQC)
{
  char *probname = "Rcplex";
  int numcols    = INTEGER(numcols_p)[0]; 
  int numrows    = INTEGER(numrows_p)[0];
  int objsen     = INTEGER(objsen_p)[0];
  double *lb     = REAL(lb_p); 
  double *ub     = REAL(ub_p);

  char sense[numrows];
  char vtype[numcols];
  //  double dj[numrows];
  int isQP       = INTEGER(isQP_p)[0];
  int isMIP      = INTEGER(isMIP_p)[0];
  int isQCP       = INTEGER(isQCP_p)[0];
  int trace      = INTEGER(getListElement(control,"trace"))[0];

  SEXP res = NULL; /* set to avoid uninitialized warning */
  SEXP xopt;
  SEXP obj;
  SEXP solstat;
  SEXP extra;
  SEXP lambda;
  SEXP slack;
  SEXP nodecnt;

  int status;
  int i, j;
  int cur_numrows, cur_numcols;

  /*
   * solution pools not supprted until cplex 11.0
   */
  #if CPX_VERSION >= 1100 
  int num_sol = 1;
  SEXP tmp;
  #endif

  /* set maxnumcalls before init */
  max_numcalls = INTEGER(maxcalls)[0];
  /* initialize CPLEX environment */
  Rcplex_init();

  //if(trace) Rprintf("Rcplex: num variables=%d num constraints=%d\n",numcols,numrows);

  //if(trace) Rprintf("Rcplex: isQP=%d isQCP=%d isMIP=%d\n", isQP, isQCP, isMIP);

  /*
   * solution pools not supported until cplex 11.0
   */
  #if CPX_VERSION < 1100
  if (INTEGER(num_poplim)[0] > 1) {
    warning("Multiple solutions not supported in CPLEX version");
    INTEGER(num_poplim)[0] = 1;
  }
  #endif 

  /* lb and ub */
  for (j = 0; j < numcols; ++j) {
    lb[j] = R_finite(lb[j]) ? lb[j] : -CPX_INFBOUND;
    ub[j] = R_finite(ub[j]) ? ub[j] : CPX_INFBOUND;
  }

  /* set constraint inequality directions */
  for (i = 0; i < numrows; ++i) {
    sense[i] = CHAR(STRING_ELT(Rsense, i))[0];
  }

  /* set variable types */
  if (isMIP) {
    for (i = 0; i < numcols; ++i) {
      vtype[i] = CHAR(STRING_ELT(Rvtype, i))[0];
    }
  }

  /* set parameters given in control list */
  setparams(env, control, isQP, isMIP);
  lp = CPXcreateprob(env, &status, probname);

  /* check memory problems */
  if (lp == NULL) {
    my_error(("Failed to create LP.\n"));
  }

  /* copy problem data into lp */
  /* Linear part of objective function c^t x; Ax ~ b*/
  /* Linear part of constraints */
  status = CPXcopylp(env, lp, numcols, numrows, objsen, REAL(cvec),
		     REAL(bvec), sense,
		     INTEGER(VECTOR_ELT(Amat, 0)), 
		     INTEGER(VECTOR_ELT(Amat, 1)),
		     INTEGER(VECTOR_ELT(Amat, 2)),
		     REAL(VECTOR_ELT(Amat, 3)),
		     lb, ub, NULL);
  if (status) {
    my_error(("Failed to copy problem data.\n"));
  }
  //if(trace) Rprintf("Rcplex: done copying linear part\n");

  /* Quadratic part of objective function 1/2 x^t Q x */
  /*int CPXcopyquad(CPXCENVptr env, CPXLPptr lp, 
                    const int * qmatbeg, 
                    const int * qmatcnt, 
                    const int * qmatind, 
                    const double * qmatval) */ 
  /* The arrays qmatbeg and qmatcnt should be of length at least
     CPXgetnumcols(env,lp). The arrays qmatind and qmatval should be
     of length at least qmatbeg[numcols-1]+qmatcnt[numcols-1]. CPLEX
     requires only the nonzero coefficients grouped by column in the
     array qmatval. The nonzero elements of every column must be
     stored in sequential locations in this array with qmatbeg[j]
     containing the index of the beginning of column j and qmatcnt[j]
     containing the number of entries in column j. Note that the
     components of qmatbeg must be in ascending order. For each k,
     qmatind[k] indicates the column number of the corresponding
     coefficient, qmatval[k]. These arrays are accessed as explained
     above. */

  if (isQP) {
    status = CPXcopyquad(env, lp,
			 INTEGER(VECTOR_ELT(Qmat, 0)),
			 INTEGER(VECTOR_ELT(Qmat, 1)),
			 INTEGER(VECTOR_ELT(Qmat, 2)),
			 REAL(VECTOR_ELT(Qmat, 3)));
    if (status) {
      my_error(("Failed to copy quadratic term of problem data.\n"));
    }
  //  if(trace) Rprintf("Rcplex: done copying quadratic part\n");
  }
  
  /* Quadratic part of constraints a_i^tx + 1/2 x^t Q_i x <= r_i*/
  /* refcallablelibrary v12: p. 492 */
  /* int CPXaddqconstr(CPXCENVptr env, CPXLPptr lp, int linnzcnt, 
                       int quadnzcnt, double rhs, int sense, 
                       const int * linind, const double * linval, 
                       const int * quadrow, const int * quadcol, 
                       const double * quadval, const char *
                       lname_str) */
  /* The nonzero coefficients of the quadratic terms must be stored in
     sequential locations in the arrays quadrow, quadcol and quadval
     from positions 0 to quadnzcnt-1. Each pair, quadrow[i],
     quadcol[i], indicates the variable indices of the quadratic term,
     and quadval[i] the corresponding coefficient. */

  if (isQCP) {
    /* Q_i will be checked if positive semi definite by the callable
       lib */
    for(i = 0; i < INTEGER(nQC)[0]; i++){
      
      status = CPXaddqconstr (env, lp, 
			      INTEGER( VECTOR_ELT(VECTOR_ELT(QC, 0), 2) )[i],
			      INTEGER( VECTOR_ELT(VECTOR_ELT(QC, 0), 3) )[i],
			      REAL(    VECTOR_ELT(QC, 2) )[i],
			      CHAR(    STRING_ELT(VECTOR_ELT(QC, 1), i) )[0],
			      INTEGER( VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(QC, 0), 1), i), 0) ),
			      REAL(    VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(QC, 0), 1), i), 1) ),
			      INTEGER( VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(QC, 0), 0), i), 0) ),
                              INTEGER( VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(QC, 0), 0), i), 1) ),
                              REAL(    VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(QC, 0), 0), i), 2) ),
			      NULL);
      if (status) {
	my_error(("Failed to add quadratic constraints to problem data.\n"));
      }
      //if(trace) Rprintf("Rcplex: done copying quadratic constraints\n");
    }
  }
  
  if (isMIP) {
    //if(trace) Rprintf("Rcplex: try setting MIP type\n");
    status = CPXcopyctype(env, lp, vtype);
    if (status) {
      my_error(("Failed to copy vtype.\n"));
    }
    //if(trace) Rprintf("Rcplex: done setting MIP type\n");
  }
  
  /* solve problem */
  if(isMIP) {
    //if(trace) Rprintf("Rcplex: try optimizing MIP type\n");
    status = CPXmipopt(env, lp);
    //if(trace) Rprintf("Rcplex: done optimizing MIP type\n");
    
    /*
     * solutions pool not supported for versions of cplex < 11
     * e.g., CPX_PARAM_POPULATELIM is not defined
     */

    #if CPX_VERSION >= 1100
    if(INTEGER(num_poplim)[0] > 1){
    /* in MIPs it is possible to have more than 1 solution. If
       num_poplim > 1 we now try to find these solutions.
       For populating the solution pool a couple of parameters are
       relevant:
       1. the 'absolute gap' solution pool parameter
          CPX_PARAM_SOLNPOOLAGAP (see Rcplex C control parameters), 
       2. the 'solution pool intensity' parameter
          CPX_PARAM_SOLNPOOLINTENSITY (see Rcplex C control parameters),
       3. the 'limit on number of solutions generated' for solution
          pool CPX_PARAM_POPULATELIM */
      status = CPXsetintparam(env, CPX_PARAM_POPULATELIM, 
			      INTEGER(num_poplim)[0] - 1);
      if (status) {
	my_error(("Failed to set 'populate limit' parameter.\n"));
      }
      /* now populate the solutions pool */
      status = CPXpopulate(env, lp);
    }
    #endif
  } 
  else if (isQCP) {
   status = CPXbaropt (env, lp); 
  }
  else if (isQP) {
    status = CPXqpopt(env, lp);
  }
  else {
    status = CPXlpopt(env, lp);
  }
  
  /* a status value of zero does not necessarily mean that an optimal
     solution exists. Examples of an exit status of non-zero here
     include CPXERR_NO_PROBLEM, CPXERR_NO_MEMORY, ... */
  if (status) {
    my_error(("Failed to optimize problem."));
  }
  //if(trace) Rprintf("Rcplex: optimizer finished work.\n");
  PROTECT(solstat = allocVector(INTSXP,  1));

  /* retrieve status of optimization */
  *INTEGER(solstat) = CPXgetstat(env, lp);

  if (isMIP) {
    /* MIP optimization */
    /* in MIPs it is possible to have more than 1 solution
       FIXME: use solution pool by default? */
    
    cur_numrows = CPXgetnumrows(env, lp);
    cur_numcols = CPXgetnumcols(env, lp);

    /* if the 'n' parameter is set in R we always return a list with
       a number of elements equal to the number of solutions */
    if(INTEGER(num_poplim)[0] > 1){
      /*
       * solution pools not supported until cplex 11.0
       */
      #if CPX_VERSION >= 1100 
      num_sol = CPXgetsolnpoolnumsolns(env, lp);

      /* MIP optimization results -> more than 1 solution */
      PROTECT(res = allocVector(VECSXP, num_sol));

      /* now retrieve multiple solutions if any */
      for( i = 0; i < num_sol; i++){
	PROTECT(tmp  = allocVector(VECSXP, 4));
    	PROTECT(xopt = allocVector(REALSXP,  cur_numcols));
    	status = CPXgetsolnpoolx(env, lp, i, REAL(xopt), 0, cur_numcols - 1);
    	
	SET_VECTOR_ELT(tmp, 0, xopt);

	PROTECT(obj   = allocVector(REALSXP, 1));
	status = CPXgetsolnpoolobjval(env, lp, i, REAL(obj));
	/* if no integer solution exists, return NA */
	if (status) {
	  *REAL(obj) = NA_REAL;
	}
	SET_VECTOR_ELT(tmp, 1, obj);
	
	/* extra info */
    	PROTECT(slack  = allocVector(REALSXP,  cur_numrows));
    	status = CPXgetmipslack(env, lp, REAL(slack), 0, cur_numrows - 1);
    	PROTECT(nodecnt = allocVector(INTSXP, 1));
	*INTEGER(nodecnt) = CPXgetnodecnt(env, lp);

	SET_VECTOR_ELT(tmp, 2, solstat); 

	PROTECT(extra   = allocVector(VECSXP,  2));
	SET_VECTOR_ELT(extra, 0, nodecnt);
	SET_VECTOR_ELT(extra, 1, slack);
	SET_VECTOR_ELT(tmp, 3, extra);
	
	/* add solution to return vector */
	SET_VECTOR_ELT(res, i, tmp);

    	UNPROTECT(6);
      } /* END FOR */
      #endif /* end #if CPX_VERSION >= 1100 */

    } /* END multiple solutions */
    else {
      /* MIP optimization 1 solution */
      PROTECT(res   = allocVector(VECSXP,  4));
      PROTECT(obj   = allocVector(REALSXP, 1));
      PROTECT(xopt  = allocVector(REALSXP, numcols));
      PROTECT(extra = allocVector(VECSXP,  2));
      PROTECT(slack = allocVector(REALSXP, numrows));

      status = CPXgetmipobjval(env, lp, REAL(obj));
      /* if no integer solution exists, return NA */
      if (status) {
	*REAL(obj) = NA_REAL;
      }

      status = CPXgetmipx(env, lp, REAL(xopt), 0, cur_numcols - 1);
      if (status) {
	for(i = 0; i < cur_numcols; i++)
	  REAL(xopt)[i] = NA_REAL;
      }
      
      status = CPXgetmipslack(env, lp, REAL(slack), 0, cur_numrows - 1);
      if (status) {
	for(i = 0; i < cur_numrows; i++)
	  REAL(slack)[i] = NA_REAL;
      }

      /* Provide some little extra information */
      PROTECT(nodecnt = allocVector(INTSXP, 1));
      *INTEGER(nodecnt) = CPXgetnodecnt(env, lp);
      SET_VECTOR_ELT(extra, 0, nodecnt);
      SET_VECTOR_ELT(extra, 1, slack);
    
      /* Create return vector for MIP 1 solution*/
      SET_VECTOR_ELT(res, 0, xopt); 
      SET_VECTOR_ELT(res, 1, obj); 
      SET_VECTOR_ELT(res, 2, solstat); 
      SET_VECTOR_ELT(res, 3, extra);
      
      UNPROTECT(5);
    } /* END MIP optimization 1 solution */
  } /* END MIP optimization */
  else {
    /* continuous optimization */
    PROTECT(obj   = allocVector(REALSXP, 1));
    PROTECT(xopt  = allocVector(REALSXP, numcols));
    PROTECT(extra = allocVector(VECSXP,  2));
    PROTECT(slack = allocVector(REALSXP, numrows));

    status = CPXgetobjval(env, lp, REAL(obj));
    if (status) {
      *REAL(obj) = NA_REAL;
    }

    cur_numrows = CPXgetnumrows(env, lp);
    cur_numcols = CPXgetnumcols(env, lp);
    status = CPXgetx(env, lp, REAL(xopt), 0, cur_numcols - 1);
    if (status) {
      for(i = 0; i < cur_numcols; i++)
	REAL(xopt)[i] = NA_REAL;
    }

    status = CPXgetslack(env, lp, REAL(slack), 0, cur_numrows - 1);
    if (status) {
      for(i = 0; i < cur_numrows; i++)
	REAL(slack)[i] = NA_REAL;
    }

    /* Provide some little extra information */
    PROTECT(lambda = allocVector(REALSXP, numrows));
    status = CPXgetpi(env, lp, REAL(lambda), 0, cur_numrows - 1);
    if (status) {
      for(i = 0; i < cur_numrows; i++)
	REAL(lambda)[i] = NA_REAL;
    }
    SET_VECTOR_ELT(extra, 0, lambda);
    SET_VECTOR_ELT(extra, 1, slack);

    /* Create return vector for continuous solution */
    PROTECT(res = allocVector(VECSXP,4));
    SET_VECTOR_ELT(res, 0, xopt); 
    SET_VECTOR_ELT(res, 1, obj); 
    SET_VECTOR_ELT(res, 2, solstat); 
    SET_VECTOR_ELT(res, 3, extra);
    
    UNPROTECT(5);
  } /* END IF continuous optimization */
  

  UNPROTECT(2); /* unprotect return-vector res and solstat */
  
  /* Reset all CPLEX parameters to default values */
  status = CPXsetdefaults(env);
  if (status) {
    Rprintf("Status: %d", status);
    my_error(("Failed to set parameters to default.\n"));
  }
  return(res);
}

