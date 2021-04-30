// The actual solving procedure is called using the function Rcplex 
#include "Rcplex2.h"

int max_numcalls;

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
	     SEXP maxcalls,
	     SEXP tuning)
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
  int trace      = INTEGER(getListElement(control,"trace"))[0];

  SEXP res = NULL; /* set to avoid uninitialized warning */
  SEXP xopt;
  SEXP epgap;
  SEXP obj;
  SEXP solstat;
  SEXP extra;
  SEXP lambda;
  SEXP slack;
  SEXP nodecnt;

  int status;
  int i, j;
  int cur_numrows, cur_numcols;

  int tstat;
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

  if(trace) Rprintf("Rcplex: num variables=%d num constraints=%d\n",numcols,numrows);

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
  status = CPXcopylp(env, lp, numcols, numrows, objsen,REAL(cvec),
		     REAL(bvec), sense,
		     INTEGER(VECTOR_ELT(Amat, 0)), 
		     INTEGER(VECTOR_ELT(Amat, 1)),
		     INTEGER(VECTOR_ELT(Amat, 2)),
		     REAL(VECTOR_ELT(Amat, 3)),
		     lb, ub, NULL);
  if (status) {
    my_error(("Failed to copy problem data.\n"));
  }

  if (isQP) {
    status = CPXcopyquad(env, lp,
			 INTEGER(VECTOR_ELT(Qmat, 0)),
			 INTEGER(VECTOR_ELT(Qmat, 1)),
			 INTEGER(VECTOR_ELT(Qmat, 2)),
			 REAL(VECTOR_ELT(Qmat, 3)));
    if (status) {
      my_error(("Failed to copy quadratic term of problem data.\n"));
    }
  }

  if (isMIP) {
    status = CPXcopyctype(env, lp, vtype);
    if (status) {
      my_error(("Failed to copy vtype.\n"));
    }
  }

  /* solve problem */
  if(isMIP) {

    /* MARCIN ADDED set starts */
    setstarts(env, lp, control, isMIP);

    int do_tuning = asLogical(tuning);
    if (do_tuning){
      /* set the tuning total tilim to hard-setted value 600s */
      SEXP tuning_control = setListElement(control, "tilim", 600);
      /* temporarily change the tilim*/
      setparams(env, tuning_control, isQP, isMIP);
      
      /* tune the hidden parameters */
      status = CPXtuneparam(env, lp,
			    0, NULL, NULL, // zero fixed int par
			    0, NULL, NULL,
			    0, NULL, NULL, &tstat);
      if (status) {
	my_error(("Failed to tune params."));
      }
      /* Write the optimized parameters to file */
      status = CPXwriteparam (env, "tuned.prm");

      setparams(env, control, isQP, isMIP);
      /* MARCIN ADDED set starts */
      setstarts(env, lp, control, isMIP);
    } else {
      /* Write the optimized parameters to file */
      status = CPXwriteparam (env, "raw.prm");
    }

    status = CPXmipopt(env, lp);
    
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
      PROTECT(res   = allocVector(VECSXP,  5));
      PROTECT(obj   = allocVector(REALSXP, 1));
      PROTECT(xopt  = allocVector(REALSXP, numcols));
      PROTECT(epgap  = allocVector(REALSXP, 1)); /* added by Marcin */
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

      status = CPXgetmiprelgap(env, lp, REAL(epgap)); /* added by Marcin */

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
      SET_VECTOR_ELT(res, 4, epgap); /* added by Marcin */
      
      UNPROTECT(6);
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

