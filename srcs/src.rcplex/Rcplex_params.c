#include "Rcplex2.h"

void setparams(CPXENVptr env, SEXP control, int isQP, int isMIP) {
  int i, status, value;
  const char *cur_parm;
  SEXP names;
  
  /* get list names */
  PROTECT(names = getAttrib(control, R_NamesSymbol));

  status = 1; /* avoid warning */

  /* for each element in 'control' try to set the corresponding
     parameter */
  for (i = 0; i < length(control); i++) {

    cur_parm = CHAR(STRING_ELT(names, i));
    
    /* trace - CPX_PARAM_SCRIND */
    if(strcmp(cur_parm, "trace") == 0) {
      status = CPXsetintparam(env, CPX_PARAM_SCRIND, 
		  *INTEGER(VECTOR_ELT(control, i)) ? CPX_ON : CPX_OFF);
    }
    /* method */
    else if(strcmp(cur_parm, "method") == 0) {
      switch((value = *INTEGER(VECTOR_ELT(control, i)))) {
      case 0:
	value = CPX_ALG_AUTOMATIC;
	break;
      case 1:
	value = CPX_ALG_PRIMAL;
	break;
      case 2:
	value = CPX_ALG_DUAL;
	break;
      case 3:
	value = CPX_ALG_NET;
	break;
      case 4:
	value = CPX_ALG_BARRIER;
	break;
      default:
	warning("Unknown optimization method %d, using default\n", value);
	value = CPX_ALG_AUTOMATIC;
      }
      /* Do we have a QP or a LP?*/
      if (isQP)
	status = CPXsetintparam(env, CPX_PARAM_QPMETHOD, value);
      else
	status = CPXsetintparam(env, CPX_PARAM_LPMETHOD, value);
    }
    else if(strcmp(cur_parm, "preind") == 0) {
      status = CPXsetintparam(env, CPX_PARAM_PREIND,
		  *INTEGER(VECTOR_ELT(control, i)) ? CPX_ON : CPX_OFF);
    }
    else if(strcmp(cur_parm, "aggind") == 0) {
      status = CPXsetintparam(env, CPX_PARAM_AGGIND,
			      *INTEGER(VECTOR_ELT(control, i)));
    }
    else if(strcmp(cur_parm, "itlim") == 0) {
      status = CPXsetintparam(env, CPX_PARAM_ITLIM, 
			      *INTEGER(VECTOR_ELT(control, i)));
    }
    else if(strcmp(cur_parm, "epgap") == 0) {
      status = CPXsetdblparam(env, CPX_PARAM_EPGAP,
			      *REAL(VECTOR_ELT(control, i)));
    }
    else if(strcmp(cur_parm, "epagap") == 0) {
      status = CPXsetdblparam(env, CPX_PARAM_EPAGAP, 
    			      *REAL(VECTOR_ELT(control, i)));
    }
    else if(strcmp(cur_parm,"tilim") == 0) {
      status = CPXsetdblparam(env, CPX_PARAM_TILIM, 
			      *REAL(VECTOR_ELT(control, i)));
    }
    else if(strcmp(cur_parm, "mipemphasis") == 0) {
      switch((value = *INTEGER(VECTOR_ELT(control, i)))) {
      case 0:
	value = CPX_MIPEMPHASIS_BALANCED;
	break;
      case 1:
	value = CPX_MIPEMPHASIS_FEASIBILITY;
	break;
      case 2:
	value = CPX_MIPEMPHASIS_OPTIMALITY;
	break;
      case 3:
	value = CPX_MIPEMPHASIS_BESTBOUND;
	break;
      case 4:
	value = CPX_MIPEMPHASIS_HIDDENFEAS;
	break;
      default:
	warning("Unknown mip emphasis setting %d, using default\n", value);
	value = CPX_MIPEMPHASIS_BALANCED;
      }
      status = CPXsetintparam(env, CPX_PARAM_MIPEMPHASIS, value);
    }
    else if(strcmp(cur_parm, "disjcuts") == 0) {
      status = CPXsetintparam(env, CPX_PARAM_DISJCUTS, 
			      *INTEGER(VECTOR_ELT(control,i)));
    }
    else if(strcmp(cur_parm, "cliques") == 0) {
      status = CPXsetintparam(env, CPX_PARAM_CLIQUES, 
			      *INTEGER(VECTOR_ELT(control, i)));
    }
    else if(strcmp(cur_parm,"nodesel") == 0) {
      switch((value = *INTEGER(VECTOR_ELT(control, i)))) {
      case 0:
	value = CPX_NODESEL_DFS;
	break;
      case 1:
	value = CPX_NODESEL_BESTBOUND;
	break;
      case 2:
	value = CPX_NODESEL_BESTEST;
	break;
      case 3:
	value = CPX_NODESEL_BESTEST_ALT;
	break;
      default:
	warning("Unknown node selection strategy %d, using default\n", value);
	value = CPX_NODESEL_BESTBOUND;
      }
      status = CPXsetintparam(env, CPX_PARAM_NODESEL, value);
    }
    else if(strcmp(cur_parm, "probe") == 0) {
      status = CPXsetintparam(env, CPX_PARAM_PROBE,
			      *INTEGER(VECTOR_ELT(control, i)));
    }
    else if(strcmp(cur_parm, "varsel") == 0) {
      switch((value = *INTEGER(VECTOR_ELT(control, i)))) {
      case -1:
	value = CPX_VARSEL_MININFEAS;
	break;
      case 0:
	value = CPX_VARSEL_DEFAULT;
	break;
      case 1:
	value = CPX_VARSEL_MAXINFEAS;
	break;
      case 2:
	value = CPX_VARSEL_PSEUDO;
	break;
      case 3:
	value = CPX_VARSEL_STRONG;
	break;
      case 4:
	value = CPX_VARSEL_PSEUDOREDUCED;
	break;
      default:
	warning("Unknown variable selection strategy %d, using default\n", 
		value);
	value = CPX_VARSEL_DEFAULT;
      }
      status = CPXsetintparam(env, CPX_PARAM_VARSEL, value);
    }
    else if(strcmp(cur_parm, "flowcovers") == 0) {
      status = CPXsetintparam(env, CPX_PARAM_FLOWCOVERS,
			      *INTEGER(VECTOR_ELT(control, i)));
    }
    else if(strcmp(cur_parm, "mipstart") == 0) {
    }

    else if(strcmp(cur_parm, "solnpoolagap") == 0){
      /* solution pool parameters */
      #if CPX_VERSION >= 1100 
      status = CPXsetdblparam(env, CPX_PARAM_SOLNPOOLAGAP,
      			      *REAL(VECTOR_ELT(control, i)));
      #endif
    }
    else if(strcmp(cur_parm, "solnpoolgap") == 0){
      /* solution pool parameters */
      #if CPX_VERSION >= 1100
      status = CPXsetdblparam(env, CPX_PARAM_SOLNPOOLGAP,
			      *REAL(VECTOR_ELT(control, i)));
      #endif
    }
    else if(strcmp(cur_parm, "solnpoolintensity") == 0){
      #if CPX_VERSION >= 1100
      status = CPXsetintparam(env, CPX_PARAM_SOLNPOOLINTENSITY, 
    			      *INTEGER(VECTOR_ELT(control, i)));
      #endif
    }
    else {
      /* If parameter not known print a warning */
      warning("Unknown CPLEX parameter %s. Ignoring it.\n", cur_parm);
    }

    if (status)
      my_error(("Failure to set parameter %s, error %d.\n", cur_parm, status));
  }
  
  UNPROTECT(1);
}

/* ADDED BY MARCIN!!
 * 
 */
void setstarts(CPXENVptr env, CPXLPptr lp, SEXP control, int isMIP) {
  int i, status, value;
  const char *cur_parm;
  SEXP names;
  
  /* get list names */
  PROTECT(names = getAttrib(control, R_NamesSymbol));

  status = 1; /* avoid warning */

  /* for each element in 'control' try to set the corresponding
     parameter */
  if (isMIP) {
    for (i = 0; i < length(control); i++) {

      cur_parm = CHAR(STRING_ELT(names, i));

      if(strcmp(cur_parm, "mipstart") == 0) {
        SEXP values_vec = VECTOR_ELT(control, i);
        double* values = REAL(values_vec);
        int mcnt = 1;
        int beg = 0;
        int nzcnt = LENGTH(values_vec);
        SEXP tmp = allocVector(INTSXP, nzcnt);
        int *varindices = INTEGER(tmp);
        int effortlevel = CPX_MIPSTART_AUTO;
        char *mipstartname = "start";

        for (int i = 0; i < nzcnt; i++) {
          varindices[i] = i;
        }

        status = CPXaddmipstarts(env, lp, mcnt, nzcnt, &beg, varindices, values, &effortlevel, &mipstartname);

        if (status)
          my_error(("Failure to set parameter %s in setstarts, error %d.\n", cur_parm, status));
      }
    }
  }
  
  UNPROTECT(1);
}
