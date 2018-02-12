// Utility functions like CPLEX initialization/closing routines and
// list accessors

#include "Rcplex.h"
#include <time.h>

void Rcplex_wait (int seconds) {
  clock_t endwait;
  endwait = clock() + seconds * CLOCKS_PER_SEC;
  while (clock() < endwait) {}
}

void Rcplex_init(void) {
  int numtries = 10;
  int status;
  char errmsg[1024];

  /* Initialize CPLEX environment */
  if (env == NULL) {
    env = CPXopenCPLEX (&status);
    while(env == NULL && numtries > 0) {
      Rcplex_wait(30);
      numtries--;
    }
    if (env == NULL) {
      CPXgeterrorstring(env, status, errmsg);
      error("Could not open CPLEX environment.\n%s\n",errmsg);
    }
    REprintf("CPLEX environment opened\n");
    numcalls = max_numcalls;
  } 
  else {
    numcalls--;
  }
  forceCplxClose = 0;
}

void Rcplex_close(void) {
  forceCplxClose = 1;
  Rcplex_free();
}

void Rcplex_free(void) {
  int status1, status2;
  char errmsg[1024];

  status1 = status2 = 0;
  if (lp != NULL) {
    status1 = CPXfreeprob(env,&lp);
    /*    REprintf("Freed CPLEX problem\n");*/
    lp = NULL;
  }

  if (env != NULL && (numcalls == 0 || forceCplxClose)) {
    status2 = CPXcloseCPLEX(&env);
    REprintf("Closed CPLEX environment\n");
    env = NULL;
  }

  if (status1 || status2) {
    status2 ? strcpy(errmsg,"env close ok") : CPXgeterrorstring(env,status2,errmsg);
    error("Rcplex_free failed: free problem code: %d\nClose environment msg: %s\n",status1,errmsg);
  }
  return;
}

SEXP getListElement(SEXP list, char *str) {
  SEXP element = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;

  for (i=0; i < length(list); i++) {
    if (strcmp(CHAR(STRING_ELT(names,i)), str) == 0) {
      element = VECTOR_ELT(list,i);
      break;
    }
  }
  return element;
}    
