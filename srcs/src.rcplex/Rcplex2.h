#ifndef _RCPLEX_H
#define _RCPLEX_H

#include <ilcplex/cplex.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <unistd.h>

#define my_error(x) { forceCplxClose = 1; error x; }

/* Global Variables */
CPXENVptr env;
CPXLPptr lp;
int numcalls;
int max_numcalls;
int forceCplxClose;

/* Function definitions */
SEXP getListElement( SEXP, char* );
SEXP setListElement( SEXP, char*, double);
void setparams( CPXENVptr, SEXP, int, int );
void setstarts( CPXENVptr, CPXLPptr, SEXP, int );
void Rcplex_init( void );
void Rcplex_close( void );
void Rcplex_free( void );

#endif
