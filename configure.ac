dnl R interface to CPLEX solvers
dnl
dnl A working installation of CPLEX is required
dnl
dnl	USAGE: autoconf
dnl
dnl Global variables:
dnl	   PKG_CFLAGS, PKG_CPPFLAGS and PKG_LIBS
dnl
dnl PKG_CFALGS is no longer required (is optional now)
dnl
dnl
dnl The resulting configure script implements the following:
dnl
dnl 1. If PKG_CFLAGS, PKG_CPPFLAGS and PKG_LIBS are defined, use them and exit
dnl
dnl 2. If user specifies locations for the CPLEX include file and library as well as C flags we set the global variables accordingly
dnl    This can be specified through the following arguments (all three MUST be given)
dnl	    --with-cplex-include=<include-DIR>
dnl	    --with-cplex-cflags=<required C flags>
dnl	    --with-cplex-lib=<library-DIR> <lib flags>
dnl
dnl 3. If user specifies the location of the CPLEX installation, we will guess at the location of the example file
dnl    i.e. guess the system and libformat and set global variables using the found Makefile. This can be set with
dnl	    --with-cplex-dir=<cplex-DIR>
dnl
dnl 4. If none of the above are given, the script tries to find the location of the CPLEX installation and then does the same as 3.


AC_INIT([Rcplex], 0.3.2)

FOUND_ALL="no"

dnl Check if global variables are already set, if so we are done
if test [-n "$PKG_CPPFLAGS" -a -n "$PKG_LIBS"] ; then
   FOUND_ALL="yes"
else
   if test [-n "$PKG_CPPFLAGS" -o -n "$PKG_LIBS"] ; then
      AC_MSG_ERROR([All three of PKG_CFLAGS, PKG_CPPFLAGS and PKG_LIBS must be given together])
   fi
fi

dnl get arguments
AC_ARG_WITH([cplex-include],
	[AS_HELP_STRING([--with-cplex-include],[specify location of cplex include directory])],
	CPLEX_INCLUDE=$withval,
	CPLEX_INCLUDE="")

AC_ARG_WITH([cplex-cflags],
	[AS_HELP_STRING([--with-cplex-cflags],[specify C flags for compilation of interface])],
	CPLEX_CFLAGS=$withval,
	CPLEX_CFLAGS="")

AC_ARG_WITH([cplex-lib],
	[AS_HELP_STRING([--with-cplex-lib],[specify location of CPLEX callable library and other libraries])],
	CPLEX_LIB=$withval,
	CPLEX_LIB="")

AC_ARG_WITH([cplex-dir],
	[AS_HELP_STRING([--with-cplex-dir],[specify location of CPLEX installation directory])],
	CPLEX_DIR=$withval,
	CPLEX_DIR="")

dnl Check if include flags and libs given using --with arguments
if test ["$FOUND_ALL"="no"] ; then
   if test [-n "$CPLEX_INCLUDE" -a -n "$CPLEX_LIB" ] ; then
      if test [! -d "$CPLEX_INCLUDE"] ; then
	 AC_MSG_ERROR([CPLEX include directory $CPLEX_INCLUDE does not exist])
      fi

dnl      if test [! -d "$CPLEX_LIB"] ; then
dnl	 AC_MSG_ERROR([CPLEX library directory $CPLEX_LIB does not exist])
dnl      fi

      PKG_CPPFLAGS="-I$CPLEX_INCLUDE"
      if test [-n "$CPLEX_CFLAGS"] ; then
            PKG_CFLAGS="$CPLEX_CFLAGS"   
      fi
      PKG_LIBS="$CPLEX_LIB"
      FOUND_ALL="yes"
   else
      if test [-n "$CPLEX_INCLUDE" -o -n "$CPLEX_LIBS" ] ; then
	 AC_MSG_ERROR([All options --with-cplex-include, --with-cplex-cflags and --with-cplex-lib must be given together])
      fi
   fi
fi

dnl No information given try to guess it from location of interactive optimizer
if test ["$FOUND_ALL" = "no" -a -z "$CPLEX_DIR"] ; then
   AC_PATH_PROG(CPLEX_BIN, cplex)
   CPLEX_DIR=`AS_DIRNAME(["$CPLEX_BIN"])`
   ## Autoconf docs say that basename is not portable and one should use
   ## expr instead.  Let us use AS_BASENAME, even though this is not in
   ## the docs ...
   SYSTEM=`AS_BASENAME(["$CPLEX_DIR"])`
   CPLEX_DIR=`AS_DIRNAME(["$CPLEX_DIR"])`
   CPLEX_DIR=`AS_DIRNAME(["$CPLEX_DIR"])`
   if test [-n "$CPLEX_DIR"] ; then
      echo "Trying to use $CPLEX_DIR as CPLEX installation location"
   fi
fi

dnl Use the given installation directory --with-cplex-dir 
if test ["$FOUND_ALL" = "no"] ; then
   if test [! -d "$CPLEX_DIR/include"] ; then
      AC_MSG_ERROR([CPLEX include directory ${CPLEX_DIR}/include does not exist])
   fi

   if test [-z "$SYSTEM" ] ; then
      ## Set system as the first found
      SYSTEM=`ls -l ${CPLEX_DIR}/bin | awk '($1 ~ /d/) {print $NF}' | head -n1 -`
      if test [! -d "${CPLEX_DIR}/bin/${SYSTEM}" ] ; then
            AC_MSG_ERROR([CPLEX library location could not be found])
      fi
   fi

   ## Guess libformat: just grab the first that comes up.
   LIBFORMAT=`ls -l ${CPLEX_DIR}/lib/${SYSTEM} | awk '($1 ~ /d/) {print $NF}' | head -n1 -`
   if test [-z "$LIBFORMAT" -a ! -d "${CPLEX_DIR}/lib/${SYSTEM}/${LIBFORMAT}"] ; then
      AC_MSG_ERROR([CPLEX library location could not be found])
   fi

   CPLEX_MAKEFILE="${CPLEX_DIR}/examples/${SYSTEM}/${LIBFORMAT}/Makefile"
   if test [! -f "$CPLEX_MAKEFILE"] ; then
      AC_MSG_ERROR([CPLEX example Makefile could not be found])
   fi

dnl check for awk
   AC_PROG_AWK

   CPLEXLIBDIR="${CPLEX_DIR}/lib/${SYSTEM}/${LIBFORMAT}"
   CPLEX_LIBS="-L${CPLEXLIBDIR} `${AWK} 'BEGIN {FS = " = "} /^CLNFLAGS/ {print $2}' ${CPLEX_MAKEFILE}`"
   if test [-z CPLEX_LIBS] ; then
      AC_MSG_ERROR([Could not get library info from CPLEX example Makefile])
   fi

   CPLEX_COPT=`${AWK} 'BEGIN {FS = " = "} /^COPT/ {print $2}' ${CPLEX_MAKEFILE}`
   if test [-z CPLEX_COPT] ; then
      AC_MSG_ERROR([Could not get compiler options from CPLEX example Makefile])
   fi

   PKG_CFLAGS="$CPLEX_COPT"
   PKG_CPPFLAGS="-I$CPLEX_DIR/include"
   PKG_LIBS="$CPLEX_LIBS"
   FOUND_ALL="yes"
fi


if test "$FOUND_ALL" = "no" ; then	
   AC_MSG_ERROR([No cplex installation found])
fi

dnl Find compiler and compiler flags
: ${R_HOME=`R RHOME`}	  

if test -z "${R_HOME}" ; then
   echo "could not determine R_HOME"
   exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CPP=`"${R_HOME}/bin/R" CMD config CPP`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_PROG_CC
AC_PROG_CPP

dnl Check for header and library
LIBS="${PKG_LIBS}"
CFLAGS="${CFLAGS} ${PKG_CFLAGS}"
CPPFLAGS="${CPPFLAGS} ${PKG_CPPFLAGS}"

AC_CHECK_HEADER(ilcplex/cplex.h)
if test "${ac_cv_header_ilcplex_cplex_h}" == no; then
   AC_MSG_ERROR([CPLEX header cplex.h not found]);
fi

AC_SEARCH_LIBS(CPXversion,cplex, ,AC_MSG_ERROR([Could not link CPLEX library using $PKG_LIBS]))

dnl substitute in src/Makevars
AC_SUBST(PKG_CFLAGS)
AC_SUBST(PKG_LIBS)
AC_SUBST(PKG_CPPFLAGS)
AC_OUTPUT(src/Makevars)

exit 0

