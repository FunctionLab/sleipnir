# This file is provided under the Creative Commons Attribution 3.0 license.
#
# You are free to share, copy, distribute, transmit, or adapt this work
# PROVIDED THAT you attribute the work to the authors listed below.
# For more information, please see the following web page:
# http://creativecommons.org/licenses/by/3.0/
#
# This file is a component of the Sleipnir library for functional genomics,
# authored by:
# Curtis Huttenhower (chuttenh@princeton.edu)
# Mark Schroeder
# Maria Chikina
# Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
#
# If you use this library, the included executable tools, or any related
# code in your work, please cite the following publication:
# Curtis Huttenhower, Mark Schroeder, Maria Chikina, and
# Olga G. Troyanskaya.
# "The Sleipnir library for computational functional genomics"

#! TODO:
#! b. Warn for pthread.
#! d. Gengetopt checks.
#! e. Improve boost logic?

AC_PREREQ(2.59)
AC_INIT([libSleipnir], [1.0.0], [ogt@princeton.edu])

AC_COPYRIGHT([This software is copyrighted.])

# Any file under src will do.
AC_CONFIG_SRCDIR([src/stdafx.h])
# doesn't work with automake?
#AC_CONFIG_HEADER([config.h])
#AM_CONFIG_HEADER([config.h])
m4_pattern_allow([AM_INIT_AUTOMAKE])
AM_INIT_AUTOMAKE

#SLEIPNIR_HOST_TYPE
AC_LANG([C++])

# why oh why does this break things?
#AC_PROG_EGREP
EGREP=egrep

# some don't have this??
#AC_PROG_GREP
GREP=$EGREP

topleveldir=$PWD

# make the path absolute
#!! this may not be portable
qualify_path () {
  d=${!1}
  if ! echo $d | $EGREP -q '^/'; then
    eval "$1=$topleveldir/$d"
  fi
}


# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

# Checks for header files.
AC_HEADER_DIRENT
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h netinet/in.h sys/socket.h sys/time.h pthread.h])

AC_CHECK_HEADERS([arpa/inet.h])   # tools/BNServer

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_FUNC_STRTOD
AC_FUNC_VPRINTF
AC_FUNC_STAT   # tools/Answerer, tools/MEFIT

AC_CHECK_FUNCS([memmove memset munmap pow socket sqrt strchr strerror strrchr strstr])

AC_CHECK_FUNCS([inet_ntoa])   # tools/BNServer
AC_CHECK_FUNCS([strtol])  # tools/Data2Features, tools/Data2Bnt


# Common AC functions (see www.gnu.org/autoconf/manual)
#
# AC_DEFINE(var[,val[,desc]])
# Specify a C preprocessor definition. Val defaults to 1.
#
# AC_SUBST(var)
# Output variable (to automake).
#
# AC_ARG_WITH(n,help[,given[,not-given]])
# For "with" package n, use this help string (usually formed with AS_HELP_STRING), perform "given" if given, "not-given" if not.
#
# AC_PATH_PROG(v,p[,nf[,path]])
# Set $v to absolute name of p if p is found. Set $v to "nf" if not found. Use PATH=$path if specified.
#
# AC_LINK_IFELSE(in[,true[,false]])
# Run the complier and linker on "in", and execute "true" or "false" depending on success.
#
# AC_LANG_PROGRAM(prologue,body)
# Create source file using "prologue" and then "body" with main function. Usually used in an AC_LINK_IFELSE or similar.
#
# AC_CHECK_HEADER(h[,found[,not-found]])
# Check for header. Do "found" if found, "not-found" if not.


## LOCAL_PATH_PROG(prog name, local name,
##   action if found, extra action if found installed, extra action if found in with path,
##   action if not found, action if user requested without)
##
## Set prog_(local name)_dir to dir containing it, prog_(local name) to the path, prog_(local name)_info.
## Also set prog_(local name)_state to "ok" or "no".
define([LOCAL_PATH_PROG],
[
if test "x$with_$2" = "xno"; then  # user said no
  prog_$2_info=NONE
  prog_$2_state=no
  $7
elif test "x$with_$2" != "xcheck"; then
#  if echo with_$2 | $EGREP -q '^extlib/'; then extlib=true; else extlib=false; fi
  prog_$2_info=$with_$2
  prog_$2_state=ok
#  qualify_path with_$2
  prog_$2=$with_$2
  $3
  $5
else
  AC_PATH_PROG([prog_$2], [$1])
  if test "x$prog_$2" = "x"; then
    prog_$2_info="not found"
    prog_$2_state=no
    $6
  else
#    prog_$2_info="found installed"
    prog_$2_info=$prog_$2
    prog_$2_state=ok
    $3
    $4
  fi
fi
if test "x$prog_$2" != "x"; then
  prog_$2_dir=`AS_DIRNAME([$prog_$2])`
  prog_$2_prefix=$prog_$2_dir/..     # assumes bin and include are at same level
#  $4
fi
])


## LOCAL_CHECK_LIB(1 lib name, 2 local name, 3 function,
##   4 action if found, 5 extra action if found installed, 6 extra action if found in 'with' path,
##   7 action if not found, 8 action if not found installed, 9 action if not found installed in 'with' path,
##   10 action if requested without,
##   11 headers, 12 extra libs
##   13 extra cflags for testing 'with' path)
## Set $(local_name)_info to "NO", location, or "found installed".
define([LOCAL_CHECK_LIB],
[
if test "x${$2_state}" = "xno"; then
  $2_info=NO
  $10
elif test "x${$2_state}" = "xwarn"; then
  $2_info="INVALID PATH"
  $7
  $9
else
  found=true
  oLIBS=$LIBS
  LIBS="$LIBS -l$1 $12"
  if test "x${$2_state}" = "xtry"; then
# CFLAGS not respected here - how are we supposed to do this?
#    oCFLAGS=$CFLAGS
#    export CFLAGS="$CFLAGS $13"
    LIBS="$13 $LIBS"
  fi
  AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[$11]], [[$3]])],
    [
     if test "x${$2_state}" = "xtry"; then
       $2_info=${with_$2}
#       qualify_path with_$2
       $4
       $6
     else
       $2_info="found installed"
       $4
       $5
     fi
    ],
    [found=false])
  LIBS=$oLIBS
#  if test "x${$2_state}" = "xtry"; then
#    CFLAGS=$oCFLAGS
#  fi
  if ! $found; then
    $2_info="not found"
    if test "x${$2_state}" = "xtry"; then
      $7
      $9
      true   # in case 7 and 9 are empty
    else
      $7
      $8
      true   # in case 7 and 8 are empty
    fi
  fi
fi
])


## LOCAL_CHECK_APPEND_PATHS(base path, paths to append, common path to look for, var to set,
##   action if not found)
## Look for the common path in one of the appendable paths under the base path.
## If found, set the variable to the concatenation of the prefix and the appended path,
## not including the common path. If not found use the base dir.
define([LOCAL_CHECK_APPEND_PATHS],
[
#result=$1
result=
found=false
echo -n "looking for $3... "
for p in $2; do
  if test -e "$1/$p/$3"; then
    result=$1/$p
    echo found $1/$p/$3
    found=true
    break
  fi
done
$4=$result
if ! $found; then
#  if ! test -e $result/$3; then
    echo not found
#    true   # in case the following is empty
    $5
#  fi
fi
])


## GENGETOPT
# If the user specifies a path, use it.
# If the user doesn't specify anything, check for it.
# Gengetopt is required to update cmdline.c and cmdline.h (in tools) if out of date.
#! CFLAGS assumes headers in same directory as binaries for installed
AC_ARG_WITH([gengetopt],
        [AS_HELP_STRING([--with-gengetopt=PATH], [prefix of gengetopt installation])],
        [
          if test "x$with_gengetopt" = "xyes"; then
            with_gengetopt=check
          elif test "x$with_gengetopt" = "xno"; then
            with_gengetopt = no
          else
            qualify_path with_gengetopt
          fi
        ],
        [with_gengetopt=check])
LOCAL_PATH_PROG([gengetopt], [gengetopt],
        [GENGETOPT_PREFIX=$prog_gengetopt_prefix
         GENGETOPT_INCLUDES="-I${GENGETOPT_PREFIX}/include"],  dnl found
        [],                                                    dnl found installed
        [],                                                    dnl found in with path
        [],                                                    dnl not found
        [])                                                    dnl requested without
## set CFLAGS and LIBS
if test "x$GENGETOPT_PREFIX" != "x"; then
  #! check (in LOCAL_PATH_PROG?)
  GENGETOPT_CFLAGS=$GENGETOPT_INCLUDES
  GENGETOPT_LIBS=
fi
AC_SUBST(prog_gengetopt)
AC_SUBST(prog_gengetopt_dir)
#AC_SUBST(GENGETOPT_PREFIX)
AC_SUBST(GENGETOPT_CFLAGS)
AC_SUBST(GENGETOPT_LIBS)


## READLINE
# If the user specifies a path, use it.
# If not specified, check for it.
# If no readline, do not build OntoShell.
AC_ARG_WITH([readline],
        [AS_HELP_STRING([--with-readline=PATH], [prefix of readline installation])],
        [
         if test "x$with_readline" = "xyes"; then
           readline_state=check
         elif test "x$with_readline" != "xno"; then
           readline_state=try
           qualify_path with_readline
           LOCAL_CHECK_APPEND_PATHS([$with_readline], [. include], [readline/readline.h], [READLINE_INCLUDE_DIR], [readline_state=warn])
           LOCAL_CHECK_APPEND_PATHS([$with_readline], [. lib], [libreadline.a], [READLINE_LIB_DIR], [readline_state=warn])
         else
           readline_state=no
         fi
        ],
        [readline_state=check])
AC_LANG_PUSH([C])
LOCAL_CHECK_LIB([readline], [readline], [readline("foo")],
        [readline_state=ok],                                       dnl found
        [READLINE_LIBS="-lreadline -lncurses"],                    dnl and found installed
        [
          READLINE_CFLAGS="-I $READLINE_INCLUDE_DIR"
          READLINE_LIBS="-L $READLINE_LIB_DIR -lreadline -lncurses"
        ],                                                         dnl and found in specified path
        [],                                                        dnl not found
        [readline_state=no],                                       dnl and not found installed
        [readline_state=warn],                                     dnl and not found in specified path
        [],                                                        dnl requested without
        [#include <stdio.h>
         #include "readline/readline.h"],
        [-lncurses],
        [-L $READLINE_LIB_DIR -I $READLINE_INCLUDE_DIR])
AC_LANG_POP
if test "x$readline_state" != "xok"; then
  AC_DEFINE([NO_READLINE], [1])
fi
#  # apparently on Macs, rl_completion_matches() is named completion_matches()
#  # hopefully this is an unusual exception
#  AC_LINK_IFELSE(
#    [AC_LANG_PROGRAM([#include <stdio.h>
#                     #include "readline/readline.h"],
#                     [rl_completion_matches])],
#    [],
#    [AC_DEFINE([rl_completion_matches], [completion_matches])])
#AC_SUBST(READLINE_PREFIX)
AC_SUBST(READLINE_CFLAGS)
AC_SUBST(READLINE_LIBS)


## PTHREAD
# If the user specifies a path, use it.
# If not specified, check for it.
# If no readline, do not build OntoShell.
AC_ARG_WITH([pthread],
        [AS_HELP_STRING([--with-pthread=PATH], [prefix of pthread installation])],
        [
         if test "x$with_pthread" = "xyes"; then
           pthread_state=check
         elif test "x$with_pthread" != "xno"; then
           pthread_state=try
           qualify_path with_pthread
           LOCAL_CHECK_APPEND_PATHS([$with_pthread], [. include], [pthread.h], [PTHREAD_INCLUDE_DIR], [pthread_state=warn])
           LOCAL_CHECK_APPEND_PATHS([$with_pthread], [. lib], [libpthread.a], [PTHREAD_LIB_DIR], [pthread_state=warn])
         else
           pthread_state=no
         fi
        ],
        [pthread_state=check])
AC_LANG_PUSH([C])
LOCAL_CHECK_LIB([pthread], [pthread], [pthread_create],
        [pthread_state=ok],                                        dnl found
        [PTHREAD_LIBS="-lpthread"],                                dnl and found installed
        [
          PTHREAD_CFLAGS="-I $PTHREAD_INCLUDE_DIR"
          PTHREAD_LIBS="-L $PTHREAD_LIB_DIR -lpthread"
        ],                                                         dnl and found in specified path
        [],                                                        dnl not found
        [pthread_state=no],                                        dnl and not found installed
        [pthread_state=warn],                                      dnl and not found in specified path
        [],                                                        dnl requested without
        [#include <stdio.h>
         #include "pthread.h"],
        [],
        [-L $PTHREAD_LIB_DIR -I $PTHREAD_INCLUDE_DIR])
AC_LANG_POP
if test "x$pthread_state" != "xok"; then
  AC_DEFINE([NO_PTHREAD], [1])
fi
#AC_SUBST(PTHREAD_PREFIX)
AC_SUBST(PTHREAD_CFLAGS)
AC_SUBST(PTHREAD_LIBS)


## LOG4CPP
# If the user specifies a path, use it.
# If the user says not to use it, use the stub.
# If the user doesn't say anything, check for it. Use the stub if not found.
# While checking, look for a source distribution in extlib/.
# If using stub, set USE_LOG4CPP_STUB=1.
if test "x$pthread_state" != "xok"; then
  log4cpp_info=stub
  log4cpp_state=no
else
AC_ARG_WITH([log4cpp],
        [AS_HELP_STRING([--with-log4cpp=PATH], [prefix of log4cpp installation or build directory])],
        [
         if test "x$with_log4cpp" = "xyes"; then
           log4cpp_state=check
         elif test "x$with_log4cpp" != "xno"; then
           log4cpp_state=try
           qualify_path with_log4cpp
           LOCAL_CHECK_APPEND_PATHS([$with_log4cpp], [. include], [log4cpp/Category.hh], [LOG4CPP_INCLUDE_DIR], [log4cpp_state=warn])
           LOCAL_CHECK_APPEND_PATHS([$with_log4cpp], [. lib src/.libs], [liblog4cpp.a], [LOG4CPP_LIB_DIR], [log4cpp_state=warn])
         else
           log4cpp_state=no
         fi
        ],
        [log4cpp_state=check])
LOCAL_CHECK_LIB([log4cpp], [log4cpp], [log4cpp::Category::getInstance("foo")],
        [log4cpp_state=ok],                                       dnl found
        [LOG4CPP_LIBS="-llog4cpp"],                               dnl and found installed
        [
         LOG4CPP_CFLAGS="-I $LOG4CPP_INCLUDE_DIR"
         LOG4CPP_LIBS="-L $LOG4CPP_LIB_DIR -llog4cpp"
        ],                                                        dnl and found in specified path
        [],                                                       dnl not found
        [log4cpp_state=no],                                       dnl and not found installed
        [log4cpp_state=warn],                                     dnl and not found in specified path
        [],                                                       dnl requested without
        [#include "log4cpp/Category.hh"],
        [-lpthread],
        [-L $LOG4CPP_LIB_DIR -I $LOG4CPP_INCLUDE_DIR])
fi
if test "x$log4cpp_state" != "xok"; then
  log4cpp_info=stub
  AC_DEFINE([USE_LOG4CPP_STUB], [1])
fi
#AC_SUBST(LOG4CPP_PREFIX)
AC_SUBST(LOG4CPP_CFLAGS)
AC_SUBST(LOG4CPP_LIBS)

dnl QIAN ==============================
USE_OPENMP=true
sleipnir_cv_openmp=yes
AC_ARG_ENABLE(openmp,
        [AS_HELP_STRING([--disable-openmp],
        [disable OpenMP SMP multi-threading])],
        [if test "$enableval" = no; then sleipnir_cv_openmp=no; USE_OPENMP=false; fi], [])
if $USE_OPENMP; then
        dnl case "$canonical_host_type" in
        dnl	*-*-mingw* | *-*-cygwin* | *-*-gnu*)
                AC_MSG_CHECKING([for support of OpenMP])
                XCFLAGS="$CFLAGS"
                CFLAGS="$CFLAGS -fopenmp"
                AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
                        #include <omp.h>
                        #include <stdio.h>
                ]], [[
                          _Pragma("omp parallel")
                          printf("Hello, world.\n");
                          return 0;
                        ]])],sleipnir_cv_openmp=yes, sleipnir_cv_openmp=no)
                dnl AC_MSG_FAILURE([Error has occurred])
                AC_MSG_RESULT([$sleipnir_cv_openmp])
                        if test "$sleipnir_cv_openmp" = yes; then
                          AC_DEFINE(HAVE_OPENMP,1,[Define if compiler supports OpenMP])
                          CXXFLAGS="$CXXFLAGS -fopenmp"
                          USE_OPENMP=true
                        else
                          CFLAGS="$XCFLAGS"
                          USE_OPENMP=false
                        fi
                        dnl SLEIPNIR_CHECK_OPENMP(-fopenmp)
        dnl	;;
        dnl	*-*-msdosmsvc)
        dnl	## FIXME is this the right flag for MSVC?
        dnl	SLEIPNIR_CHECK_OPENMP(-openmp)
        dnl	;;
        dnl	## Add other compilers supporting OpenMP here
        dnl #esac
fi

#
#if $USE_OPENMP; then
#  AC_MSG_WARN([])
#  AC_MSG_WARN([You used the EXPERIMENTAL --enable-openmp option.])
#  warn_msg_printed=true
#fi




#AC_ARG_WITH([omp],
#	[AS_HELP_STRING([--with-omp], [Enable OpenMP directives])],
#	[
#		if test "x$with_omp" == "xno"; then
#			omp_state=no
#		else
#			omp_state=check
#		fi
#	],
#	[omp_state=check])


## SMILE
# If the user specifies a path, use it.
# If the user says not to use it, then don't.
# If the user doesn't specify anything, check for it.
# Also check for libsmilexml.a and define SMILEXML_LIB if found (not used?)
# If none, set NO_SMILE=1.
AC_ARG_WITH([smile],
        [AS_HELP_STRING([--with-smile=PATH], [prefix of SMILE installation])],
        [
        if test "x$with_smile" = "xyes"; then
            smile_state=check
        elif test "x$with_smile" != "xno"; then
            smile_state=try
                qualify_path with_smile
            LOCAL_CHECK_APPEND_PATHS([$with_smile], [. include smile include/smile], [smile.h], [SMILE_INCLUDE_DIR], [smile_state=warn])
            LOCAL_CHECK_APPEND_PATHS([$with_smile], [. lib smile], [libsmile.a], [SMILE_LIB_DIR], [smile_state=warn])
        else
            smile_state=no
        fi
    ],
        [smile_state=check])
        LOCAL_CHECK_LIB([smile], [smile], [DSL_isFinite(1.0)],

    [smile_state=ok],                                         dnl found
    [
            SMILE_LIBS="-lsmile"
        #SMILEXML_LIB
        ],                                                       dnl and found installed
        [
#    if test -f "$SMILE_LIB_DIR/libsmilexml.a"; then
#      AC_DEFINE([SMILEXML_LIB], [1])
#      SMILEXML_LIB="-lsmilexml"
#    fi
         SMILE_CFLAGS="-I $SMILE_INCLUDE_DIR"
#         SMILE_LIBS="-L $SMILE_LIB_DIR $SMILEXML_LIB -lsmile"
         SMILE_LIBS="-L $SMILE_LIB_DIR -lsmile"
        ],                                                        dnl and found in specified path
        [],                                                       dnl not found
    [smile_state=no],                                         dnl and not found installed
    [smile_state=warn],                                       dnl and not found in specified path
    [],                                                       dnl requested without
        [#include "smile.h"],
    [],
    [-L $SMILE_LIB_DIR -I $SMILE_INCLUDE_DIR])
if test "x$smile_state" != "xok"; then
    AC_DEFINE([NO_SMILE], [1])
fi
#AC_SUBST(SMILE_PREFIX)
AC_SUBST(SMILE_CFLAGS)
AC_SUBST(SMILE_LIBS)


## LibSVM
# If the user specifies a path, use it.
# If the user says not to use it, then don't.
# Does not check for it.
# If none, set NO_LIBSVM = 1.
AC_ARG_WITH([libsvm],
        [AS_HELP_STRING([--with-libsvm=PATH], [prefix of LibSVM installation])],
        [
         if test "x$with_libsvm" = "xyes"; then
           libsvm_state=check
         elif test "x$with_libsvm" != "xno"; then
           libsvm_state=try
           qualify_path with_libsvm
           LOCAL_CHECK_APPEND_PATHS([$with_libsvm], [. include], [libsvm/svm.h], [LIBSVM_INCLUDE_DIR], [libsvm_state=warn])
           LOCAL_CHECK_APPEND_PATHS([$with_libsvm], [. lib], [liblibsvm.a], [LIBSVM_LIB_DIR], [libsvm_state=warn])
         else
           libsvm_state=no
         fi
        ],
        [libsvm_state=check])
#AC_LANG_PUSH([C])
LOCAL_CHECK_LIB([libsvm], [libsvm], [ struct svm_node *x_space ],
        [libsvm_state=ok],                                      dnl found
        [LIBSVM_LIBS="-llibsvm"],                              dnl and found installed
        [
         LIBSVM_CFLAGS="-I $LIBSVM_INCLUDE_DIR"
         LIBSVM_LIBS="-L $LIBSVM_LIB_DIR -llibsvm"
        ],                                                        dnl and found in specified path
        [],                                                       dnl not found
        [libsvm_state=no],                                      dnl and not found installed
        [libsvm_state=warn],                                    dnl and not found in specified path
        [],                                                       dnl requested without
        [
         #include <sys/types.h>
         #include "libsvm/svm.h"],
        [],
        [-L $LIBSVM_LIB_DIR -I $LIBSVM_INCLUDE_DIR])
## LOCAL_CHECK_LIB(1 lib name, 2 local name, 3 function,
##   4 action if found,
##   5 extra action if found installed,
##   6 extra action if found in 'with' path,
##   7 action if not found,
##   8 action if not found installed,
##   9 action if not found installed in 'with' path,
##   10 action if requested without,
##   11 headers, 12 extra libs
##   13 extra cflags for testing 'with' path)
## Set $(local_name)_info to "NO", location, or "found installed".

#AC_LANG_POP
if test "x$libsvm_state" != "xok"; then
  AC_DEFINE([NO_LIBSVM], [1])
fi

AC_SUBST(LIBSVM_CFLAGS)
AC_SUBST(LIBSVM_LIBS)


## SVM PERF
# If the user specifies a path, use it.
# If the user says not to use it, then don't.
# If the user doesn't specify anything, check for it.
# While checking, look for a source distribution in extlib/.
# If none, set NO_SVM_PERF=1.
AC_ARG_WITH([svm_perf],
        [AS_HELP_STRING([--with-svm-perf=PATH], [prefix of SVM Perf installation])],
        [
         if test "x$with_svm_perf" = "xyes"; then
           svm_perf_state=check
         elif test "x$with_svm_perf" != "xno"; then
           svm_perf_state=try
           qualify_path with_svm_perf
           LOCAL_CHECK_APPEND_PATHS([$with_svm_perf], [. include], [svm_struct_api.h], [SVM_PERF_INCLUDE_DIR], [svm_perf_state=warn])
           LOCAL_CHECK_APPEND_PATHS([$with_svm_perf], [. lib], [libsvmperf.a], [SVM_PERF_LIB_DIR], [svm_perf_state=warn])
         else
           svm_perf_state=no
         fi
        ],
        [svm_perf_state=check])
AC_LANG_PUSH([C])
LOCAL_CHECK_LIB([svmperf], [svm_perf], [optimize_svm],
        [svm_perf_state=ok],                                      dnl found
        [SVM_PERF_LIBS="-lsvmperf"],                              dnl and found installed
        [
         SVM_PERF_CFLAGS="-I $SVM_PERF_INCLUDE_DIR"
         SVM_PERF_LIBS="-L $SVM_PERF_LIB_DIR -lsvmperf"
        ],                                                        dnl and found in specified path
        [],                                                       dnl not found
        [svm_perf_state=no],                                      dnl and not found installed
        [svm_perf_state=warn],                                    dnl and not found in specified path
        [],                                                       dnl requested without
        [#include <sys/types.h>
         #include "svm_light/svm_common.h",
         #include "svm_light/svm_learn.h"],
        [],
        [-L $SVM_PERF_LIB_DIR -I $SVM_PERF_INCLUDE_DIR])
AC_LANG_POP
if test "x$svm_perf_state" != "xok"; then
  AC_DEFINE([NO_SVM_PERF], [1])
fi
#AC_SUBST(SVM_PERF_PREFIX)
AC_SUBST(SVM_PERF_CFLAGS)
AC_SUBST(SVM_PERF_LIBS)
#        #include "svm/svm_struct_api.h"
#        #include "svm/svm_struct_api_types.h"],

## SVM MULTICLASS
# If the user specifies a path, use it.
# If the user says not to use it, then don't.
# If the user doesn't specify anything, check for it.
# While checking, look for a source distribution in extlib/.
# If none, set NO_SVM_MULTICLASS=1.
AC_ARG_WITH([svm_multiclass],
  [AS_HELP_STRING([--with-svm-multiclass=PATH], [prefix of SVM Multiclass installation])],
  [
         if test "x$with_svm_multiclass" = "xyes"; then
     svm_multiclass_state=check
         elif test "x$with_svm_multiclass" != "xno"; then
           svm_multiclass_state=try
     qualify_path with_svm_multiclass
           LOCAL_CHECK_APPEND_PATHS([$with_svm_multiclass], [. include], [svm_multiclass/svm_struct_api.h], [SVM_MULTICLASS_INCLUDE_DIR], [svm_multiclass_state=warn])
           LOCAL_CHECK_APPEND_PATHS([$with_svm_multiclass], [. lib], [svm_multiclass/libsvmmulticlass.a], [SVM_MULTICLASS_LIB_DIR], [svm_multiclass_state=warn])
         else
           svm_multiclass_state=no
         fi
        ],
  [svm_multiclass_state=check])
AC_LANG_PUSH([C])
LOCAL_CHECK_LIB([svmmulticlass], [svm_multiclass], [optimize_svm],
        [svm_multiclass_state=ok],                                      dnl found
        [SVM_MULTICLASS_LIBS="-lsvmmulticlass"],                              dnl and found installed
        [
         SVM_MULTICLASS_CFLAGS="-I $SVM_MULTICLASS_INCLUDE_DIR"
         SVM_MULTICLASS_LIBS="-L $SVM_MULTICLASS_LIB_DIR/svm_multiclass -lsvmmulticlass"
        ],                                                        dnl and found in specified path
  [],                                                       dnl not found
        [svm_multiclass_state=no],                                      dnl and not found installed
        [svm_multiclass_state=warn],                                    dnl and not found in specified path
        [],                                                       dnl requested without
  [#include <sys/types.h>
   #include "svm_light/svm_common.h",
   #include "svm_light/svm_learn.h"],
        [],
        [-L $SVM_MULTICLASS_LIB_DIR/svm_multiclass -I $SVM_MULTICLASS_INCLUDE_DIR/svm_multiclass])
AC_LANG_POP
if test "x$svm_multiclass_state" != "xok"; then
  AC_DEFINE([NO_SVM_MULTICLASS], [1])
fi
#AC_SUBST(SVM_MULTICLASS_PREFIX)
AC_SUBST(SVM_MULTICLASS_CFLAGS)
AC_SUBST(SVM_MULTICLASS_LIBS)
#  #include "svm/svm_struct_api.h"
#  #include "svm/svm_struct_api_types.h"],


## SVM HIERARCHY
# If the user specifies a path, use it.
# If the user says not to use it, then don't.
# If the user doesn't specify anything, check for it.
# While checking, look for a source distribution in extlib/.
# If none, set NO_SVM_HIERARCHY=1.
AC_ARG_WITH([svm_hierarchy],
  [AS_HELP_STRING([--with-svm-hierarchy=PATH], [prefix of SVM Hierarchy installation])],
  [
         if test "x$with_svm_hierarchy" = "xyes"; then
     svm_hierarchy_state=check
         elif test "x$with_svm_hierarchy" != "xno"; then
           svm_hierarchy_state=try
     qualify_path with_svm_hierarchy
           LOCAL_CHECK_APPEND_PATHS([$with_svm_hierarchy], [. include], [svm_hierarchy/svm_struct_api.h], [SVM_HIERARCHY_INCLUDE_DIR], [svm_hierarchy_state=warn])
           LOCAL_CHECK_APPEND_PATHS([$with_svm_hierarchy], [. lib], [svm_hierarchy/libsvmhierarchy.a], [SVM_HIERARCHY_LIB_DIR], [svm_hierarchy_state=warn])
         else
           svm_hierarchy_state=no
         fi
        ],
  [svm_hierarchy_state=check])
AC_LANG_PUSH([C])
LOCAL_CHECK_LIB([svmhierarchy], [svm_hierarchy], [optimize_svm],
        [svm_hierarchy_state=ok],                                      dnl found
        [SVM_HIERARCHY_LIBS="-lsvmhierarchy"],                              dnl and found installed
        [
         SVM_HIERARCHY_CFLAGS="-I $SVM_HIERARCHY_INCLUDE_DIR"
         SVM_HIERARCHY_LIBS="-L $SVM_HIERARCHY_LIB_DIR/svm_hierarchy -lsvmhierarchy"
        ],                                                        dnl and found in specified path
  [],                                                       dnl not found
        [svm_hierarchy_state=no],                                      dnl and not found installed
        [svm_hierarchy_state=warn],                                    dnl and not found in specified path
        [],                                                       dnl requested without
  [#include <sys/types.h>
   #include "svm_light/svm_common.h",
   #include "svm_light/svm_learn.h"],
        [],
        [-L $SVM_HIERARCHY_LIB_DIR/svm_hierarchy -I $SVM_HIERARCHY_INCLUDE_DIR/svm_hierarchy])
AC_LANG_POP
if test "x$svm_hierarchy_state" != "xok"; then
  AC_DEFINE([NO_SVM_HIERARCHY], [1])
fi
#AC_SUBST(SVM_HIERARCHY_PREFIX)
AC_SUBST(SVM_HIERARCHY_CFLAGS)
AC_SUBST(SVM_HIERARCHY_LIBS)
#  #include "svm/svm_struct_api.h"
#  #include "svm/svm_struct_api_types.h"],


## BOOST
# If the user specifies a path, use it.
# If not specified, check for it.
# BOOST graph lib is required for BNServer.
AC_ARG_WITH([boost],
        [AS_HELP_STRING([--with-boost=PATH], [prefix of Boost installation])],
        [
         if test "x$with_boost" = "xyes"; then
           boost_state=check
         elif test "x$with_boost" != "xno"; then
           boost_state=try
           qualify_path with_boost
           LOCAL_CHECK_APPEND_PATHS([$with_boost], [. .. include], [boost/graph/graphviz.hpp], [BOOST_INCLUDE_DIR], [boost_state=warn])
           LOCAL_CHECK_APPEND_PATHS([$with_boost], [. lib stage/lib], [libboost_graph.a], [BOOST_GRAPH_LIB_DIR], [boost_state=warn])
           LOCAL_CHECK_APPEND_PATHS([$with_boost], [. lib stage/lib], [libboost_regex.a], [BOOST_REGEX_LIB_DIR], [boost_state=warn])
           if test "x$BOOST_GRAPH_LIB_DIR" != "x$BOOST_REGEX_LIB_DIR"; then
             boost_state=warn
           else
             BOOST_LIB_DIR=$BOOST_GRAPH_LIB_DIR
           fi
         else
           boost_state=no
         fi
        ],
        [boost_state=check])
LOCAL_CHECK_LIB([boost_graph], [boost], [],
        [boost_state=ok],                                         dnl found
        [],                                                       dnl and found installed
        [],                                                       dnl and found in specified path
        [],                                                       dnl not found
        [boost_state=no],                                         dnl and not found installed
        [boost_state=warn],                                       dnl and not found in specified path
        [],                                                       dnl requested without
        [#include "boost/graph/graphviz.hpp"],
        [],
        [-L $BOOST_LIB_DIR -I $BOOST_INCLUDE_DIR])
if test "x$boost_state" = "xok"; then
LOCAL_CHECK_LIB([boost_regex], [boost], [],
        [boost_state=ok],                                         dnl found
        [BOOST_LIBS="-lboost_graph -lboost_regex"],               dnl and found installed
        [
         BOOST_CFLAGS="-I $BOOST_INCLUDE_DIR"
         BOOST_LIBS="-L $BOOST_LIB_DIR -lboost_graph -lboost_regex"
        ],                                                        dnl and found in specified path
        [],                                                       dnl not found
        [boost_state=no],                                         dnl and not found installed
        [boost_state=warn],                                       dnl and not found in specified path
        [],                                                       dnl requested without
        [#include "boost/graph/graphviz.hpp"],
        [],
        [-L $BOOST_LIB_DIR -I $BOOST_INCLUDE_DIR])
fi
if test "x$boost_state" != "xok"; then
  AC_DEFINE([NO_BOOST], [1])
fi
#AC_SUBST(BOOST_PREFIX)
AC_SUBST(BOOST_CFLAGS)
AC_SUBST(BOOST_LIBS)

#AC_CHECK_LIB([m],[cos])
#AC_CHECK_LIB([gslcblas],[cblas_dgemm])
#AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

## GSL library

AC_ARG_WITH([gsl],
        [AS_HELP_STRING([--with-gsl=PATH], [prefix of GSL installation or build dir])],
        [
         if test "x$with_gsl" = "xyes"; then
                        gsl_state=check
         elif test "x$with_gsl" != "xno"; then
                gsl_state=try
           qualify_path with_gsl
           LOCAL_CHECK_APPEND_PATHS([$with_gsl], [. include], [gsl/gsl_blas.h], [GSL_INCLUDE_DIR], [gsl_state=warn])
           LOCAL_CHECK_APPEND_PATHS([$with_gsl], [. lib], [libgsl.a], [GSL_LIB_DIR], [gsl_state=warn])
         else
           gsl_state=no
         fi
        ],
        [gsl_state=check])
LOCAL_CHECK_LIB([gsl], [gsl], [gsl_rng_env_setup()],
        [gsl_state=ok],                                 dnl found
        [GSL_LIBS="-lgsl -lgslcblas"],                              dnl and found installed
        [
         GSL_CFLAGS="-I $GSL_INCLUDE_DIR"
         GSL_LIBS="-L $GSL_LIB_DIR -lgsl -lgslcblas"
        ],                                                        dnl and found in specified path
        [],                                                       dnl not found
        [gsl_state=no],                                 dnl and not found installed
        [gsl_state=warn],                               dnl and not found in specified path
        [],                                                       dnl requested without
        [#include <stdio.h>
         #include <gsl/gsl_rng.h>],
        [-lgsl -lgslcblas],
        [-L $GSL_LIB_DIR -I $GSL_INCLUDE_DIR])
if test "x$gsl_state" != "xok"; then
  AC_DEFINE([NO_GSL], [1])
fi
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)



## Vowpal Wabbit
# If the user specifies a path, use it.
# If the user says not to use it, then don't.
# If the user doesn't specify anything, check for it.
# While checking, look for a source distribution in extlib/.
# If none, set NO_VOWPAL_WABBIT=1.
if test "x$pthread_state" != "xok" -o "x$boost_state" != "xok"; then
  vowpal_wabbit_info=no
  vowpal_wabbit_state=no
else
AC_ARG_WITH([vowpal_wabbit],
        [AS_HELP_STRING([--with-vowpal-wabbit=PATH], [prefix of Vowpal Wabbit installation or build dir])],
        [
         if test "x$with_vowpal_wabbit" = "xyes"; then
           vowpal_wabbit_state=check
         elif test "x$with_vowpal_wabbit" != "xno"; then
           vowpal_wabbit_state=try
           qualify_path with_vowpal_wabbit
           LOCAL_CHECK_APPEND_PATHS([$with_vowpal_wabbit], [. include], [vw.h], [VOWPAL_WABBIT_INCLUDE_DIR], [vowpal_wabbit_state=warn])
           LOCAL_CHECK_APPEND_PATHS([$with_vowpal_wabbit], [. lib], [libvw.a], [VOWPAL_WABBIT_LIB_DIR], [vowpal_wabbit_state=warn])
         else
           vowpal_wabbit_state=no
         fi
        ],
        [vowpal_wabbit_state=check])
LOCAL_CHECK_LIB([vw], [vowpal_wabbit], [start_noop()],
        [vowpal_wabbit_state=ok],                                 dnl found
        [VOWPAL_WABBIT_LIBS="-lvw"],                              dnl and found installed
        [
         VOWPAL_WABBIT_CFLAGS="-I $VOWPAL_WABBIT_INCLUDE_DIR"
         VOWPAL_WABBIT_LIBS="-L $VOWPAL_WABBIT_LIB_DIR -lvw"
        ],                                                        dnl and found in specified path
        [],                                                       dnl not found
        [vowpal_wabbit_state=no],                                 dnl and not found installed
        [vowpal_wabbit_state=warn],                               dnl and not found in specified path
        [],                                                       dnl requested without
        [#include "boost/program_options.hpp"
         #include "zlib.h"
         #include "noop.h"],
        [-lboost_program_options -lpthread -lz],
        [-L $VOWPAL_WABBIT_LIB_DIR -I $VOWPAL_WABBIT_INCLUDE_DIR])
fi
if test "x$vowpal_wabbit_state" != "xok"; then
  AC_DEFINE([NO_VOWPAL_WABBIT], [1])
fi
#AC_SUBST(VOWPAL_WABBIT_PREFIX)
AC_SUBST(VOWPAL_WABBIT_CFLAGS)
AC_SUBST(VOWPAL_WABBIT_LIBS)


## Now define variables to be used to reference this library.
sleipnirdocdir=$datadir/doc/
AC_SUBST(sleipnirdocdir)
SLEIPNIR_PREFIX=src
qualify_path SLEIPNIR_PREFIX
SLEIPNIR_CFLAGS="-I${SLEIPNIR_PREFIX}"
SLEIPNIR_LIBS="-L${SLEIPNIR_PREFIX} -lSleipnir"
AC_SUBST(SLEIPNIR_PREFIX)
AC_SUBST(SLEIPNIR_CFLAGS)
AC_SUBST(SLEIPNIR_LIBS)


# uncomment to print out all variables
#set


AC_ARG_WITH([tools],
        [AS_HELP_STRING([--without-tools], [exclude tools and build lib only])],
        [if test "x$with_tools" = "xno"; then with_tools=false; else with_tools=true; fi],
        [with_tools=true])


AM_CONDITIONAL([WITH_TOOLS], [$with_tools])
AM_CONDITIONAL([WITH_SMILE_TOOLS], [test "x$smile_state" = "xok"])
AM_CONDITIONAL([WITH_SVM_TOOLS], [test "x$svm_perf_state" = "xok"])
AM_CONDITIONAL([WITH_SVM_MULTICLASS_TOOLS], [test "x$svm_multiclass_state" = "xok"])
AM_CONDITIONAL([WITH_SVM_HIERARCHY_TOOLS], [test "x$svm_hierarchy_state" = "xok"])
AM_CONDITIONAL([WITH_SVM_LIBSRC], [test "x$svm_perf_state" = "xok"])
AM_CONDITIONAL([WITH_VW_TOOLS], [test "x$vowpal_wabbit_state" = "xok"])
AM_CONDITIONAL([WITH_READLINE_TOOLS], [test "x$readline_state" = "xok"])
AM_CONDITIONAL([WITH_BOOST_TOOLS], [test "x$boost_state" = "xok"])
AM_CONDITIONAL([WITH_GSL_TOOLS], [test "x$gsl_state" = "xok"])
AM_CONDITIONAL([WITH_GENGETOPT], [test "x$prog_gengetopt_state" = "xok"])
AM_CONDITIONAL([WITH_LIBSVM_TOOLS], [test "x$libsvm_state" = "xok"])#
AM_CONDITIONAL([WITH_LIBSVM_LIBSRC], [test "x$libsvm_state" = "xok"])#

AC_CONFIG_FILES([Makefile \
                 src/Makefile \
                 tools/Makefile \
                 tools/Answerer/Makefile \
                 tools/Cliquer/Makefile \
                 tools/Clinician/Makefile \
                 tools/Clusterer/Makefile \
                 tools/Clusters2Dab/Makefile \
                 tools/COALESCE/Makefile \
                 tools/Combiner/Makefile \
                 tools/NetworkCombiner/Makefile \
                 tools/DChecker/Makefile \
                 tools/Dat2Dab/Makefile \
                 tools/Dat2Graph/Makefile \
                 tools/Data2Bnt/Makefile \
                 tools/Data2Features/Makefile \
                 tools/Data2Sql/Makefile \
                 tools/DataDumper/Makefile \
                 tools/Distancer/Makefile \
                 tools/Explainer/Makefile \
                 tools/Filterer/Makefile \
                 tools/Funcaeologist/Makefile \
                 tools/Funcifier/Makefile \
                 tools/Funcographer/Makefile \
                 tools/Hubber/Makefile \
                 tools/KNNImputer/Makefile \
                 tools/Mat2Txt/Makefile \
                 tools/MCluster/Makefile \
                 tools/Matcher/Makefile \
                 tools/MIer/Makefile \
                 tools/MIed/Makefile \
                 tools/NetMiner/Makefile \
                 tools/Normalizer/Makefile \
                 tools/Orthologer/Makefile \
                 tools/Overlapper/Makefile \
                 tools/PCLPlotter/Makefile \
                 tools/PCLServer/Makefile \
                 tools/PCL2Bin/Makefile \
                 tools/Randomizer/Makefile \
                 tools/Seqs2Ngrams/Makefile \
                 tools/SMRF/Makefile \
                 tools/SVDer/Makefile \
                 tools/Synthesizer/Makefile \
                 tools/Txt2Bin/Makefile \
                 tools/BNConverter/Makefile \
                 tools/BNCreator/Makefile \
                 tools/SparseNetCombiner/Makefile \
                 tools/SpeciesConnector/Makefile \
                 tools/BNEvaluator/Makefile \
                 tools/BNFunc/Makefile \
                 tools/BNTester/Makefile \
                 tools/BNTruster/Makefile \
                 tools/BNs2Txt/Makefile \
                 tools/BNUnraveler/Makefile \
                 tools/BNWeaver/Makefile \
                 tools/Contexter/Makefile \
                 tools/Counter/Makefile \
                 tools/Data2DB/Makefile \
                 tools/SeekReader/Makefile \
                 tools/SeekTest/Makefile \
                 tools/SeekMiner/Makefile \
                 tools/SeekPrep/Makefile \
                 tools/SeekServer/Makefile \
                 tools/SeekPValue/Makefile \
                 tools/SeekEvaluator/Makefile \
                 tools/SeekGeneRecommender/Makefile \
                 tools/SeekAggregatedDataset/Makefile \
         tools/DBCombiner/Makefile \
                 tools/DSLConverter/Makefile \
                 tools/Dab2Dad/Makefile \
                 tools/Edges2Posteriors/Makefile \
                 tools/MEFIT/Makefile \
                 tools/Data2Svm/Makefile \
                 tools/SVMer/Makefile \
                 tools/SVMperfer/Makefile \
                 tools/SVMmulticlass/Makefile \
                 tools/SVMhierarchy/Makefile \
                 tools/SVMperfing/Makefile \
                 tools/LibSVMer/Makefile \
                 tools/VWer/Makefile \
                 tools/OntoShell/Makefile \
                 tools/BNServer/Makefile])


if test "x$pthread_state" != "xok"; then
   LOCAL_STATE=warn
else
   LOCAL_STATE=ok
fi


if test "x$LOCAL_STATE" = "xok"; then
AC_OUTPUT
fi


# Print out a summary.
echo ""
if $with_tools; then
echo "Used for library and tools:"
else
echo "Dependencies:"
fi
echo "  log4cpp                 = $log4cpp_info"
echo "  SMILE                   = $smile_info"
echo "  SVM Perf                = $svm_perf_info"
echo "  SVM Multiclass          = $svm_multiclass_info"
echo "  SVM Hierarchy          = $svm_hierarchy_info"
echo "  Vowpal Wabbit           = $vowpal_wabbit_info"
echo "  pthread                 = $pthread_info"
echo "  gsl                     = $gsl_info"
echo "  OpenMP                  = $sleipnir_cv_openmp"
echo "  LibSVM                  = $libsvm_info"
if $with_tools; then
echo "Used for tools:"
echo "  boost graph, regex      = $boost_info"
echo "  readline                = $readline_info"
echo "  gengetopt               = $prog_gengetopt_info"
fi

echo ""
echo "Installation prefix = $prefix"

# Print these warnings here, so that they are at the end.

if test "x$log4cpp_info" = "xstub"; then
  cat << EOF

** WARNING: Using a simple replacement stub for log4cpp.
If this is insufficient rerun configure with --with-log4cpp=(path).
Sometimes there is no installed static library for log4cpp.
The full log4cpp library also requires pthreads.
EOF
fi

if test "x$log4cpp_state" = "xwarn"; then
  cat << EOF

The path to log4cpp may be incorrect. I looked for Category.hh
and liblog4cpp.a and did not find both.
EOF

fi

if test "x$sleipnir_cv_openmp" != "xyes"; then
        cat << EOF
** BUILDING WITHOUT OpenMP SUPPORT
OpenMP is essential for parallization in SEEK toolset.
SeekServer, SeekGeneRecommender, SeekMiner will run on single thread.
This is either because you have disabled OpenMP using --disable-openmp, or
you have an outdated version of GCC that does not support OpenMP.
In the latter case, please obtain the latest GCC compiler from http://gcc.gnu.org/.
EOF
fi

if test "x$gsl_state" != "xok"; then
        cat << EOF
** BUILDING WITHOUT GSL
GNU GSL is required for building SEEK toolset.
SeekServer, SeekGeneRecommender, SeekMiner would not be built.
GNU GSL is available at http://www.gnu.org/software/gsl/.
EOF
fi

if test "x$smile_state" != "xok"; then
  cat << EOF

** BUILDING WITHOUT SMILE
SMILE is strongly recommended.
Bayes net functionality will be missing.
SMILE is available as a binary distribution from http://genie.sis.pitt.edu/.
EOF

fi

if test "x$smile_state" = "xwarn"; then
  cat << EOF
** WARNING: The path to SMILE may be incorrect.
I looked for smile.h and libsmile.a and did not find both.
EOF

fi

if test "x$svm_perf_state" != "xok"; then
  cat << EOF

** BUILDING WITHOUT SVM PERF
SVM Perf is strongly recommended.
It is available from http://svmlight.joachims.org/.
EOF

fi

if test "x$svm_perf_state" = "xwarn"; then
  cat << EOF
** WARNING: The path to SVM Perf may be incorrect.
I looked for svm_struct_api.h and libsvmperf.a and did not
find both. NOTE: You may need to make libsvmperf.a. See the README.
EOF

fi



if test "x$svm_multiclass_state" != "xok"; then
  cat << EOF

** BUILDING WITHOUT SVM MULTICLASS
SVM Multiclass is strongly recommended.
It is available from http://svmlight.joachims.org/.
EOF

fi

if test "x$svm_multiclass_state" = "xwarn"; then
  cat << EOF
** WARNING: The path to SVM Multiclass may be incorrect.
I looked for svm_struct_api.h and libsvmmulticlass.a and did not
find both. NOTE: You may need to make libsvmmulticlass.a. See the README.
EOF

fi


if test "x$svm_hierarchy_state" != "xok"; then
  cat << EOF

** BUILDING WITHOUT SVM HIERARCHY
SVM Hierarchy is strongly recommended.
It is available from http://svmlight.joachims.org/.
EOF

fi

if test "x$svm_hierarchy_state" = "xwarn"; then
  cat << EOF
** WARNING: The path to SVM Hierarchy may be incorrect.
I looked for svm_struct_api.h and libsvmhierarchy.a and did not
find both. NOTE: You may need to make libsvmhierarchy.a. See the README.
EOF

fi



if test "x$smile_state" = "xwarn"; then
  cat << EOF
** WARNING: The path to SMILE may be incorrect.
I looked for smile.h and libsmile.a and did not find both.
EOF

fi

##
if test "x$libsvm_state" != "xok"; then
  cat << EOF

** BUILDING WITHOUT LIBSVM
LIBSVM is strongly recommended.
It is available from http://www.csie.ntu.edu.tw/~cjlin/libsvm/.
EOF

fi

if test "x$libsvm_state" = "xwarn"; then
  cat << EOF
** WARNING: The path to LibSVM may be incorrect.
I looked for  libsvm.a and did not
find. NOTE: You may need to make liblibsvm.a.

ar rcs liblibsvm.a *.o

See the README.
EOF

fi
##


if $with_tools; then

  if test "x$prog_gengetopt_state" = "xno"; then
  cat << EOF

** WARNING: gengetopt was not found. This is probably fine, but if any
of the cmdline.c or cmdline.h files are missing or out of date, they
will not be rebuilt.
EOF

  fi

  if test "x$readline_state" != "xok"; then
  cat << EOF

** BUILDING WITHOUT READLINE
OntoShell will not be made. It's possible the problem could be with
your ncurses library.
EOF

  fi

  if test "x$readline_state" = "xwarn"; then
  cat << EOF
** WARNING: The path to readline may be incorrect.
I looked for readline/readline.h and libreadline.a and did not
find both.
EOF

  fi

  if test "x$boost_state" != "xok"; then
  cat << EOF

** BUILDING WITHOUT BOOST
BNServer will not be made (requires boost_graph and boost_regex libraries).
Boost is available from http://boost.org/.
EOF

  fi

  if test "x$boost_state" = "xwarn"; then
  cat << EOF
** WARNING: The Boost path(s) may be incorrect.
I looked for boost/graph/graphviz.hpp, libboost_graph.a, libboost_regex.a
and did not find them all.
EOF

  fi

#  if test "x$with_boost" != "x" -a "x$with_boost" != "xno"; then
#  cat << EOF
#
#** A WARNING ABOUT BOOST
#Boost installs header files under a version-dependent directory,
#for example in /usr/local/include/boost/boost-1_34_1/boost/. I have
#no reliable way of determining the Boost version, so you may run into
#a problem with missing Boost includes. If this is the case,
#simply use --with-boost-includes and --with-boost-graph-libs to direct
#me to their actual locations, for example:
#  --with-boost-includes=/usr/local/include/boost/boost-1_34_1
#  --with-boost-graph-libs=/usr/local/lib/libboost_graph-mt-gcc34_1_34_1.dylib
#EOF
#
#  fi

else

  cat <<EOF

** Tools will not be made.
EOF

fi

if test "x$vowpal_wabbit_state" != "xok"; then
  cat << EOF

** BUILDING WITHOUT VOWPAL WABBIT
Vowpal Wabbit is optional.
It is available from http://hunch.net/~vw/.
Vowpal Wabbit also requires pthreads, boost, and zlib.
EOF
# boost::program_options

fi

if test "x$vowpal_wabbit_state" = "xwarn"; then
  cat << EOF
** WARNING: The path to Vowpal Wabbit may be incorrect.
I looked for vw.h and libvw.a and did not find both.
EOF

fi

echo ""

if test "x$LOCAL_STATE" = "xok"; then
  echo "Now run make."
else
  echo "** Some requirements were not met. Please review and try again. **"
fi

echo ""

#set
