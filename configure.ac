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
# Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
#
# If you use this library, the included executable tools, or any related
# code in your work, please cite the following publication:
# Curtis Huttenhower, Mark Schroeder, and Olga G. Troyanskaya.
# "The Sleipnir library for computational functional genomics"

#! TODO:
#! b. Warn for pthread.
#! d. Gengetopt checks.
#! e. Fix boost logic?

AC_PREREQ(2.59)
AC_INIT([libSleipnir], [0.9.0], [ogt@princeton.edu])

AC_COPYRIGHT([This software is copyrighted.])

# Any file under src will do.
AC_CONFIG_SRCDIR([src/stdafx.h])
# doesn't work with automake?
#AC_CONFIG_HEADER([config.h])
#AM_CONFIG_HEADER([config.h])
m4_pattern_allow([AM_INIT_AUTOMAKE])
AM_INIT_AUTOMAKE

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


## LOCAL_PATH_PROG(prog name, local name,
##   action if found, extra action if found installed, extra action if found in with path,
##   action if not found, action if user requested without)
##
## Set prog_(local name)_dir to dir containing it, prog_(local name) to the path, prog_(local name)_info.
## Also set prog_(local name)_state to "ok" or "no".
define([LOCAL_PATH_PROG],
[
if test "x$with_$2" = "xno"; then
  prog_$2_info=NONE
  prog_$2_state=no
#  true   # in case $7 isn't defined or is empty
  $7
elif test "x$with_$2" != "xcheck"; then
#  if echo with_$2 | $EGREP -q '^extlib/'; then extlib=true; else extlib=false; fi
  prog_$2_info=$with_$2
  prog_$2_state=ok
  qualify_path with_$2
  prog_$2=$with_$2
  $3
  $5
else
  AC_PATH_PROG([prog_$2], [$1])
  if test "x$prog_$2" = "x"; then
#    true   # in case is not defined or is empty
    prog_$2_info="not found"
    prog_$2_state=no
    $6
  else
    prog_$2_info="found installed"
    prog_$2_state=ok
    $3
    $4
  fi
fi
if test "x$prog_$2" != "x"; then
  prog_$2_dir=`AS_DIRNAME([$prog_$2])`
  prog_$2_prefix=$prog_$2_dir/..     # assumes bin and include are at same level
  $4
fi
])


## LOCAL_CHECK_LIB(lib name, local name, function,
##   action if found, extra action if found installed, extra action if found in with path,
##   action if not found, action if user requested without,
##   includes)
## Set $(local_name)_info to "NONE", location, or "installed" (can be overridden with above actions).
## Also set $(local_name)_state to "ok" or "no".
define([LOCAL_CHECK_LIB],
[
if test "x$with_$2" = "xno"; then
  $2_info=NONE
  $2_state=no
  true   # in case not defined or is empty
  $8
elif test "x$with_$2" != "xcheck"; then
# We rely on the AC_ARG_WITH call to set $with_(local name) to "check"
# if the user does not specify it. That way it is never undefined.
# That may not be necessary.
#  if echo with_$2 | $EGREP -q '^extlib/'; then extlib=true; else extlib=false; fi
  $2_info=$with_$2
  $2_state=ok
  qualify_path with_$2
  $4
  $6
else
  # In AC_LINK_IFELSE, the last argument is optional. If it is specified but is
  # empty, however, we get a syntax error. But, we have to be careful,
  # because arguments to this macro may contain quotes, semicolons, etc.
  # So we work around it by always supplying a last argument.
  flag=false
  AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[$9]], [[$3]])],
    [$2_info="found installed"
     $2_state=ok
     $4
     $5],
    [flag=true])
  if $flag; then
    # we also can't have an empty if block...
#    true
    $2_info="not found"
    $2_state=no
    $7
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
result=$1
found=false
for p in $2; do
  if test -e "$1/$p/$3"; then
    result=$1/$p
    found=true
    break
  fi
done
$4=$result
if ! $found; then
  if ! test -e $result; then
    true   # in case the following is empty
    $5
  fi
fi
])


## GENGETOPT
# If the user specifies a path, use it.
# If the user doesn't specify anything, check for it.
# While checking, look for a source distribution in extlib/.
# Gengetopt is required (for tools).
#! CFLAGS assumes headers in same directory as binaries for installed
AC_ARG_WITH([gengetopt],
	[AS_HELP_STRING([--with-gengetopt=PATH], [prefix of gengetopt installation])],
	[],
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
AC_SUBST(GENGETOPT_PREFIX)
AC_SUBST(GENGETOPT_CFLAGS)
AC_SUBST(GENGETOPT_LIBS)


## LOG4CPP
# If the user specifies a path, use it.
# If the user says not to use it, use the stub.
# If the user doesn't say anything, check for it. Use the stub if not found.
# While checking, look for a source distribution in extlib/.
# If using stub, set USE_LOG4CPP_STUB=1.
AC_ARG_WITH([log4cpp],
	[AS_HELP_STRING([--with-log4cpp=PATH], [prefix of log4cpp installation or build directory])],
	[],   # $with_log4cpp will either be a path or "no"
	[with_log4cpp=check])
LOCAL_CHECK_LIB([log4cpp], [log4cpp], [log4cpp::Category::getInstance],
	[],                                                        dnl found
	[],                                                        dnl and found installed
	[LOG4CPP_PREFIX=$with_log4cpp],                            dnl and found in specified path
	[log4cpp_info="stub (library not found)"
	 log4cpp_state=stub
	 AC_DEFINE([USE_LOG4CPP_STUB], [1])],                      dnl not found
	[log4cpp_info=stub
         log4cpp_state=stub
         AC_DEFINE([USE_LOG4CPP_STUB], [1])],                      dnl requested without
	[#include "log4cpp/Category.hh"])
# if stub, there might be additional info, so use grep
if test "x$log4cpp_state" != "xstub"; then
  ## set LOG4CPP_CFLAGS and LOG4CPP_LIBS
  if test "x$LOG4CPP_PREFIX" != "x"; then
    LOCAL_CHECK_APPEND_PATHS([$LOG4CPP_PREFIX], [include], [log4cpp/Category.hh], [LOG4CPP_INCLUDE_DIR], [log4cpp_state=warn])
    LOCAL_CHECK_APPEND_PATHS([$LOG4CPP_PREFIX], [lib src/.libs], [liblog4cpp.a], [LOG4CPP_LIB_DIR], [log4cpp_state=warn])
    LOG4CPP_CFLAGS="-I$LOG4CPP_INCLUDE_DIR"
    LOG4CPP_LIBS="-L$LOG4CPP_LIB_DIR"
#    if ! test -f "$LOG4CPP_INCLUDE_DIR/log4cpp/Category.hh"; then log4cpp_state=warn; fi
#    if ! test -f "$LOG4CPP_LIB_DIR/liblog4cpp.a"; then log4cpp_state=warn; fi
  fi
  LOG4CPP_LIBS="$LOG4CPP_LIBS -llog4cpp"
fi
AC_SUBST(LOG4CPP_PREFIX)
AC_SUBST(LOG4CPP_CFLAGS)
AC_SUBST(LOG4CPP_LIBS)


## SMILE
# If the user specifies a path, use it.
# If the user says not to use it, then don't.
# If the user doesn't specify anything, check for it.
# While checking, look for a source distribution in extlib/.
# If none, set NO_SMILE=1.
AC_ARG_WITH([smile],
	[AS_HELP_STRING([--with-smile=PATH], [prefix of SMILE installation])],
	[],
	[with_smile=check])
# it's hard to find a decent static function to use...
LOCAL_CHECK_LIB([smile], [smile], [RoundProbability],
        [],                                                        dnl found
        [],                                                        dnl and found installed
        [SMILE_PREFIX=$with_smile],                                dnl and found in specified path
	[AC_DEFINE([NO_SMILE], [1])],                              dnl not found
        [AC_DEFINE([NO_SMILE], [1])],                              dnl requested without
	[#include "smile.h"])
if test "x$smile_state" != "xno"; then
  ## set SMILE_CFLAGS and SMILE_LIBS
  if test "x$SMILE_PREFIX" != "x"; then
    LOCAL_CHECK_APPEND_PATHS([$SMILE_PREFIX], [include], [smile.h], [SMILE_INCLUDE_DIR], [smile_state=warn])
    LOCAL_CHECK_APPEND_PATHS([$SMILE_PREFIX], [lib], [libsmile.a], [SMILE_LIB_DIR], [smile_state=warn])
#    if ! test -f "$SMILE_INCLUDE_DIR/smile.h"; then smile_state=warn; fi
#    if ! test -f "$SMILE_LIB_DIR/libsmile.a"; then smile_state=warn; fi
    if test -f "$SMILE_LIB_DIR/libsmilexml.a"; then
      AC_DEFINE([SMILEXML_LIB], [1])
      SMILEXML_LIB="-lsmilexml"
    fi
    SMILE_CFLAGS="-I$SMILE_INCLUDE_DIR"
    SMILE_LIBS="-L$SMILE_LIB_DIR"
  fi
  SMILE_LIBS="$SMILE_LIBS $SMILEXML_LIB -lsmile"
fi
AC_SUBST(SMILE_PREFIX)
AC_SUBST(SMILE_CFLAGS)
AC_SUBST(SMILE_LIBS)


## SVM LIGHT
# If the user specifies a path, use it.
# If the user says not to use it, then don't.
# If the user doesn't specify anything, check for it.
# While checking, look for a source distribution in extlib/.
# If none, set NO_SVM_LIGHT=1.
AC_ARG_WITH([svm_light],
	[AS_HELP_STRING([--with-svm-light=PATH], [prefix of SVM Light installation])],
	[],
	[with_svm_light=check])
AC_LANG_PUSH([C])
LOCAL_CHECK_LIB([svmlight], [svm_light], [optimize_svm],
        [],                                                        dnl found
        [],                                                        dnl and found installed
        [SVM_LIGHT_PREFIX=$with_svm_light],                        dnl and found in specified path
	[AC_DEFINE([NO_SVM_LIGHT], [1])],                          dnl not found
        [AC_DEFINE([NO_SVM_LIGHT], [1])],                          dnl requested without
	[#include <sys/types.h>
	 #include "svm_common.h"
	 #include "svm_learn.h"])
AC_LANG_POP
if test "x$svm_light_state" != "xno"; then
  ## set SVM_LIGHT_CFLAGS and SVM_LIGHT_LIBS
  if test "x$SVM_LIGHT_PREFIX" != "x"; then
    LOCAL_CHECK_APPEND_PATHS([$SVM_LIGHT_PREFIX], [include], [svm_common.h], [SVM_LIGHT_INCLUDE_DIR], [svm_light_state=warn])
    LOCAL_CHECK_APPEND_PATHS([$SVM_LIGHT_PREFIX], [lib], [libsvmlight.a], [SVM_LIGHT_LIB_DIR], [svm_light_state=warn])
#    if ! test -f "$SVM_LIGHT_INCLUDE_DIR/svm_common.h"; then svm_light_state=warn; fi
#    if ! test -f "$SVM_LIGHT_LIB_DIR/libsvmlight.a"; then svm_light_state=warn; fi
    SVM_LIGHT_CFLAGS="-I$SVM_LIGHT_INCLUDE_DIR"
    SVM_LIGHT_LIBS="-L$SVM_LIGHT_LIB_DIR"
  fi
  SVM_LIGHT_LIBS="$SVM_LIGHT_LIBS -lsvmlight"
fi
AC_SUBST(SVM_LIGHT_PREFIX)
AC_SUBST(SVM_LIGHT_CFLAGS)
AC_SUBST(SVM_LIGHT_LIBS)


## BOOST
# If the user specifies a path, use it.
# If not specified, check for it.
# BOOST graph-mt lib is required.
AC_ARG_WITH([boost_includes],
	[AS_HELP_STRING([--with-boost-includes=PATH], [location of boost header files])],
	[],   # $with_boost_includes will either be a path or "no"
	[with_boost_includes=check])
AC_ARG_WITH([boost_graph_lib],
	[AS_HELP_STRING([--with-boost-graph-lib=PATH], [path to boost_graph-mt library file])],
	[],   # $with_boost_graph_lib will either be a path or "no"
        [with_boost_graph_lib=check])
# ensure consistency in options
if test "x$with_boost_includes" = "xno" -o "x$with_boost_graph_lib" = "xno"; then
  with_boost_includes=no
  with_boost_graph_lib=no
  boost_includes_info=NONE
#  boost_includes_state=no
  boost_graph_lib_info=NONE
  boost_state=no
elif test "x$with_boost_includes" = "xcheck" -o "x$with_boost_graph_lib" = "xcheck"; then
  with_boost_includes=no
  with_boost_graph_lib=no
  boost_includes_info="not specified (required if using with-boost-graph-lib)"
#  boost_includes_state=no
  boost_graph_lib_info="not specified (required if using with-boost-includes)"
  boost_state=no
else
  LOCAL_CHECK_LIB([boost_graph-mt], [boost_graph_lib], [boost::detail::graph::read_graphviz],
	[],                                                         dnl found
	[boost_includes_info="found installed"],                    dnl and found installed
	[BOOST_GRAPH_LIB=$with_boost_graph_lib
	 boost_includes_info=$with_boost_includes],                 dnl and found in specified path
	[boost_includes_info=$boost_graph_lib_info],                dnl not found
	[boost_includes_info=$boost_graph_lib_info],                dnl requested without
	[#include "boost/graph/graphviz.hpp"])
fi
if test "x$boost_state" != "xno"; then
  ## set BOOST_CFLAGS and BOOST_LIBS
  if test "x$BOOST_GRAPH_LIB" != "x"; then
    # boost_graph_lib_info and boost_includes_info are paths
    # now look for the headers
    BOOST_LIBS=$BOOST_GRAPH_LIB
    qualify_path with_boost_includes
    LOCAL_CHECK_APPEND_PATHS([$with_boost_includes], [..], [boost/config.hpp], [BOOST_INCLUDES_DIR], [boost_state=warn])
    if ! test -f $BOOST_GRAPH_LIB; then boost_state=warn; fi
    if ! echo $BOOST_GRAPH_LIB | $GREP -q "libboost_graph-mt"; then boost_state=warn; fi
    BOOST_CFLAGS="-I$BOOST_INCLUDES_DIR"
  fi
fi
#AC_SUBST(BOOST_PREFIX)
AC_SUBST(BOOST_CFLAGS)
AC_SUBST(BOOST_LIBS)


## READLINE
# If the user specifies a path, use it.
# If not specified, check for it.
# If no readline, do not build OntoShell.
AC_ARG_WITH([readline],
	[AS_HELP_STRING([--with-readline=PATH], [prefix of readline installation])],
	[],   # $with_readline will either be a path or "no"
	[with_readline=check])
AC_LANG_PUSH([C])
LOCAL_CHECK_LIB([readline], [readline], [readline],
	[],                                                        dnl found
	[],                                                        dnl and found installed
	[READLINE_CFLAGS=$with_readline],                          dnl and found in specified path
	[],                                                        dnl not found
	[],                                                        dnl requested without
	[#include <stdio.h>
	 #include "readline/readline.h"])
AC_LANG_POP
if test "x$readline_state" != "xno"; then
  ## set READLINE_CFLAGS and READLINE_LIBS
  if test "x$READLINE_PREFIX" != "x"; then
    LOCAL_CHECK_APPEND_PATHS([$READLINE_PREFIX], [include], [readline/readline.h], [READLINE_INCLUDE_DIR], [readline_state=warn])
    LOCAL_CHECK_APPEND_PATHS([$READLINE_PREFIX], [lib], [libreadline.a], [READLINE_LIB_DIR], [readline_state=warn])
#    if ! test -f "$READLINE_INCLUDE_DIR/readline/readline.h"; then readline_state=warn; fi
#    if ! test -f "$READLINE_LIB_DIR/libreadline.a"; then readline_state=warn; fi
    READLINE_CFLAGS="-I$READLINE_INCLUDE_DIR"
    READLINE_LIBS="-L$READLINE_LIB_DIR"
  fi
  # it should be very rare indeed that ncurses is not there if readline is there
  READLINE_LIBS="$READLINE_LIBS -lreadline -lncurses"
  # apparently on Macs, rl_completion_matches() is named completion_matches()
  # hopefully this is an unusual exception
  AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([#include <stdio.h>
		      #include "readline/readline.h"],
		      [rl_completion_matches])],
    [],
    [AC_DEFINE([rl_completion_matches], [completion_matches])])
fi
AC_SUBST(READLINE_PREFIX)
AC_SUBST(READLINE_CFLAGS)
AC_SUBST(READLINE_LIBS)


## READLINE
# If the user specifies a path, use it.
# If not specified, check for it.
# If no readline, do not build OntoShell.
AC_ARG_WITH([pthread],
	[AS_HELP_STRING([--with-pthread=PATH], [prefix of pthread installation])],
	[],   # $with_pthread will either be a path or "no"
	[with_pthread=check])
AC_LANG_PUSH([C])
LOCAL_CHECK_LIB([pthread], [pthread], [pthread_create],
	[],                                                        dnl found
	[],                                                        dnl and found installed
	[PTHREAD_CFLAGS=$with_pthread],                            dnl and found in specified path
	[],                                                        dnl not found
	[],                                                        dnl requested without
	[#include <stdio.h>
	 #include "pthread.h"])
AC_LANG_POP
if test "x$pthread_state" != "xno"; then
  ## set PTHREAD_CFLAGS and PTHREAD_LIBS
  if test "x$PTHREAD_PREFIX" != "x"; then
    LOCAL_CHECK_APPEND_PATHS([$PTHREAD_PREFIX], [include], [pthread.h], [PTHREAD_INCLUDE_DIR], [pthread_state=warn])
    LOCAL_CHECK_APPEND_PATHS([$PTHREAD_PREFIX], [lib], [libpthread.a], [PTHREAD_LIB_DIR], [pthread_state=warn])
    PTHREAD_CFLAGS="-I$PTHREAD_INCLUDE_DIR"
    PTHREAD_LIBS="-L$PTHREAD_LIB_DIR"
  fi
  PTHREAD_LIBS="$PTHREAD_LIBS -lpthread"
fi
AC_SUBST(PTHREAD_PREFIX)
AC_SUBST(PTHREAD_CFLAGS)
AC_SUBST(PTHREAD_LIBS)


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
AM_CONDITIONAL([WITH_SMILE_TOOLS], [test "x$smile_state" != "xno"])
AM_CONDITIONAL([WITH_SVM_TOOLS], [test "x$svm_light_state" != "xno"])
AM_CONDITIONAL([WITH_READLINE_TOOLS], [test "x$readline_state" != "xno"])
AM_CONDITIONAL([WITH_BOOST_TOOLS], [test "x$boost_state" != "xno"])
AM_CONDITIONAL([WITH_GENGETOPT], [test "x$prog_gengetopt_state" != "xno"])

AC_CONFIG_FILES([Makefile \
                 src/Makefile \
                 tools/Makefile \
		 tools/Answerer/Makefile \
		 tools/Cliquer/Makefile \
		 tools/Clusterer/Makefile \
		 tools/Clusters2Dab/Makefile \
		 tools/Combiner/Makefile \
		 tools/DChecker/Makefile \
		 tools/Dat2Dab/Makefile \
		 tools/Dat2Graph/Makefile \
		 tools/Data2Bnt/Makefile \
		 tools/Data2Features/Makefile \
		 tools/Data2Sql/Makefile \
		 tools/DataDumper/Makefile \
		 tools/Distancer/Makefile \
		 tools/Explainer/Makefile \
		 tools/Funcaeologist/Makefile \
		 tools/Funcifier/Makefile \
		 tools/Funcographer/Makefile \
		 tools/Hubber/Makefile \
		 tools/KNNImputer/Makefile \
		 tools/MCluster/Makefile \
		 tools/MIer/Makefile \
		 tools/Normalizer/Makefile \
		 tools/Orthologer/Makefile \
		 tools/Overlapper/Makefile \
		 tools/Randomizer/Makefile \
		 tools/Seqs2Ngrams/Makefile \
		 tools/Txt2Bin/Makefile \
                 tools/BNConverter/Makefile \
		 tools/BNCreator/Makefile \
		 tools/BNEvaluator/Makefile \
		 tools/BNFunc/Makefile \
		 tools/BNTester/Makefile \
		 tools/BNTruster/Makefile \
		 tools/BNUnraveler/Makefile \
		 tools/BNWeaver/Makefile \
		 tools/Contexter/Makefile \
		 tools/Data2DB/Makefile \
		 tools/DSLConverter/Makefile \
		 tools/Dab2Dad/Makefile \
		 tools/Edges2Posteriors/Makefile \
		 tools/MEFIT/Makefile \
                 tools/Data2Svm/Makefile \
		 tools/SVMer/Makefile \
                 tools/OntoShell/Makefile \
                 tools/BNServer/Makefile])


AC_OUTPUT


# Print out a summary.
echo ""
if $with_tools; then
echo "Used for library and tools:"
else
echo "Dependencies:"
fi
echo "  log4cpp         = $log4cpp_info"
echo "  SMILE           = $smile_info"
echo "  SVM light       = $svm_light_info"
echo "  pthread         = $pthread_info"
if $with_tools; then
echo "Used for tools:"
echo "  boost graph lib = $boost_graph_lib_info"
echo "  boost includes  = $boost_includes_info"
echo "  readline        = $readline_info"
echo "  gengetopt       = $prog_gengetopt_info"
fi

# Print these warnings here, so that they are at the end.

# if stub there might be additional info, so use grep
if test "x$log4cpp_state" = "xstub"; then
  cat << EOF

** Using a simple replacement stub for log4cpp.
If this is insufficient rerun configure with --with-log4cpp=(path).
EOF

elif test "x$log4cpp_state" = "xwarn"; then
  cat << EOF

** WARNING: The path to log4cpp may be incorrect.
I looked for Category.hh and liblog4cpp.a ane didn't
find one (or both).
EOF

fi

if test "x$smile_state" = "xno"; then
  cat << EOF

** BUILDING WITHOUT SMILE
SMILE is strongly recommended.
Bayes net functionality will be missing.
SMILE is available as a binary distribution from http://genie.sis.pitt.edu/.
EOF

elif test "x$smile_state" = "xwarn"; then
  cat << EOF

** WARNING: The path to SMILE may be incorrect.
I looked for smile.h and libsmile.a and didn't
find one (or both).
EOF

fi

if test "x$svm_light_state" = "xno"; then
  cat << EOF

** BUILDING WITHOUT SVM LIGHT
SVM light is strongly recommended.
It is available from http://svmlight.joachims.org/.
EOF

elif test "x$svm_light_state" = "xwarn"; then
  cat << EOF

** WARNING: The path to SVM light may be incorrect.
I looked for svm_common.h and libsvmlight.a and didn't
find one (or both).
EOF

fi

if $with_tools; then

  if test "x$prog_gengetopt_state" = "xno"; then
  cat << EOF

** Gengetopt was not found. This is probably fine, but if any of the
cmdline.c or cmdline.h files are missing or out of date, they
will not be rebuilt.
EOF

  fi

  if test "x$readline_state" = "xno"; then
  cat << EOF

** BUILDING WITHOUT READLINE
OntoShell will not be made.
EOF

  elif test "x$readline_state" = "xwarn"; then
  cat << EOF

** WARNING: The path to readline may be incorrect.
I looked for readline/readline.h and libreadline.a and didn't
find one (or both).
EOF

  fi

  if test "x$boost_state" = "xno"; then
  cat << EOF

** BUILDING WITHOUT BOOST
BNServer will not be made (requires boost_graph-mt library).
Both --with-boost-includes and --with-boost-graph-lib are required for Boost.
Boost is available from http://boost.org/.
EOF

  elif test "x$boost_state" = "xwarn"; then
  cat << EOF

** WARNING: The Boost path(s) may be incorrect.
I looked for boost/config.hpp and a library whose name
contains libboost_graph-mt and didn't find one (or both).
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

echo ""
echo "Now run make."

#set
