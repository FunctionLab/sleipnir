# TODO
# a. What is the deal with AM_CONFIG_HEADER and AC_CONFIG_HEADER?
# b. Installation.
# c. Ranlib SMILE
# d. Check SMILE version!
#
# Potential portability concerns:
# a. programs (ls, egrep, dirname, true, false)

AC_PREREQ(2.59)
AC_INIT([libSleipnir], [0.9.0], [ogt@princeton.edu])

AC_COPYRIGHT([This software is copyrighted.])

# Any file under src will do.
AC_CONFIG_SRCDIR([src/stdafx.h])
# doesn't work with automake
#AC_CONFIG_HEADER([config.h])
#AM_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE


topleveldir=$PWD

# make the path absolute
#!! this may not be portable
qualify_path () {
  d=${!1}
  if ! echo $d | egrep -q '^/'; then
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
AC_CHECK_HEADERS([float.h netinet/in.h sys/socket.h sys/time.h])

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

#!! AC_CHECK_FUNCS: stop if not found?
AC_CHECK_FUNCS([memmove memset munmap pow socket sqrt strchr strerror strrchr strstr])

AC_CHECK_FUNCS([inet_ntoa])   # tools/BNServer
AC_CHECK_FUNCS([mkdir])  # tools/BNWrangler
AC_CHECK_FUNCS([strtol])  # tools/Data2Features, tools/Data2Bnt


## LOCAL_PATH_PROG(prog name, local name, dir in extlib,
##   action if found, extra action if found in extlib, action if not found)
##
## Set prog_(local name)_dir to dir containing it, prog_(local name) to the path, prog_(local name)_info.
define([LOCAL_PATH_PROG],
[
AC_PATH_PROG([prog_$2], [$1])
if test "x$prog_$2" = "x"; then
  extpath=extlib/$3/$1
  if test -x $extpath; then
    prog_$2_info=$extpath
    prog_$2=$extpath
    qualify_path prog_$2
    prog_$2_dir=extlib/$3
    qualify_path prog_$2_dir
    prog_$2_prefix=$prog_$2_dir
    $4
    $5
  fi
else
  # prog_$2 already set
  prog_$2_info=installed
  prog_$2_dir=`dirname $prog_$2`
  # assumes bin and include are at same level
  prog_$2_prefix=$prog_$2_dir/..
  $4
fi
if test "x$prog_$2" = "x"; then
  $6
fi
])


## LOCAL_CHECK_LIB(lib name, local name, function,
##   action if found, extra action if found installed, extra action if found in with path,
##   action if not found, action if user requested without)
## Set $(local_name)_info to "NONE", location, or "installed" (can be overridden with above actions).
define([LOCAL_CHECK_LIB],
[
if test "x$with_$2" = "xno"; then
  $2_info=NONE
  $8
elif test "x$with_$2" != "xcheck"; then
# We rely on the AC_ARG_WITH call to set $with_(local name) to "check"
# if the user does not specify it. That way it is never undefined.
# That may not be necessary.
  if echo with_$2 | egrep -q '^extlib/'; then extlib=true; else extlib=false; fi
  $2_info=$with_$2
  qualify_path with_$2
  $4
  $6
else
  # In AC_CHECK_LIB, the last argument is optional. If it is specified but is
  # empty, however, we get a syntax error. But, we have to be careful,
  # because arguments to this macro may contain quotes, semicolons, etc.
  # So we work around it by always supplying a last argument.
  flag=false
  AC_CHECK_LIB([$1], [$3],
    [$2_info=installed
     $4
     $5],
    [flag=true])
  if $flag; then
    # we also can't have an empty if block...
    true
    $7
  fi
fi
])


## LOCAL_CHECK_APPEND_PATHS(base dir, paths to append, prefix to set, var to set)
## Look for one of the paths under the base dir. If found, set the variable
## to the concatenation of the prefix, the discovered path.
define([LOCAL_CHECK_APPEND_PATHS],
[
result=$3${$1}
for p in $2; do
  if test -d "${$1}/$p"; then
    result=$3${$1}/$p
    break
  fi
done
$4=$result
])


## GENGETOPT
# If the user specifies a path, use it.
# If the user doesn't specify anything, check for it.
# While checking, look for a source distribution in extlib/.
# Gengetopt is required (for tools).
#! CFLAGS assumes headers in same directory as binaries for installed
AC_ARG_WITH([gengetopt],
	[AS_HELP_STRING([--with-gengetopt=<path>], [prefix of gengetopt installation])],
	[],
	[with_gengetopt=check])
LOCAL_PATH_PROG([gengetopt], [gengetopt], [gengetopt-2.22/src],
        [GENGETOPT_PREFIX=$prog_gengetopt_prefix; GENGETOPT_INCLUDES="-I${GENGETOPT_PREFIX}/include"],
        [GENGETOPT_PREFIX=$prog_gengetopt_dir; GENGETOPT_INCLUDES="-I${prog_gengetopt_dir}"],
        [echo "GENGETOPT IS REQUIRED" 1>&2; exit 1])
## set CFLAGS and LIBS
if test "x$GENGETOPT_PREFIX" != "x"; then
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
	[AS_HELP_STRING([--with-log4cpp=<path>], [prefix of log4cpp installation or build directory])],
	[],   # $with_log4cpp will either be a path or "no"
	[with_log4cpp=check])
LOCAL_CHECK_LIB([log4cpp], [log4cpp], [log4cpp::Category::getInstance],
	[],                                                        dnl found (sets with_log4cpp_info=installed)
	[],                                                        dnl and found installed
	[LOG4CPP_PREFIX=$with_log4cpp],                            dnl and found in specified path
	[log4cpp_info="stub (not found)";                          dnl not found
	 AC_DEFINE([USE_LOG4CPP_STUB], [1])],
	[log4cpp_info=stub; AC_DEFINE([USE_LOG4CPP_STUB], [1])])   dnl requested without
## set LOG4CPP_CFLAGS and LOG4CPP_LIBS
if test "x$LOG4CPP_PREFIX" != "x"; then
  LOCAL_CHECK_APPEND_PATHS([LOG4CPP_PREFIX], [include], [-I], [LOG4CPP_CFLAGS])
  LOCAL_CHECK_APPEND_PATHS([LOG4CPP_PREFIX], [lib src/.libs], [-L], [LOG4CPP_LIBS])
  # It may be useful to allow a built (but not installed) path for log4cpp,
  # that is, a source distribution. However the log4cpp distribution does not
  # have a lib/, it has a src/.libs/.
#  if test -d ${LOG4CPP_PREFIX}/lib; then
#    LOG4CPP_LIBS="-L${LOG4CPP_PREFIX}/lib -llog4cpp"
#  elif test -d ${LOG4CPP_PREFIX}/src/.libs; then
#    LOG4CPP_LIBS="-L${LOG4CPP_PREFIX}/src/.libs -llog4cpp"
#  else
#    log4cpp_info=
#    LOG4CPP_PREFIX=
#    #! revert to using the stub?
#    echo "COULD NOT LOCATE LOG4CPP LIBRARY AT SPECIFIED PATH" 1>&2
#    exit 1
#  fi
fi
LOG4CPP_LIBS="$LOG4CPP_LIBS -llog4cpp"
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
	[AS_HELP_STRING([--with-smile=<path>], [prefix of SMILE installation])],
	[],
	[with_smile=check])
LOCAL_CHECK_LIB([smile], [smile], [DSL_header::GetID],
        [],                                                        dnl found
        [],                                                        dnl and found installed (sets smile_info=installed)
        [SMILE_PREFIX=$with_smile],                                dnl and found in specified path
	[],                                                        dnl not found (smile_info not set)
        [AC_DEFINE([NO_SMILE], [1])])                              dnl requested without
## set SMILE_CFLAGS and SMILE_LIBS
if test "x$SMILE_PREFIX" != "x"; then
  LOCAL_CHECK_APPEND_PATHS([SMILE_PREFIX], [include], [-I], [SMILE_CFLAGS])
  LOCAL_CHECK_APPEND_PATHS([SMILE_PREFIX], [lib], [-L], [SMILE_LIBS])
fi
SMILE_LIBS="$SMILE_LIBS -lsmile"
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
	[AS_HELP_STRING([--with-svm-light=<path>], [prefix of SVM Light installation])],
	[],
	[with_svm_light=check])
LOCAL_CHECK_LIB([svmlight], [svm_light], [kernel],
        [],                                                        dnl found (sets svm_light_info=installed)
        [],                                                        dnl and found installed
        [SVM_LIGHT_PREFIX=$with_svm_light],                        dnl and found in specified path
	[],                                                        dnl not found (svm_light_info not set)
        [AC_DEFINE([NO_SVM_LIGHT], [1])])                          dnl requested without
## set SVM_LIGHT_CFLAGS and SVM_LIGHT_LIBS
if test "x$SVM_LIGHT_PREFIX" != "x"; then
  LOCAL_CHECK_APPEND_PATHS([SVM_LIGHT_PREFIX], [include], [-I], [SVM_LIGHT_CFLAGS])
  LOCAL_CHECK_APPEND_PATHS([SVM_LIGHT_PREFIX], [lib], [-L], [SVM_LIGHT_LIBS])
fi
SVM_LIGHT_LIBS="$SVM_LIGHT_LIBS -lsvm"
AC_SUBST(SVM_LIGHT_PREFIX)
AC_SUBST(SVM_LIGHT_CFLAGS)
AC_SUBST(SVM_LIGHT_LIBS)


# Now define variables to be used to reference this library.
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


AC_CONFIG_FILES([Makefile src/Makefile \
		 tools/Makefile \
		 tools/Answerer/Makefile \
		 tools/BNConverter/Makefile \
		 tools/BNCreator/Makefile \
		 tools/BNEvaluator/Makefile \
		 tools/BNFunc/Makefile \
		 tools/BNServer/Makefile \
		 tools/BNTester/Makefile \
		 tools/BNTruster/Makefile \
		 tools/BNUnraveler/Makefile \
		 tools/BNWeaver/Makefile \
		 tools/BNWrangler/Makefile \
		 tools/Cliquer/Makefile \
		 tools/Clusterer/Makefile \
		 tools/Clusters2Dab/Makefile \
		 tools/Combiner/Makefile \
		 tools/Contexter/Makefile \
		 tools/DChecker/Makefile \
		 tools/DSLConverter/Makefile \
		 tools/Dab2Dad/Makefile \
		 tools/Dat2Dab/Makefile \
		 tools/Dat2Graph/Makefile \
		 tools/Data2Bnt/Makefile \
		 tools/Data2DB/Makefile \
		 tools/Data2Features/Makefile \
		 tools/Data2Sql/Makefile \
		 tools/Data2Svm/Makefile \
		 tools/DataDumper/Makefile \
		 tools/Distancer/Makefile \
		 tools/Edges2Posteriors/Makefile \
		 tools/Explainer/Makefile \
		 tools/Funcaeologist/Makefile \
		 tools/Funcifier/Makefile \
		 tools/Funcographer/Makefile \
		 tools/Hubber/Makefile \
		 tools/KNNImputer/Makefile \
		 tools/MCluster/Makefile \
		 tools/MEFIT/Makefile \
		 tools/MIer/Makefile \
		 tools/Normalizer/Makefile \
		 tools/OntoShell/Makefile \
		 tools/Orthologer/Makefile \
		 tools/Overlapper/Makefile \
		 tools/Randomizer/Makefile \
		 tools/SVMer/Makefile \
		 tools/Seqs2Ngrams/Makefile \
		 tools/Txt2Bin/Makefile ])

AC_OUTPUT


# Print out a summary.
echo ""
echo "  log4cpp     = $log4cpp_info"
echo "  SMILE       = $smile_info"
echo "  SVM light   = $svm_light_info"
echo ""
echo "  gengetopt   = $prog_gengetopt_info"

ok=true

if test "x$log4cpp_info" = "xstub"; then
  cat << EOF

Using a simple replacement stub for log4cpp.
If this is insufficient rerun configure with --with-log4cpp.
EOF
fi

if test "x$smile_info" = "x"; then
  ok=false
cat << EOF

SMILE NOT FOUND
To compile without SMILE, use --without-smile.
Otherwise use --with-smile=(path) to specify the location.
It is available as a binary distribution from http://genie.sis.pitt.edu/.
EOF

elif echo $smile_info | grep -q "NONE"; then
cat << EOF

BUILDING WITHOUT SMILE
SMILE is strongly recommended.
Bayes net functionality will be missing.
It is available as a binary distribution from http://genie.sis.pitt.edu/.
EOF

fi

if test "x$svm_light_info" = "x"; then
  ok=false
cat << EOF

SVM LIGHT NOT FOUND
To compile without SVM Light, use --without-svm-light.
Otherwise use --with-svm-light=(path) to specify the location.
It is available from http://svmlight.joachims.org/.
EOF

elif echo $svm_light_info | grep -q "NONE"; then
cat << EOF

BUILDING WITHOUT SVM LIGHT
SVM light is strongly recommended.
It is available from http://svmlight.joachims.org/.
EOF

fi

if $ok; then
  echo ""
  echo "Now run make."
else
  echo ""
  echo "Please correct the above errors before continuing."
fi

#set
