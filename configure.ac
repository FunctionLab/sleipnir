# TODO
# a. What is the deal with AM_CONFIG_HEADER and AC_CONFIG_HEADER?
# c. Check for SMILE.
# d. Check for SVM light.
# f. Build tools.
# g. Installation.

AC_PREREQ(2.59)
AC_INIT([libSleipnir], [0.9.0], [ogt@princeton.edu])

#!! update
AC_COPYRIGHT([This software is copyrighted.])

# Any file under src will do.
AC_CONFIG_SRCDIR([src/stdafx.h])
# doesn't work with automake
#AC_CONFIG_HEADER([config.h])
#AM_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE

topleveldir=$PWD

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

# Checks for libraries.
#!! how does this work with --with?
#AC_CHECK_LIB([log4cpp], [log4cpp::Category::getInstance],
#	[],
#	[echo "Error: log4cpp not found"; exit -1])

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


## LOG4CPP
# If the user specifies a path for log4cpp, use it.
# If the user says not to use it, or if he doesn't say anything, check for it.
AC_ARG_WITH([log4cpp],
	[AS_HELP_STRING([--with-log4cpp=<path>], [prefix of log4cpp installation or build directory])],
	[],   # $with_log4cpp will either be a path or "no"
	[with_log4cpp=check])
## set LOG4CPP_PREFIX and log4cpp_info
if test "x$with_log4cpp" = "xcheck"; then
  LOG4CPP_PREFIX=
  AC_CHECK_LIB([log4cpp], [log4cpp::Category::getInstance],
        [log4cpp_info="-llog4cpp"],
        [log4cpp_info="stub"])
elif test "x$with_log4cpp" = "xno"; then
  LOG4CPP_PREFIX=
  log4cpp_info=stub
else
  #! possibly check path
  qualify_path with_log4cpp
  LOG4CPP_PREFIX=$with_log4cpp
  log4cpp_info=$with_log4cpp
fi
## set LOG4CPP_CFLAGS and LOG4CPP_LIBS (and -DUSE_LOG4CPP_STUB=1 if stub)
if test "x$log4cpp_info" = "xstub"; then
  AC_DEFINE([USE_LOG4CPP_STUB], [1])
  LOG4CPP_CFLAGS=
  LOG4CPP_LIBS=
elif test "x$LOG4CPP_PREFIX" = "x"; then
  LOG4CPP_CFLAGS=
  LOG4CPP_LIBS=
else
  LOG4CPP_CFLAGS="-I${LOG4CPP_PREFIX}/include"
  # It may be useful to allow a built (but not installed) path for log4cpp,
  # that is, a source distribution. However the log4cpp distribution does not
  # have a lib/, it has a src/.libs/.
  if test -d ${LOG4CPP_PREFIX}/lib; then
    LOG4CPP_LIBS="-L${LOG4CPP_PREFIX}/lib"
  elif test -d ${LOG4CPP_PREFIX}/src/.libs; then
    LOG4CPP_LIBS="-L${LOG4CPP_PREFIX}/src/.libs"
  else
    # revert to using the stub
    echo "COULD NOT LOCATE LOG4CPP LIBRARY AT SPECIFIED PATH" 1>&2
    AC_DEFINE([USE_LOG4CPP_STUB], [1])
    LOG4CPP_PREFIX=
    LOG4CPP_CFLAGS=
    LOG4CPP_LIBS=
    log4cpp_info="stub (lib not found)"
  fi
fi
## include substitutions in output files
AC_SUBST(LOG4CPP_PREFIX)
AC_SUBST(LOG4CPP_CFLAGS)
AC_SUBST(LOG4CPP_LIBS)


## SMILE
AC_ARG_WITH([smile],
	[AS_HELP_STRING([--with-smile=<path>], [prefix of SMILE installation])],
	[],
	[with_smile=check])
## set SMILE_PREFIX and smile_info (and -DNO_SMILE if excluded)
if test "x$with_smile" = "xno"; then
  AC_DEFINE([NO_SMILE], [1])
  SMILE_PREFIX=
  smile_info="NONE"
elif test "x$with_smile" = "xcheck"; then
  SMILE_PREFIX=
  AC_CHECK_HEADER([smile.h],
        [smile_info=installed],
        [smile_info="NONE (not found)"; AC_DEFINE([NO_SMILE], [1])])
#  AC_CHECK_LIB([smile], [????],
#        [smile_info="-lsmile"],
#        [smile_info="not found"])
else
  #! possibly check path
  qualify_path with_smile
  SMILE_PREFIX=$with_smile
  smile_info=$with_smile
fi
## set SMILE_CFLAGS and SMILE_LIBS
if test "x$SMILE_PREFIX" != "x"; then
  #! shouldn't it be $SMILE_PREFIX/include?
  SMILE_CFLAGS="-I${SMILE_PREFIX}"
  #SMILE_LIBS="-L${SMILE_PREFIX}/lib"
fi
## include substitutions in output files
AC_SUBST(SMILE_PREFIX)
AC_SUBST(SMILE_CFLAGS)
#AC_SUBST(SMILE_LIBS)


## SVM_LIGHT
AC_ARG_WITH([svm_light],
	[AS_HELP_STRING([--with-svm-light=<path>], [prefix of SVM Light installation])],
	[],
	[with_svm_light=check])
## set SVM_LIGHT_PREFIX and svm_light_info (and -DNO_SVM_LIGHT if excluded)
if test "x$with_svm_light" = "xno"; then
  AC_DEFINE([NO_SVM_LIGHT], [1])
  svm_light_info="NONE"
elif test "x$with_svm_light" = "xcheck"; then
  SVM_LIGHT_PREFIX=
  AC_CHECK_LIB([svmlight], [kernel],
        [svm_light_info=installed],
        [svm_light_info="NONE (not found)"; AC_DEFINE([NO_SVM_LIGHT], [1])])
else
  #! possibly check path
  qualify_path with_svm_light
  SVM_LIGHT_PREFIX=$with_svm_light
  svm_light_info=$with_svm_light
fi
## set SVM_LIGHT_CFLAGS and SVM_LIGHT_LIBS
if test "x$SVM_LIGHT_PREFIX" != "x"; then
  #! shouldn't it be $SVM_LIGHT_PREFIX/include?
  SVM_LIGHT_CFLAGS="-I${SVM_LIGHT_PREFIX}"
  #SVM_LIGHT_LIBS="-L${SVM_LIGHT_PREFIX}/lib"
fi
AC_SUBST(SVM_LIGHT_PREFIX)
AC_SUBST(SVM_LIGHT_CFLAGS)
#AC_SUBST(SVM_LIGHT_LIBS)


sleipnirdocdir=$datadir/doc/
AC_SUBST(sleipnirdocdir)


# uncomment to print out all variables
#set


AC_CONFIG_FILES([Makefile src/Makefile])

AC_OUTPUT


#AC_CONFIG_SUBDIRS()


# Print out a summary.
echo ""
echo "  log4cpp    = $log4cpp_info"
echo "  SMILE      = $smile_info"
echo "  SVM light  = $svm_light_info"
echo ""

if echo $smile_info | grep -q "NONE"; then
  cat << EOF
BUILDING WITHOUT SMILE
SMILE is strongly recommended.
Bayes net functionality will be missing.
It is available as a binary distribution from http://genie.sis.pitt.edu/.

EOF
fi

if echo $svm_light_info | grep -q "NONE"; then
  cat << EOF
BUILDING WITHOUT SVM LIGHT
SVM light is strongly recommended.
It is available from http://svmlight.joachims.org/.

EOF
fi

echo "Now run make."

#set
