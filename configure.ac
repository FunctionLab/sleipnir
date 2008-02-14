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
  echo d $d
  if ! echo $d | egrep -q '^/'; then
    eval "$1=$topleveldir/$d"
  fi
}


## LOG4CPP
# If the user specifies a path for log4cpp, use it.
# If the user says not to use it, or if he doesn't say anything, check for it.
AC_ARG_WITH([log4cpp],
	[AS_HELP_STRING([--with-log4cpp=<path>], [prefix of log4cpp installation or build directory])],
	[if test "x$with_log4cpp" = "xno"; then with_log4cpp=check; fi],
	[with_log4cpp=check])
if test "x$with_log4cpp" = "xcheck"; then
  #!! need to look for it
  LOG4CPP_PREFIX=
  log4cpp_info="?"
else
#  if ! echo $with_log4cpp | egrep -q '^/'; then
#    with_log4cpp=$topleveldir/$with_log4cpp
#  fi
  qualify_path with_log4cpp
  LOG4CPP_PREFIX=$with_log4cpp
  log4cpp_info=$with_log4cpp
fi
if test "x$LOG4CPP_PREFIX" = "x"; then
  AC_DEFINE([USE_LOG4CPP_STUB], [1])
  log4cpp_info=stub
else
  # It may be useful to allow a built (but not installed) path for log4cpp,
  # however the log4cpp distribution does not have a lib/, it has a src/.libs/.
  if test -d ${LOG4CPP_PREFIX}/lib; then
    LOG4CPP_LIBS="-L${LOG4CPP_PREFIX}/lib"
  elif test -d ${LOG4CPP_PREFIX}/src/.libs; then
    LOG4CPP_LIBS="-L${LOG4CPP_PREFIX}/src/.libs"
  else
    #!! should use something besides just echo?
#    echo "LOG4CPP LIB NOT FOUND - USING STUB"
    AC_DEFINE([USE_LOG4CPP_STUB], [1])
    log4cpp_info="stub (lib not found)"
  fi
fi
AC_SUBST(LOG4CPP_PREFIX)
LOG4CPP_CFLAGS="-I${LOG4CPP_PREFIX}/include"
AC_SUBST(LOG4CPP_LIBS)
AC_SUBST(LOG4CPP_CFLAGS)


## SMILE
AC_ARG_WITH([smile],
	[AS_HELP_STRING([--with-smile=<path>], [prefix of SMILE installation])],
	[],
	[with_smile=check])
if test "x$with_smile" = "xno"; then
  smile_info="NONE"
  AC_DEFINE([NO_SMILE], [1])
elif test "x$with_smile" = "xcheck"; then
  #!! check for it
  SMILE_PREFIX=
  smile_info="?"
else
#  if ! echo $with_smile | egrep -q '^/'; then
#    with_smile=$topleveldir/$with_smile
#  fi
  qualify_path with_smile
  SMILE_PREFIX=$with_smile
  smile_info=$with_smile
fi
AC_SUBST(SMILE_PREFIX)
if test "x$SMILE_PREFIX" != "x"; then
  #!! $SMILE_PREFIX/include?
  SMILE_CFLAGS="-I${SMILE_PREFIX}"
  #SMILE_LIBS="-L${SMILE_PREFIX}/lib/"
fi
#AC_SUBST(SMILE_LIBS)
AC_SUBST(SMILE_CFLAGS)


## SVM_LIGHT
AC_ARG_WITH([svm_light],
	[AS_HELP_STRING([--with-svm-light=<path>], [prefix of SVM Light installation])],
	[],
	[with_svm_light=check])
if test "x$with_svm_light" = "xno"; then
  svm_light_info="NONE"
  AC_DEFINE([NO_SVM_LIGHT], [1])
elif test "x$with_svm_light" = "xcheck"; then
  #!! check for it
  SVM_LIGHT_PREFIX=
  svm_light_info="?"
else
#  if ! echo $with_svm_light | egrep -q '^/'; then
#    with_svm_light=$topleveldir/$with_svm_light
#  fi
  qualify_path with_svm_light
  SVM_LIGHT_PREFIX=$with_svm_light
  svm_light_info=$with_svm_light
fi
AC_SUBST(SVM_LIGHT_PREFIX)
if test "x$SVM_LIGHT_PREFIX" != "x"; then
  #!! $SVM_LIGHT_PREFIX/include?
  SVM_LIGHT_CFLAGS="-I${SVM_LIGHT_PREFIX}"
  #SVM_LIGHT_LIBS="-L${SVM_LIGHT_PREFIX}/lib/"
fi
#AC_SUBST(SVM_LIGHT_LIBS)
AC_SUBST(SVM_LIGHT_CFLAGS)


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

sleipnirdocdir=$datadir/doc/
AC_SUBST(sleipnirdocdir)

#set

AC_CONFIG_FILES([Makefile src/Makefile])

# deprecated
#AC_OUTPUT([./src/Makefile])
AC_OUTPUT

#AC_CONFIG_SUBDIRS()

# Print out a summary.
echo ""
echo "  log4cpp    = $log4cpp_info"
echo "  SMILE      = $smile_info"
echo "  SVM light  = $svm_light_info"
echo ""

if test "x$smile_info" = "xNONE"; then
  cat << EOF
BUILDING WITHOUT SMILE
SMILE is strongly recommended.
Bayes net functionality will be missing.
It is available as a binary distribution from http://genie.sis.pitt.edu.

EOF
fi

if test "x$svm_light_info" = "xNONE"; then
  cat << EOF
BUILDING WITHOUT SVM LIGHT
SVM light is strongly recommended.
It is available from http://svmlight.joachims.org.

EOF
fi

echo "Now run make."

#set
