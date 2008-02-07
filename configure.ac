AC_INIT(libSleipnir, 0.9.0, ogt@princeton.edu)
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/stdafx.h])
# AC_CONFIG_HEADER([config.h])

AC_ARG_WITH(log4cpp,
	[  --with-log4cpp=<path>   prefix of log4cpp installation],
	[LOG4CPP_PREFIX=$with_log4cpp],
	[LOG4CPP_PREFIX="/usr/"])
AC_SUBST(LOG4CPP_PREFIX)
LOG4CPP_LIBS="-L${LOG4CPP_PREFIX}/lib/"
LOG4CPP_CFLAGS="-I${LOG4CPP_PREFIX}/include/"
AC_SUBST(LOG4CPP_LIBS)
AC_SUBST(LOG4CPP_CFLAGS)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

# Checks for libraries.
AC_CHECK_LIB(log4cpp, log4cpp::Category::getInstance, [],
	[echo "Error: log4cpp not found"; exit -1])

# Checks for header files.
AC_HEADER_DIRENT
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h netinet/in.h sys/socket.h sys/time.h])

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
AC_CHECK_FUNCS([memmove memset munmap pow socket sqrt strchr strerror strrchr strstr])

AC_SUBST(sleipnirdocdir)
sleipnirdocdir=$datadir/doc/

AC_CONFIG_FILES([Makefile])
AC_OUTPUT([./src/Makefile])
