bin_PROGRAMS = LibSVMer
AM_CPPFLAGS = $(SLEIPNIR_CFLAGS) $(LIBSVM_CFLAGS) $(LOG4CPP_CFLAGS) $(PTHREAD_CFLAGS)
LDADD = $(SLEIPNIR_LIBS) $(LIBSVM_LIBS) $(LOG4CPP_LIBS) $(PTHREAD_LIBS)

if WITH_GENGETOPT
gengetopt = $(prog_gengetopt)
BUILT_SOURCES = cmdline.c cmdline.h
cmdline.c cmdline.h: LibSVMer.ggo
	$(gengetopt) -iLibSVMer.ggo --default-optional -u -N -e
endif
LibSVMer_SOURCES = LibSVMer.cpp stdafx.cpp cmdline.c \
	stdafx.h cmdline.h
