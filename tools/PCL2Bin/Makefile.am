bin_PROGRAMS = PCL2Bin
AM_CPPFLAGS = $(SLEIPNIR_CFLAGS) $(LOG4CPP_CFLAGS) $(PTHREAD_CFLAGS)
LDADD = $(SLEIPNIR_LIBS) $(LOG4CPP_LIBS) $(PTHREAD_LIBS)

if WITH_GENGETOPT
gengetopt = $(prog_gengetopt)
BUILT_SOURCES = cmdline.c cmdline.h
cmdline.c cmdline.h: PCL2Bin.ggo
	$(gengetopt) -iPCL2Bin.ggo --default-optional -u -N -e
endif
PCL2Bin_SOURCES = PCL2Bin.cpp cmdline.c \
	stdafx.h cmdline.h
