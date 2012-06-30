bin_PROGRAMS = Liblinearer
AM_CPPFLAGS = $(SLEIPNIR_CFLAGS) $(LOG4CPP_CFLAGS) $(PTHREAD_CFLAGS)
LDADD = $(SLEIPNIR_LIBS) $(LOG4CPP_LIBS) $(PTHREAD_LIBS)

if WITH_GENGETOPT
gengetopt = $(prog_gengetopt)
BUILT_SOURCES = cmdline.c cmdline.h
cmdline.c cmdline.h: SVMperfer.ggo
	$(gengetopt) -iSVMperfer.ggo --default-optional -u -N -e
endif
Liblinearer_SOURCES = stdafx.cpp SVMperfer.cpp cmdline.c \
	stdafx.h cmdline.h
