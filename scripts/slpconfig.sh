#!/bin/sh
#
CONFIGDIR=/shared/hg/sleipnir/trunk/
PREFIXDIR=/shared/proj/sleipnir/
EXTDIR=/shared/hg/sleipnir/extlib/
#
#
SMILEDIR=${EXTDIR}smile_1_1_linux64_gcc_4_1_2/
#
SVMDIR=${EXTDIR}svm_perf/
#
GENGETDIR=${PREFIXDIR}
#
BOOSTINCLDIR=${PREFIXDIR}include/boost/
BOOSTGRPHDIR=${PREFIXDIR}lib/
BOOSTGRPHPGM=libboost_graph.a
#
#
WITHSMILE="--with-smile=${SMILEDIR}"
WITHSVM="--with-svm-perf=${SVMDIR}"
WITHGENGET="--with-gengetopt=${GENGETDIR}"
WITHBOOSTINCL="--with-boost-includes=${BOOSTINCLDIR}"
WITHBOOSTGRPH="--with-boost-graph-lib=${BOOSTGRPHDIR}${BOOSTGRPHPGM}"
#
#
DATESTAMP=`date "+%Y_%m_%d_%H%M"`
LOGDIR=/shared/proj/build/
LOGFILE="${LOGDIR}sleipnir.config.${DATESTAMP}.log"
#
# ---
#
cd ${CONFIGDIR}
./configure  --prefix=${PREFIXDIR} LDFLAGS=-static ${WITHSMILE} ${WITHSVM} ${WITHGENGET} ${WITHBOOSTINCL} ${WITHBOOSTGRPH} 2>&1 >> ${LOGFILE}
#
# ---
#notes:
# - gengetopt installed to prefixdir /rm.00323
# - current versions used:
#   boost:     boost_1_42_0
#   gengetopt: gengetopt-2.22
#   log4cpp:   log4cpp-1.0
#   smile:     smile_1_1_linux64_gcc_4_1_2
#   svm:
#   /rm.00401
