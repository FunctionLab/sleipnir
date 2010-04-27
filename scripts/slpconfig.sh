#!/bin/sh
#
CONFIGDIR=/shared/hg/sleipnir/trunk/
PREFIXDIR=/shared/proj/sleipnir/
EXTDIR=/shared/hg/sleipnir/extlib/
#
#
SMILEDIR=${EXTDIR}smile_1_1_linux64_gcc_4_1_2/
#
VWDIR=${EXTDIR}vowpal_wabbit_v4.1/
#
LOG4CPPDIR=${EXTDIR}log4cpp-1.0/
#
SVMDIR=${EXTDIR}svm_perf/
#
GENGETOPT=${EXTDIR}gengetopt-2.22/bin/gengetopt
#
BOOSTINCLDIR=${EXTDIR}boost_1_42_0/include/boost/
BOOSTGRPHDIR=${EXTDIR}boost_1_42_0/lib/
BOOSTGRPHPGM=libboost_graph.a
#
#
WITHSMILE="--with-smile=${SMILEDIR}"
WITHSVM="--with-svm-perf=${SVMDIR}"
WITHGENGET="--with-gengetopt=${GENGETOPT}"
WITHBOOSTINCL="--with-boost-includes=${BOOSTINCLDIR}"
WITHBOOSTGRPH="--with-boost-graph-lib=${BOOSTGRPHDIR}${BOOSTGRPHPGM}"
WITHVW="--with-vowpal-wabbit=${VWDIR}"
WITHLOG4CPP="--with-log4cpp=${LOG4CPPDIR}"
#
#
DATESTAMP=`date "+%Y_%m_%d_%H%M"`
LOGDIR=/shared/proj/build/
LOGFILE="${LOGDIR}sleipnir.config.${DATESTAMP}.log"
#
# ---
#
cd ${CONFIGDIR}
./configure  --prefix=${PREFIXDIR} LDFLAGS=-static ${WITHSMILE} ${WITHSVM} ${WITHGENGET} ${WITHBOOSTINCL} ${WITHBOOSTGRPH} ${WITHVW} ${WITHLOG4CPP} 2>&1 >> ${LOGFILE}
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
