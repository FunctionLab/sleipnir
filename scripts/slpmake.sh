#!/bin/sh
#
HGDIR=/shared/hg/sleipnir/
MAKEDIR=/shared/hg/sleipnir/trunk
SCRIPTDIR=/shared/hg/sleipnir/scripts/
#
DATESTAMP=`date "+%Y_%m_%d_%H%M"`
LOGDIR=/shared/proj/build/
LOGFILE="${LOGDIR}sleipnir.build.${DATESTAMP}.log"
#
# ---
#
cd ${HGDIR}
INDIR=`pwd`
#
echo "[inDir] " ${INDIR}    2>&1 >> ${LOGFILE}
echo "[step:hgUpdate] "     2>&1 >> ${LOGFILE}
hg update                   2>&1 >> ${LOGFILE}
#
cd ${MAKEDIR}
INDIR=`pwd`
#
echo "[inDir] " ${INDIR}    2>&1 >> ${LOGFILE}
echo "[step:clean] "        2>&1 >> ${LOGFILE}
make -j18 clean             2>&1 >> ${LOGFILE}
rm -f tools/*/cmdline.[ch]  2>&1 >> ${LOGFILE}
#
cd ${SCRIPTDIR}
INDIR=`pwd`
#
echo "[inDir] " ${INDIR}    2>&1 >> ${LOGFILE}
echo "[step:config] "       2>&1 >> ${LOGFILE}
./slpconfig.sh              2>&1 >> ${LOGFILE}
#
cd ${MAKEDIR}
INDIR=`pwd`
#
echo "[inDir] " ${INDIR}    2>&1 >> ${LOGFILE}
echo "[step:make] "         2>&1 >> ${LOGFILE}
make -j18                   2>&1 >> ${LOGFILE}
#
echo "[inDir] " ${INDIR}    2>&1 >> ${LOGFILE}
echo "[step:install] "      2>&1 >> ${LOGFILE}
make install                2>&1 >> ${LOGFILE}
#
#
#mail -s "Sleipnir: Out-of-CycleBuild" hut-dev@hsphsun3.harvard.edu < ${LOGFILE};
mail -s "Sleipnir: Nightly Build" hut-dev@hsphsun3.harvard.edu < ${LOGFILE};
#mail -s "Sleipnir: Test Build" rmcgover < ${LOGFILE};
