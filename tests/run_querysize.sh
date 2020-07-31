#!/bin/bash

# options: [-v] -s <seek_path> -b <seek_bin_dir>"

testname="seek_querysize"

source core_run.sh "$@" -t $testname

# Run the test
if [ ! -z $verbose ]; then
  echo "python $testdir/pytools/test_seek_querysizes.py -s $seekdir -b $seekbin -g /tmp/$testname -o $tmpdir $verbose"
fi
python $testdir/pytools/test_seek_querysizes.py -s $seekdir -b $seekbin -g /tmp/$testname -o $tmpdir $verbose
exit $?
