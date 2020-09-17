#!/bin/bash

# options: [-v] -s <seek_path> -b <seek_bin_dir>"

testname="seek_paramtest"

source core_run.sh "$@" -t $testname

if [ ! -z $verbose ]; then
  echo "python $testdir/pytools/test_seek_params.py -s $seekdir -b $seekbin -g /tmp/$testname -o $tmpdir"
fi

# Run the test
python $testdir/pytools/test_seek_params.py -s $seekdir -b $seekbin -g /tmp/$testname -o $tmpdir
exit $?
