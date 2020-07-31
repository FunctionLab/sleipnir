#!/bin/bash

# usage:  echo "run_all.sh [-v] -s <seek_path> -b <seek_bin_dir>"

testdir=$(dirname $(realpath -s $0))  # the path of this script

failed_tests=()

# Run the tests
bash $testdir/run_paramtest.sh "$@"
if [ $? -ne 0 ]; then
  failed_tests+=('paramtest')  
fi

bash $testdir/run_querysize.sh "$@"
if [ $? -ne 0 ]; then
  failed_tests+=('querysize')  
fi

bash $testdir/run_bioinform.sh "$@" -t short
if [ $? -ne 0 ]; then
  failed_tests+=('bioinform -t short')  
fi

if [ ${#failed_tests[*]} -ne 0 ]; then
  echo "Tests Failed: " ${failed_tests[*]}
else
  echo "All tests passed"
fi
