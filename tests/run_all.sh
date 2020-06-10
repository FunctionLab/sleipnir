#!/bin/bash

usage() {
  echo "$(basename $0): -s <seek_path> -b <seek_bin_dir>"
}

testdir=$(dirname $0)  # the path of this script

seekdir="/data/gwallace/seek/Seek"

while getopts s:b: option; do
  case "${option}" in
    s) seekdir=${OPTARG};;
    b) seekbin=${OPTARG};;
    *) usage; exit -1;;
  esac
done

#if [ $OPTIND -eq 1 ]; then
#  echo "$0: No options were supplied";
#  usage
#  exit -1
#fi

if [ -z $seekbin ]; then
  seekbin="$seekdir/bin"
fi

failed_tests=()

# Run the tests
bash $testdir/run_paramtest.sh -s $seekdir -b $seekbin
if [ $? -ne 0 ]; then
  failed_tests+=('paramtest')  
fi

bash $testdir/run_querysize.sh -s $seekdir -b $seekbin
if [ $? -ne 0 ]; then
  failed_tests+=('querysize')  
fi

bash $testdir/run_bioinform.sh -s $seekdir -b $seekbin -t short
if [ $? -ne 0 ]; then
  failed_tests+=('bioinform -t short')  
fi

if [ ${#failed_tests[*]} -ne 0 ]; then
  echo "Tests Failed: " ${failed_tests[*]}
else
  echo "All tests passed"
fi
