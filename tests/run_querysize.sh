#!/bin/bash

testname="seek_querysize"

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

# Requires a seek_env file that sets LD_LIBRARY_PATH for seek binaries
source $seekdir/seek_env

# Requires the conda environment is setup
# use "conda env create --file conda_environment.yml" in the sleipnir/tests/pytools dir
source ~/.bashrc
if [ -z $CONDA_DEFAULT_ENV ] || [ $CONDA_DEFAULT_ENV != "genomics" ]; then
  conda activate genomics
fi

# Make a temp output directory for the results
tmpdir=$(mktemp -d -t test_output_$testname-$(date +%Y-%m-%d-%H-%M-%S)-XXXX)

# Check if the gold standard result files are already unpacked
if [ ! -d /tmp/$testname ]; then
  # Unpack the gold standard result files
  tar xzvf $testdir/gold_standard_results/$testname"_goldstandard.tgz" --directory /tmp
fi

# Run the test
python $testdir/pytools/test_seek_querysizes.py -s $seekdir -b $seekbin -g /tmp/$testname -o $tmpdir
exit $?
