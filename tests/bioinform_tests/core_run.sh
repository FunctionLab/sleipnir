#!/bin/bash

usage() {
  echo "$(basename $0): [-v] -s <seek_path> -b <seek_bin_dir> -t <testname>"
}

testdir=$(dirname $(realpath -s $0))  # the path of this script

seekdir="/data/gwallace/seek/Seek"

while getopts s:b:t:v option; do
  case "${option}" in
    s) seekdir=${OPTARG};;
    b) seekbin=${OPTARG};;
    t) testname=${OPTARG};;
    v) verbose=" -v";;
    *) usage; exit -1;;
  esac
done

if [ -z $testname ]; then
  echo "Inner script error: $0: must provide testname";
  usage
  exit -1
fi

if [ -z $seekbin ]; then
  seekbin="$seekdir/bin"
fi

# Requires a seek_env file that sets LD_LIBRARY_PATH for seek binaries
if [ ! -f $seekdir/seek_env ]; then
  echo "Error: File seek_env not found in dir $seekdir."
  echo "Must be run where the SEEK DB is intalled and seek_env is present"
  usage
  exit -1
fi
echo "source $seekdir/seek_env"
source $seekdir/seek_env

# Requires the conda environment is setup
# use "conda env create --file conda_environment.yml" in the sleipnir/tests/pytools dir
source ~/.bashrc
if [ -z $CONDA_DEFAULT_ENV ] || [ $CONDA_DEFAULT_ENV != "genomics" ]; then
  conda activate genomics
fi

# Make a temp output directory for the results
tmpdir=$(mktemp -d -t test_output_$testname-$(date +%Y-%m-%d-%H-%M-%S)-XXXX)

if [ ! -d /tmp/$USER ]; then
    mkdir /tmp/$USER
fi

goldStdDir=/tmp/$USER/$testname

# location of the tar packed gold standard files in the repositoy
goldTarball=$testdir/gold_standard_results/$testname"_goldstandard.tgz"
# size of the gold standard tar files
tarBallSize=$(du -s -B 1 $goldTarball | cut -f1)
# size of tmp unpacked test directory
workSize=$(du -s -B 1 $goldStdDir | cut -f1)

# Check if the gold standard result files are already unpacked and complete
if [ ! -d $goldStdDir ] || [ $workSize -lt $tarBallSize ]; then
  # Unpack the gold standard result files
  tar xzvf $goldTarball --directory /tmp/$USER
fi
