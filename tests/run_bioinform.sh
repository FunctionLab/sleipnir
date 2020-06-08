#!/bin/bash

usage() {
  echo "$(basename $0): -t [tiny, short, medium, long]"
}

seekdir="/data/gwallace/seek/Seek"

while getopts t: option; do
  case "${option}" in
    t) test_type=${OPTARG};;
    *) usage; exit -1;;
  esac
done

if [ $OPTIND -eq 1 ]; then
  echo "$0: No options were supplied";
  usage
  exit -1
fi

testname="bioinform_"$test_type

# Get the path of this script
testdir=$(dirname $0)

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

# Make the config file
bash $testdir/make_bioinform_config.sh -t $testname -s $seekdir -d $testdir -o $tmpdir

# Run the test
python $testdir/pytools/test_seek_bioinform.py -c $tmpdir/$testname"_config.toml" -o $tmpdir

