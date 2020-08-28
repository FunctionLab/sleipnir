#!/bin/bash

# options: [-v] -s <seek_path> -b <seek_bin_dir> -t [long, medium, short, tiny]"

# get args into array starting with first param
# and pull out the -t for test_type option
args=("${@:1}")
for idx in ${!args[@]}; do
  if [ ${args[$idx]} == '-t' ]; then
    test_type=${args[$(($idx+1))]}
  fi
done

if [ -z $test_type ]; then
  echo "Error: must supply -t test_type [long, medium, short, tiny]"
  exit -1
fi

testname="bioinform_"$test_type

source core_run.sh "$@" -t $testname

# Make the config file
bash $testdir/make_bioinform_config.sh -t $testname -s $seekdir -b $seekbin -d $testdir -o $tmpdir

# Run the test
if [ ! -z $verbose ]; then
  echo "python $testdir/pytools/test_seek_bioinform.py -s $seekdir -b $seekbin -c $tmpdir/$testname"_config.toml" -o $tmpdir $verbose"
fi
python $testdir/pytools/test_seek_bioinform.py -s $seekdir -b $seekbin -c $tmpdir/$testname"_config.toml" -o $tmpdir $verbose
exit $?

