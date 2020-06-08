#!/bin/bash
# Script to make the configuration file necessary to run specific tests in the supplied directories.

usage() {
  echo "$(basename $0): -t <test_name> -d <testsrc_dir> -o <output_dir>"
}

# default values
seek_path="/data/gwallace/seek/Seek"
output_dir="/tmp/test_output"
test_name="bioinform_tiny"

while getopts t:s:d:o: option; do
  case "${option}" in
    t) test_name=${OPTARG};;
    # seek_path is the directory with the seek database and binaries
    s) seek_path=${OPTARG};;
    # test_src_dir is location of source code for running the tests
    d) test_src_dir=${OPTARG};;
    # output_dir is where the test output will be written
    o) output_dir=${OPTARG};;
    *) usage; exit -1;;
  esac
done

if [ $OPTIND -eq 1 ]; then 
  echo "$0: No options were supplied";
  usage
  exit -1
fi

# goldstd_path is where the gold standard expected results will be unpacked for reference
goldstd_path=/tmp/$test_name
output_file=$output_dir/$test_name"_config.toml"

if [ ! -d $output_dir ]; then
  mkdir $output_dir
fi

echo "# Configurations for running test_seek_bioinform" > $output_file

echo "# full path to seek database installation to test" >> $output_file
echo "seekPath = \"$seek_path\"" >> $output_file

echo "# full path to seek binaries to test" >> $output_file
echo "binPath = \"$seek_path/bin\"" >> $output_file

echo "# full path to the query input directory" >> $output_file
echo "queryPath = \"$goldstd_path/queries\"" >> $output_file

echo "# full path to known-good output files and gold standard files for the test" >> $output_file
echo "goldStdPath = \"$goldstd_path/results\"" >> $output_file

echo "# file with list of genes to include when evaluating results" >> $output_file
echo "geneFile = \"$test_src_dir/ref/genes.txt\"" >> $output_file

echo "# file with gene map" >> $output_file
echo "geneMap = \"$test_src_dir/ref/gene_map.txt\"" >> $output_file

echo "# output directory for the test results" >> $output_file
echo "outputPath = \"$output_dir\"" >> $output_file

echo "# Recall depth (in percent) for calculating precision" >> $output_file
echo "recallPct = 0.1" >> $output_file

echo "# Max diff between the test and gold standard results" >> $output_file
echo "maxPctDiff = 10" >> $output_file

