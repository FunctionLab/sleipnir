#!/bin/bash
# Script to make the configuration file necessary to run specific tests in the supplied directories.

usage() {
  echo "$(basename $0): -t <test_name> -d <testsrc_dir> -o <output_dir>"
}

# default values
seek_path="/data/gwallace/seek/Seek"
output_dir="/tmp/test_output"
test_name="bioinform_tiny"

while getopts t:s:b:d:o:g: option; do
  case "${option}" in
    t) test_name=${OPTARG};;
    # seek_path is the directory with the seek database
    s) seek_path=${OPTARG};;
    # seek_bin is the directory with the seek binaries
    b) seek_bin=${OPTARG};;
    # test_src_dir is location of source code for running the tests
    d) test_src_dir=${OPTARG};;
    # path to the gold standard result files
    g) goldstd_path=${OPTARG};;
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

source $seek_path/seek_env

if [ -z $seek_bin ]; then
  seek_bin="$seek_path/bin"
fi

if [ -z $goldstd_path ]; then
  # goldstd_path is where the gold standard expected results will be unpacked for reference
  goldstd_path=/tmp/$USER/$test_name
fi

if [ -z $GENE_MAP_FILE ]; then
  # gene_map_file env var should be set in seek_env
  GENE_MAP_FILE="gene_map.txt"
fi

gene_map_full_path="$seek_path/$GENE_MAP_FILE"

output_file=$output_dir/$test_name"_config.toml"

if [ ! -d $output_dir ]; then
  mkdir $output_dir
fi

echo "# Configurations for running test_seek_bioinform" > $output_file

echo "# full path to seek database installation to test" >> $output_file
echo "seekPath = \"$seek_path\"" >> $output_file

echo "# full path to seek binaries to test" >> $output_file
echo "binPath = \"$seek_bin\"" >> $output_file

echo "# full path to the query input directory" >> $output_file
echo "queryPath = \"$goldstd_path/queries\"" >> $output_file

echo "# full path to known-good output files and gold standard files for the test" >> $output_file
echo "goldStdPath = \"$goldstd_path/results\"" >> $output_file

echo "# file with gene map" >> $output_file
echo "geneMap = \"$gene_map_full_path\"" >> $output_file

echo "# file with list of genes to include when evaluating results" >> $output_file
echo "geneFile = \"$output_dir/include_genes.txt\"" >> $output_file

echo "# output directory for the test results" >> $output_file
echo "outputPath = \"$output_dir\"" >> $output_file

echo "# Recall depth (in percent) for calculating precision" >> $output_file
echo "recallPct = 0.1" >> $output_file

echo "# Max diff between the test and gold standard results" >> $output_file
echo "maxPctDiff = 10" >> $output_file

if [ ! -f $gene_map_full_path ]; then
  echo "Error: gene_map.txt file not found: $gene_map_full_path"
  exit -1
fi
# create the include genes by keeping the second column of gene_map.txt and 
#  changing newlines to spaces.
cut -f2 $gene_map_full_path | tr '\n' ' ' > $output_dir/include_genes.txt

