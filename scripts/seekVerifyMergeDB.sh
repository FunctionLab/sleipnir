#!/usr/bin/bash

usage() {
  echo "$(basename $0): [-v] -b <seek_bin_dir> -d <existing_db_dir> -n <new_db_dir> -c <combined_db_dir>"
}

while getopts b:c:d:n:v option; do
  case "${option}" in
    b) seekBinDir=${OPTARG};;  # directory with Sleipnir binaries
    c) combined=${OPTARG};;  # directory to put merged results
    d) dir1=${OPTARG};;  # main DB
    n) dir2=${OPTARG};;  # DB that was merged into it main DB
    v) verbose=" -v";;
    *) usage; exit -1;;
  esac
done

db_count=$(ls -1 $dir1/*.db | wc -l)
echo "Checking $db_count db files"
for i in $(seq -f "%08g" 0 $db_count); do 
  echo "verifying file $i.db"
  $seekBinDir/verifyMergedDBFiles -1 $dir1/$i.db -2 $dir2/$i.db -c $combined/$i.db
  if [ $? -ne 0 ]; then
    echo "ERROR: verify file $i.db failed"
    exit -1
  fi
done

echo "PASS: All files verified"