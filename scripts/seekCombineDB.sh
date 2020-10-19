#!/usr/bin/bash

# This script will combine to Seek databases together into a new combined database
#  - provide the path to the existing database -d
#  - provide the path to the new database to merge in -n
#  - provide the path to output the combined database -o

usage() {
  echo "$(basename $0): [-v] -b <seek_bin_dir> -d <existing_db_dir> -n <new_db_to_merge> -o <output_dir>"
}

while getopts b:o:d:n:v option; do
  case "${option}" in
    b) seekBinDir=${OPTARG};;  # directory with Sleipnir binaries
    o) outDir=${OPTARG};;  # directory to put merged results
    d) dir1=${OPTARG};;  # main DB
    n) dir2=${OPTARG};;  # DB to merge into it main DB
    v) verbose=" -v";;
    *) usage; exit -1;;
  esac
done

# Step 0: create the dirs and copy static files over
mkdir -p $outDir

echo "Copying in common files over..."
read -p "Press enter to continue"
# check gene_maps are equal
diff $dir1/gene_map.txt $dir2/gene_map.txt
if [ $? -ne 0 ]; then
  echo "ERROR: gene_map files differ"
  exit -1
fi
cp $dir1/gene_map.txt $outDir/

# check quant files are equal
diff $dir1/quant2 $dir2/quant2
if [ $? -ne 0 ]; then
  echo "ERROR: quant files differ"
  exit -1
fi
cp $dir1/quant2 $outDir/

# check db_lists are equal
diff $dir1/db_list $dir2/db_list
if [ $? -ne 0 ]; then
  echo "ERROR: db_lists differ"
  exit -1
fi
cp $dir1/db_list $outDir/db_list

# Check that there are no overlap in datasets between the two DBs
# Note: `unexpand -a` converts series of spaces to tab
# Note: `uniq -d` reports only duplicate lines
dups=$(cat $dir1/dataset_platform $dir2/dataset_platform | unexpand -a | sort -b | uniq -d)
if [ ! -z "$dups" ]; then 
  echo "ERROR: duplicate datasets detected: $dups"
  exit -1 
fi

echo "Combine dataset_platform text files..."
read -p "Press enter to continue"
# Step 1: Combine the dataset_platform files (sometimes called dataset.map)
# Previously had it check for and remove duplicate datasets.
# BUT this can't be done this way, can't change the order of the datasets because it is the order of entries in the DB files
# cat $dir1/dataset_platform $dir2/dataset_platform | unexpand -a | sort -b --unique > $outDir/dataset_platform
# Note: check for dataset overlap must be done separtely, such as in dataset_incremental_merge.py
cat $dir1/dataset_platform $dir2/dataset_platform > $outDir/dataset_platform

echo "Combine dataset_size text files..."
read -p "Press enter to continue"
# Step 2: Combine the dataset_size files
# Check for and remove duplicate datasets
cat $dir1/dataset_size $dir2/dataset_size > $outDir/dataset_size 

echo "Combine prep directories..."
read -p "Press enter to continue"
# step 3: Combine prep directories
mkdir -p $outDir/prep
cp $dir2/prep/* $outDir/prep/
cp $dir1/prep/* $outDir/prep/

echo "Combine sinfo directories..."
read -p "Press enter to continue"
# step 4: Combine sinfo directories
mkdir -p $outDir/sinfo
cp $dir2/sinfo/* $outDir/sinfo/
cp $dir1/sinfo/* $outDir/sinfo/

echo "Combine pcl directories..."
read -p "Press enter to continue"
# step x: Combine pcl files
mkdir -p $outDir/pclbin
cp $dir2/pclbin/* $outDir/pclbin/
cp $dir1/pclbin/* $outDir/pclbin/

echo "Run DBCombiner to combine db files..."
read -p "Press enter to continue"
# Step 5: Combine the DB files
mkdir -p $outDir/db
for f1 in $(ls -1 $dir1/db/*.db); do
    # Step through each db file and combine them one at a time
    # Substitue dir2 for dir1
    f2=$(echo $f1 | sed 's,'"$dir1"','"$dir2"',')
    echo $f2
    echo $f1 > tmp_dbcomb.txt
    echo $f2 >> tmp_dbcomb.txt
    $seekBinDir/DBCombiner -C -i $outDir/gene_map.txt --db ./tmp_dbcomb.txt -D $outDir/db
done

echo "Run SeekPrep to combine platform statistics files..."
read -p "Press enter to continue"
# Step 6: Calculate the platform averages
mkdir -p $outDir/plat
$seekBinDir/SeekPrep -f -m -i $outDir/gene_map.txt -1 $dir1/plat -2 $dir2/plat -D $outDir/plat

# Previous method of SeekPrep before the platform merge option was added
#find $outDir/db -name "*.db" > $outDir/db_list
#./bin/SeekPrep -f -P -i $outDir/gene_map.txt -b $outDir/db_list -I $outDir/prep -A $outDir/dataset.map -Q $outDir/quant2 -D $outDir/plat

