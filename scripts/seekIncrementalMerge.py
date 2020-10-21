import os
import argparse
import glob
from datetime import datetime
from seekCreateDB import createSeekDB

# This script will create a new database from PCL files that matches the parameters
#   of an existing small and large database.
# It will then prepare the data for combining into the small database.
# At the end it will print the command to run to merge the new database 
#   into the small database.

# Notes:
# If you want to create smallDB, then specify the same referenace DB for large and small db inputs
# If you want to create a newDB, then specify both small and large DB for reference to make sure no overlap in datasets
# If you want to merge smallDB to largeDB - just seekDBCombine.sh script


def loadDatasetPlatformMap(filename):
    if not os.path.exists(filename):
      return None
    dsetDict = {}
    with open(filename, 'r') as fh:
      for line in fh:
        vals = line.strip().split()
        if len(vals) == 3:
          # Likely of type (pcl_file, name, platform)
          assert 'pcl' in vals[0].lower(), "First of 3 columns should specify pcl file name"
          assert 'pcl' not in vals[1].lower(), "Reference to pcl found in 2nd column"
          assert 'pcl' not in vals[2].lower(), "Reference to pcl found in 3rd column"
          dsetDict[vals[1]] = vals[2]
        elif len(vals) == 2:
          # Likely of type (pcl_file, name, platform)
          assert 'pcl' not in line.lower(), "check dataset file, seems to reference pcl files"
          dsetDict[vals[0]] = vals[1]
        else:
          assert len(vals) == 2 or len(vals)==3, "Incorrect number of columns in dataset file, expect 2 or 3" 
    return dsetDict


def main(args):
  # Step 1: Load the dataset lists for the new, small and large DBs and make sure the are all disjoint
  # check existence of the dataset map files
  if not os.path.exists(args.smallDsetFile):
    raise FileNotFoundError("Small DB dataset_platform map not found: " + args.smallDsetFile)
  if not os.path.exists(args.largeDsetFile):
    raise FileNotFoundError("Large DB dataset_platform map not found: " + args.largeDsetFile)
  if not os.path.exists(args.newDsetFile):
    raise FileNotFoundError("New dataset_platform map not found: " + args.newDsetFile)

  smallDsetDict = loadDatasetPlatformMap(args.smallDsetFile)
  largeDsetDict = loadDatasetPlatformMap(args.largeDsetFile)
  newDsetDict = loadDatasetPlatformMap(args.newDsetFile)

  smallDsets = set(smallDsetDict.keys())
  largeDsets = set(largeDsetDict.keys())
  newDsets = set(newDsetDict.keys())
  
  large_small_overlap = largeDsets.intersection(smallDsets)
  large_new_overlap = largeDsets.intersection(newDsets)
  small_new_overlap = smallDsets.intersection(newDsets)

  # Datasets must be disjoint between all the databases
  assert not large_small_overlap, "large and small datasets overlap " + str(large_small_overlap)
  assert not large_new_overlap, "large and new datasets overlap " + str(large_new_overlap)
  assert not small_new_overlap, "small and new datasets overlap " + str(small_new_overlap)

  # Step 2:
  # Make a new DB from the dataset PCL files
  # check how many db files large and small db have
  assert os.path.isdir(args.dirLargeDB)
  assert os.path.isdir(args.dirSmallDB)
  assert os.path.isdir(args.dirLargeDB+"/db")
  assert os.path.isdir(args.dirSmallDB+"/db")

  numDBFiles = len(glob.glob1(args.dirLargeDB+"/db", "*.db"))
  numSmallDBFiles = len(glob.glob1(args.dirSmallDB+"/db", "*.db"))
  assert numDBFiles == numSmallDBFiles, "numDBFiles mismatch between large and small db"
  # assert numDBFiles == 1000

  dateStr = datetime.now().strftime("%Y%m%d_%H%M%S")
  incrDBDirName = os.path.join(args.outDir, "incr_dset_" + dateStr)
  res = createSeekDB(sleipnirBinDir=args.sleipnirBinDir,
                     inputDatasetFile=args.newDsetFile,
                     pclDir=args.dirNewPCL,
                     refDir=args.dirLargeDB,
                     output_dir=incrDBDirName,
                     numDBFiles=numDBFiles)
  assert res == True, "createSeekDB failed"
  print(f'Incremental database created in {incrDBDirName}')

  # Step 3:
  # Merge the new DB into the small DB using bash script
  # Output directory must be specified
  # Recommend after completion of merge do the following by hand:
  #   Rename small DB directory to prev.num
  #   Rename combined DB directory to small DB name
  scriptPath = os.path.dirname(os.path.realpath(__file__))

  print("Incremental Dataset Build Complete. Ready for merge!")
  print('*** NEXT STEP ***')
  print("Run the following command to perform the merge")
  print('### COMMAND ###')
  print(f'bash {scriptPath}/seekCombineDB.sh -o {args.outDir} -b {args.sleipnirBinDir} -d {args.dirSmallDB} -n {incrDBDirName}')

  print('\nOptional final step: verify the merged db contests')
  print('### COMMAND ###')
  print(f'bash {scriptPath}/seekVerifyMergeDB.sh -b {args.sleipnirBinDir} -d {args.dirSmallDB}/db -n {incrDBDirName}/db -c {args.outDir}/db')
  # Step 4: 
  # Check size of small DB relative to large DB, and recommend combining at some size/percentage threshold



if __name__ == "__main__":
    argParser = argparse.ArgumentParser()
    argParser.add_argument('--dirNewPCL', '-p', type=str, required=True,
                           help='directory containing the PCL files for the new datasets')
    argParser.add_argument('--newDsetFile', '-dn', type=str, required=True,
                           help='text file listing the new datasets and corresponding platforms, one per line')
    argParser.add_argument('--largeDsetFile', '-dl', type=str, required=True,
                           help='text file listing the large DB datasets and corresponding platforms, one per line')
    argParser.add_argument('--smallDsetFile', '-ds', type=str, required=True,
                           help='text file listing the small DB datasets and corresponding platforms, one per line')
    argParser.add_argument('--dirSmallDB', '-s', type=str, required=True,
                           help='directory of existing small DB')
    argParser.add_argument('--dirLargeDB', '-l', type=str, required=True,
                           help='directory of existing large DB')
    argParser.add_argument('--sleipnirBinDir', '-b', type=str, required=True,
                           help='Directory where the Sleipnir binaries are installed')
    argParser.add_argument('--outDir', '-o', type=str, required=True,
                           help='Output directory to write new database into')
    args = argParser.parse_args()
    main(args)


## Previous
  # smallDsetDict = loadDatasetPlatformMap(os.path.join(args.dirSmallDB, "dataset.map"))
  # if smallDsetDict is None:
  #   smallDsetDict = loadDatasetPlatformMap(os.path.join(args.dirSmallDB, "dataset_platform.txt"))

  # smallDsetFile = os.path.join(args.dirSmallDB, "dataset_platform.txt")
  # largeDsetFile = os.path.join(args.dirLargeDB, "dataset_platform.txt")
