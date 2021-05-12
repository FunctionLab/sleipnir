"""
This script implements an incremental merge strategy, where there is a largeDB,
and a smallDB. In general new PCL dataset files are merged into the smallDB (in
order to speed up the merge). When the smallDB grows to some threshold, say 10%
of the size of the largeDB, then the small and largeDB are merged together.

This script handles the common path of merging a new set of PCL dataset files
into a smallDB.

Inputs:
dirLargeDB: directory of the largeDB - for reference only and to validate
  that the datasets are disjoint and gene maps match.

dirSmallDB: directory of the smallDB - the db files and metadata that the
  new PCL data will be combined into.

dirNewPCL: the directory containing the new PCL files to merge in.

newDsetFile: the list of the new datasets to merge in, the first column is
  the name of the pcl files found in dirNewPCL

largeDsetFile: list of datasets in the largeDB

smallDsetFile: list of datasets in the smallDB

sleipnirBinDir: sleipnir binaries

outDir: where the final merge of smallDB and newPCL will be written (somewhere other
  than the smallDB dir so that verification can happen)

Example usage:
conda activate genomics
python seekIncrementalMerge.py -p <newPclDir> -dn <newDatasetList> \
  -dl <largeDatasetList> -ds <smallDatasetList> -l <largeDBDir> \
  -s <smallDBDir> -b <sleipnirBinaries> -o <mergedDir> \
"""

import os
import sys
import argparse
import glob
import subprocess
from datetime import datetime
currPath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(currPath)
import seekUtils as sutils
from seekCreateDB import createSeekDB

# This script will create a new database from PCL files
# It will use the parameters of an existing small and large database (which must agree).
# It will then prepare the data for combining into the small database.
# At the end it will print the command to run to merge the new database into the small database.

# Notes:
# If you just want to create a smallDB, use the seekCreateDB.py script

# TODO - make it optional to specify a newPCL diretory. If none is specified
# then merge the small and large db.


def copyFile(fileName, srcDir, destDir):
  src = os.path.join(srcDir, fileName)
  dst = os.path.join(destDir, fileName)
  ret = subprocess.run(f'cp {src} {dst}', shell=True)
  assert ret.returncode == 0

def checkFilesMatch(fileName, dir1, dir2):
  f1 = os.path.join(dir1, fileName)
  f2 = os.path.join(dir2, fileName)
  ret = subprocess.run(f'diff {f1} {f2}', shell=True)
  assert ret.returncode == 0

def concatenateFiles(fileName, dir1, dir2, outDir):
  f1 = os.path.join(dir1, fileName)
  f2 = os.path.join(dir2, fileName)
  dst = os.path.join(outDir, fileName)
  ret = subprocess.run(f'cat {f1} {f2} > {dst}', shell=True)
  assert ret.returncode == 0

def combineDirs(subDirName, dir1, dir2, outDir):
  d1 = os.path.join(dir1, subDirName)
  d2 = os.path.join(dir2, subDirName)
  dst = os.path.join(outDir, subDirName)
  os.makedirs(dst, exist_ok=True)
  ret = subprocess.run(f'cp -a {d1}/* {dst}', shell=True)
  assert ret.returncode == 0
  ret = subprocess.run(f'cp -a {d2}/* {dst}', shell=True)
  assert ret.returncode == 0

def main(args):
    refCfg = sutils.getDefaultConfig()
    refCfg.datasetsFile = os.path.basename(args.smallDsetFile)

    os.makedirs(args.outDir, exist_ok=True)

    # STEP 01: Load the dataset lists for the new, small and large DBs and make sure
    #   they are all disjoint
    # check existence of the dataset map files
    if not os.path.exists(args.smallDsetFile):
        raise FileNotFoundError("Small DB dataset_platform map not found: " + args.smallDsetFile)
    if not os.path.exists(args.largeDsetFile):
        raise FileNotFoundError("Large DB dataset_platform map not found: " + args.largeDsetFile)
    if not os.path.exists(args.newDsetFile):
        raise FileNotFoundError("New dataset_platform map not found: " + args.newDsetFile)

    smallDsets = sutils.readDatasetList(args.smallDsetFile)
    largeDsets = sutils.readDatasetList(args.largeDsetFile)
    newDsets = sutils.readDatasetList(args.newDsetFile)

    smallDsets = set([dset[1] for dset in smallDsets])
    largeDsets = set([dset[1] for dset in largeDsets])
    newDsets = set([dset[1] for dset in newDsets])

    large_small_overlap = largeDsets.intersection(smallDsets)
    large_new_overlap = largeDsets.intersection(newDsets)
    small_new_overlap = smallDsets.intersection(newDsets)

    # Datasets must be disjoint between all the databases
    assert not large_small_overlap, "large and small datasets overlap " + str(large_small_overlap)
    assert not large_new_overlap, "large and new datasets overlap " + str(large_new_overlap)
    assert not small_new_overlap, "small and new datasets overlap " + str(small_new_overlap)

    # STEP 02: Do some checks
    # check that large and small DBs have the same number/name of db files
    assert os.path.isdir(args.dirLargeDB)
    assert os.path.isdir(args.dirSmallDB)
    largeDBFileDir = os.path.join(args.dirLargeDB, refCfg.dbDir)
    smallDBFileDir = os.path.join(args.dirSmallDB, refCfg.dbDir)
    assert os.path.isdir(largeDBFileDir)
    assert os.path.isdir(smallDBFileDir)
    largeDBFiles = glob.glob1(largeDBFileDir, "*.db")
    largeDBFiles.sort()
    smallDBFiles = glob.glob1(smallDBFileDir, "*.db")
    smallDBFiles.sort()
    assert largeDBFiles == smallDBFiles
    # assert numDBFiles == numSmallDBFiles, "numDBFiles mismatch between large and small db"

    checkFilesMatch(refCfg.geneMapFile, args.dirLargeDB, args.dirSmallDB)
    checkFilesMatch(refCfg.quantFile, args.dirLargeDB, args.dirSmallDB)

    # STEP 03: Create a new database from the new incremental pcl files
    #   make a temp directory to hold the incremental database
    dateStr = datetime.now().strftime("%Y%m%d_%H%M%S")
    incrDBDirName = os.path.join(args.outDir, "incr_dset_" + dateStr)
    os.makedirs(incrDBDirName)
    # copy over geneMap and quant files needed from refDB
    copyFile(refCfg.geneMapFile, args.dirLargeDB, incrDBDirName)
    copyFile(refCfg.quantFile, args.dirLargeDB, incrDBDirName)
    #   set the configs
    newCfg = sutils.getDefaultConfig()
    newCfg.binDir = args.sleipnirBinDir
    newCfg.pclDir = args.dirNewPCL
    newCfg.datasetsFile = args.newDsetFile
    newCfg.inDir = incrDBDirName
    newCfg.outDir = incrDBDirName
    newCfg.numDbFiles = len(largeDBFiles)
    sutils.checkConfig(newCfg)
    #  create the db
    res = createSeekDB(newCfg, None, runAll=True, concurrency=8)
    assert res == True, "createSeekDB failed"
    print(f'Incremental database created in {incrDBDirName}')

    # STEP 04: Combine metadata
    copyFile(refCfg.geneMapFile, args.dirLargeDB, args.outDir)
    copyFile(refCfg.quantFile, args.dirLargeDB, args.outDir)
    dsetFileBaseName = os.path.basename(args.smallDsetFile)
    concatenateFiles(dsetFileBaseName, args.dirSmallDB, incrDBDirName, args.outDir)
    concatenateFiles(refCfg.dsetSizeFile, args.dirSmallDB, incrDBDirName, args.outDir)
    combineDirs('prep', args.dirSmallDB, incrDBDirName, args.outDir)
    combineDirs('sinfo', args.dirSmallDB, incrDBDirName, args.outDir)
    combineDirs('pclbin', args.dirSmallDB, incrDBDirName, args.outDir)


    # STEP 05: Run the db combiner
    if args.yesToPrompts is False:
      reply = input('Proceed with dbCombiner command? ' + '(y/n): ')
      reply.lower().strip()
      if reply[0] != 'y':
        return -1
    mergedCfg = sutils.getDefaultConfig()
    mergedCfg.binDir = args.sleipnirBinDir
    mergedCfg.inDir = args.dirSmallDB
    mergedCfg.outDir = args.outDir
    mergedCfg.datasetsFile = dsetFileBaseName
    mergedCfg.numDbFiles = len(largeDBFiles)
    sutils.checkConfig(mergedCfg)
    dbDirsToCombine = [smallDBFileDir, newCfg.dbDir]
    sutils.parallelCombineThreadDBs(mergedCfg, dbDirsToCombine, concurrency=8)

    # STEP 06: Run verification of the combined DB files
    if args.yesToPrompts is False:
      reply = input('Proceed with db verification command? ' + '(y/n): ')
      reply.lower().strip()
      if reply[0] != 'y':
        return -1
    sutils.verifyCombinedDBs(mergedCfg, dbDirsToCombine, concurrency=8)

    print("Run SeekPrep to combine platform statistics files...")
    # STEP 07: Calculate the platform averages
    smallPlatDir = os.path.join(args.dirSmallDB, 'plat')
    incrPlatDir = os.path.join(incrDBDirName, 'plat')
    outPlatDir = os.path.join(args.outDir, 'plat')
    os.makedirs(outPlatDir, exist_ok=True)
    cmd = f'{args.sleipnirBinDir}/SeekPrep -f -m -i {mergedCfg.geneMapFile} ' \
          f'-1 {smallPlatDir} -2 {incrPlatDir} -D {outPlatDir}'
    ret = subprocess.run(cmd, shell=True)
    assert ret.returncode == 0

    # FINALIZATION STEPS
    # Recommend after completion of merge do the following by hand:
    #   Rename small DB directory to prev.num
    #   Rename combined DB directory to small DB name
    # Check size of small DB relative to large DB, and recommend combining
    #   at some size/percentage threshold.
    return 0

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
    argParser.add_argument('--yesToPrompts', '-y', default=False, action='store_true',
                           help='Answer yes to all prompts')
    args = argParser.parse_args()
    res = main(args)
    sys.exit(res)
