import os
import sys
import glob
import argparse
import seekUtils as utils

"""
Program to copy pcl files for datasets which appear in the dataset_platform_map
from a source diretory into an output pcl directory and name them in the
format experiment.platform.pcl.
Note:
1. Pcl files in the source directory that aren't in the dset_platform_map will
be copied into the excludedPcls subdirectory.
2. Datasets in the dset_platform_map with no corresponding pcl file in the source
directory will be listed in the missingPclFiles.txt file
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--dsetMap', '-m', required=True,
                        help="dset_plat.map file with datasets to copy")
    parser.add_argument('--inDir', '-i', required=True,
                        help="directory with original pcl files")
    parser.add_argument('--outDir', '-o', required=True,
                        help="directory to output renamed pcl files")
    args = parser.parse_args()

    cmdType = 'cp'
    if args.inDir == args.outDir:
        reply = input('input and output directory are the same, rename files inplace? (y/n): ')
        reply.lower().strip()
        if reply[0] != 'y':
            sys.exit(-1)
        cmdType = 'mv'

    # Get the dataset platform map
    dsetPlatMap = utils.readPlatMap(args.dsetMap)
    dsetSet = set(dsetPlatMap.keys())

    # Get list of pcl files in input directory
    pclList = glob.glob1(args.inDir, '*.pcl')
    pclSet = set([pclFile.split('.')[0] for pclFile in pclList])

    # Intersection of sets from pcl files and dataset list
    dsetNames = pclSet.intersection(dsetSet)

    # pcl files not found in dsetSet list
    excludedPclFiles = pclSet - dsetSet

    # datasets with no corresponding Pcl file
    missingPclFiles = dsetSet - pclSet

    if len(excludedPclFiles) > 0:
        # put excluded Pcl files in a sub-directory
        excludedDir = os.path.join(args.outDir, 'excludedPcls')
        os.makedirs(excludedDir, exist_ok=True)
        for pcl in excludedPclFiles:
            fileName = pcl + '.pcl'
            inFile = os.path.join(args.inDir, fileName)
            outFile = os.path.join(excludedDir, fileName)
            cmd = f'{cmdType} {inFile} {outFile}'
            os.system(cmd)

    if len(missingPclFiles) > 0:
        # write out a list of datasets with no corresponding pcl file
        missingFile = os.path.join(args.outDir, 'missingPclFiles.txt')
        with open(missingFile, 'w') as fp:
            fp.write("\n".join(missingPclFiles))

    datasetList = []
    mapList = []
    # copy or move the pcl files to the output directory
    for dset in dsetNames:
        platform = dsetPlatMap[dset]
        inFile = os.path.join(args.inDir, dset + '.pcl')
        outFile = os.path.join(args.outDir, f'{dset}.{platform}.pcl')
        datasetList.append(f'{dset}.{platform}.pcl')
        mapList.append(f'{dset}.{platform}\t{platform}')
        cmd = f'{cmdType} {inFile} {outFile}'
        os.system(cmd)

    # write out the dataset list
    with open(os.path.join(args.outDir, 'datasets.txt'), 'w') as fp:
        fp.write("\n".join(datasetList))

    # write out revised dset_platform_map
    platMapFilename = os.path.basename(args.dsetMap)
    with open(os.path.join(args.outDir, platMapFilename), 'w') as fp:
        fp.write("\n".join(mapList))
