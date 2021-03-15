"""
Create a new Seek DB from pcl files. The inputs are pcl files, a gene_map.txt file
which is the ordered list of genes to include, a quant file which is the binning
values to map floats to 8 bits values (bins), and a dataset list file which
includes the pcl filename for the dataset and the platform type.

When specifying the --all option, the output is a full Seek DB, including DB files,
prep, plat, sinfo, dab and pclbin directories.

Note: requires the python conda environment specified in the conda_environment.yml
To create the env 'conda env create --file conda_environment.yml'

Example usage:
conda activate genomics
python seekCreateDB -b <path_sleipnir_binaries> -i <input_dir> -o <output_dir> --all

Where the gene_map.txt and quant2 and pcl files can be found in the input_dir.
"""
# This script is an adaptation of the prepare_seek.py script that
# comes with the breast cancer example dataset on the Seek website
import os
import re
import sys
import resource
import argparse
import subprocess
currPath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(currPath)
import seekUtils as sutils
from structDict import StructDict

# TODO - command line args: --all --sinfo --prep --plat --db --dab --dsetsize --pclbin

def createSeekDB(cfg, tasksToRun, runAll=False, concurrency=8):
    sutils.checkConfig(cfg)

    if tasksToRun is None:
        tasksToRun = StructDict()

    # copy the geneMapFile and quant file to output dir
    os.system(f"cp {cfg.geneMapFile} {cfg.outDir}")
    os.system(f"cp {cfg.quantFile} {cfg.outDir}")
    # output a dset list file
    dsetFileName = os.path.join(cfg.outDir, os.path.basename(cfg.datasetsFile))
    os.system(f"ls -1 {cfg.pclDir} > {dsetFileName}")

    if tasksToRun.all is True or runAll is True:
        tasksToRun.pclbin = True
        tasksToRun.dab = True
        tasksToRun.prep = True
        tasksToRun.makeDB = True
        tasksToRun.plat = True
        tasksToRun.sinfo = True
        tasksToRun.dsetSize = True

    if tasksToRun.pclbin:
        # Convert pcl to pclbin
        sutils.Pcl2Pclbin(cfg, concurrency)
    if tasksToRun.dab:
        # Create the DAB files from pcl files
        sutils.PclToDabFiles(cfg, concurrency)
        # link the quant files for the DAB files
        sutils.linkQuantFilesForDabs(cfg)
    if tasksToRun.prep:
        # Calculate the gene zscore averages and stdev
        sutils.makeGenePrepFiles(cfg, concurrency)
    if tasksToRun.makeDB:
        # Create the DB files
        sutils.makeDBFiles(cfg, concurrency)
    if tasksToRun.plat:
        # Calculate the platform wide gene averages
        sutils.makePlatFiles(cfg, concurrency)
    if tasksToRun.sinfo:
        # Calculate the avg and stddev gene correleations per dataset (sinfo)
        sutils.sinfoCreate(cfg, concurrency)
    if tasksToRun.dsetSize:
        # Create the dataset size file
        sutils.makeDsetSizeFile(cfg)
    return True


#set up script for SEEK
if __name__=="__main__":
    cfg = sutils.getDefaultConfig()

    #input directory (containing the PCL's)
    #input file: dataset list (3-column format: file, name, platform)
    #input directory (containing the setting file: quant2, gene_map.txt)
    #path to sleipnir build
    argParser = argparse.ArgumentParser()
    argParser.add_argument('--all', default=False, action='store_true',
                           help='run all creation processes')
    argParser.add_argument('--pclbin', default=False, action='store_true',
                           help='run pcl2pclbin')
    argParser.add_argument('--dab', default=False, action='store_true',
                           help='run pcl2Dab')
    argParser.add_argument('--prep', default=False, action='store_true',
                           help='create prep directory files')
    argParser.add_argument('--makeDB', default=False, action='store_true',
                           help='run dab2DB')
    argParser.add_argument('--plat', default=False, action='store_true',
                           help='create platform directory files')
    argParser.add_argument('--sinfo', default=False, action='store_true',
                           help='create sinfo directory files')
    argParser.add_argument('--dsetSize', default=False, action='store_true',
                           help='create dataset size file')

    argParser.add_argument('--sleipnirBinDir', '-b', type=str, required=True,
                           help='Directory where the Sleipnir binaries are installed')
    argParser.add_argument('--inDir', '-i', type=str, required=True,
                           help='Directory of existing reference files (i.e. for reference gene_list and quant files)')
    argParser.add_argument('--outDir', '-o', type=str, required=True,
                           help='Output directory to write new database into')
    argParser.add_argument('--pclDir', '-p', type=str, required=False, default=cfg.pclDir,
                           help='Directory containing the PCL files for the new datasets')
    argParser.add_argument('--datasetFile', '-d', type=str, required=False, default=cfg.datasetsFile,
                           help='Text file listing the datasets, one dataset per line, three columns (pcl_file, name, platform)')
    argParser.add_argument('--numDBFiles', '-n', type=int, required=False, default=cfg.numDbFiles,
                           help='Number of output DB files to spread gene data across (should match refDB number)')
    argParser.add_argument('--concurrency', '-m', type=int, required=False, default=4,
                           help='Number of parallel processes to run for each task')
    args = argParser.parse_args()

    tasksToRun = StructDict()
    tasksToRun.all = args.all
    tasksToRun.pclbin = args.pclbin
    tasksToRun.dab = args.dab
    tasksToRun.prep = args.prep
    tasksToRun.makeDB = args.makeDB
    tasksToRun.plat = args.plat
    tasksToRun.sinfo = args.sinfo
    tasksToRun.dsetSize = args.dsetSize

    if not any(tasksToRun.values()):
        print("No task types specified: specify one or more, for example --dab or --all etc.")
        sys.exit(-1)

    cfg.binDir = args.sleipnirBinDir
    cfg.inDir = args.inDir
    cfg.outDir = args.outDir
    cfg.pclDir = args.pclDir
    cfg.datasetsFile = args.datasetFile
    cfg.numDbFiles = args.numDBFiles
    sutils.checkConfig(cfg)

    # check max open files setting is sufficient, i.e. ulimit -n
    softFileLimit, _ = resource.getrlimit(resource.RLIMIT_NOFILE)
    if softFileLimit < args.numDBFiles:
        print("Max open files limit is insufficient: {}".format(softFileLimit))
        sys.exit(-1)

    createSeekDB(cfg, tasksToRun, concurrency=args.concurrency)
