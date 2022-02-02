"""
This module has a set of functions for making Seek DB files and metadata.
The functions have been parallelized where possible.
"""
import os
import re
import sys
import math
import glob
import time
import toml
import shutil
import random
import fnmatch
import resource
import tempfile
import subprocess
import numpy as np
currPath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(currPath)
from runParallelJobs import runParallelJobs
from structDict import StructDict

# TODO - parallelize makeDsetSizeFile() and makePlatFiles()

# Each funtion operates with a cfg struct that defines file and binary locations
defaultConfig = StructDict({
    'binDir': None,
    'inDir': None,
    'outDir': None,
    'pclDir': 'pcl',
    'dabDir': 'dab',
    'dbDir': 'db',
    'sinfoDir': 'sinfo',
    'prepDir': 'prep',
    'platformDir': 'plat',
    'geneMapFile': 'gene_map.txt',
    'quantFile': 'quant2',
    'datasetsFile': 'datasets.txt',
    'datasetPlatMapFile': 'dsetPlatMap.txt',
    'dsetSizeFile': 'dset_size.txt',
    'numDbFiles': 1000,
    'configVerified': False,
})

def getDefaultConfig():
    return defaultConfig.copy()

def to_camel_case(snake_str):
    components = snake_str.split('_')
    # Capitalize the first letter of each component except the first one
    # Use the 'title' method and join them together.
    return components[0].lower() + ''.join(x.title() for x in components[1:])

def loadConfig(configFilename):
    tomlCfg = toml.load(configFilename)
    DB1 = tomlCfg['Database']['DB1']
    cfg = getDefaultConfig()
    for k, v in DB1.items():
        newKey = to_camel_case(k)
        cfg[newKey] = v
    # Translate for a couple of legacy 'DB' names to new config names
    cfg.numDbFiles = cfg.numberOfDb
    cfg.datasetPlatMapFile = cfg.dsetMapFile
    return cfg

def checkConfig(cfg):
    """
    Does checks on the cfg and standardizes the file locations to absolute paths
    """
    cfg.outDir = os.path.abspath(cfg.outDir)
    cfg.inDir = os.path.abspath(cfg.inDir)
    if not os.path.exists(cfg.outDir):
        os.makedirs(cfg.outDir)
    # set some paths relative to input directory if path missing
    for elem in ('pclDir', 'geneMapFile', 'quantFile', 'datasetsFile'):
        # check if the elem contains a path, i.e. '/'
        if not os.sep in cfg[elem]:
            # it doesn't contain a path, so set one
            cfg[elem] = os.path.join(cfg.inDir, cfg[elem])
        if not os.path.exists(cfg[elem]):
            raise FileNotFoundError(f'{elem}: {cfg[elem]}')
    # set some paths relative to output directory if path missing
    for elem in ('dabDir', 'dbDir', 'dsetSizeFile', 'sinfoDir', 'prepDir', 'platformDir', 'datasetPlatMapFile'):
        # check if the elem contains a path, i.e. '/'
        if not os.sep in cfg[elem]:
            # it doesn't contain a path, so set one
            cfg[elem] = os.path.join(cfg.outDir, cfg[elem])
        if not os.path.exists(cfg[elem]):
            if elem.endswith('Dir'):
                os.makedirs(cfg[elem])
    cfg.configVerified = True


def prepCmd(executableName, tagName, cfg):
    """
    Does some checks and command prep:
    - Verifies that checkConfig has been run
    - Builds out the full path to the command binary
    - Verifies the the command binary exists
    """
    if cfg.configVerified is False:
        raise ValueError(f'Please run checkConfig(cfg) prior to calling {tagName}')
    cmdName = os.path.join(cfg.binDir, executableName)
    if not os.path.isfile(cmdName):
        raise FileNotFoundError(f'Binary {cmdName} not found')
    return cmdName

def readGeneMapFile(geneFile):
    genes = []
    with open(geneFile) as fp:
        for line in fp:
            id, gene = line.split()
            gene = gene.rstrip("\n")
            genes.append(gene)
    return genes

def readSeekBinaryResultFile(dataFile):
    vals = []
    with open(dataFile, 'rb') as fp:
        # The first 8 byte (long int) is the number of elements stored
        headerVals = np.fromfile(fp, count=1, dtype=np.ulonglong)
        numVals = headerVals[0]
        # The remaining are 4 byte float values, numVal of them
        vals = np.fromfile(fp, dtype=np.float32)
        assert len(vals) == numVals
    return vals

def readDatasetList(dsetFile):
    """
    Read in datasets file and try to parse the various possible layouts, including
    one, two or three column variations
    """
    fp = open(dsetFile)
    count = 0
    dset_list = []
    for line in fp:
        count += 1
        # check if line has spaces in it
        if " " in line:
            raise ValueError(f"seekUtils: readDatasetList: line({count}): spaces dectected in line")
        cols = line.rstrip("\n").split("\t")
        # expecting first column like 'dataset.platform.pcl'
        parts = cols[0].split(".")
        if len(parts) < 3:
            msg = f"seekUtils: readDatasetList: line({count}): {cols[0]}: expecting at first column of form datasetName.platformName.pcl"
            raise ValueError(msg)
        if parts[-1] != 'pcl':
            print(f'Skipping non-pcl file {cols[0]}')
        if len(cols) == 1:
            # next column is 'dataset.platform'
            cols.append(f'{parts[0]}.{parts[1]}')
            # last column is 'platform'
            cols.append(parts[1])
        if len(cols) == 2:
            # for two columns, expecting 'datasets.platform.pcl platform'
            if cols[1] != parts[1]:
                 print(f"Error line({count}): For two columns, expecting 'dataset.platform.pcl platform'")
            cols.insert(1, f'{parts[0]}.{parts[1]}')
        if len(cols) > 3:
            msg = f"seekUtils: readDatasetList: line({count}): contains too many columns!"
            raise ValueError(msg)
        else:
            dset_list.append(tuple(cols))
    fp.close()
    return dset_list

def readPlatMap(platMapFile):
    """
    Read in platform mapping file and return the map from dataset name to platform name
    """
    dsetPlatMap = {}
    with open(platMapFile, 'r') as fp:
        for line in fp:
            cols = line.rstrip("\n").split("\t")
            parts = cols[0].split(".")
            dsetName = parts[0]
            if len(cols) > 1:
                platName = cols[-1]
            elif len(parts) > 1:
                platName = parts[1]
            else:
                msg = f"seekUtils: readPlatMap: found only one column with no platform information"
                raise ValueError(msg)
            dsetPlatMap[dsetName] = platName
    return dsetPlatMap

def writeDatasetList(dset_list, dsetFile):
    fw = open(dsetFile, "w")
    for (file_name, name, platform) in dset_list:
        m = re.match("(.*).pcl", file_name)
        datasetName = m.group(1)
        fw.write("%s\n" % datasetName)
    fw.close()


def writeDatasetPlatformMap(dset_list, dsetFile):
    fw = open(dsetFile, "w")
    for (file_name, name, platform) in dset_list:
        # m = re.match("(.*).pcl", file_name)
        # datasetName = m.group(1)
        fw.write("%s\t%s\n" % (name, platform))
    fw.close()

def splitFile(origFile, numParts):
    """
    Split a text file into N different parts
    Returns list of filenames of the partial files
    """
    outputDir = os.path.dirname(origFile)
    filename = os.path.basename(origFile)
    # Read in all the file lines
    lines = []
    with open(origFile, 'r') as fp:
        lines = fp.readlines()
    numPerFile = math.ceil(len(lines) / numParts)
    partialFiles = []
    for idx in range(numParts):
        pFile = os.path.join(outputDir, f'{idx:02d}_{filename}')
        partialFiles.append(pFile)
        with open(pFile, 'w') as fp:
            startPoint = numPerFile * idx
            endPoint = startPoint + numPerFile
            fp.writelines(lines[startPoint:endPoint])
    return partialFiles

def makeRandomQueryFile(cfg, numQueries, outFile):
    """
    Make random queries between 1 and 100 genes
    The number of queries made is always increment of 100
    """
    genes = readGeneMapFile(cfg.geneMapFile)
    # Will create queries of size between 1 and 100
    # A set of queries will have one of each size
    numSets = math.ceil(numQueries / 100)
    with open(outFile, 'w') as fp:
        for _ in range(numSets):
            for querySize in range(1, 101):
                random.shuffle(genes)
                fp.write(" ".join(genes[:querySize]) + "\n")
    return

def Pcl2Pclbin(cfg, concurrency=8):
    """
    Convert PCL files to PCL.bin files
    """
    cmdName = prepCmd('PCL2Bin', 'Pcl2Pclbin', cfg)
    datasets = readDatasetList(cfg.datasetsFile)
    pclBinDir = os.path.join(cfg.outDir, "pclbin")
    if not os.path.exists(pclBinDir):
        os.makedirs(pclBinDir)
    taskList = []
    for (pclFile, dsetPlat, platform) in datasets:
        pclFileFullPath = os.path.join(cfg.pclDir, pclFile)
        pclBinFullPath = os.path.join(pclBinDir, pclFile + '.bin')
        cmd = f"{cmdName} -i {pclFileFullPath} -o {pclBinFullPath}"
        taskList.append(cmd)
    startTime = time.time()
    runParallelJobs(taskList, concurrency=concurrency)
    print('Pcl2Pclbin: completion time {}s'.format(time.time() - startTime))


def PclToDabFiles(cfg, limit_genes=False, concurrency=8):
    """
    Convert PCL files to DAB files, which computes gene correlation values across all
    samples in a PCL file.
    """
    print("DAB File Creation: Calculating correlation matrices...\n")
    cmdName = prepCmd('Distancer', 'PclToDabFiles', cfg)
    if not os.path.exists(cfg.dabDir):
        os.makedirs(cfg.dabDir, exist_ok=True)
    datasets = readDatasetList(cfg.datasetsFile)
    # Create the list of genes to include in the DABs from the gene_map.txt
    # Distancer expects the input gene file to only contain one column of gene names
    # so we cut and keep the second column from the gene_map file
    if limit_genes is True:
        tmpFile = tempfile.NamedTemporaryFile(delete=False)
        tmpGeneFile = tmpFile.name
        os.system(f'cut -f 2 {cfg.geneMapFile} > {tmpGeneFile}')
    # create the dab files
    taskList = []
    for (pclFile, dsetPlat, platform) in datasets:
        m = re.match("(.*).pcl", pclFile)
        datasetName = m.group(1)
        pclFileFullPath = os.path.join(cfg.pclDir, pclFile)
        dabFileFullPath = os.path.join(cfg.dabDir, f'{datasetName}.dab')
        cmd = f"{cmdName} -i {pclFileFullPath} -o {dabFileFullPath} -s 0 -t 1"
        if limit_genes is True:
            cmd += f" -g {tmpGeneFile}"
        taskList.append(cmd)
    startTime = time.time()
    runParallelJobs(taskList, concurrency=concurrency)
    if limit_genes is True:
        os.remove(tmpGeneFile)
    print('PclToDabFiles: completion time {}s'.format(time.time() - startTime))


def linkQuantFilesForDabs(cfg):
    """
    Create a per-DAB file link to the quant file. This is needed later by Data2DB.
    """
    # link to the new quant file
    if cfg.configVerified is False:
        raise ValueError('Please run checkConfig(cfg) prior to calling linkQuantFiles')
    print("Linking quant files for each DAB file\n")
    # copy quant file to the output dir
    baseName = os.path.basename(cfg.quantFile)
    targetQuantFile = os.path.join(cfg.outDir, baseName)
    os.system(f'cp {cfg.quantFile} {targetQuantFile}')
    datasets = readDatasetList(cfg.datasetsFile)
    for (pclFile, dsetPlat, platform) in datasets:
        m = re.match("(.*).pcl", pclFile)
        datasetName = m.group(1)
        quantFile = os.path.join(cfg.dabDir, f'{datasetName}.quant')
        if os.path.exists(quantFile):
            os.system(f"rm -rf {quantFile}")
        cmd = f"ln -s {targetQuantFile} {quantFile}"
        # print(cmd)
        os.system(cmd)


def makeGenePrepFiles(cfg, concurrency=8):
    """
    Make the Prep directory files, which contains the gene average values for each dataset
    and the gene presence for each dataset.
    """
    cmdName = prepCmd('SeekPrep', 'makeGenePrepFiles', cfg)
    prepDir = os.path.join(cfg.outDir, "prep")
    if not os.path.exists(prepDir):
        os.makedirs(prepDir)
    datasets = readDatasetList(cfg.datasetsFile)
    taskList = []
    for (pclFile, dsetPlat, platform) in datasets:
        m = re.match("(.*).pcl", pclFile)
        datasetName = m.group(1)
        dabFileFullPath = os.path.join(cfg.dabDir, f'{datasetName}.dab')
        print("SeekPrep Processing %s...\n" % datasetName)
        cmd = f"{cmdName} -d -a -B {dabFileFullPath} -i {cfg.geneMapFile} -D {prepDir}"
        taskList.append(cmd)
        cmd = f"{cmdName} -d -p -B {dabFileFullPath} -i {cfg.geneMapFile} -D {prepDir}"
        taskList.append(cmd)
    startTime = time.time()
    runParallelJobs(taskList, concurrency=concurrency)
    print('makeGenePrepFiles: completion time {}s'.format(time.time() - startTime))


def sinfoCreate(cfg, concurrency=8):
    """
    Make sinfo directory contents, which contains the gene average and stddev values for
    each dataset.
    """
    cmdName = prepCmd('SeekPrep', 'sinfoCreate', cfg)
    pclBinDir = os.path.join(cfg.outDir, "pclbin")
    if not os.path.exists(pclBinDir):
        raise FileNotFoundError(f"Pclbin directory not found {pclBinDir}")
    sinfoDir = os.path.join(cfg.outDir, "sinfo")
    if not os.path.exists(sinfoDir):
        os.makedirs(sinfoDir)
    datasets = readDatasetList(cfg.datasetsFile)
    taskList = []
    for (pclFile, dsetPlat, platform) in datasets:
        inputPclbinFile = os.path.join(pclBinDir, f'{pclFile}.bin')
        cmd = f"{cmdName} -i {cfg.geneMapFile} -D {sinfoDir} -e -V {inputPclbinFile} -s"
        taskList.append(cmd)
    startTime = time.time()
    runParallelJobs(taskList, concurrency=concurrency)
    print('sinfoCreate: completion time {}s'.format(time.time() - startTime))


def gvarCreate(cfg, concurrency=8):
    """
    Make gvar directory contents, which contains the gene average and var values for
    each dataset.
    """
    cmdName = prepCmd('SeekPrep', 'gvarCreate', cfg)
    pclBinDir = os.path.join(cfg.outDir, "pclbin")
    if not os.path.exists(pclBinDir):
        raise FileNotFoundError(f"Pclbin directory not found {pclBinDir}")
    gvarDir = os.path.join(cfg.outDir, "gvar")
    if not os.path.exists(gvarDir):
        os.makedirs(gvarDir)
    datasets = readDatasetList(cfg.datasetsFile)
    taskList = []
    for (pclFile, dsetPlat, platform) in datasets:
        inputPclbinFile = os.path.join(pclBinDir, f'{pclFile}.bin')
        cmd = f"{cmdName} -i {cfg.geneMapFile} -D {gvarDir} -e -V {inputPclbinFile} -v"
        taskList.append(cmd)
    startTime = time.time()
    runParallelJobs(taskList, concurrency=concurrency)
    print('gvarCreate: completion time {}s'.format(time.time() - startTime))


def makePlatFiles(cfg, concurrency=8):
    """
    Calculate platform-wide gene average and stddev
    This is parallelized within SeekPrep using OpenMP (omp)
    Requires:
        - DB files already made
        - gene prep files already made
    """
    cmdName = prepCmd('SeekPrep', 'makePlatFiles', cfg)
    platDir = os.path.join(cfg.outDir, "plat")
    if not os.path.exists(platDir):
        os.makedirs(platDir)
    prepDir = os.path.join(cfg.outDir, "prep")
    datasets = readDatasetList(cfg.datasetsFile)
    # write out dsetPlatMap if needed
    if not os.path.exists(cfg.datasetPlatMapFile):
        writeDatasetPlatformMap(datasets, cfg.datasetPlatMapFile)
    # write out temporary files for the SeekPrep call
    # write out the db file list
    dbList = glob.glob(os.path.join(cfg.dbDir, '*.db'))
    tmpFile1 = tempfile.NamedTemporaryFile(delete=False)
    dbListFile = tmpFile1.name
    with open(dbListFile, 'w') as fp:
        fp.write(os.linesep.join(dbList))
    cmd = f"{cmdName} -i {cfg.geneMapFile} -D {platDir} -f -P -b {dbListFile} -I {prepDir} " \
          f"-A {cfg.datasetPlatMapFile} -Q {cfg.quantFile}"
    os.system(cmd)
    os.remove(dbListFile)


def makeDsetSizeFile(cfg, concurrency=8):
    """Create dataset size file, which contains the number of samples per dataset
    TODO - add parallelization
    """
    # TODO - this could be parallelized if multiprocess returns a value for each task
    #   have it return the dset size (in a result queue) and collect them all.
    cmdName = prepCmd('SeekPrep', 'makeDsetSizeFile', cfg)
    pclBinDir = os.path.join(cfg.outDir, 'pclbin')
    dsetSizeFile = os.path.join(cfg.outDir, 'dset_size.txt')
    # create and truncate the dsetSize file
    with open(dsetSizeFile, 'w'):
        pass  # opening for write automatically truncates if it exists
    datasets = readDatasetList(cfg.datasetsFile)
    for (pclFile, dsetPlat, platform) in datasets:
        if pclFile.endswith('.pcl'):
            pclFile = pclFile + '.bin'
        pclFileFullPath = os.path.join(pclBinDir, pclFile)
        cmd = f"{cmdName} -e -S -i /tmp -D {cfg.outDir} -V {pclFileFullPath} >> {dsetSizeFile}"
        os.system(cmd)


def dsetFileName(dbDirName):
    return dbDirName + "_dataset_list.txt"


def dsetPlatMapFileName(dbDirName, cfg):
    filename = os.path.basename(cfg.datasetsFile)
    return dbDirName + "_" + filename


def Task_DABToDBFiles(cfg, threadDbDir, threadDatasets: list):
    """
    A task that creates a DB directory with a set of .db files from a subset of the datasets.
    This is a parallel task which will be called from runParallelJobs, it will be invoked once
    per task in the task queue, thus there will be N out DB directories, one per thread.
    """
    cmdName = os.path.join(cfg.binDir, 'Data2DB')
    if not os.path.isfile(cmdName):
        raise FileNotFoundError(f'Binary {cmdName} not found')
    if not os.path.exists(threadDbDir):
        os.makedirs(threadDbDir, exist_ok=True)
    threadDatasetFile = dsetFileName(threadDbDir)
    writeDatasetList(threadDatasets, threadDatasetFile)
    writeDatasetPlatformMap(threadDatasets, dsetPlatMapFileName(threadDbDir, cfg))
    cmd = f"{cmdName} -i {cfg.geneMapFile} -d {cfg.dabDir} -D {threadDbDir} -u -B 50 -f {cfg.numDbFiles} -x {threadDatasetFile}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    return True


def parallelMakePerThreadDB(cfg, concurrency=8):
    """
    Divide datasets between the threads and build N sets of DB files
    Return:
        The list of db directories which were created (one per thread).
        These will need to be combined in a next step
    """
    cmdName = prepCmd('Data2DB', 'parallelMakePerThreadDB', cfg)
    datasets = readDatasetList(cfg.datasetsFile)
    dbDir = os.path.join(cfg.outDir, 'thread_work', 'db')
    dbDir = os.path.abspath(dbDir)
    # Number of datasets per thread
    numDsPerThread = math.ceil(len(datasets) / concurrency)
    taskList = []
    dbDirsToCombine = []
    for i in range(concurrency):
        threadDbDir = dbDir + f'{i:03d}'
        beg_idx = i * numDsPerThread
        end_idx = beg_idx + numDsPerThread
        threadDatasetPart = datasets[beg_idx:end_idx]
        if len(threadDatasetPart) == 0:
            break
        args = (cfg, threadDbDir, threadDatasetPart)
        cmd = {'func': Task_DABToDBFiles, 'args': args}
        taskList.append(cmd)
        dbDirsToCombine.append(threadDbDir)
    startTime = time.time()
    runParallelJobs(taskList, concurrency=concurrency, isPyFunction=True)
    print('parallelMakePerThreadDB: completion time {}s'.format(time.time() - startTime))
    # Also concatenate the dataset_map files together
    combinedListFile = os.path.join(cfg.outDir, 'combined_dataset_map.txt')
    if os.path.isfile(combinedListFile):
        os.remove(combinedListFile)
    for dbDir in dbDirsToCombine:
        partDsetFile =  dsetPlatMapFileName(dbDir, cfg)
        if not os.path.isfile(partDsetFile):
            raise FileNotFoundError(f"Error combining dataset lists, missing file {partDsetFile}")
        os.system(f"cat {partDsetFile} >> {combinedListFile}")
    return dbDirsToCombine


def Task_CombineDBFiles(cfg, dbDirsToCombine, dbFileName):
    """
    This task combines all instances of a db file (of the same name) across all DB directories.
    For example combine all 000000.db files from dbDirs [db000, db001, db002]. This is the
    second step of making DB files in parallel, the first creates many DB dirs with part
    of the datasets represented in each. This step combines them all together into one DB dir.
    """
    cmdName = os.path.join(cfg.binDir, 'DBCombiner')
    if not os.path.isfile(cmdName):
        raise FileNotFoundError(f'Binary {cmdName} not found')
    # Create a temp file to write the dbFile paths into
    #  This will be a list of, for example, all 00000.db files from each thread directory
    tmpFile = tempfile.NamedTemporaryFile(delete=False)
    with open(tmpFile.name, 'w') as fp:
        for dirName in dbDirsToCombine:
            dbDirFile = os.path.join(dirName, dbFileName)
            fp.write(dbDirFile + os.linesep)
    cmd = f"{cmdName} -C -i {cfg.geneMapFile} --db {tmpFile.name} -D {cfg.dbDir}"
    print(f'{cmd}')
    subprocess.run(cmd, shell=True)
    os.remove(tmpFile.name)
    # TODO - perhaps remove the threadDB files after combined
    return True


def parallelCombineThreadDBs(cfg, dbDirsToCombine, concurrency=8):
    """
    Merge the N sets of DB files created by the threads together into the
    final single DB files
    """
    # get the listing of db files
    cmdName = prepCmd('DBCombiner', 'parallelCombineThreadDBs', cfg)
    assert len(dbDirsToCombine) > 0
    # The list of db files should be identical for each db directory (same count and names)
    dbFileList = fnmatch.filter(os.listdir(dbDirsToCombine[0]), '*.db')
    cfg.numDbFiles = len(dbFileList)
    # We will create one task per db file name (i.e. combine all 00000001.db files together)
    taskList = []
    for dbFileName in dbFileList:
        args = (cfg, dbDirsToCombine, dbFileName)
        task = {'func': Task_CombineDBFiles, 'args': args}
        taskList.append(task)
    startTime = time.time()
    runParallelJobs(taskList, concurrency=concurrency, isPyFunction=True)
    print('parallelCombineThreadDBs: completion time {}s'.format(time.time() - startTime))


def verifyCombinedDBs(cfg, dbDirsToCombine, concurrency=8):
    """
    Verifies that the per-thread sets of DB files were correctly merged into the
    final combined set fo DB files.
    """
    cmdName = os.path.join(cfg.binDir, 'verifyMergedDBFiles')
    assert len(dbDirsToCombine) > 0
    dbFileList = fnmatch.filter(os.listdir(dbDirsToCombine[0]), '*.db')
    # write out dbDirsToCombine
    tmpFile = tempfile.NamedTemporaryFile(delete=False)
    dbDirsFile = tmpFile.name
    with open(dbDirsFile, 'w') as fp:
        fp.write(os.linesep.join(dbDirsToCombine))
    taskList = []
    for dbFileName in dbFileList:
        cmd = f'{cmdName} -i {dbDirsFile} -c {cfg.dbDir} -f {dbFileName}'
        taskList.append(cmd)
    startTime = time.time()
    runParallelJobs(taskList, concurrency=concurrency)
    print('verifyCombinedDBs: completion time {}s'.format(time.time() - startTime))
    os.remove(tmpFile.name)

def makeDBFiles(cfg, concurrency=8):
    """
    Function to make the database .db files in parallel, let N = concurrency
        - First makes N sets of DB files (threads each create one)
        - Second, combines (merges) the N set into the final single set of DB files
        - Finally, verifies that the merged result matches the N inputs
    """
    checkConfig(cfg)
    dbDirsToCombine = parallelMakePerThreadDB(cfg, concurrency)
    parallelCombineThreadDBs(cfg, dbDirsToCombine, concurrency)
    verifyCombinedDBs(cfg, dbDirsToCombine, concurrency)
    # Remove the temp dbDirs created by the threads
    for tmpDir in dbDirsToCombine:
        shutil.rmtree(tmpDir)

def runSeekMiner(cfg, queryFile, outputDir, concurrency=8):
    """Run SeekMiner on a set of queries specified in a file"""
    # Divide the queries into N files for the concurrent tasks
    # Alternate: Unix command 'split --number=l/5 inputfile outputprefix'
    tmpFiles = splitFile(queryFile, concurrency)
    # Create the list of commands for to run concurrently
    seekMinerBin = os.path.join(cfg.binDir, 'SeekMiner')
    seekMinerCmd = \
        f'time {seekMinerBin} ' \
        f'-x {cfg.datasetPlatMapFile} -i {cfg.geneMapFile} ' \
        f'-d {cfg.dbDir} -p {cfg.prepDir} -n {cfg.numDbFiles} ' \
        f'-P {cfg.platformDir} -Q {cfg.quantFile} ' \
        f'-u {cfg.sinfoDir} -R {cfg.dsetSizeFile} ' \
        f'-b 200 -V CV -I LOI -z z_score -m -M -O '
    # Add to a list of tasks
    resDirs = []
    taskList = []
    for idx, qFile in enumerate(tmpFiles):
        resDir = os.path.join(outputDir, f'{idx:02d}')
        os.makedirs(resDir, exist_ok=True)
        resDirs.append(resDir)
        task = seekMinerCmd + f'-q {qFile} -o {resDir}'
        taskList.append(task)
    # Run the tasks
    startTime = time.time()
    runParallelJobs(taskList, concurrency=concurrency)
    print('runSeekMiner: completion time {}s'.format(time.time() - startTime))
    # Move all the result files to the top-level directory
    renumberMoveFiles(resDirs, outputDir)
    # remove all the tmp result directories
    for resDir in resDirs:
        shutil.rmtree(resDir)
    # remove the tmp queryFiles
    for tfile in tmpFiles:
        os.remove(tfile)
    return

def renumberMoveFiles(dirs, outputDir):
    """
    Given a set of directories, each containing files with numeric
    names (e.g. 1.gscore, 10.gscore, etc.), move the files to
    the output directory and renumber such that the first dir
    will be 0 to N, the second dir files will become N+1 to M, etc.
    """
    totalIdx = -1
    for fromDir in dirs:
        # list all files starting with a number (in numeric order)
        pattern = os.path.join(fromDir, "[0-9]*.*")
        files = glob.glob(pattern)
        # Sort the files by number
        files.sort(key=lambda k: int(os.path.basename(k).split('.')[0]))
        # There may be multiple files with the same numbered name, like
        #   0.results.txt, 0.gscores, 0.dweights etc.
        # So track the prevFileNum and only change when we encounter a new num
        prevFileNum = None
        # rename and move them to output directory
        for file in files:
            filename = os.path.basename(file)
            filenum, rest = filename.split('.', 1)
            assert filenum is not None
            if int(filenum) != prevFileNum:
                prevFileNum = int(filenum)
                totalIdx += 1
            newName = os.path.join(outputDir, f'{totalIdx}.{rest}')
            os.rename(file, newName)
