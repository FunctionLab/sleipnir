"""
This module has a set of functions for make Seek DB files and metadata.
The functions have been parallelized where possible.
TODO - parallelize makeDsetSizeFile() and makePlatFiles()
"""
import os
import re
import sys
import math
import glob
import time
import fnmatch
import resource
import tempfile
import subprocess
currPath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(currPath)
from runParallelJobs import runParallelJobs
from structDict import StructDict

# Each funtion operates with a cfg struct that defines file and binary locations
defaultConfig = StructDict({    
    'binDir': None,
    'inDir': None,
    'outDir': None,
    'pclDir': 'pcl',
    'dabDir': 'dab',
    'dbDir': 'db',
    'geneMapFile': 'gene_map.txt',
    'quantFile': 'quant2',
    'datasetsFile': 'datasets.txt',
    'numDbFiles': 1000,
    'configVerified': False,
})

def checkConfig(cfg):
    """
    Does checks on the cf and standardizes the file locations to absolute paths
    """
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
    for elem in ('dabDir', 'dbDir'):
        # check if the elem contains a path, i.e. '/'
        if not os.sep in cfg[elem]:
            # it doesn't contain a path, so set one
            cfg[elem] = os.path.join(cfg.outDir, cfg[elem])
        if not os.path.exists(cfg[elem]):
            os.makedirs(cfg[elem])
    cfg.configVerified = True
   

def read_dataset_list(dsetFile):
    fp = open(dsetFile)
    count = 0
    dset_list = []
    for line in fp:
        count += 1
        cols = line.rstrip("\n").split("\t")
        # expecting first column like 'dataset.platform.pcl'
        parts = cols[0].split(".")
        if len(parts) < 2:
            print(f"Error line({count}): expecting at least two columns per line")
            break
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
            print(f"Error line({count}): contains too many columns!\n");
            break
        else:
            dset_list.append(tuple(cols))
    fp.close()
    return dset_list


def write_dataset_list(dset_list, dsetFile):
    fw = open(dsetFile, "w")
    for (file_name, name, platform) in dset_list:
        m = re.match("(.*).pcl", file_name)
        datasetName = m.group(1)
        fw.write("%s\n" % datasetName)
    fw.close()


def write_dataset_platform_map(dset_list, dsetFile):
    fw = open(dsetFile, "w")
    for (file_name, name, platform) in dset_list:
        # m = re.match("(.*).pcl", file_name)
        # datasetName = m.group(1)
        fw.write("%s\t%s\n" % (name, platform))
    fw.close()


def prepCmd(funcName, executableName, cfg):
    if cfg.configVerified is False:
        raise ValueError(f'Please run checkConfig(cfg) prior to calling {funcName}')
    cmdName = os.path.join(cfg.binDir, executableName)
    if not os.path.isfile(cmdName):
        raise FileNotFoundError(f'Binary {cmdName} not found')
    return cmdName


def Pcl2Pclbin(cfg, concurrency=8):
    cmdName = prepCmd('Pcl2Pclbin', 'PCL2Bin', cfg)
    datasets = read_dataset_list(cfg.datasetsFile)
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


def PclToDabFiles(cfg, concurrency=8):
    print("DAB File Creation: Calculating correlation matrices...\n")
    cmdName = prepCmd('PclToDabFiles', 'Distancer', cfg)
    if not os.path.exists(cfg.dabDir):
        os.makedirs(cfg.dabDir, exist_ok=True)
    datasets = read_dataset_list(cfg.datasetsFile)
    # create the dab files
    taskList = []
    for (pclFile, dsetPlat, platform) in datasets:
        m = re.match("(.*).pcl", pclFile)
        datasetName = m.group(1)
        pclFileFullPath = os.path.join(cfg.pclDir, pclFile)
        dabFileFullPath = os.path.join(cfg.dabDir, f'{datasetName}.dab')
        cmd = f"{cmdName} -i {pclFileFullPath} -o {dabFileFullPath} -s 0 -t 1"
        taskList.append(cmd)
    startTime = time.time()
    runParallelJobs(taskList, concurrency=concurrency)
    print('PclToDabFiles: completion time {}s'.format(time.time() - startTime))


def linkQuantFilesForDabs(cfg):
    # link to the new quant file
    if cfg.configVerified is False:
        raise ValueError('Please run checkConfig(cfg) prior to calling linkQuantFiles')
    print("Linking quant files for each DAB file\n")
    # copy quant file to the output dir
    baseName = os.path.basename(cfg.quantFile)
    targetQuantFile = os.path.join(cfg.outDir, baseName)
    os.system(f'cp {cfg.quantFile} {targetQuantFile}')
    datasets = read_dataset_list(cfg.datasetsFile)
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
    cmdName = prepCmd('makeGenePrepFiles', 'SeekPrep', cfg)
    prepDir = os.path.join(cfg.outDir, "prep")
    if not os.path.exists(prepDir):
        os.makedirs(prepDir)
    datasets = read_dataset_list(cfg.datasetsFile)
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
    cmdName = prepCmd('sinfoCreate', 'SeekPrep', cfg)
    pclBinDir = os.path.join(cfg.outDir, "pclbin")
    if not os.path.exists(pclBinDir):
        raise FileNotFoundError(f"Pclbin directory not found {pclBinDir}")
    sinfoDir = os.path.join(cfg.outDir, "sinfo")
    if not os.path.exists(sinfoDir):
        os.makedirs(sinfoDir)
    datasets = read_dataset_list(cfg.datasetsFile)
    taskList = []
    for (pclFile, dsetPlat, platform) in datasets:
        inputPclbinFile = os.path.join(pclBinDir, f'{pclFile}.bin')
        cmd = f"{cmdName} -i {cfg.geneMapFile} -D {sinfoDir} -e -V {inputPclbinFile} -s"
        taskList.append(cmd)
    startTime = time.time()
    runParallelJobs(taskList, concurrency=concurrency)
    print('sinfoCreate: completion time {}s'.format(time.time() - startTime))


def makePlatFiles(cfg, concurrency=8):
    '''
    Calculate platform-wide gene score averages - not parallelized
    TODO - this could be parallelized if SeekPrep combineplat is updated to take a list of input dirs
    Requires: 
        - DB files already made
        - gene prep files already made
    '''
    cmdName = prepCmd('makePlatFiles', 'SeekPrep', cfg)
    platDir = os.path.join(cfg.outDir, "plat")
    if not os.path.exists(platDir):
        os.makedirs(platDir)
    prepDir = os.path.join(cfg.outDir, "prep")
    # write out temporary files for the SeekPrep call
    datasets = read_dataset_list(cfg.datasetsFile)
    dsetPlatMapFile = os.path.join(cfg.outDir, 'dset_plat_map.txt')
    write_dataset_platform_map(datasets, dsetPlatMapFile)
    # write out the db file list
    dbList = glob.glob(os.path.join(cfg.dbDir, '*.db'))
    tmpFile1 = tempfile.NamedTemporaryFile(delete=False)
    dbListFile = tmpFile1.name
    with open(dbListFile, 'w') as fp:
        fp.write(os.linesep.join(dbList))
    cmd = f"{cmdName} -i {cfg.geneMapFile} -D {platDir} -f -P -b {dbListFile} -I {prepDir} -A {dsetPlatMapFile} -Q {cfg.quantFile}"
    os.system(cmd)
    os.remove(dbListFile)


def makeDsetSizeFile(cfg, concurrency=8):
    '''Create dataset size file - not parallelized'''
    # TODO - this could be parallelized if multiprocess returns a value for each task
    #   have it return the dset size (in a result queue) and collect them all.
    cmdName = prepCmd('makeDsetSizeFile', 'SeekPrep', cfg)
    pclBinDir = os.path.join(cfg.outDir, 'pclbin')
    dsetSizeFile = os.path.join(cfg.outDir, 'dset_size.txt')
    # create and truncate the dsetSize file
    with open(dsetSizeFile, 'w'):
        pass  # opening for write automatically truncates if it exists
    datasets = read_dataset_list(cfg.datasetsFile)
    for (pclFile, dsetPlat, platform) in datasets:
        if pclFile.endswith('.pcl'):
            pclFile = pclFile + '.bin'
        pclFileFullPath = os.path.join(pclBinDir, pclFile)
        cmd = f"{cmdName} -e -S -i /tmp -D {cfg.outDir} -V {pclFileFullPath} >> {dsetSizeFile}"
        os.system(cmd)


def dsetFileName(dbDirName):
    return os.path.join(dbDirName + "_dataset_list.txt")


def dsetPlatMapFileName(dbDirName):
    return os.path.join(dbDirName + "_dataset_map.txt")


def Task_DABToDBFiles(cfg, threadDbDir, threadDatasets: list):
    '''
    A task that creates a DB directory with a set of .db files from a subset of the datasets.
    This is a parallel task which will be called from runParallelJobs, it will be invoked once 
    per task in the task queue.
    '''
    cmdName = os.path.join(cfg.binDir, 'Data2DB')
    if not os.path.isfile(cmdName):
        raise FileNotFoundError(f'Binary {cmdName} not found')
    if not os.path.exists(threadDbDir):
        os.makedirs(threadDbDir, exist_ok=True)
    threadDatasetFile = dsetFileName(threadDbDir)
    write_dataset_list(threadDatasets, threadDatasetFile)
    write_dataset_platform_map(threadDatasets, dsetPlatMapFileName(threadDbDir))
    cmd = f"{cmdName} -i {cfg.geneMapFile} -d {cfg.dabDir} -D {threadDbDir} -u -B 50 -f {cfg.numDbFiles} -x {threadDatasetFile}"
    print(cmd)
    subprocess.run(cmd, shell=True)


def parallelMakePerThreadDB(cfg, concurrency=8):
    ''' Divide datasets between the threads and build sets of DB files
        Return:
            The list of db directories which were created (one per thread).
            These will need to be combined in a next step
    '''
    cmdName = prepCmd('parallelMakePerThreadDB', 'Data2DB', cfg)
    datasets = read_dataset_list(cfg.datasetsFile)
    dbDir = os.path.join(cfg.outDir, 'thread_work', 'db')
    dbDir = os.path.abspath(dbDir)
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
    return dbDirsToCombine


def Task_CombineDBFiles(cfg, dbDirsToCombine, dbFileName):
    '''
    This task combines all instances of a db file (of the same name) across all DB directories.
    For example combine all 000000.db files from dbDirs [db000, db001, db002]. This is the
    second step of making DB files in parallel, the first creates many DB dirs with part 
    of the datasets represented in each. This step combines them all together into one DB dir.
    '''
    cmdName = os.path.join(cfg.binDir, 'DBCombiner')
    if not os.path.isfile(cmdName):
        raise FileNotFoundError(f'Binary {cmdName} not found')
    # Create a temp file to write the dbFile paths into
    tmpFile = tempfile.NamedTemporaryFile(delete=False)
    fp = open(tmpFile.name, 'w')
    for dirName in dbDirsToCombine:
        dbDirFile = os.path.join(dirName, dbFileName)
        fp.write(dbDirFile + os.linesep)
    fp.close()
    cmd = f"{cmdName} -C -i {cfg.geneMapFile} --db {tmpFile.name} -D {cfg.dbDir}"
    print(f'{cmd}')
    subprocess.run(cmd, shell=True)
    os.remove(tmpFile.name)
    # TODO - perhaps remove the threadDB files after combined


def parallelCombineThreadDBs(cfg, dbDirsToCombine, concurrency=8):
    # merge the sets of DB files together
    # get the listing of db files
    cmdName = prepCmd('parallelCombineThreadDBs', 'DBCombiner', cfg)
    assert len(dbDirsToCombine) > 0
    # The list of db files should be identical for each db directory (same count and names)
    dbFileList = fnmatch.filter(os.listdir(dbDirsToCombine[0]), '*.db')
    # We will create one task per db file name (i.e. combine all 00000001.db files together)
    taskList = []
    for dbFileName in dbFileList:
        args = (cfg, dbDirsToCombine, dbFileName)
        task = {'func': Task_CombineDBFiles, 'args': args}
        taskList.append(task)
    startTime = time.time()
    runParallelJobs(taskList, concurrency=concurrency, isPyFunction=True)
    # Also concatenate the dataset_map files together
    combinedListFile = os.path.join(cfg.outDir, 'combined_dataset_map.txt')
    if os.path.isfile(combinedListFile):
        os.remove(combinedListFile)
    for dbDir in dbDirsToCombine:
        partDsetFile =  dsetPlatMapFileName(dbDir)
        if not os.path.isfile(partDsetFile):
            raise FileNotFoundError("Error combining dataset lists, missing file {partDsetFile}")
        os.system(f"cat {partDsetFile} >> {combinedListFile}")
    print('parallelCombineThreadDBs: completion time {}s'.format(time.time() - startTime))


def verifyCombinedDBs(cfg, dbDirsToCombine, concurrency=8):
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
    checkConfig(cfg)
    dbDirsToCombine = parallelMakePerThreadDB(cfg, concurrency)
    parallelCombineThreadDBs(cfg, dbDirsToCombine, concurrency)
    verifyCombinedDBs(cfg, dbDirsToCombine, concurrency)


