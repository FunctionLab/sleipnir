import pytest
import os
import subprocess
import glob
import fnmatch
import filecmp
import tempfile
import seekUtils as sutils

testDir = os.path.dirname(__file__)
seekScriptsDir = os.path.dirname(testDir)
sleipnirDir = os.path.dirname(os.path.dirname(seekScriptsDir))
sleipnirBin = os.path.join(sleipnirDir, 'Debug')
use_tempfile = False

class TestIncrDB:
    temp_dir = None

    def setup_class(cls):
        # Make a tmp directory and copy input mockDB files to it
        # TODO uncomment
        if use_tempfile:
            cls.temp_dir = tempfile.TemporaryDirectory()
            tmpDirName = cls.temp_dir.name
        else:
            username = os.environ.get('USER')
            tmpDirName = os.path.join('/tmp', username, 'testSeekIncrMerge')
            if os.path.exists(tmpDirName):
                os.system(f'rm -rf {tmpDirName}/*')
            else:
                os.makedirs(tmpDirName)
        assert os.path.exists(tmpDirName)
        cls.incrDbDir = os.path.join(tmpDirName, 'incrDb')
        cls.smallDbDir = os.path.join(tmpDirName, 'smallDb')
        cls.mockDbDir = os.path.join(tmpDirName, 'mockDb')
        cls.mergeDbDir = os.path.join(tmpDirName, 'mergeDb')
        cls.verifyDbDir = os.path.join(tmpDirName, 'verifyDb')
        os.makedirs(cls.incrDbDir)
        os.makedirs(cls.smallDbDir)
        os.makedirs(cls.mockDbDir)
        os.makedirs(cls.mergeDbDir)
        os.makedirs(cls.verifyDbDir)
        print(f'### Using temp directory {tmpDirName}')
        # Setup the temp test directory
        testInputsDir = os.path.join(testDir, 'inputs', 'incr_merge_db')
        # copy over incrDB files
        os.system(f'cp {testInputsDir}/gene_map.txt {cls.incrDbDir}')
        os.system(f'cp {testInputsDir}/quant2 {cls.incrDbDir}')
        os.makedirs(f'{cls.incrDbDir}/pcl')
        os.system(f'cp -a {testInputsDir}/incr_pcl/* {cls.incrDbDir}/pcl/')
        os.system(f'ls -1 {cls.incrDbDir}/pcl > {cls.incrDbDir}/dset_list.txt')
        # copy over smallDB files
        os.system(f'cp {testInputsDir}/gene_map.txt {cls.smallDbDir}')
        os.system(f'cp {testInputsDir}/quant2 {cls.smallDbDir}')
        os.makedirs(f'{cls.smallDbDir}/pcl')
        os.system(f'cp -a {testInputsDir}/smalldb_pcl/* {cls.smallDbDir}/pcl/')
        os.system(f'ls -1 {cls.smallDbDir}/pcl > {cls.smallDbDir}/dset_list.txt')
        # create a verifyDb from scratch from small and incr db
        os.system(f'cp {testInputsDir}/gene_map.txt {cls.verifyDbDir}')
        os.system(f'cp {testInputsDir}/quant2 {cls.verifyDbDir}')
        os.makedirs(f'{cls.verifyDbDir}/pcl')
        os.system(f'cp -a {testInputsDir}/smalldb_pcl/* {cls.verifyDbDir}/pcl/')
        os.system(f'cp -a {testInputsDir}/incr_pcl/* {cls.verifyDbDir}/pcl/')
        os.system(f'ls -1 {cls.verifyDbDir}/pcl > {cls.verifyDbDir}/dset_list.txt')
        # copy over mockDB files
        mockInputs = os.path.join(testDir, 'inputs', 'mockDb')
        os.system(f'cp -a {mockInputs}/* {cls.mockDbDir}/')

    def teardown_class(cls):
        if use_tempfile:
            if cls.temp_dir is not None:
                cls.temp_dir.cleanup()

    def test_makeMockDB(self):
        # setup the config file
        mockDir = TestIncrDB.mockDbDir
        # cfg = sutils.getDefaultConfig()
        # cfg.inDir = f'{mockDir}'
        # cfg.outDir = f'{mockDir}'
        # cfg.binDir = os.path.join(sleipnirDir, 'Debug')
        # cfg.datasetsFile = 'dset_list.txt'
        # sutils.checkConfig(cfg)
        cmd = f'python {seekScriptsDir}/seekCreateDB.py --all -d pcl_list.txt ' \
              f'-i {mockDir} -o {mockDir} -b {sleipnirBin}'
        ret = subprocess.run(cmd, shell=True)
        assert ret.returncode == 0

    def test_makeSmallDB(self):
        # setup the config file
        smallDir = TestIncrDB.smallDbDir
        # cfg = sutils.getDefaultConfig()
        # cfg.inDir = f'{smallDir}'
        # cfg.outDir = f'{smallDir}'
        # cfg.binDir = os.path.join(sleipnirDir, 'Debug')
        # cfg.datasetsFile = 'dset_list.txt'
        # sutils.checkConfig(cfg)
        cmd = f'python {seekScriptsDir}/seekCreateDB.py --all -d dset_list.txt ' \
              f'-i {smallDir} -o {smallDir} -b {sleipnirBin}'
        ret = subprocess.run(cmd, shell=True)
        assert ret.returncode == 0

    def test_incrMerge(self):
        largeDB = TestIncrDB.mockDbDir
        smallDB = TestIncrDB.smallDbDir
        mergeDB = TestIncrDB.mergeDbDir
        incrDB = TestIncrDB.incrDbDir
        verifyDB = TestIncrDB.verifyDbDir
        largeDset = os.path.join(largeDB, 'pcl_list.txt')
        smallDset = os.path.join(smallDB, 'dset_list.txt')
        incrDset = os.path.join(incrDB, 'dset_list.txt')
        incrPcl = os.path.join(incrDB, 'pcl')
        cmd = f'python seekIncrementalMerge.py -l {largeDB} -s {smallDB} ' \
              f'-o {mergeDB} -dl {largeDset} -ds {smallDset} -dn {incrDset} ' \
              f'-p {incrPcl} -b {sleipnirBin} -y'
        ret = subprocess.run(cmd, shell=True)
        assert ret.returncode == 0

        # Combine all the PCL files and build the full database from scratch
        # Create the verifyDB from scratch from the combined small, incr dbs
        cmd = f'python {seekScriptsDir}/seekCreateDB.py --all -d dset_list.txt ' \
              f'-i {verifyDB} -o {verifyDB} -b {sleipnirBin}'
        print(cmd)
        ret = subprocess.run(cmd, shell=True)
        assert ret.returncode == 0

        # Compare the full build DB files to the incremental build DB files
        # dbFiles = glob.glob1(os.path.join(verifyDB , 'db'),  '*.db')
        dirs = ['db', 'plat', 'prep', 'sinfo']
        for d in dirs:
            files = glob.glob1(os.path.join(verifyDB, d),  '*')
            files.sort()
            for f in files:
                vfile = os.path.join(verifyDB, d, f)
                ifile = os.path.join(mergeDB, d, f)
                assert filecmp.cmp(ifile, vfile, shallow=False)

