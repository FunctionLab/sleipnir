import pytest
import os
import filecmp
import subprocess
import tempfile
import seekUtils as sutils

testDir = os.path.dirname(__file__)
seekScriptsDir = os.path.dirname(testDir)
sleipnirDir = os.path.dirname(os.path.dirname(seekScriptsDir))
sleipnirBinDir = os.path.join(sleipnirDir, 'Debug')
use_tempfile = False

class TestSeekCreateDB:
    temp_dir = None
    mockDir = ''
    cfg = None

    def setup_class(cls):
        # Make a tmp directory and copy input mockDB files to it
        if use_tempfile:
            cls.temp_dir = tempfile.TemporaryDirectory()
            tmpDirName = cls.temp_dir.name
        else:
            tmpDirName = '/tmp/testSeekCreateDb'
            if os.path.exists(tmpDirName):
                os.system(f'rm -rf {tmpDirName}/*')
            else:
                os.makedirs(tmpDirName)
        assert os.path.exists(tmpDirName)
        cls.mockDir = os.path.join(tmpDirName, 'mockDb')
        print(f'### Using temp directory {tmpDirName}')
        # copy the test mock directory to the temp directory
        mockDbDir = os.path.join(testDir, 'inputs/mockDb')
        os.system(f'cp -a {mockDbDir} {tmpDirName}')


    def teardown_class(cls):
        if use_tempfile:
            if cls.temp_dir is not None:
                cls.temp_dir.cleanup()


    def test_createDB(self):
        mockDbDir = TestSeekCreateDB.mockDir
        # Create the DB
        cmd = f'python {seekScriptsDir}/seekCreateDB.py --all ' \
              f'-b {sleipnirBinDir} -i {mockDbDir} -o {mockDbDir} ' \
              f'-d pcl_list.txt -m 2'
        subprocess.run(cmd, shell=True)

        cfg = sutils.getDefaultConfig()
        cfg.inDir = mockDbDir
        cfg.outDir = mockDbDir
        cfg.binDir = os.path.join(sleipnirDir, 'Debug')
        cfg.datasetsFile = 'dset_plat_map.txt'
        sutils.checkConfig(cfg)
        # Run SeekMiner query
        os.makedirs(os.path.join(mockDbDir, 'results'))
        cmd = f'{sleipnirBinDir}/SeekMiner -x {cfg.datasetsFile} -i ' \
              f'{cfg.geneMapFile} -d {cfg.dbDir} -p {mockDbDir}/prep ' \
              f'-P {mockDbDir}/plat -Q {cfg.quantFile} -u {mockDbDir}/sinfo ' \
              f'-n 6 -b 20  -V CV -I LOI -z z_score -m -M -O ' \
              f'-q {mockDbDir}/query.txt -o {mockDbDir}/results -Y -T 2 -t 1'
        print(cmd)
        ret = subprocess.run(cmd, shell=True)
        assert ret.returncode == 0
        expected_results = os.path.join(mockDbDir, 'query_result.txt')
        seekminer_results = os.path.join(mockDbDir, 'results', '0.results.txt')
        assert filecmp.cmp(expected_results, seekminer_results, shallow=False)

    def test_seekRPC(self):
        # TODO - run seekRPCServer and seekRPCClient with the same query agains mockDB
        pass
