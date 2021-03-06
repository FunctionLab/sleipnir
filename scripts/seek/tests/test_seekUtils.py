import pytest
import os
import sys
import glob
import fnmatch
import filecmp
import tempfile
testDir = os.path.dirname(__file__)
seekScriptsDir = os.path.dirname(testDir)
sys.path.append(seekScriptsDir)
import seekUtils as sutils

sleipnirDir = os.path.dirname(os.path.dirname(seekScriptsDir))
use_tempfile = False

class TestSeekUtils:
    temp_dir = None
    mockDir = ''
    cfg = None

    def setup_class(cls):
        # Make a tmp directory and copy input mockDB files to it
        if use_tempfile:
            cls.temp_dir = tempfile.TemporaryDirectory()
            tmpDirName = cls.temp_dir.name
        else:
            username = os.environ.get('USER')
            tmpDirName = os.path.join('/tmp', username, 'testSeekUtils')
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
        # setup the config file
        cfg = sutils.getDefaultConfig()
        cfg.inDir = f'{cls.mockDir}'
        cfg.outDir = os.path.join(cls.mockDir, 'outdir')
        cfg.binDir = os.path.join(sleipnirDir, 'Debug')
        cfg.datasetsFile = 'pcl_list.txt'
        cls.cfg = cfg

    def teardown_class(cls):
        if use_tempfile:
            if cls.temp_dir is not None:
                cls.temp_dir.cleanup()

    def test_config(self):
        sutils.checkConfig(TestSeekUtils.cfg)

    def test_readDatasetPclList(self):
        # make several versions of the pcl list
        cfg = TestSeekUtils.cfg
        sutils.checkConfig(cfg)
        datasetPclFile = os.path.join(cfg.inDir, cfg.datasetsFile)
        dsetOneColFile = os.path.join(cfg.inDir, 'dset_one_col')
        dsetTwoColFile = os.path.join(cfg.inDir, 'dset_two_col')
        os.system(f"cut -f 1 {datasetPclFile} > {dsetOneColFile}")
        os.system(f"cut -f 1,3 {datasetPclFile} > {dsetTwoColFile}")
        d1 = sutils.readDatasetList(dsetOneColFile)
        d2 = sutils.readDatasetList(dsetTwoColFile)
        d3 = sutils.readDatasetList(datasetPclFile)
        assert d1 == d2
        assert d2 == d3

    @pytest.mark.dependency(name='pclbin')
    def test_pclbin(self):
        cfg = TestSeekUtils.cfg
        sutils.checkConfig(cfg)
        sutils.Pcl2Pclbin(cfg, concurrency=4)
        pclBinDir = os.path.join(cfg.outDir, 'pclbin')
        pclFileList = fnmatch.filter(os.listdir(cfg.pclDir), '*.pcl')
        pclbinFileList = fnmatch.filter(os.listdir(pclBinDir), '*.pcl.bin')
        assert len(pclFileList) == len(pclbinFileList)

    @pytest.mark.dependency(depends=['pclbin'])
    def test_sinfo(self):
        cfg = TestSeekUtils.cfg
        sutils.checkConfig(cfg)
        # next test creating the sinfo files from pclbin files
        sutils.sinfoCreate(cfg, concurrency=6)
        sinfoDir = os.path.join(cfg.outDir, 'sinfo')
        pclBinDir = os.path.join(cfg.outDir, 'pclbin')
        sinfoFileList = fnmatch.filter(os.listdir(sinfoDir), '*.sinfo')
        pclbinFileList = fnmatch.filter(os.listdir(pclBinDir), '*.pcl.bin')
        assert len(sinfoFileList) == len(pclbinFileList)

    @pytest.mark.dependency(name='pclToDab')
    def test_PclToDabFiles(self):
        cfg = TestSeekUtils.cfg
        sutils.checkConfig(cfg)
        sutils.PclToDabFiles(cfg)
        dabFileList = fnmatch.filter(os.listdir(cfg.dabDir), '*.dab')
        pclFileList = fnmatch.filter(os.listdir(cfg.pclDir), '*.pcl')
        assert len(dabFileList) == len(pclFileList)

    @pytest.mark.dependency(depends=['pclToDab'])
    def test_linkQuantFiles(self):
        cfg = TestSeekUtils.cfg
        sutils.checkConfig(cfg)
        sutils.linkQuantFilesForDabs(cfg)
        quantFileList = glob.glob(os.path.join(cfg.dabDir, '*.quant'))
        dabFileList = fnmatch.filter(os.listdir(cfg.dabDir), '*.dab')
        assert len(dabFileList) == len(quantFileList)
        assert filecmp.cmp(quantFileList[0], quantFileList[1], shallow=False)
        assert filecmp.cmp(quantFileList[1], quantFileList[2], shallow=False)

    @pytest.mark.dependency(name='genePrep', depends=['pclToDab'])
    def test_makeGenePrepFiles(self):
        cfg = TestSeekUtils.cfg
        sutils.checkConfig(cfg)
        sutils.makeGenePrepFiles(cfg, concurrency=5)
        prepDir = os.path.join(cfg.outDir, "prep")
        dabFileList = fnmatch.filter(os.listdir(cfg.dabDir), '*.dab')
        gavgFileList = fnmatch.filter(os.listdir(prepDir), '*.gavg')
        gpresFileList = fnmatch.filter(os.listdir(prepDir), '*.gpres')
        assert len(dabFileList) == len(gavgFileList)
        assert len(dabFileList) == len(gpresFileList)

    @pytest.mark.dependency(name='makeDB')
    def test_parallelMakeDB(self):
        cfg = TestSeekUtils.cfg
        sutils.checkConfig(cfg)
        # first create the per-thread db files based on partial dataset lists
        dbDirsToCombine = sutils.parallelMakePerThreadDB(cfg, concurrency=4)
        # check results
        threadWorkDir = os.path.join(cfg.outDir, 'thread_work')
        filePattern = os.path.join(threadWorkDir, 'db*')
        threadDbList = [d for d in glob.glob(filePattern) if os.path.isdir(d)]
        numFilesPerCollection = None
        for tdb in threadDbList:
            dbFiles = fnmatch.filter(os.listdir(tdb), '*.db')
            if numFilesPerCollection is None:
                assert len(dbFiles) > 0
                numFilesPerCollection = len(dbFiles)
            else:
                assert len(dbFiles) == numFilesPerCollection
        dsetFilename = os.path.basename(cfg.datasetsFile)
        threadDsetFiles = glob.glob(os.path.join(threadWorkDir, '*' + dsetFilename))
        threadDsetFiles.sort()
        assert len(threadDsetFiles) == len(threadDbList)

        # Now combine the DB files together
        sutils.parallelCombineThreadDBs(cfg, dbDirsToCombine, concurrency=8)
        dbFiles = fnmatch.filter(os.listdir(cfg.dbDir), '*.db')
        assert len(dbFiles) == numFilesPerCollection
        combinedListFile = os.path.join(cfg.outDir, 'combined_dataset_map.txt')
        dsetList = open(combinedListFile).read().splitlines()
        threadList = []
        for threadDsetFile in threadDsetFiles:
            threadDsets = open(threadDsetFile).read().splitlines()
            threadList.extend(threadDsets)
        assert dsetList == threadList

        # Verify the db files
        sutils.verifyCombinedDBs(cfg, dbDirsToCombine, concurrency=4)

    @pytest.mark.dependency(depends=['genePrep', 'makeDB'])
    def test_makePlatFiles(self):
        cfg = TestSeekUtils.cfg
        sutils.checkConfig(cfg)
        sutils.makePlatFiles(cfg)
        platDir = os.path.join(cfg.outDir, "plat")
        platFiles = fnmatch.filter(os.listdir(platDir), 'all_platforms.*')
        assert len(platFiles) == 4

    @pytest.mark.dependency(depends=['pclbin'])
    def teste_makeDsetSizeFile(self):
        cfg = TestSeekUtils.cfg
        sutils.checkConfig(cfg)
        sutils.makeDsetSizeFile(cfg)
        dsetSizeFile = os.path.join(cfg.outDir, 'dset_size.txt')
        lineCount = len(open(dsetSizeFile).readlines())
        pclFileList = fnmatch.filter(os.listdir(cfg.pclDir), '*.pcl')
        assert lineCount == len(pclFileList)
