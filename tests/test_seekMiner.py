import pytest
import os
import sys
import time
import glob
import filecmp
import hashlib
import tempfile
import subprocess

testDir = os.path.abspath(os.path.dirname(__file__))
sleipnirDir = os.path.dirname(testDir)
sleipnirBin = os.path.join(sleipnirDir, 'Debug')
seekScriptsDir = os.path.join(sleipnirDir, 'scripts', 'seek')
sys.path.append(seekScriptsDir)
pytoolsDir = os.path.join(testDir, 'bioinform_tests', 'pytools')
sys.path.append(pytoolsDir)
import seekUtils as sutils
from pytestHelper import createSampleDatabase
from rank_correlation import files_rank_correlation
import compare_result_files as cmpResFiles

min_result_correlation = 0.95

class TestSeekMiner:
    cfg = None

    def setup_class(cls):
        # Step 01: Make the breast cancer example DB if needed
        sampleBcDir = createSampleDatabase()
        sampleConfigFile = os.path.join(sampleBcDir, 'sampleBC-config.toml')
        # modify config file paths, sub '/path' with path to sampleBcDir
        sampleBcDirEscaped = sampleBcDir.replace('/', '\\/')
        cmd = f"sed -i '' -e 's/\\/path/{sampleBcDirEscaped}/' {sampleConfigFile}"
        subprocess.run(cmd, shell=True)
        cfg = sutils.loadConfig(sampleConfigFile)
        cfg.inDir = sampleBcDir
        cfg.outDir = sampleBcDir
        cfg.binDir = sleipnirBin
        sutils.checkConfig(cfg)
        cls.sampleBcDir = sampleBcDir
        cls.cfg = cfg


    def test_queries(self):
        # Run the various example tests againt the seekRPC server and verify
        #   we get the expected results.
        sampleBcDir = TestSeekMiner.sampleBcDir
        cfg = TestSeekMiner.cfg

        inputQueries = f'{sampleBcDir}/queries/input_queries.txt'
        outputResults = f'{sampleBcDir}/results'
        os.makedirs(outputResults, exist_ok=True)

        print("## Run SeekMiner Queries ##")
        cmd = f'{sleipnirBin}/SeekMiner -x {cfg.datasetPlatMapFile} ' \
              f'-i {cfg.geneMapFile} -d {cfg.dbDir} -p {cfg.prepDir} ' \
              f'-P {cfg.platformDir} -Q {cfg.quantFile} -u {cfg.sinfoDir} ' \
              f'-U {cfg.gvarDir} -n {cfg.numDbFiles} -b 100 -V CV -I LOI -z z_score ' \
              f'-m -M -O -q {inputQueries} -o {outputResults} -Y -T 2 -t 1'
        print(f'### {cmd}')
        ret = subprocess.run(cmd, shell=True)
        assert ret.returncode == 0

        # combine all the query result files together
        filePattern = os.path.join(outputResults, "*.results.txt")
        resultFiles = sorted(glob.glob(filePattern), key=os.path.getmtime)
        resultFiles = " ".join(resultFiles)
        combinedResults = os.path.join(outputResults, 'all_res.txt')
        cmd = f'cat {resultFiles} > {combinedResults}'
        ret = subprocess.run(cmd, shell=True)
        assert ret.returncode == 0

        expectedResults = f'{sampleBcDir}/queries/seekminer_results.txt'
        corrs = files_rank_correlation(expectedResults, combinedResults)
        print('Result Correlations: {}'.format(corrs))
        correlation_errors = 0
        for corr in corrs:
            if corr < min_result_correlation:
                correlation_errors += 1
                print('ERROR: Result correlation too low')
        assert correlation_errors == 0

    def test_utilsRunSeekMiner(self):
        sampleBcDir = TestSeekMiner.sampleBcDir
        tmpDirObj = tempfile.TemporaryDirectory()
        tmpDir = tmpDirObj.name
        # grab 10 queries to run
        inputQFile = os.path.join(sampleBcDir, 'randQueries.txt')
        randQFile = os.path.join(tmpDir, 'shortQFile.txt')
        os.system(f'head -4 {inputQFile} > {randQFile}')
        sutils.runSeekMiner(TestSeekMiner.cfg, randQFile, tmpDir, concurrency=3)
        # See if 3.result.txt matches
        file1 = os.path.join(tmpDir, '3.results.txt')
        file2 = os.path.join(sampleBcDir, 'randTestInputs/3.results.txt')
        filesMatch = filecmp.cmp(file1, file2, shallow=False)
        if filesMatch:
            pass
        else:
            # check if they are very close
            pctDiff = cmpResFiles.get_pct_error(file1, file2, skiprows=1)
            assert pctDiff < 1  # less than 1% difference
        pass

    def test_buildSeekPrepStats(self):
        # Test the multi-threaded version of calculating the plat stats
        sampleBcDir = TestSeekMiner.sampleBcDir
        cfg = TestSeekMiner.cfg
        # make tmp dir for plat stats results
        tmpDirObj = tempfile.TemporaryDirectory()
        tmpDir = tmpDirObj.name
        tmpPlatDir = os.path.join(tmpDir,  'plat')
        os.makedirs(tmpPlatDir)
        # make list of db files
        dbList = glob.glob(os.path.join(sampleBcDir, 'db', '*.db'))
        dbListFile = os.path.join(tmpDir, 'dblist.txt')
        with open(dbListFile, 'w') as fp:
            fp.write(os.linesep.join(dbList))

        # run command
        cmd = f'{sleipnirBin}/SeekPrep -i {cfg.geneMapFile} -Q {cfg.quantFile} ' \
              f'-A {cfg.datasetPlatMapFile} -I {cfg.prepDir} -b {dbListFile} ' \
              f'-D {tmpPlatDir} -f -P'
        print(f'### {cmd}')
        ret = subprocess.run(cmd, shell=True)
        assert ret.returncode == 0

        # check against expected md5sums
        md5Expected = {
            "all_platforms.gplatavg" : "4cbd7c548934ef02d74f3f54906dc010",
            "all_platforms.gplatcount" : "deb0c431dd924497efb33c0bb610a9ed",
            "all_platforms.gplatorder" : "c20abb609c9de9c469e341bee3a182ad",
            "all_platforms.gplatstdev" : "87f50ee50fdadbec22d07c1600226661",
        }
        for fileName in md5Expected.keys():
            filePath = os.path.join(tmpPlatDir, fileName)
            with open(filePath, 'rb') as fp:
                data = fp.read()
            md5Result = hashlib.md5(data).hexdigest()
            assert md5Result == md5Expected[fileName]
        pass
