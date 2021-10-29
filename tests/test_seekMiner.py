import pytest
import os
import sys
import time
import glob
import filecmp
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

min_result_correlation = 0.95

class TestSeekMiner:
    cfg = None

    def setup_class(cls):
        # Step 01: Make the breast cancer example DB if needed
        sampleBcDir = createSampleDatabase()
        sampleConfigFile = os.path.join(sampleBcDir, 'sampleBC-config.toml')
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
        assert filecmp.cmp(file1, file2, shallow=False)
