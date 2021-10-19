import pytest
import os
import sys
import time
import glob
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
        cls.sampleBcDir = createSampleDatabase()


    def test_queries(self):
        # Run the various example tests againt the seekRPC server and verify
        #   we get the expected results.
        sampleBcDir = TestSeekMiner.sampleBcDir

        inputQueries = f'{sampleBcDir}/queries/input_queries.txt'
        outputResults = f'{sampleBcDir}/results'
        os.makedirs(outputResults, exist_ok=True)

        cfg = sutils.getDefaultConfig()
        cfg.inDir = sampleBcDir
        cfg.outDir = sampleBcDir
        cfg.binDir = sleipnirBin
        cfg.geneMapFile = 'bc_gene_map.txt'
        sutils.checkConfig(cfg)

        print("## Run SeekMiner Queries ##")
        cmd = f'{sleipnirBin}/SeekMiner -x {cfg.datasetPlatMapFile} -i ' \
              f'{cfg.geneMapFile} -d {cfg.dbDir} -p {sampleBcDir}/prep ' \
              f'-P {sampleBcDir}/plat -Q {cfg.quantFile} -u {sampleBcDir}/sinfo ' \
              f'-U {sampleBcDir}/gvar -n 100 -b 100  -V CV -I LOI -z z_score ' \
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
