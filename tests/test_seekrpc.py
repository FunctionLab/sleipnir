import pytest
import os
import sys
import time
import glob
import subprocess
import filecmp
import tempfile

testDir = os.path.abspath(os.path.dirname(__file__))
testInputsDir = os.path.join(testDir, 'test_inputs')
testOutputsDir = os.path.join(testDir, 'test_outputs')
sleipnirDir = os.path.dirname(testDir)
sleipnirBin = os.path.join(sleipnirDir, 'Debug')
sampleBcDir = os.path.join(testOutputsDir, 'sampleBC')
seekScriptsDir = os.path.join(sleipnirDir, 'scripts', 'seek')
sys.path.append(seekScriptsDir)
pytoolsDir = os.path.join(testDir, 'bioinform_tests', 'pytools')
sys.path.append(pytoolsDir)
import seekUtils as sutils
from rank_correlation import files_rank_correlation

use_tempfile = False
min_result_correlation = 0.95

class TestSeekUtils:
    cfg = None

    def setup_class(cls):
        # Step 01: Make the breast cancer example DB if needed
        dbDir = os.path.join(sampleBcDir, 'db')
        dbFiles = glob.glob1(dbDir, '*.db')
        if len(dbFiles) < 100:
            assert False  # don't do this for now
            os.makedirs(sampleBcDir, exist_ok=True)
            inputBCFiles = os.path.join(testInputsDir, 'breast-cancer-sample')
            os.system(f'cp -a {inputBCFiles}/* {sampleBcDir}')
            # Create the sample DB
            cmd = f'python {seekScriptsDir}/seekCreateDB.py --all -g bc_gene_map.txt ' \
                f'-n 100 -m 4 -i {sampleBcDir} -o {sampBcDir} -b {sleipnirBin} --dab-use-gene-set'
            ret = subprocess.run(cmd, shell=True)
            assert ret.returncode == 0

        # Step 02: Start the SeekRPC server running
        cmd = f'pushd {sleipnirBin}; ./SeekRPC -c {sampleBcDir}/sampleBC-config.toml'
        cls.SeekServerProc = subprocess.Popen(cmd, shell=True)


    def teardown_class(cls):
        if cls.SeekServerProc:
            cls.SeekServerProc.kill()
        if use_tempfile:
            if cls.temp_dir is not None:
                cls.temp_dir.cleanup()

    def test_queries(self):
        # Run the various example tests againt the seekRPC server and verify
        #   we get the expected results.
        print("## Run Client ##")
        inputQueries = f'{sampleBcDir}/queries/input_queries.txt'
        outputResults = f'{sampleBcDir}/queries/seekrpc_results.txt'
        cmd = f'pushd {sleipnirBin}; ./SeekRPCClient -s sampleBC ' \
              f'-q {inputQueries} -o {outputResults}'
        clientProc = subprocess.Popen(cmd, shell=True)
        clientProc.wait()
        expectedResults = f'{sampleBcDir}/queries/seekminer_results.txt'
        corrs = files_rank_correlation(expectedResults, outputResults)
        print('Result Correlations: {}'.format(corrs))
        correlation_errors = 0
        for corr in corrs:
            if corr < min_result_correlation:
                correlation_errors += 1
                print('ERROR: Result correlation too low')
        assert correlation_errors == 0
