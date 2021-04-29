import pytest
import os
import sys
import time
import subprocess
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
from pytestHelper import createSampleDatabase
from rank_correlation import files_rank_correlation

use_tempfile = False
min_result_correlation = 0.95

class TestSeekRPC:
    cfg = None

    def setup_class(cls):
        # Step 01: Make the breast cancer example DB if needed
        sampleBcDir = createSampleDatabase()

        # Step 02: Start the SeekRPC server running
        seekrpcConfigFile = os.path.join(sampleBcDir, 'sampleBC-config.toml')
        # modify config file paths, sub '/path' with path to sampleBcDir
        sampleBcDirEscaped = sampleBcDir.replace('/', '\\/')
        cmd = f"sed -i '' -e 's/\\/path/{sampleBcDirEscaped}/' {seekrpcConfigFile}"
        subprocess.run(cmd, shell=True)
        # Run the server
        cmd = f'{sleipnirBin}/SeekRPC -c {seekrpcConfigFile}'
        cls.SeekServerProc = subprocess.Popen(cmd, shell=True)
        print(f'### {cmd}')
        time.sleep(.5)


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
        cmd = f'{sleipnirBin}/SeekRPCClient -s sampleBC ' \
              f'-q {inputQueries} -o {outputResults}'
        print(f'### {cmd}')
        clientProc = subprocess.Popen(cmd, shell=True)
        clientProc.wait()
        assert clientProc.returncode == 0
        expectedResults = f'{sampleBcDir}/queries/seekminer_results.txt'
        corrs = files_rank_correlation(expectedResults, outputResults)
        print('Result Correlations: {}'.format(corrs))
        correlation_errors = 0
        for corr in corrs:
            if corr < min_result_correlation:
                correlation_errors += 1
                print('ERROR: Result correlation too low')
        assert correlation_errors == 0
