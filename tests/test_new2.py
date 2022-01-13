import pytest
import os
import sys
import time
import random
import subprocess
import tempfile
import numpy as np
import scipy.stats as stats

testDir = os.path.abspath(os.path.dirname(__file__))
sleipnirDir = os.path.dirname(testDir)
testInputsDir = os.path.join(testDir, 'test_inputs')
testOutputsDir = os.path.join(testDir, 'test_outputs')
sampleBcDir = os.path.join(testOutputsDir, 'sampleBC')
sleipnirBin = os.path.join(sleipnirDir, 'Debug')
seekScriptsDir = os.path.join(sleipnirDir, 'scripts', 'seek')
sys.path.append(seekScriptsDir)
pytoolsDir = os.path.join(testDir, 'bioinform_tests', 'pytools')
sys.path.append(pytoolsDir)
seekRpcDir = os.path.join(sleipnirDir, 'tools', 'SeekRPC')
seekRpcPyDir = os.path.join(seekRpcDir, 'gen-py')
sys.path.append(seekRpcPyDir)
import seekUtils as sutils
from pytestHelper import createSampleDatabase
from rank_correlation import files_rank_correlation
from seek_rpc import SeekRPC
from seek_rpc.ttypes import SeekQueryArgs, SeekQueryParams, QueryStatus
from seek_rpc.ttypes import SearchMethod, DistanceMeasure

use_tempfile = False
min_result_correlation = 0.95
testPort = 9123

class TestSeekRPC:
    cfgFile = None

    def setup_class(cls):
        # Step 01: Make the breast cancer example DB if needed
        sampleBcDir = createSampleDatabase()

        # Step 02: Start the SeekRPC server running
        seekrpcConfigFile = os.path.join(sampleBcDir, 'sampleBC-config.toml')
        TestSeekRPC.cfgFile = seekrpcConfigFile
        # modify config file paths, sub '/path' with path to sampleBcDir
        sampleBcDirEscaped = sampleBcDir.replace('/', '\\/')
        
        print(f"### CMD: sed -i '' -e 's/\\/path/{sampleBcDirEscaped}/' {seekrpcConfigFile}")
        if os.path.exists(seekrpcConfigFile):
            with open(seekrpcConfigFile, 'r') as fp:
                print(fp.read())
        else:
            print(f"### FILE doesn't exist {seekrpcConfigFile}")

        cmd = f"sed -i '' -e 's/\\/path/{sampleBcDirEscaped}/' {seekrpcConfigFile}"
        subprocess.run(cmd, shell=True)
        print("===================")
        with open(seekrpcConfigFile, 'r') as fp:
            print(fp.read())
        # Run the server
        cmd = f'{sleipnirBin}/SeekRPC -c {seekrpcConfigFile} -p {testPort}'
        print(f'### {cmd}')
        cls.SeekServerProc = subprocess.Popen(cmd, shell=True)
        # sleep for 5 secs to accomodate initial run of a new db which builds
        #   pvalue bins from the random queries
        time.sleep(5)
        print(f'### Sleep completed, server started')

    def teardown_class(cls):
        if cls.SeekServerProc:
            cls.SeekServerProc.kill()
        if use_tempfile:
            if cls.temp_dir is not None:
                cls.temp_dir.cleanup()


    def test_queries(self):
        pass