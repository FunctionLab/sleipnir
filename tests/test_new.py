import pytest
import os
import sys
import time
import random
import subprocess
import numpy as np

testDir = os.path.abspath(os.path.dirname(__file__))
sleipnirDir = os.path.dirname(testDir)
testInputsDir = os.path.join(testDir, 'test_inputs')

class TestSeekRPC:
    cfgFile = None

    def setup_class(cls):
        # Step 01: Make the breast cancer example DB if needed
        sampleBcDir = os.path.join(testInputsDir, 'breast-cancer-sample')

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

    def test_one(self):
        assert True

