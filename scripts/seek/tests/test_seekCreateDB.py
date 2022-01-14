import pytest
import os
import sys
import time
import filecmp
import resource
import subprocess
import tempfile
testDir = os.path.dirname(__file__)
seekScriptsDir = os.path.dirname(testDir)
sys.path.append(seekScriptsDir)
import seekUtils as sutils

sleipnirDir = os.path.dirname(os.path.dirname(seekScriptsDir))
sleipnirBinDir = os.path.join(sleipnirDir, 'Debug')
seekRpcDir = os.path.join(sleipnirDir, 'tools/SeekRPC')
use_tempfile = False
testPort = 9123

class TestSeekCreateDB:
    temp_dir = None
    mockDir = ''
    cfg = None

    def setup_class(cls):
        # Check that os max open files is at least 1024
        if resource.getrlimit(resource.RLIMIT_NOFILE)[0] < 1024:
            print("Please increase max open files limit using bash command "
                  "'ulimit -n 1024' in order to run tests")
            assert(False)
        # Make a tmp directory and copy input mockDB files to it
        if use_tempfile:
            cls.temp_dir = tempfile.TemporaryDirectory()
            tmpDirName = cls.temp_dir.name
        else:
            username = os.environ.get('USER')
            tmpDirName = os.path.join('/tmp', username, 'testSeekCreateDb')
            if os.path.exists(tmpDirName):
                os.system(f'rm -rf {tmpDirName}/*')
            else:
                os.makedirs(tmpDirName)
        assert os.path.exists(tmpDirName)
        cls.mockDir = os.path.join(tmpDirName, 'mockDb')
        print(f'### Using temp directory {tmpDirName}')
        # copy the test mock directory to the temp directory
        inputMockDbDir = os.path.join(testDir, 'inputs/mockDb')
        os.system(f'cp -a {inputMockDbDir} {tmpDirName}')


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
        ret = subprocess.run(cmd, shell=True)
        assert ret.returncode == 0

        # fix up the config file paths
        seekrpcConfigFile = os.path.join(mockDbDir, 'mockDb-config.toml')
        # modify paths in config files, sub '/path' with path to mockDbDir
        mockDbDirEscaped = mockDbDir.replace('/', '\\/')
        cmd = f"sed -i '' -e 's/\\/path/{mockDbDirEscaped}/' {seekrpcConfigFile}"
        subprocess.run(cmd, shell=True)

        # Test created db using SeekMiner
        cfg = sutils.loadConfig(seekrpcConfigFile)
        cfg.inDir = mockDbDir
        cfg.outDir = mockDbDir
        cfg.binDir = os.path.join(sleipnirDir, 'Debug')
        cfg.datasetsFile = 'pcl_list.txt'
        sutils.checkConfig(cfg)
        # Run SeekMiner query
        queryFile = os.path.join(mockDbDir, 'query.txt')
        resultsDir = os.path.join(mockDbDir, 'results')
        os.makedirs(os.path.join(mockDbDir, 'results'))
        cmd = f'{sleipnirBinDir}/SeekMiner -x {cfg.datasetPlatMapFile} -i ' \
              f'{cfg.geneMapFile} -d {cfg.dbDir} -p {cfg.prepDir} ' \
              f'-P {cfg.platformDir} -Q {cfg.quantFile} -u {cfg.sinfoDir} ' \
              f'-U {cfg.gvarDir} -n 6 -b 20  -V CV -I LOI -z z_score ' \
              f'-m -M -O -q {queryFile} -o {resultsDir} -Y -T 2 -t 1'
        ret = subprocess.run(cmd, shell=True)
        assert ret.returncode == 0
        expected_results = os.path.join(mockDbDir, 'query_result.txt')
        seekminer_results = os.path.join(resultsDir, '0.results.txt')
        assert filecmp.cmp(expected_results, seekminer_results, shallow=False)

        # Test created db using SeekRPC
        # Run the server
        cmd = f'{sleipnirBinDir}/SeekRPC -p {testPort} -c {seekrpcConfigFile}'
        # Note: for Popen we need to use shell=False for Popen or the .kill()
        #  method on the process will kill the shell but not the child process
        #  which in this case is SeekRPC. The other option is to use shell=True
        #  but prepend the command with 'exec ' + cmd. This will exec SeekRPC
        #  into the shell process.
        cmdlist = [f'{sleipnirBinDir}/SeekRPC', '-p', f'{testPort}', '-c', f'{seekrpcConfigFile}']
        print(f'Run: {cmd}')
        SeekServerProc = subprocess.Popen(cmdlist, shell=False)
        time.sleep(.5)

        # Run the SeekRPC clients
        seekrpcResultsFile = os.path.join(resultsDir, 'seekRPC_results.txt')
        cmd = f'{sleipnirBinDir}/SeekRPCClient -p {testPort} -s mockDB ' \
              f'-q {queryFile} -o {seekrpcResultsFile}'
        clientProc = subprocess.Popen(cmd, shell=True)
        clientProc.wait()
        assert filecmp.cmp(expected_results, seekrpcResultsFile, shallow=False)

        # run the pcl python client
        pclResultsFile = os.path.join(resultsDir, 'pclClientResults.txt')
        cmd = f'python {seekRpcDir}/PclRpcClient.py -p {testPort} -s mockDB ' \
              f'-g 10998,90634 -q 23658,23659 -d GSE13494.GPL570,GSE17215.GPL3921 -o {pclResultsFile}'
        print(f"Run: {cmd}")
        clientProc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        clientProc.wait()
        output1, err = clientProc.communicate()
        print(f"## OUTPUT: {output1}")
        assert clientProc.returncode == 0

        # run the pvalue python client
        cmd = f'python {seekRpcDir}/PvalueRpcClient.py -p {testPort} -s mockDB ' \
              f'-g 10998,90634 -v 0.45,0.95'
        clientProc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        clientProc.wait()
        output1, err = clientProc.communicate()
        print(f"## OUTPUT: {output1}")
        assert clientProc.returncode == 0

        # Stop the RPC server
        SeekServerProc.kill()

