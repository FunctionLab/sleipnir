import pytest
import os
import sys
import time
import subprocess
import tempfile

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
from pytestHelper import createSampleDatabase
from rank_correlation import files_rank_correlation
from seek_rpc import SeekRPC
from seek_rpc.ttypes import SeekQuery, QueryParams, QueryStatus 
from seek_rpc.ttypes import SearchMethod, DistanceMeasure

use_tempfile = False
min_result_correlation = 0.95
testPort = 9123

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
        cmd = f'{sleipnirBin}/SeekRPC -c {seekrpcConfigFile} -p {testPort}'
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
        # Test the cpp client
        inputQueries = f'{sampleBcDir}/queries/input_queries.txt'
        outputResults = f'{sampleBcDir}/queries/seekrpc_results.txt'
        cmd = f'{sleipnirBin}/SeekRPCClient -p {testPort} -s sampleBC ' \
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

        # Test the python client
        cmd = f'python {seekRpcDir}/SeekRPCClient.py -p {testPort} -s sampleBC ' \
              f'-g 8091,6154,5810,9183'
        print(f'### {cmd}')
        clientProc = subprocess.Popen(cmd, shell=True)
        clientProc.wait()
        assert clientProc.returncode == 0

    def test_failed_query(self):
        # Test failed request through cpp client
        cmd = f'{sleipnirBin}/SeekRPCClient -p {testPort} -s invalidSpecies ' \
              f'-g 8091,6154,5810,9183'
        print(f'### {cmd}')
        clientProc = subprocess.Popen(cmd, shell=True)
        clientProc.wait()
        assert clientProc.returncode == 255  # retval of -1

        # Test failed request through python client
        cmd = f'python {seekRpcDir}/SeekRPCClient.py -p {testPort} -s invalidSpecies ' \
              f'-g 8091,6154,5810,9183'
        print(f'### {cmd}')
        clientProc = subprocess.Popen(cmd, shell=True)
        clientProc.wait()
        assert clientProc.returncode == 255  # retval of -1


    def test_async_client(self):
        # Run the query through the python rpc client
        from thrift.transport import TTransport, TSocket
        from thrift.protocol.TBinaryProtocol import TBinaryProtocol
        socket = TSocket.TSocket('localhost', testPort)
        transport = TTransport.TBufferedTransport(socket)
        protocol = TBinaryProtocol(transport)
        client = SeekRPC.Client(protocol)
        transport.open()

        params = QueryParams(distanceMeasure=DistanceMeasure.ZScoreHubbinessCorrected,
                            minQueryGenesFraction=0.5,
                            minGenomeFraction=0.5,
                            useGeneSymbols=False)
        genes = ['55755', '64859', '348654', '79791', '7756', '8555', '835', '5347']
        datasets = ['GSE45584.GPL6480', 'GSE24468.GPL570', 'GSE3744.GPL570']
        query = SeekQuery(species='sampleBC', genes=genes, datasets=datasets, parameters=params)

        # Do an async query using isQueryComplete to test
        task_id = client.seekQueryAsync(query)
        while True:
            complete = client.isQueryComplete(task_id)
            if complete:
                break;
            print("### Waiting for query to complete ...")
            time.sleep(.1)
        result = client.getQueryResult(task_id, block=True)
        assert result.success is True
        assert len(result.geneScores) > 0
        assert len(result.datasetWeights) == 3  # because query specified 3 datasets

        # Do an async query using non-blocking getQueryResult()
        genes = ['5884', '9575', '51343', '57805', '29980', '8091', '6154', '51776']
        query = SeekQuery(species='sampleBC', genes=genes, parameters=params)
        task_id = client.seekQueryAsync(query)
        while True:
            result = client.getQueryResult(task_id, block=False)
            if result.status is not QueryStatus.Incomplete:
                break;
            print("### Waiting for query to complete ...")
            time.sleep(.1)
        assert result.success is True
        assert len(result.geneScores) > 0
        assert len(result.datasetWeights) > 0

        transport.close()

    def test_simulate_weight(self):
        # Run the query through the python rpc client
        from thrift.transport import TTransport, TSocket
        from thrift.protocol.TBinaryProtocol import TBinaryProtocol
        socket = TSocket.TSocket('localhost', testPort)
        transport = TTransport.TBufferedTransport(socket)
        transport.open()
        protocol = TBinaryProtocol(transport)
        client = SeekRPC.Client(protocol)

        genes = ['55755', '64859', '348654', '79791', '7756', '8555', '835', '5347']

        # Run EqualWeighting without simulateWeights
        params = QueryParams(searchMethod=SearchMethod.EqualWeighting,
                            simulateWeights=False)
        query = SeekQuery(species='sampleBC', genes=genes, parameters=params)
        result = client.seekQuery(query)
        assert result.success is True
        assert len(result.geneScores) > 0
        assert len(result.datasetWeights) == 0  # because no simulate weights

        # Run EqualWeighting with simulateWeights
        params = QueryParams(searchMethod=SearchMethod.EqualWeighting,
                            simulateWeights=True)
        query = SeekQuery(species='sampleBC', genes=genes, parameters=params)
        result = client.seekQuery(query)
        assert result.success is True
        assert len(result.geneScores) > 0
        assert len(result.datasetWeights) > 0  # because simulate weights set

        # Run OrderStatistics without simulateWeights
        params = QueryParams(searchMethod=SearchMethod.OrderStatistics,
                            simulateWeights=False)
        query = SeekQuery(species='sampleBC', genes=genes, parameters=params)
        result = client.seekQuery(query)
        assert result.success is True
        assert len(result.geneScores) > 0
        assert len(result.datasetWeights) == 0  # because no simulate weights

        # Run OrderStatistics with simulateWeights
        params = QueryParams(searchMethod=SearchMethod.OrderStatistics,
                            simulateWeights=True)
        query = SeekQuery(species='sampleBC', genes=genes, parameters=params)
        result = client.seekQuery(query)
        assert result.success is True
        assert len(result.geneScores) > 0
        assert len(result.datasetWeights) > 0  # because simulate weights set