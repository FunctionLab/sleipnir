import pytest
import os
import sys
import time
import subprocess
import tempfile
import numpy as np

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
from seek_rpc.ttypes import SeekQueryArgs, SeekQueryParams, QueryStatus
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
        time.sleep(3)

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

        params = SeekQueryParams(distanceMeasure=DistanceMeasure.ZScoreHubbinessCorrected,
                            minQueryGenesFraction=0.5,
                            minGenomeFraction=0.5,
                            useGeneSymbols=False)
        genes = ['55755', '64859', '348654', '79791', '7756', '8555', '835', '5347']
        datasets = ['GSE45584.GPL6480', 'GSE24468.GPL570', 'GSE3744.GPL570']
        queryArgs = SeekQueryArgs(species='sampleBC', genes=genes, datasets=datasets, parameters=params)

        # Do an async query using isQueryComplete to test
        task_id = client.seekQueryAsync(queryArgs)
        while True:
            complete = client.isQueryComplete(task_id)
            if complete:
                break;
            print("### Waiting for query to complete ...")
            time.sleep(.1)
        result = client.getSeekResult(task_id, block=True)
        assert result.success is True
        assert len(result.geneScores) > 0
        assert len(result.datasetWeights) == 3  # because query specified 3 datasets

        # Do an async query using non-blocking getSeekResult()
        genes = ['5884', '9575', '51343', '57805', '29980', '8091', '6154', '51776']
        queryArgs = SeekQueryArgs(species='sampleBC', genes=genes, parameters=params)
        task_id = client.seekQueryAsync(queryArgs)
        while True:
            result = client.getSeekResult(task_id, block=False)
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
        params = SeekQueryParams(searchMethod=SearchMethod.EqualWeighting,
                            simulateWeights=False)
        queryArgs = SeekQueryArgs(species='sampleBC', genes=genes, parameters=params)
        result = client.seekQuery(queryArgs)
        assert result.success is True
        assert len(result.geneScores) > 0
        assert len(result.datasetWeights) == 0  # because no simulate weights

        # Run EqualWeighting with simulateWeights
        params = SeekQueryParams(searchMethod=SearchMethod.EqualWeighting,
                            simulateWeights=True)
        queryArgs = SeekQueryArgs(species='sampleBC', genes=genes, parameters=params)
        result = client.seekQuery(queryArgs)
        assert result.success is True
        assert len(result.geneScores) > 0
        assert len(result.datasetWeights) > 0  # because simulate weights set

        # Run OrderStatistics without simulateWeights
        params = SeekQueryParams(searchMethod=SearchMethod.OrderStatistics,
                            simulateWeights=False)
        queryArgs = SeekQueryArgs(species='sampleBC', genes=genes, parameters=params)
        result = client.seekQuery(queryArgs)
        assert result.success is True
        assert len(result.geneScores) > 0
        assert len(result.datasetWeights) == 0  # because no simulate weights

        # Run OrderStatistics with simulateWeights
        params = SeekQueryParams(searchMethod=SearchMethod.OrderStatistics,
                            simulateWeights=True)
        queryArgs = SeekQueryArgs(species='sampleBC', genes=genes, parameters=params)
        result = client.seekQuery(queryArgs)
        assert result.success is True
        assert len(result.geneScores) > 0
        assert len(result.datasetWeights) > 0  # because simulate weights set

    def runPclQuery(client, pclArgs):
        result = client.pclQuery(pclArgs)
        assert result.success is True
        assert len(result.datasetSizes) > 0
        return result

    def checkPclVals(vals, expectedOutputFile):
        # Read in expected result from test_input file
        with open(os.path.join(testInputsDir, expectedOutputFile)) as fp:
            lines = fp.readlines()
        expectedVals = [float(item) for item in lines]
        assert np.allclose(vals, expectedVals, atol=1e-6)
        print("Success!")


    def test_pclQueries(self):
        # Run the various example tests againt the PclRPC server and verify
        #   we get the expected results.
        print("## Run Client ##")
        # Run the query through the python rpc client
        from thrift.transport import TTransport, TSocket
        from thrift.protocol.TBinaryProtocol import TBinaryProtocol
        socket = TSocket.TSocket('localhost', testPort)
        transport = TTransport.TBufferedTransport(socket)
        transport.open()
        protocol = TBinaryProtocol(transport)
        client = SeekRPC.Client(protocol)

        datasets = ['GSE13494.GPL570.pcl', 'GSE17215.GPL3921.pcl']
        genes = ['10998', '10994']
        queryGenes = ['23658', '23659']

        # Test gene expression
        settings = SeekRPC.PclSettings(
            rbp = -1,
            outputNormalized = True,
            outputGeneExpression = True,
        )
        pclArgs = SeekRPC.PclQueryArgs(
            species='sampleBC',
            genes=genes,
            datasets=datasets,
            settings=settings)
        result = TestSeekRPC.runPclQuery(client, pclArgs)
        TestSeekRPC.checkPclVals(result.geneExpressions, "pclTestGeneExpr.txt")

        # Test gene coexpression
        settings = SeekRPC.PclSettings(
            rbp = -1,
            outputNormalized = True,
            outputGeneCoexpression = True,
        )
        pclArgs = SeekRPC.PclQueryArgs(
            species='sampleBC',
            genes=genes,
            queryGenes=queryGenes,  # queryGenes must be present for geneCoexpression calc
            datasets=datasets,
            settings=settings)
        result = TestSeekRPC.runPclQuery(client, pclArgs)
        TestSeekRPC.checkPclVals(result.geneCoexpressions, "pclTestGeneCoExpr.txt")

        # test query expression
        settings = SeekRPC.PclSettings(
            rbp = -1,
            outputNormalized = True,
            outputQueryExpression = True,
        )
        pclArgs = SeekRPC.PclQueryArgs(
            species='sampleBC',
            queryGenes=queryGenes,
            datasets=datasets,
            settings=settings)
        result = TestSeekRPC.runPclQuery(client, pclArgs)
        TestSeekRPC.checkPclVals(result.queryExpressions, "pclTestQueryExpr.txt")

        # test query Coexpression
        settings = SeekRPC.PclSettings(
            rbp = -1,
            outputNormalized = True,
            outputQueryCoexpression = True,
        )
        pclArgs = SeekRPC.PclQueryArgs(
            species='sampleBC',
            queryGenes=queryGenes,
            datasets=datasets,
            settings=settings)
        result = TestSeekRPC.runPclQuery(client, pclArgs)
        TestSeekRPC.checkPclVals(result.queryCoexpressions, "pclTestQueryCoExpr.txt")

    def test_pclSettings(self):
    # Test that uses different dataset names, like with .pcl, .pcl.bin and no pcl.
    # Also within settings file specifying pcl dir as /pcl and /pclbin
        # Test using the python client
        cmd = f'python {seekRpcDir}/PclRpcClient.py -p {testPort} -s sampleBC ' \
              f'-d GSE13494.GPL570.pcl,GSE17215.GPL3921.pcl ' \
              f'-g 10998,10994 -q 23658,23659 --zexp --zcoexp'
        print(f'### {cmd}')
        clientProc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        clientProc.wait()
        assert clientProc.returncode == 0
        output1, err = clientProc.communicate()

        cmd = f'python {seekRpcDir}/PclRpcClient.py -p {testPort} -s sampleBC ' \
              f'-d GSE13494.GPL570,GSE17215.GPL3921 ' \
              f'-g 10998,10994 -q 23658,23659 --zexp --zcoexp'
        print(f'### {cmd}')
        clientProc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        clientProc.wait()
        assert clientProc.returncode == 0
        output2, err = clientProc.communicate()
        out2 = output2.decode("utf-8")
        assert out2.find('Gene Expressions') > 0
        assert output1 == output2

    def test_pclAsyncClient(self):
        # Run the query through the python rpc client
        from thrift.transport import TTransport, TSocket
        from thrift.protocol.TBinaryProtocol import TBinaryProtocol
        socket = TSocket.TSocket('localhost', testPort)
        transport = TTransport.TBufferedTransport(socket)
        protocol = TBinaryProtocol(transport)
        client = SeekRPC.Client(protocol)
        transport.open()

        datasets = ['GSE13494.GPL570.pcl', 'GSE17215.GPL3921.pcl']
        genes = ['10998', '10994']
        queryGenes = ['23658', '23659']

        # Test gene expression
        settings = SeekRPC.PclSettings(
            rbp = -1,
            outputNormalized = True,
            outputGeneExpression = True,
        )
        pclArgs = SeekRPC.PclQueryArgs(
            species='sampleBC',
            genes=genes,
            datasets=datasets,
            settings=settings
        )

        # Do an async query using isQueryComplete to test
        task_id = client.pclQueryAsync(pclArgs)
        while True:
            complete = client.isQueryComplete(task_id)
            if complete:
                break;
            print("### Waiting for query to complete ...")
            time.sleep(.1)
        result = client.getPclResult(task_id, block=True)
        assert result.success is True
        assert len(result.datasetSizes) > 0
        TestSeekRPC.checkPclVals(result.geneExpressions, "pclTestGeneExpr.txt")

        # Do an async query using non-blocking getPclResult()
        pclArgs.settings.outputGeneExpression = False
        pclArgs.settings.outputGeneCoexpression = True
        pclArgs.queryGenes = queryGenes
        task_id = client.pclQueryAsync(pclArgs)
        while True:
            result = client.getPclResult(task_id, block=False)
            if result.status is not QueryStatus.Incomplete:
                break;
            print("### Waiting for query to complete ...")
            time.sleep(.1)
        assert result.success is True
        assert len(result.datasetSizes) > 0
        TestSeekRPC.checkPclVals(result.geneCoexpressions, "pclTestGeneCoExpr.txt")

        transport.close()

    def readSeekBinaryResultFile(dataFile):
        vals = []
        with open(dataFile, 'rb') as f:
            # The first 8 byte (long int) is the number of elements stored
            headerVals = np.fromfile(f, count=1, dtype=np.ulonglong)
            numVals = headerVals[0]
            # The remaining are 4 byte float values, numVal of them
            vals = np.fromfile(f, dtype=np.float32)
            assert len(vals) == numVals
        return vals

    def test_pvalueQuery(self):
        # x - Need a consistent queryFile to generate the random query results
        # x - Need a query and query results file (use one from random dir)
        # x - Take the specific query and get the gscores and ranks (take one from random dir)
        # x - Using the legacy pvalue server get the pvalues for scores and ranks for that query
        # Run the new pvalue server and query all and partial

        # Query gene score pvalues with all gscores
        # Run the queries through the python rpc client
        from thrift.transport import TTransport, TSocket
        from thrift.protocol.TBinaryProtocol import TBinaryProtocol
        socket = TSocket.TSocket('localhost', testPort)
        transport = TTransport.TBufferedTransport(socket)
        transport.open()
        protocol = TBinaryProtocol(transport)
        client = SeekRPC.Client(protocol)

        # Load the query results from a pre-selected canned query.
        randInputsDir = os.path.join(sampleBcDir, 'randTestInputs')
        # The results include the list of genes and gene scores.
        # Read in the gene scores from a seek .gscore binary file:
        gscoreFile = os.path.join(randInputsDir, '3.gscore')
        gscores = TestSeekRPC.readSeekBinaryResultFile(gscoreFile)

        # Given 1000 random queries used to build the pvalue tables, accuracy
        #  should only be to about .001.
        scorePvalueFile = os.path.join(randInputsDir, '3.score_pvalues')
        expectedScorePvalues = TestSeekRPC.readSeekBinaryResultFile(scorePvalueFile)
        rankPvalueFile = os.path.join(randInputsDir, '3.rank_pvalues')
        expectedRankPvalues = TestSeekRPC.readSeekBinaryResultFile(rankPvalueFile)

        # Do the score based pvalue query for all genes
        pvalueArgs = SeekRPC.PValueGeneArgs(
            species = 'sampleBC',
            # don't specify the geneIDs if getting scores for all genes, in which
            #  case it is assumed the gscores are in the gene_map order
            geneScores = gscores,
            useRank = False
        )
        result = client.pvalueGenes(pvalueArgs)
        assert result.success is True
        resPvalues = np.array(result.pvalues, dtype=np.float32)
        isEquivalent = np.allclose(resPvalues, expectedScorePvalues, rtol=.05, atol=.001)
        assert isEquivalent == True

        # Do the rank based pvalue query for all genes
        pvalueArgs = SeekRPC.PValueGeneArgs(
            species = 'sampleBC',
            # don't specify the geneIDs if getting scores for all genes, in which
            #  case it is assumed the gscores are in the gene_map order
            geneScores = gscores,
            useRank = True
        )
        result = client.pvalueGenes(pvalueArgs)
        assert result.success is True
        resPvalues = np.array(result.pvalues, dtype=np.float32)
        isEquivalent = np.allclose(resPvalues, expectedRankPvalues, rtol=.05, atol=.001)
        import pdb; pdb.set_trace()
        assert isEquivalent == True
        pass