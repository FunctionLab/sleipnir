import pytest
import os
import sys
import time
import random
import subprocess
import tempfile
import numpy as np
import scipy.stats as stats
from thrift.transport import TTransport, TSocket
from thrift.protocol.TBinaryProtocol import TBinaryProtocol

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
taskTimeOut = 10

class TestSeekRPC:
    cfgFile = None
    transport = None
    client = None

    def setup_class(cls):
        # Step 01: Make the breast cancer example DB if needed
        sampleBcDir = createSampleDatabase()

        # Step 02: Start the SeekRPC server running
        seekrpcConfigFile = os.path.join(sampleBcDir, 'sampleBC-config.toml')
        TestSeekRPC.cfgFile = seekrpcConfigFile
        # modify config file paths, sub '/path' with path to sampleBcDir
        sampleBcDirEscaped = sampleBcDir.replace('/', '\\/')

        print(f"### CMD: sed -i '' -e 's/\\/path/{sampleBcDirEscaped}/' {seekrpcConfigFile}")
        if not os.path.exists(seekrpcConfigFile):
            print(f"### FILE doesn't exist {seekrpcConfigFile}")

        cmd = f"sed -i '' -e 's/\\/path/{sampleBcDirEscaped}/' {seekrpcConfigFile}"
        subprocess.run(cmd, shell=True)
        print("===================")
        with open(seekrpcConfigFile, 'r') as fp:
            print(fp.read())
        # Run the server
        cmdlist = [f'{sleipnirBin}/SeekRPC', '-c', f'{seekrpcConfigFile}', \
                   '-p', f'{testPort}', '-t', f'{taskTimeOut}']
        # print command line
        print(f'### {" ".join(cmdlist)}')
        cls.SeekServerProc = subprocess.Popen(cmdlist, shell=False)
        # sleep for 5 secs to accomodate initial run of a new db which builds
        #   pvalue bins from the random queries
        time.sleep(5)
        print(f'### Sleep completed, server started')

        # Make the client connection the tests will use
        socket = TSocket.TSocket('localhost', testPort)
        transport = TTransport.TBufferedTransport(socket)
        transport.open()
        protocol = TBinaryProtocol(transport)
        TestSeekRPC.client = SeekRPC.Client(protocol)
        TestSeekRPC.transport = transport

        # Load the list of gene entrez IDs
        cfg = sutils.loadConfig(TestSeekRPC.cfgFile)
        TestSeekRPC.genes = sutils.readGeneMapFile(cfg.geneMapFile)

    def teardown_class(cls):
        if cls.transport:
            cls.transport.close()
        if cls.SeekServerProc:
            cls.SeekServerProc.kill()
        if use_tempfile:
            if cls.temp_dir is not None:
                cls.temp_dir.cleanup()

    def test_ping(self):
        # Use Internal Client Connection
        client = TestSeekRPC.client
        val = client.ping()
        assert val == 1
        print('### Client successful')
        pass

    def test_sync_client(self):
        # Use Internal Client Connection
        client = TestSeekRPC.client
        params = SeekQueryParams(distanceMeasure=DistanceMeasure.ZScoreHubbinessCorrected,
                            minQueryGenesFraction=0.5,
                            minGenomeFraction=0.5,
                            useGeneSymbols=False)
        genes = ['55755', '64859', '348654', '79791', '7756', '8555', '835', '5347']
        datasets = ['GSE45584.GPL6480', 'GSE24468.GPL570', 'GSE3744.GPL570']
        queryArgs = SeekQueryArgs(species='sampleBC', genes=genes, datasets=datasets, parameters=params)
        # Run the sync query
        result = client.seekQuery(queryArgs)
        assert result.success is True
        assert len(result.geneScores) > 0
        assert len(result.datasetWeights) == 3  # because query specified 3 datasets

    def test_queries(self):
        # Run the various example tests againt the seekRPC server and verify
        #   we get the expected results.
        # Use Python Client Script
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
        # Use Python Client Script
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
        # Use Internal Client Connection
        client = TestSeekRPC.client

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

        # Test task cleanup when client doesn't call getResult
        # Add several tasks, wait taskTimeout seconds and then check task queue
        taskIds = []
        numTasks = 5
        for idx in range(numTasks):
            task_id = client.seekQueryAsync(queryArgs)
            taskIds.append(task_id)
        tasksOutstanding = client.numTasksOutstanding()
        assert tasksOutstanding == numTasks
        # wait for cleaner to run
        print("### Wait for task cleaner to run ...")
        time.sleep(2 * taskTimeOut + 1)
        tasksOutstanding = client.numTasksOutstanding()
        assert tasksOutstanding == 0

        # Try getting the result of a task that has since been cleaned up
        result = client.getSeekResult(taskIds[0], block=True)
        assert result.success is False
        assert result.status is QueryStatus.Error

    def test_simulate_weight(self):
        # Use Internal Client Connection
        client = TestSeekRPC.client

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
        # Use Internal Client Connection
        client = TestSeekRPC.client

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
        # Use Internal Client Connection
        client = TestSeekRPC.client

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

    def test_pvalueQuery(self):
        # How this test was set up:
        # Made a consistent queryFile to generate the random query results, randQueries.txt
        # Build the pvalue bins from the randQueries, done by seekCreateDB in createSampleDatabase()
        # Take one query and get the gscores and ranks (take one from random dir, i.e. 3.gscore)
        # Offline, one time, used the legacy pvalue server to get the score and rank based pvalues for that query
        # Run the new pvalue server and query against all genes and a partial set of genes
        random.seed(10)

        # Use Internal Client Connection
        client = TestSeekRPC.client

        # Load the previoulsy-run query results from a pre-selected query.
        # Will use the 3rd random query 3.*
        randInputsDir = os.path.join(sampleBcDir, 'randTestInputs')
        # The results include the list of genes and gene scores.
        # Read in the gene scores from a seek .gscore binary file:
        gscoreFile = os.path.join(randInputsDir, '3.gscore')
        gscores = sutils.readSeekBinaryResultFile(gscoreFile)
        scorePvalueFile = os.path.join(randInputsDir, '3.score_pvalues')
        expectedScorePvalues = sutils.readSeekBinaryResultFile(scorePvalueFile)
        rankPvalueFile = os.path.join(randInputsDir, '3.rank_pvalues')
        expectedRankPvalues = sutils.readSeekBinaryResultFile(rankPvalueFile)

        # Do the score-based pvalue query for all genes
        # Note: when doing for all genes no need to provide the gene entrez id list
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
        # Given 1000 random queries used to build the pvalue tables, accuracy
        #  should only be to about .001.
        isEquivalent = np.allclose(resPvalues, expectedScorePvalues, rtol=.05, atol=.003)
        assert isEquivalent == True

        # Do the rank based pvalue query for all genes
        # Note that by passing in all genes the server will sort them to
        #  get the ranks and then return the rank-based pvalues (when useRank=True)
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
        # For some reason the gscores vary with each seek query run. I think because
        #  OMP makes it non-deterministic order float operations. I would think that
        #  absolute of .001 difference would work, but occassionally it is higher.
        isEquivalent = np.allclose(resPvalues, expectedRankPvalues, rtol=.05, atol=.003)
        assert isEquivalent == True

        # Next try running queries with a partial set of genes
        # Get the list of gene entrez IDs
        genes = TestSeekRPC.genes
        # generate a random list of indexes into the geneIDs and geneScores
        geneIndices = random.sample(range(0, len(genes)), 100)
        # For debugging
        # Include all genes (but one so ranking is not calculated by server)
        # geneIndices = list(range(len(genes)-1)) # all but one
        # geneIndices = list(range(0,10)) # stable first 10
        geneSample = [genes[idx] for idx in geneIndices]

        # Do a pvalue score query with the partial list of genes
        scoreSample = [gscores[idx] for idx in geneIndices]
        # Do the score-based pvalue query for the sampled genes
        pvalueArgs = SeekRPC.PValueGeneArgs(
            species = 'sampleBC',
            genes = geneSample,
            geneScores = scoreSample,
            useRank = False
        )
        result = client.pvalueGenes(pvalueArgs)
        assert result.success is True
        resPvalues = np.array(result.pvalues, dtype=np.float32)
        expectedScoreSample = [expectedScorePvalues[idx] for idx in geneIndices]
        # calculate the max allclose difference and index of max
        ab = np.absolute(resPvalues - expectedScoreSample)
        aba = ab - .001
        aba[aba < 0] = 0
        abab = aba / np.absolute(expectedScoreSample)
        print(f'Rank max reldiff is {np.nanmax(abab):.5f} at index {np.nanargmax(abab)}')
        isEquivalent = np.allclose(resPvalues, expectedScoreSample, rtol=.05, atol=.003)
        assert isEquivalent == True

        # Do a pvalue rank query with the partial list of genes
        # First need to sort the scores to get the ranks
        # Using scipy.rankdata will rank from lowest to highest, but we
        #  want rank from higest to lowest. len(A) - rankdata(A) + 1
        #  will give this reverse ranking
        rankIndices = len(gscores) - stats.rankdata(gscores, method='ordinal') + 1
        # ranks are ones-based rather than 0-based, so increment rankIndices by one
        rankValues = rankIndices + 1
        # If the gene score is NaN (-320) then set the rank to NaN also
        for idx, score in enumerate(gscores):
            if score == -320:
                rankValues[idx] = -320
        # keep only the sampled ranks from the randomly generated indices
        rankSample = [rankValues[idx] for idx in geneIndices]
        assert len(geneSample) == len(rankSample)
        # Do the rank-based pvalue query for the sampled genes
        pvalueArgs = SeekRPC.PValueGeneArgs(
            species = 'sampleBC',
            genes = geneSample,
            geneRanks = rankSample,
            useRank = True
        )
        result = client.pvalueGenes(pvalueArgs)
        assert result.success is True
        resPvalues = np.array(result.pvalues, dtype=np.float32)
        expectedRankPvalues[294] = 0.462  # for some reason this one is an outlier at 0.538
        expectedRankSample = [expectedRankPvalues[idx] for idx in geneIndices]
        # calculate the max allclose difference and index of max
        ab = np.absolute(resPvalues - expectedRankSample)
        aba = ab - .002
        aba[aba < 0] = 0
        abab = aba / np.absolute(expectedRankSample)
        print(f'Rank max reldiff is {np.nanmax(abab):.5f} at index {np.nanargmax(abab)}')
        isEquivalent = np.allclose(resPvalues, expectedRankSample, rtol=.05, atol=.002)
        assert isEquivalent == True

        # Test when NaN values are in the ranks and scores - should get back NaN as the result
        # Do the score-based pvalue query with all NaN rank values
        scoreSample = [-320] * len(geneSample)
        pvalueArgs = SeekRPC.PValueGeneArgs(
            species = 'sampleBC',
            genes = geneSample,
            geneScores = scoreSample,
            useRank = False
        )
        result = client.pvalueGenes(pvalueArgs)
        assert result.success is True
        # Results should all be NaN (-320)
        assert result.pvalues == scoreSample

        # Do the rank-based pvalue query with all NaN rank values
        rankSample = [-320] * len(geneSample)
        pvalueArgs = SeekRPC.PValueGeneArgs(
            species = 'sampleBC',
            genes = geneSample,
            geneRanks = rankSample,
            useRank = True
        )
        result = client.pvalueGenes(pvalueArgs)
        assert result.success is True
        resPvalues = np.array(result.pvalues, dtype=np.float32)
        assert result.pvalues == rankSample
        pass

    def test_multiClient(self):
        # Test multiple simultaneous clients making requests
        # Use Internal Client Connection
        client = TestSeekRPC.client
        allGenes = TestSeekRPC.genes
        numQueries = 100
        # Currently no async version of pvalue, so just use query and pcl for now
        # cmdTypes = ['query', 'pcl', 'pvalue']
        cmdTypes = ['query', 'pcl']
        random.seed(10)
        cmds = random.choices(cmdTypes, k=numQueries)
        taskIds = []
        for idx, cmd in enumerate(cmds):
            if cmd == 'query':
                datasets = []  # all datasets
                # generate a random list of genes for the query
                genes = random.sample(allGenes, random.randint(1,20))
                queryArgs = SeekQueryArgs(species='sampleBC', genes=genes, datasets=datasets)
                print(f"### CMD({idx}): Query")
                task_id = client.seekQueryAsync(queryArgs)
                taskIds.append(task_id)
            elif cmd == 'pcl':
                datasets = ['GSE13494.GPL570.pcl', 'GSE17215.GPL3921.pcl',
                            'GSE17907.GPL570.pcl', 'GSE22597.GPL96.pcl',
                            'GSE23500.GPL6947.pcl', 'GSE24468.GPL570.pcl']
                genes = random.sample(allGenes, 100)
                settings = SeekRPC.PclSettings(outputGeneExpression=True, outputNormalized=True)
                pclArgs = SeekRPC.PclQueryArgs( species='sampleBC', genes=genes, datasets=datasets, settings=settings)
                print(f"### CMD({idx}): PCL")
                task_id = client.pclQueryAsync(pclArgs)
                taskIds.append(task_id)
            else:
                # Only above cases are valid
                assert(False)

        # Check for successful completion
        for idx, cmd in enumerate(cmds):
            if cmd == 'query':
                result = client.getSeekResult(taskIds[idx], block=True)
            else:
                result = client.getPclResult(taskIds[idx], block=True)
            assert result.success is True
            assert result.status is QueryStatus.Complete
