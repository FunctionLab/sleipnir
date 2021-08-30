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
pclServerDir = os.path.join(sleipnirDir, 'tools', 'PCLServer')
pclRpcPyDir = os.path.join(pclServerDir, 'gen-py')
sys.path.append(pclRpcPyDir)
from pytestHelper import createSampleDatabase
# from rank_correlation import files_rank_correlation
from pcl_rpc import PclRPC
from pcl_rpc.ttypes import PclQueryArgs, PclSettings, PclStatus 

use_tempfile = False
min_result_correlation = 0.95
testPort = 9124

class TestSeekRPC:
    cfg = None

    def setup_class(cls):
        # Step 01: Make the breast cancer example DB if needed
        sampleBcDir = createSampleDatabase()

        # Step 02: Start the PclRPC server running
        pclRpcConfigFile = os.path.join(sampleBcDir, 'sampleBC-config.toml')
        # modify config file paths, sub '/path' with path to sampleBcDir
        sampleBcDirEscaped = sampleBcDir.replace('/', '\\/')
        cmd = f"sed -i '' -e 's/\\/path/{sampleBcDirEscaped}/' {pclRpcConfigFile}"
        subprocess.run(cmd, shell=True)
        # Run the server
        cmd = f'{sleipnirBin}/PclRpcServer -c {pclRpcConfigFile} -p {testPort}'
        cls.PclServerProc = subprocess.Popen(cmd, shell=True)
        print(f'### {cmd}')
        time.sleep(.5)


    def teardown_class(cls):
        if cls.PclServerProc:
            cls.PclServerProc.kill()
        if use_tempfile:
            if cls.temp_dir is not None:
                cls.temp_dir.cleanup()

    def runQuery(client, pclArgs):
        result = client.pclQuery(pclArgs)
        assert result.success is True
        assert len(result.datasetSizes) > 0
        return result

    def checkVals(vals, expectedOutputFile):
        # Read in expected result from test_input file
        with open(os.path.join(testInputsDir, expectedOutputFile)) as fp:
            lines = fp.readlines()
        expectedVals = [float(item) for item in lines]
        assert np.allclose(vals, expectedVals, atol=1e-6)
        print("Success!")


    def test_queries(self):
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
        client = PclRPC.Client(protocol)

        datasets = ['GSE13494.GPL570.pcl', 'GSE17215.GPL3921.pcl']
        genes = ['10998', '10994']
        queryGenes = ['23658', '23659']

        # Test gene expression
        settings = PclRPC.PclSettings(
            rbp = -1,
            outputNormalized = True,
            outputGeneExpression = True,
        )
        pclArgs = PclRPC.PclQueryArgs(
            species='sampleBC',
            genes=genes,
            datasets=datasets,
            settings=settings)
        result = TestSeekRPC.runQuery(client, pclArgs)
        TestSeekRPC.checkVals(result.geneExpressions, "pclTestGeneExpr.txt")

        # Test gene coexpression
        settings = PclRPC.PclSettings(
            rbp = -1,
            outputNormalized = True,
            outputGeneCoexpression = True,
        )
        pclArgs = PclRPC.PclQueryArgs(
            species='sampleBC',
            genes=genes,
            queryGenes=queryGenes,  # queryGenes must be present for geneCoexpression calc
            datasets=datasets,
            settings=settings)
        result = TestSeekRPC.runQuery(client, pclArgs)
        TestSeekRPC.checkVals(result.geneCoexpressions, "pclTestGeneCoExpr.txt")

        # test query expression
        settings = PclRPC.PclSettings(
            rbp = -1,
            outputNormalized = True,
            outputQueryExpression = True,
        )
        pclArgs = PclRPC.PclQueryArgs(
            species='sampleBC',
            queryGenes=queryGenes,
            datasets=datasets,
            settings=settings)
        result = TestSeekRPC.runQuery(client, pclArgs)
        TestSeekRPC.checkVals(result.queryExpressions, "pclTestQueryExpr.txt")

        # test query Coexpression
        settings = PclRPC.PclSettings(
            rbp = -1,
            outputNormalized = True,
            outputQueryCoexpression = True,
        )
        pclArgs = PclRPC.PclQueryArgs(
            species='sampleBC',
            queryGenes=queryGenes,
            datasets=datasets,
            settings=settings)
        result = TestSeekRPC.runQuery(client, pclArgs)
        TestSeekRPC.checkVals(result.queryCoexpressions, "pclTestQueryCoexpr.txt")

# Test that uses different dataset names, like with .pcl, .pcl.bin and no pcl.
# Also within settings file specifying pcl dir as /pcl and /pclbin
    def test_settings(self):
        # Test using the python client
        cmd = f'python {pclServerDir}/PclRpcClient.py -p {testPort} -s sampleBC ' \
              f'-d GSE13494.GPL570.pcl,GSE17215.GPL3921.pcl ' \
              f'-g 10998,10994 -q 23658,23659 --zexp --zcoexp'
        print(f'### {cmd}')
        clientProc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        clientProc.wait()
        assert clientProc.returncode == 0
        output1, err = clientProc.communicate()

        cmd = f'python {pclServerDir}/PclRpcClient.py -p {testPort} -s sampleBC ' \
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


