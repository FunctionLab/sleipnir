import os
import sys
import random
import argparse
import datetime as dt
import numpy as np
import scipy.stats as stats
from thrift.transport import TTransport, TSocket
from thrift.protocol.TBinaryProtocol import TBinaryProtocol

testDir = os.path.abspath(os.path.dirname(__file__))
sleipnirDir = os.path.dirname(testDir)
seekScriptsDir = os.path.join(sleipnirDir, 'scripts', 'seek')
sys.path.append(seekScriptsDir)
seekRpcDir = os.path.join(sleipnirDir, 'tools', 'SeekRPC', 'gen-py')
sys.path.append(seekRpcDir)

import seekUtils as sutils
from seek_rpc import SeekRPC
from seek_rpc.ttypes import SeekQueryArgs


def main(args):
    allGenes = sutils.readGeneMapFile(args.geneMapFile)
    dsetTuples = sutils.readDatasetList(args.datasetsFile)
    allDatasets = [tup[0] for tup in dsetTuples]
    cmdTypes = ['query', 'pcl']

    # Make the client connection the tests will use
    socket = TSocket.TSocket('localhost', args.port)
    transport = TTransport.TBufferedTransport(socket)
    transport.open()
    protocol = TBinaryProtocol(transport)
    client = SeekRPC.Client(protocol)
    transport = transport

    endTime = dt.datetime.min
    if args.durationHours is not None:
        endTime = dt.datetime.now() + dt.timedelta(hours=args.durationHours)
        args.numRequests = 0

    cmdIdx = 0
    # while cmdIdx < args.numRequests:
    while endTime > dt.datetime.now() or cmdIdx < args.numRequests:
        numQueries = 5
        taskIds = []
        cmds = random.choices(cmdTypes, k=numQueries)
        # Run a number of seek and pcl queries
        for idx, cmd in enumerate(cmds):
            if cmd == 'query':
                datasets = []  # all datasets
                # generate a random list of genes for the query
                genes = random.sample(allGenes, random.randint(1,20))
                queryArgs = SeekQueryArgs(species=args.species, genes=genes, datasets=datasets)
                print(f"### CMD({cmdIdx}): Query")
                task_id = client.seekQueryAsync(queryArgs)
                taskIds.append(task_id)
            elif cmd == 'pcl':
                datasets = random.sample(allDatasets, random.randint(2,5))
                genes = random.sample(allGenes, random.randint(20,100))
                settings = SeekRPC.PclSettings(outputGeneExpression=True, outputNormalized=True)
                pclArgs = SeekRPC.PclQueryArgs(species=args.species, genes=genes, datasets=datasets, settings=settings)
                print(f"### CMD({cmdIdx}): PCL")
                task_id = client.pclQueryAsync(pclArgs)
                taskIds.append(task_id)
            else:
                # Only above cases are valid
                assert(False)
            cmdIdx += 1
        # Then add a pvalue query
        genes = random.sample(allGenes, 100)
        scores = np.random.uniform(low=-1.0, high=1.0, size=100)
        pvalueArgs = SeekRPC.PValueGeneArgs(species=args.species, genes=genes, geneScores=scores)
        print(f"### CMD(***): Pvalue")
        result = client.pvalueGenes(pvalueArgs)
        # Check for successful completion
        if result.success is False:
            print(f"Pvalue Query Error: {result.status}: {result.statusMsg}")
        # Check the async tasks completion
        for idx, cmd in enumerate(cmds):
            if cmd == 'query':
                result = client.getSeekResult(taskIds[idx], block=True)
            else:
                result = client.getPclResult(taskIds[idx], block=True)
            if result.success is False:
                print(f"\tQuery Error: {result.status}: {result.statusMsg}")


if __name__=="__main__":
    argParser = argparse.ArgumentParser()
    argParser.add_argument('--port', '-p', default=9090, type=int,
                           help='server port')
    argParser.add_argument('--species', '-s', default='human', type=str,
                           help='species name')
    argParser.add_argument('--geneMapFile', '-g', type=str, required=True, default=None,
                           help='Text file containing the list of genes in the database')
    argParser.add_argument('--datasetsFile', '-d', type=str, required=True, default=None,
                           help='Text file containing the list datasets')
    argParser.add_argument('--numRequests', '-n', default=100, type=int,
                           help='total number of query requests to make')
    argParser.add_argument('--durationHours', '-t', default=None, type=float,
                           help='total number of hours to run')
    args = argParser.parse_args()
    main(args)
