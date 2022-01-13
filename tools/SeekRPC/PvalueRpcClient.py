import sys
import time
import argparse
import importlib
from thrift.transport import TTransport, TSocket
from thrift.protocol.TBinaryProtocol import TBinaryProtocol
# Since the path gen-py/seek_rpc has a hyphen in it we need to use
#   the import_module() function rather than the usual import call
SeekRPC = importlib.import_module('gen-py.seek_rpc.SeekRPC')
ttypes = importlib.import_module('gen-py.seek_rpc.ttypes')
constants = importlib.import_module('gen-py.seek_rpc.constants')

host = 'localhost'

'''
Query parameters and options:
struct PValueGeneArgs {
    1: required string species = "Unknown";
    2: optional list<string> genes;   // EntrezIDs
    3: optional list<double> geneScores;  // Seek query gene scores
    4: optional list<i32> geneRanks;  // Seek query gene ranks
    5: optional bool useRank = false;
}

struct PValueDatasetArgs {
    1: required string species = "Unknown";
    2: optional list<string> datasets;   // dataset names
    3: optional list<double> datasetWeights;  // Seek query result weights
}

struct PValueResult {
    1: required bool success;
    2: optional list<double> pvalues;
    3: optional QueryStatus status;
    4: optional string statusMsg;
}
'''

def runQuery(args):
    pvalueArgs = SeekRPC.PValueGeneArgs(species=args.species,
                                  genes=args.genes,
                                  geneRanks=args.ranks,
                                  geneScores=args.scores,
                                  useRank=args.useRank,)

    socket = TSocket.TSocket(host, args.port)
    transport = TTransport.TBufferedTransport(socket)
    protocol = TBinaryProtocol(transport)
    client = SeekRPC.Client(protocol)
    transport.open()
    version = client.getRpcVersion()
    assert version == constants.RPCVersion
    retval = -1

    result = client.pvalueGenes(pvalueArgs)
    if result.success is True:
        print(f'result: {result}')
        for idx, pval in enumerate(result.pvalues):
            print(f"{args.genes[idx]} <--> {pval:.04f}")
        retval = 0
    else:
        print(f'query error: {result.statusMsg}')
        retval = -1

    transport.close()
    return retval


if __name__ == "__main__":
    argParser = argparse.ArgumentParser()
    argParser.add_argument('--species', '-s', default='human', type=str,
                           help='species name')
    argParser.add_argument('--genes', '-g', default=None, type=str,
                           help='list of genes to query')
    argParser.add_argument('--ranks', '-r', default=None, type=str,
                           help='list of gene ranks')
    argParser.add_argument('--scores', '-v', default=None, type=str,
                           help='list of gene scores')
    argParser.add_argument('--useRank', default=False, action='store_true',
                           help='get rank based pvalues')
    argParser.add_argument('--port', '-p', default=9090, type=int,
                           help='server port')
    args = argParser.parse_args()


    if args.genes is not None:
        args.genes = args.genes.split(',')

    if args.ranks is not None:
        args.ranks = args.ranks.split(',')
        args.ranks = [int(rank) for rank in args.ranks]

    if args.scores is not None:
        args.scores = args.scores.split(',')
        args.scores = [float(score) for score in args.scores]

    print(f'Genes: {args.genes}')
    if args.ranks is not None:
        print(f'Ranks: {args.ranks}')
    if args.scores is not None:
        print(f'Scores: {args.scores}')
    retval = runQuery(args)
    exit(retval)
