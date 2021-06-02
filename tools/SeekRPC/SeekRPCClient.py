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

struct QueryParams {
    string searchMethod: SearchMethods.(default: CV, EqualWeighting, OrderStatistics, CVCustom)
    string distanceMeasure: DistanceMeasure.(default: ZScore, ZScoreHubbinessCorrected, Correlation)
    double minQueryGenesFraction = 0.0;
    double minGenomeFraction = 0.0;
    double rbpParam = 0.99;
    bool useNegativeCorrelation = false;
    bool checkDatasetSize = false;
    bool useGeneSymbols = false;
    bool simulateWeights = false;
}

struct SeekQuery {
    string species: (default: "Unknown", "human", "fly", "mouse", "worm", "yeast", "zebrafish")
    list<string> genes;
    list<string> datasets;
    QueryParams parameters;
    list<string> guideGenes;
    string outputDir = "/tmp/seek";
}
'''

def runQuery(args):
    socket = TSocket.TSocket(host, args.port)
    transport = TTransport.TBufferedTransport(socket)
    protocol = TBinaryProtocol(transport)
    client = SeekRPC.Client(protocol)
    transport.open()
    version = client.getRpcVersion()
    assert version == constants.RPCVersion

    params = SeekRPC.QueryParams(
        searchMethod = ttypes.SearchMethod.CV,
        distanceMeasure = ttypes.DistanceMeasure.ZScoreHubbinessCorrected,
        minQueryGenesFraction = 0.5,
        minGenomeFraction = 0.5,
        useGeneSymbols = args.useSymbols,
        simulateWeights = False,
    )

    genes = [gene.upper() for gene in args.genes]
    query = SeekRPC.SeekQuery(species=args.species, genes=args.genes, parameters=params)

    retval = -1
    taskId = client.seekQueryAsync(query)
    result = client.getQueryResult(taskId, block=True)
    print(f'Status: {result.statusMsg}')
    # Alternate non-blocking code, checking periodically if complete
    # taskId = client.seekQueryAsync(query)
    # while client.isQueryComplete(taskId) is False:
    #     # get status messages
    #     statusMsg = client.getProgressMessage(taskId)
    #     if len(statusMsg) > 0:
    #         print(f'Status: {statusMsg}')
    #     time.sleep(0.1)
    # result = client.getQueryResult(taskId, block=True)
    if result.success is True:
        for i, gs in enumerate(result.geneScores):
            print(f'gene: {gs.name}, {gs.value}')
            if i > 10: break

        for i, ds in enumerate(result.datasetWeights):
            print(f'dset: {ds.name}, {ds.value}')
            if i > 10: break
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
    argParser.add_argument('--useSymbols', '-S', default=False, action='store_true',
                           help='genes specified by symbol name')
    argParser.add_argument('--port', '-p', default=9090, type=int,
                           help='server port')
    args = argParser.parse_args()

    if args.genes is None:
        print('No query genes specified, use [-g gene1,gene2,...] option')
        exit(-1)
    args.genes = args.genes.split(',')

    retval = runQuery(args)
    exit(retval)


#species = 'human'
#genes = ['SMO', 'PTCH1', 'PTCH2', 'BOC'] - this one differs from web, ambiguous symbol->entrez?
#genes = ['90634', '23659']

#species = 'fly'
#genes = ['35234', '35232']
#genes = ['34930', '35234', '35232']

#species = 'yeast'
#genes = ['YGL142C', 'YHR188C']

# species = 'mock'
# genes = ['six', 'four']
# genes = ['90634', '23659']
