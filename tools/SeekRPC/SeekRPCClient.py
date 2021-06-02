import time
import argparse
from thrift.transport import TTransport, TSocket
from thrift.protocol.TBinaryProtocol import TBinaryProtocol

from pyseek import SeekRPC, constants
from pyseek.ttypes import SeekQuery, QueryParams, QueryResult, SearchMethod, DistanceMeasure

host = 'localhost'

'''
Query parameters and options:

struct QueryParams {
    string search_method: SearchMethods.(default: CV, EqualWeighting, OrderStatistics, CVCustom)
    string distance_measure: DistanceMeasure.(default: ZScore, ZScoreHubbinessCorrected, Correlation)
    double min_query_genes_fraction = 0.0;
    double min_genome_fraction = 0.0;
    double rbp_param = 0.99;
    bool useNegativeCorrelation = false;
    bool check_dataset_size = false;
    bool use_gene_symbols = false;
    bool simulate_weights = false;
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
    version = client.get_rpc_version()
    assert version == constants.RPC_Version

    params = QueryParams(search_method=SearchMethod.CV,
                         distance_measure=DistanceMeasure.ZScoreHubbinessCorrected,
                         min_query_genes_fraction=0.5,
                         min_genome_fraction=0.5,
                         use_gene_symbols=args.useSymbols,
                         simulate_weights=False)

    genes = [gene.upper() for gene in args.genes]
    query = SeekQuery(species=args.species, genes=args.genes, parameters=params)

    retval = -1
    taskId = client.seek_query_async(query)
    result = client.seek_get_result(taskId, block=True)
    print(f'Status: {result.statusMsg}')
    # Alternate non-blocking code, checking periodically if complete
    # taskId = client.seek_query_async(query)
    # while client.is_query_complete(taskId) is False:
    #     # get status messages
    #     statusMsg = client.get_progress_message(taskId)
    #     if len(statusMsg) > 0:
    #         print(f'Status: {statusMsg}')
    #     time.sleep(0.1)
    # result = client.seek_get_result(taskId, block=True)
    if result.success is True:
        for i, gs in enumerate(result.gene_scores):
            print(f'gene: {gs.name}, {gs.value}')
            if i > 10: break

        for i, ds in enumerate(result.dataset_weights):
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
