import argparse
from thrift.transport import TTransport, TSocket
from thrift.protocol.TBinaryProtocol import TBinaryProtocol

from pyseek import SeekRPC
from pyseek.ttypes import SeekQuery, QueryParams, QueryResult

host = 'localhost'
port = 9090

'''
Query parameters and options:

struct QueryParams {
    string search_method: (default: "CV", "EqualWeighting", "OrderStatistics", "CVCUSTOM")
    string distance_measure: (default: "Zscore", "ZscoreHubbinessCorrected", "Correlation")
    double min_query_genes_fraction = 0.0;
    double min_genome_fraction = 0.0;
    double rbp_param = 0.99;
    bool useNegativeCorrelation = false;
    bool check_dataset_size = false;
    bool use_gene_symbols = false;
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
    socket = TSocket.TSocket(host, port)
    transport = TTransport.TBufferedTransport(socket)
    protocol = TBinaryProtocol(transport)
    client = SeekRPC.Client(protocol)
    transport.open()

    params = QueryParams(distance_measure="ZscoreHubbinessCorrected",
                         min_query_genes_fraction=0.5,
                         min_genome_fraction=0.5,
                         use_gene_symbols=args.useSymbols)

    query = SeekQuery(species=args.species, genes=args.genes, parameters=params)

    result = client.seek_query(query)
    if result.success is True:
        for i, gs in enumerate(result.gene_scores):
            print(f'gene: {gs.name}, {gs.value}')
            if i > 100: break

        for i, ds in enumerate(result.dataset_weights):
            print(f'dset: {ds.name}, {ds.value}')
            if i > 100: break
    else:
        print(f'query error: {result.statusMsg}')

    transport.close()


if __name__ == "__main__":
    argParser = argparse.ArgumentParser()
    argParser.add_argument('--species', '-s', default='human', type=str,
                           help='species name')
    argParser.add_argument('--genes', '-g', default=None, type=str,
                           help='list of genes to query')
    argParser.add_argument('--useSymbols', '-S', default=False, action='store_true',
                           help='genes specified by symbol name')
    args = argParser.parse_args()

    if args.genes is None:
        print('No query genes specified, use [-g gene1,gene2,...] option')
        exit(-1)
    args.genes = args.genes.split(',')

    runQuery(args)
    exit(0)


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
