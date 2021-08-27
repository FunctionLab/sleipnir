import time
import argparse
import importlib
from thrift.transport import TTransport, TSocket
from thrift.protocol.TBinaryProtocol import TBinaryProtocol
# Since the path gen-py/seek_rpc has a hyphen in it we need to use
#   the import_module() function rather than the usual import call
PclRPC = importlib.import_module('gen-py.pcl_rpc.PclRPC')
ttypes = importlib.import_module('gen-py.pcl_rpc.ttypes')
constants = importlib.import_module('gen-py.pcl_rpc.constants')

host = 'localhost'

'''
Query parameters and options:
struct PclSettings {
    1: optional bool outputNormalized = false;
    2: optional bool outputGeneExpression = false;
    3: optional bool outputQueryExpression = false;
    4: optional bool outputGeneCoexpression = false;
    5: optional bool outputQueryCoexpression = false;
    6: optional double rbp = 0.99;
}

struct PclQueryArgs {
    1: required string species = "Unknown";
    2: required list<string> genes;
    3: required list<string> datasets;
    4: optional list<string> query;
    5: optional PclSettings settings;
    6: optional string outputDir = "/tmp/seek";
}

struct PclResult {
    1: required bool success;
    2: required list<i32> datasetSizes;
    3: required list<double> geneExpressions;
    4: optional list<double> geneCoexpressions;
    5: optional list<double> queryExpressions;
    6: optional list<double> queryCoexpressions;
    7: optional PclStatus status;
    8: optional string statusMsg;
}
'''

def runQuery(args):
    socket = TSocket.TSocket(host, args.port)
    transport = TTransport.TBufferedTransport(socket)
    protocol = TBinaryProtocol(transport)
    client = PclRPC.Client(protocol)
    transport.open()
    version = client.getRpcVersion()
    assert version == constants.PclRPCVersion

    settings = PclRPC.PclSettings(
        outputNormalized = True,
        outputGeneExpression = True,
        outputQueryExpression = False,
        outputGeneCoexpression = False,
        outputQueryCoexpression = False,
        rbp = -1,
    )

    pclArgs = PclRPC.PclQueryArgs(species=args.species, genes=args.genes, datasets=args.datasets, settings=settings)

    retval = -1
    result = client.pclQuery(pclArgs)
    if result.success is True:
        # Output to file if requested
        if args.outfile is not None:
            with open(args.outfile, "w") as fp:
                fp.write("\n".join(str(round(item, 6)) for item in result.geneExpressions))
        # Print to command line
        numGenes = len(args.genes)
        offset = 0;
        for i, numSamples in enumerate(result.datasetSizes):
            for j in range(numGenes):
                print(f'dset: {i}, gene: {j}')
                for k in range(numSamples):
                    print(f'val: {result.geneExpressions[offset]}')
                    offset += 1
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
    argParser.add_argument('--datasets', '-d', default=None, type=str,
                           help='list of datasets to query')
    argParser.add_argument('--outfile', '-o', default=None, type=str,
                           help='filename to output results')
    argParser.add_argument('--port', '-p', default=9010, type=int,
                           help='server port')
    args = argParser.parse_args()

    if args.genes is None or args.datasets is None:
        print('Please specify both genes and datasets, use [-g gene1,gene2,...] [-d dset1,dset2,...]')
        exit(-1)
    args.genes = args.genes.split(',')
    args.datasets = args.datasets.split(',')

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
