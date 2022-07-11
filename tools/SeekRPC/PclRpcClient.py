import sys
import time
import re
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

def writeResults(fp, args, result):
    if args.zexp is True:
        fp.write('==Gene Expressions==\n')
        fp.write("\n".join(str(round(item, 6)) for item in result.geneExpressions))
        fp.write("\n")
    if args.zcoexp is True:
        fp.write('==Gene Coexpressions==\n')
        fp.write("\n".join(str(round(item, 6)) for item in result.geneCoexpressions))
        fp.write("\n")
    if args.zqexp is True:
        fp.write('==Query Expressions==\n')
        fp.write("\n".join(str(round(item, 6)) for item in result.queryExpressions))
        fp.write("\n")
    if args.zqcoexp is True:
        fp.write('==Query Coexpressions==\n')
        fp.write("\n".join(str(round(item, 6)) for item in result.queryCoexpressions))
        fp.write("\n")


def runQuery(args):
    socket = TSocket.TSocket(host, args.port)
    transport = TTransport.TBufferedTransport(socket)
    protocol = TBinaryProtocol(transport)
    client = SeekRPC.Client(protocol)
    transport.open()
    version = client.getRpcVersion()
    assert version == constants.RPCVersion
    retval = -1

    settings = SeekRPC.PclSettings(
        outputNormalized = False,
        outputGeneExpression = args.zexp,
        outputGeneCoexpression = args.zcoexp,
        outputQueryExpression = args.zqexp,
        outputQueryCoexpression = args.zqcoexp,
        rbp = -1,
    )

    pclArgs = SeekRPC.PclQueryArgs(species=args.species,
                                  genes=args.genes,
                                  queryGenes=args.queryGenes,
                                  datasets=args.datasets,
                                  settings=settings)
    result = client.pclQuery(pclArgs)
    if result.success is True:
        print(f'result: {result}')
        # Output to file if requested
        if args.outfile is not None:
            with open(args.outfile, "w") as fp:
                writeResults(fp, args, result)

        # Print to command line
        writeResults(sys.stdout, args, result)

        # Example of accessing expressions in a gene ordered way
        if result.geneExpressions and len(result.geneExpressions) > 0:
            print('## Ordered Printout ##')
            assert len(result.datasetSizes) == len(args.datasets)
            geneRowOffset = 0
            expNameOffset = 0
            for i, numSamples in enumerate(result.datasetSizes):
                # strip any .pcl extension from the dataset name
                dsetName = re.sub('\.pcl.*', '', args.datasets[i])
                print(f'DATASET: {dsetName}')
                # Print out the experiment (sample) names for each dataset
                for j in range(numSamples):
                    print(f'{result.experimentNames[expNameOffset + j]}  ', end='')
                print()
                expNameOffset += numSamples
                for j in range(len(args.genes)):
                    print(f'GENE({args.genes[j]}): ', end='')
                    for k in range(numSamples):
                        print(f'{result.geneExpressions[geneRowOffset + k]:.3f}  ', end='')
                    print()
                    geneRowOffset += numSamples
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
    argParser.add_argument('--queryGenes', '-q', default=None, type=str,
                           help='list of qurey genes')
    argParser.add_argument('--datasets', '-d', default=None, type=str, required=True,
                           help='list of datasets to query')
    argParser.add_argument('--zexp', default=False, action='store_true',
                           help='output gene expression values')
    argParser.add_argument('--zcoexp', default=False, action='store_true',
                           help='output gene Coexpression values')
    argParser.add_argument('--zqexp', default=False, action='store_true',
                           help='output query expression values')
    argParser.add_argument('--zqcoexp', default=False, action='store_true',
                           help='output query Coexpression values')
    argParser.add_argument('--outfile', '-o', default=None, type=str,
                           help='filename to output results')
    argParser.add_argument('--port', '-p', default=9090, type=int,
                           help='server port')
    args = argParser.parse_args()


    if args.genes is not None:
        args.genes = args.genes.split(',')

    if args.queryGenes is not None:
        args.queryGenes = args.queryGenes.split(',')

    if args.datasets is not None:
        args.datasets = args.datasets.split(',')

    if args.genes is None and args.queryGenes is None:
        print('Please specify genes and/or queryGenes, use [-g gene1,gene2,...] [-q gene1,gene2,...]')
        exit(-1)

    if not any([args.zexp, args.zcoexp, args.zqexp, args.zqcoexp]):
        # no output type set, set gene expression by default
        print("No output type specified, default to gene expressions")
        args.zexp = True

    if any([args.zexp, args.zcoexp]) and args.genes is None:
        print("Must specify list of genes for gene expressions")
        exit(-1)

    if any([args.zqexp, args.zqcoexp]) and args.queryGenes is None:
        print("Must specify list of queryGenes for query expressions")
        exit(-1)

    retval = runQuery(args)
    exit(retval)
