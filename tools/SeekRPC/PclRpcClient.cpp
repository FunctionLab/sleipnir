#include <stdio.h>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <sstream>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <thrift/protocol/TBinaryProtocol.h>
#include <thrift/transport/TSocket.h>
#include <thrift/transport/TTransportUtils.h>

#include "seekerror.h"
#include "gen-cpp/SeekRPC.h"
#include "gen-cpp/seek_rpc_constants.h"

using namespace std;
using namespace SeekRPC;
using namespace apache::thrift;
using namespace apache::thrift::protocol;
using namespace apache::thrift::transport;

struct Args
{
    string species;
    vector<string> genes;
    vector<string> datasets;
    string outputFile;
    uint32_t port = 9090;
    bool verbose = false;
};

bool parseArgs(int argc, char **argv, Args &args)
{
    string usage = "Usage: PclRpcClient -s <species> -g <comma-sep-query-genes> -d <comma-sep-datasets> \n"
                   "Run a PCL query for the given species, genes and datasets.\n";

    static struct option long_options[] = {
        {"species", required_argument, 0, 's'},
        {"genes", required_argument, 0, 'g'},
        {"datasets", required_argument, 0, 'd'},
        {"output", required_argument, 0, 'o'},
        {"port", required_argument, 0, 'p'},
        {"verbose", no_argument, 0, 'v'},
        {0, 0, 0, 0}};
    int opt = 0;
    int long_index = 0;
    while ((opt = getopt_long(argc, argv, "s:g:d:o:p:v",
                              long_options, &long_index)) != -1)
    {
        switch (opt)
        {
        case 's':
            args.species = optarg;
            break;
        case 'g':
        {
            stringstream s_stream(optarg);
            while (s_stream.good())
            {
                string gene;
                getline(s_stream, gene, ',');
                args.genes.push_back(gene);
            }
        }
        break;
        case 'd':
        {
            stringstream s_stream(optarg);
            while (s_stream.good())
            {
                string dataset;
                getline(s_stream, dataset, ',');
                args.datasets.push_back(dataset);
            }
        }
        break;
        case 'o':
            args.outputFile = optarg;
            break;
        case 'p':
            args.port = atoi(optarg);
            break;
        case 'v':
            args.verbose = true;
            break;
        default:
            cerr << "Error: unrecognized options: " << opt << endl;
            cerr << usage << endl;
            return false;
        }
    }
    if (args.genes.size() == 0 || args.datasets.size() == 0)
    {
        cerr << "ERROR: Must specify both query genes and datasets" << endl;
        cerr << usage << endl;
        return false;
    }

    return true;
}

int main(int argc, char **argv)
{
    Args args;
    bool res = parseArgs(argc, argv, args);
    if (res == false)
    {
        exit(-1);
    }

    if (args.verbose == true) {
        std::cout << args.species << endl;
        for (auto gene : args.genes) {
            std::cout << gene << ",";
        }
        std::cout << endl;
    }

    shared_ptr<TTransport> socket(new TSocket("localhost", args.port));
    shared_ptr<TTransport> transport(new TBufferedTransport(socket));
    shared_ptr<TProtocol> protocol(new TBinaryProtocol(transport));
    SeekRPCClient seekClient(protocol);

    try
    {
        transport->open();
        seekClient.ping();

        // // Check for version compatibility
        double version = seekClient.getRpcVersion();
        assert(version == g_seek_rpc_constants.RPCVersion);

        PclQueryArgs pclQueryArgs;
        PclResult result;
        pclQueryArgs.species = args.species;
        pclQueryArgs.__isset.settings = true;
        pclQueryArgs.settings.__set_outputGeneExpression(true);
        pclQueryArgs.settings.__set_outputNormalized(true);
        pclQueryArgs.settings.__set_rbp(-1);
        // query.settings.__set_outputGeneCoexpression(true);
        // query.settings.__set_outputQueryExpression(true);
        // query.settings.__set_outputQueryCoexpression(true);
        pclQueryArgs.__set_genes(args.genes);
        pclQueryArgs.__set_datasets(args.datasets);

        seekClient.pclQuery(result, pclQueryArgs);
        if (result.success == false) {
            throw query_error(result.statusMsg);
        }
        assert(result.geneExpressions.size() > 0);

        int numDatasets = args.datasets.size();
        int numGenes = args.genes.size();
        int offset = 0;
        for (int i=0; i<numDatasets; i++) {
            for (int j=0; j<numGenes; j++) {
                cout << "dset: " << i << ", gene: " << j << endl;
                int numSamples = result.datasetSizes[i];
                for(int k=0; k<numSamples; k++) {
                    // Can't use offset calculation because each dataseSize is different
                    // int offset = k + j*datasetSize + i*datasetSize*numGenes; 
                    cout << "val: " << result.geneExpressions[offset] << endl;
                    offset++;
                }
            }
        }
        if (!args.outputFile.empty()) {
            // output file is specified
            using boost::algorithm::join;
            using boost::adaptors::transformed;
            std::ofstream outFile;
            outFile.open(args.outputFile, ofstream::out | ofstream::trunc);
            auto tostr = static_cast<std::string(*)(double)>(std::to_string);
            std::string joinedVals =  join(result.geneExpressions | transformed(tostr), "\n");
            outFile << joinedVals;
            outFile.close();

        }

    } catch (TException &tx) {
        std::cout << "ERROR: " << tx.what() << endl;
        exit(-1);
    } catch (exception &err) {
        std::cout << "ERROR: " << err.what() << endl;
        exit(-1);
    }

    transport->close();
}