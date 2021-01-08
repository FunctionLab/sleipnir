#include <stdio.h>
#include <iostream>
#include <getopt.h>
#include <sstream>

#include <thrift/protocol/TBinaryProtocol.h>
#include <thrift/transport/TSocket.h>
#include <thrift/transport/TTransportUtils.h>

#include "seekerror.h"
#include "gen-cpp/SeekRPC.h"

using namespace std;
using namespace SeekRPC;
using namespace apache::thrift;
using namespace apache::thrift::protocol;
using namespace apache::thrift::transport;

struct Args
{
    string species;
    vector<string> genes;
    bool useSymbols = false;
};

bool parseArgs(int argc, char **argv, Args &args)
{
    string usage = "Usage: SeekRPCClient -s <species> -g <comma-sep-query-genes> [-S]\n"
                   "Run a seek query for the given species and genes.\n";

    static struct option long_options[] = {
        {"species", required_argument, 0, 's'},
        {"genes", required_argument, 0, 'g'},
        {"useSymbols", no_argument, 0, 'S'},
        {0, 0, 0, 0}};
    int opt = 0;
    int long_index = 0;
    while ((opt = getopt_long(argc, argv, "s:g:S",
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
        case 'S':
            args.useSymbols = true;
            break;
        default:
            cerr << "Error: unrecognized options: " << opt << endl;
            cerr << usage << endl;
            return false;
        }
    }
    if (args.genes.size() == 0)
    {
        cerr << "ERROR: No query genes provided" << endl;
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
    std::cout << args.species << endl;
    for (auto gene : args.genes) {
        std::cout << gene << ",";
    }
    std::cout << endl;

    shared_ptr<TTransport> socket(new TSocket("localhost", 9090));
    shared_ptr<TTransport> transport(new TBufferedTransport(socket));
    shared_ptr<TProtocol> protocol(new TBinaryProtocol(transport));
    SeekRPCClient seekClient(protocol);

    try
    {
        transport->open();
        seekClient.ping();

        SeekQuery query;

        query.species = args.species;
        query.__isset.parameters = true;
        query.parameters.__set_search_method("CV");
        query.parameters.__set_distance_measure("ZscoreHubbinessCorrected");
        query.parameters.__set_min_genome_fraction(0.5);
        query.parameters.__set_min_query_genes_fraction(0.5);
        query.parameters.__set_use_gene_symbols(args.useSymbols);

        for (string gene : args.genes) {
            query.genes.push_back(gene);
        }

        QueryResult result;
        seekClient.seek_query(result, query);

        if (result.success == false) {
            throw query_error(result.statusMsg);
        }

        std::cout << "Gene Scores:" << endl;
        // for (auto gene_score : result.gene_scores) {
        for (int i=0; i<result.gene_scores.size() && i < 100; i++) {
            std::cout << result.gene_scores[i] << endl;
        }

        std::cout << "Dataset Weights:" << endl;
        // for (auto dset : result.dataset_weights)
        for (int i=0; i<result.dataset_weights.size() && i < 100; i++) {
            std::cout << result.dataset_weights[i] << endl;
        }
    } catch (TException &tx) {
        std::cout << "ERROR: " << tx.what() << endl;
    } catch (exception &err) {
        std::cout << "ERROR: " << err.what() << endl;
    }

    transport->close();
}