#include <stdio.h>
#include <iostream>
#include <getopt.h>
#include <sstream>

#include <thrift/protocol/TBinaryProtocol.h>
#include <thrift/transport/TSocket.h>
#include <thrift/transport/TTransportUtils.h>

#include "seekreader.h"
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
    string queryFile;
    string outputFile;
    vector<string> genes;
    uint32_t port = 9090;
    uint32_t maxResults = 2000;
    bool useSymbols = false;
    bool verbose = false;
};

bool parseArgs(int argc, char **argv, Args &args)
{
    string usage = "Usage: SeekRPCClient -s <species> -g <comma-sep-query-genes> "
                   "-q <query_file> -o <output_file> -n <max_results> [-S]\n"
                   "Run a seek query for the given species and genes.\n";

    static struct option long_options[] = {
        {"species", required_argument, 0, 's'},
        {"genes", required_argument, 0, 'g'},
        {"port", required_argument, 0, 'p'},
        {"query", required_argument, 0, 'q'},
        {"output", required_argument, 0, 'o'},
        {"max_results", required_argument, 0, 'n'},
        {"useSymbols", no_argument, 0, 'S'},
        {"verbose", no_argument, 0, 'v'},
        {0, 0, 0, 0}};
    int opt = 0;
    int long_index = 0;
    while ((opt = getopt_long(argc, argv, "s:g:o:p:q:Sv",
                              long_options, &long_index)) != -1)
    {
        switch (opt)
        {
        case 's':
            args.species = optarg;
            break;
        case 'q':
            args.queryFile = optarg;
            break;
        case 'o':
            args.outputFile = optarg;
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
        case 'p':
            args.port = atoi(optarg);
            break;
        case 'n':
            args.maxResults = atoi(optarg);
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
    if (args.genes.size() == 0 && args.queryFile.empty())
    {
        cerr << "ERROR: No query genes or query file provided" << endl;
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

    // A list of queries, each query a list of genes
    vector <vector<string>> queries;

    if (!args.queryFile.empty()) {
        // Load list of queries from a file
        if (!Sleipnir::CSeekTools::ReadMultipleQueries(args.queryFile, queries)) {
            throw init_error("Unable to load query file " + args.queryFile);
        }
    }

    if (args.genes.size() != 0) {
        // Load the single query from the command-line genes.
        vector<string> queryGenes;
        for (string gene : args.genes) {
            queryGenes.push_back(gene);
        }
        queries.push_back(queryGenes);
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
        std::ofstream outFile;
        transport->open();
        seekClient.ping();



        if (!args.outputFile.empty()) {
            // output file is specified
            outFile.open(args.outputFile, ofstream::out | ofstream::trunc);
        }

        for (int i=0; i<queries.size(); i++) {
            SeekQuery query;
            QueryResult result;
            query.species = args.species;
            query.__isset.parameters = true;
            query.parameters.__set_search_method("CV");
            query.parameters.__set_distance_measure("ZscoreHubbinessCorrected");
            // query.parameters.__set_min_genome_fraction(0.5);
            query.parameters.__set_min_genome_fraction(0);
            // query.parameters.__set_min_query_genes_fraction(0.5);
            query.parameters.__set_min_query_genes_fraction(0);
            query.parameters.__set_use_gene_symbols(args.useSymbols);

            query.genes = queries[i];
            seekClient.seek_query(result, query);

            if (result.success == false) {
                throw query_error(result.statusMsg);
            }

            if (args.verbose == true || args.outputFile.empty()) {
                std::cout << "Query " << i << ": Gene Scores:" << endl;
                // for (auto gene_score : result.gene_scores) {
                for (int i=0; i<result.gene_scores.size() && i < 10; i++) {
                    std::cout << result.gene_scores[i] << endl;
                }

                std::cout << "Query " << i << ": Dataset Weights:" << endl;
                // for (auto dset : result.dataset_weights)
                for (int i=0; i<result.dataset_weights.size() && i < 10; i++) {
                    std::cout << result.dataset_weights[i] << endl;
                }
            }

            if (outFile.is_open()) {
                // print results to a file
                for (int i=0; i<result.dataset_weights.size() && i < args.maxResults; i++) {
                    if (i > 0) outFile << " ";
                    outFile << result.dataset_weights[i].name;
                }
                outFile << endl;
                for (int i=0; i<result.gene_scores.size() && i < args.maxResults; i++) {
                    if (i > 0) outFile << " ";
                    outFile << result.gene_scores[i].name;
                }
                outFile << endl;
            }
        }
        if (outFile.is_open()) {
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
