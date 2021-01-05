#include <stdio.h>
#include <iostream>

#include <thrift/protocol/TBinaryProtocol.h>
#include <thrift/transport/TSocket.h>
#include <thrift/transport/TTransportUtils.h>

#include "gen-cpp/SeekRPC.h"

using namespace std;
using namespace SeekRPC;
using namespace apache::thrift;
using namespace apache::thrift::protocol;
using namespace apache::thrift::transport;

int main(int argc, char **argv)
{
  shared_ptr<TTransport> socket(new TSocket("localhost", 9090));
  shared_ptr<TTransport> transport(new TBufferedTransport(socket));
  shared_ptr<TProtocol> protocol(new TBinaryProtocol(transport));
  SeekRPCClient seekClient(protocol);

  try {
    transport->open();
    seekClient.ping();

    SeekQuery query;
    // QueryParams params;
    // query.parameters.search_method = "CV";
    // query.parameters.__set_distance_measure("z-score");
    // query.parameters.min_genome_fraction = 0.2;
    // query.parameters.useNegativeCorrelation = true;
    // query.genes.push_back("gene1");
    // query.datasets.push_back("dataset1");
    // query.outputDir = "/tmp";
    // query.parameters = params;
    // query.__set_parameters(params);
    // query.parameters.distance_measure = "z-score";

    query.species = "mock";
    query.genes.push_back("90634");
    query.genes.push_back("23659");

    QueryResult result;
    seekClient.seek_query(result, query);

    cout << "Gene Scores:" << endl;
    for (auto gene_score: result.gene_scores)
      cout << gene_score << endl;
    
    // for (auto score: result.gene_scores)
    //   cout << score << endl;

    cout << "Dataset Weights:" << endl;
    for (auto dset: result.dataset_weights)
      cout << dset << endl;

  } catch (TException& tx) {
    cout << "ERROR: " << tx.what() << endl;
  }

  transport->close();
}