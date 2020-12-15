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
    query.genes.push_back("gene1");
    query.datasets.push_back("dataset1");
    query.outputDir = "/tmp";
    query.parameters.distance_measure = "z-score";

    QueryResult result;
    seekClient.seek_query(result, query);

    for (auto gene: result.genes)
      cout << gene << endl;

  } catch (TException& tx) {
    cout << "ERROR: " << tx.what() << endl;
  }

  transport->close();
}