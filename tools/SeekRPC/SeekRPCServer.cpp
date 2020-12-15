#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <vector>
#include <map>
#include <cassert>
#include "seekhelper.h"
#include "seekcentral.h"
#include "SeekInterface.h"
#include "SeekRPCHandler.h"
#include <thrift/protocol/TBinaryProtocol.h>
#include <thrift/server/TSimpleServer.h>
#include <thrift/transport/TServerSocket.h>
#include <thrift/transport/TBufferTransports.h>

using namespace std;
using namespace SeekRPC;
using namespace ::apache::thrift;
using namespace ::apache::thrift::protocol;
using namespace ::apache::thrift::transport;
using namespace ::apache::thrift::server;

struct Args {
    vector<string> configFiles;
};

bool parseArgs(int argc, char **argv, Args &args)
{
    string usage = "seekService -c <species_1.toml> [-c <species_2.toml> ...]\n"
                   "Starts the SeekService to handle RPC requests. " 
                   "Initializes the species with provided config files.\n";
    static struct option long_options[] = {
        {"config", required_argument, 0, 'c'},
        {0, 0, 0, 0}};
    int opt = 0;
    int long_index = 0;
    while ((opt = getopt_long(argc, argv, "c:",
                              long_options, &long_index)) != -1)
    {
        switch (opt)
        {
        case 'c':
            args.configFiles.push_back(optarg);
            break;
        default:
            cout << "Error: unrecognized options" << endl;
            cout << usage << endl;
            return false;
        }
    }
    if (args.configFiles.size() == 0)
    {
        cout << "ERROR: No species config files provided" << endl;
        cout << usage << endl;
        return false;
    }
    return true;
}


int main(int argc, char** argv) 
{
    Args args;
    bool res = parseArgs(argc, argv, args);
    if (res == false) {
        exit(-1);
    }

    map<string, SeekSettings> speciesConfigs;
    getConfigs(args.configFiles, speciesConfigs);

    map<string, CSeekCentral> speciesSeekCentrals;
    for( auto const& [speciesName, config] : speciesConfigs ) {
        bool res;
        // speciesSeekCentrals[speciesName] = CSeekCentral(); // TODO - do we need this line?
        res = speciesSeekCentrals[speciesName].InitializeFromSeekConfig(config);
        if (res == false) {
            cerr << "Error initializing CSeekCentral for species " << speciesName << endl;
            return -1;
        }
    }

    bool startServer = true;
    if (startServer == true) {
        int port = 9090;
        SeekInterface seekInterface;
        shared_ptr<SeekRPCHandler> handler(new SeekRPCHandler(seekInterface));
        shared_ptr<TProcessor> processor(new SeekRPCProcessor(handler));
        shared_ptr<TServerTransport> serverTransport(new TServerSocket(port));
        shared_ptr<TTransportFactory> transportFactory(new TBufferedTransportFactory());
        shared_ptr<TProtocolFactory> protocolFactory(new TBinaryProtocolFactory());

        TSimpleServer server(processor, serverTransport, transportFactory, protocolFactory);
        server.serve();
    }
}
