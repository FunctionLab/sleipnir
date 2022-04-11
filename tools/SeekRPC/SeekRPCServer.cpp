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
//#include <thrift/server/TSimpleServer.h>
#include <thrift/server/TThreadedServer.h>
#include <thrift/transport/TServerSocket.h>
#include <thrift/transport/TBufferTransports.h>
#include "seekerror.h"
#include "stdafx.h"

using namespace std;
using namespace SeekRPC;
using namespace ::apache::thrift;
using namespace ::apache::thrift::protocol;
using namespace ::apache::thrift::transport;
using namespace ::apache::thrift::server;

struct Args {
    vector<string> configFiles;
    uint32_t port = 9090;
    uint32_t taskTimeoutSec = 5 * 60;  // default 5 minutes in seconds
    uint32_t maxThreads = 0;  // default will base on hardware
    bool verbose = false;
};

bool parseArgs(int argc, char **argv, Args &args)
{
    // parse options
    string usage = "Usage: seekService -c <species_1.toml> [-c <species_2.toml> ...]\n"
                   "Starts the SeekService to handle RPC requests. " 
                   "Initializes the species with provided config files.\n";
    static struct option long_options[] = {
        {"config", required_argument, 0, 'c'},
        {"port", required_argument, 0, 'p'},
        {"timeout", required_argument, 0, 't'},
        {"maxThreads", required_argument, 0, 'm'},
        {"verbose", no_argument, 0, 'v'},
        {0, 0, 0, 0}};
    int opt = 0;
    int long_index = 0;
    while ((opt = getopt_long(argc, argv, "c:p:t:m:v",
                              long_options, &long_index)) != -1)
    {
        switch (opt)
        {
        case 'c':
            args.configFiles.push_back(optarg);
            break;
        case 'p':
            args.port = atoi(optarg);
            break;
        case 't':
            args.taskTimeoutSec = atoi(optarg);
            break;
        case 'm':
            args.maxThreads = atoi(optarg);
            break;
        case 'v':
            args.verbose = true;
            break;
        default:
            cout << "Error: unrecognized options: " << opt << endl;
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

    int logLevel = Priority::INFO;
    if (args.verbose) {
        logLevel = Priority::DEBUG;
    }
    CMeta Meta(logLevel);

    try {
        if (args.maxThreads == 0) {
            args.maxThreads = thread::hardware_concurrency();
        }
        g_CatSleipnir().info("SeekRPC Setting: max concurrency: %d", args.maxThreads);
        g_CatSleipnir().info("SeekRPC Setting: max thread time(sec): %d", args.taskTimeoutSec);
        SeekInterface seekInterface(args.configFiles, args.maxThreads, args.taskTimeoutSec);

        bool startServer = true;
        if (startServer == true) {
            shared_ptr<SeekRPCHandler> handler(new SeekRPCHandler(seekInterface));
            shared_ptr<TProcessor> processor(new SeekRPCProcessor(handler));
            shared_ptr<TServerTransport> serverTransport(new TServerSocket(args.port));
            shared_ptr<TTransportFactory> transportFactory(new TBufferedTransportFactory());
            shared_ptr<TProtocolFactory> protocolFactory(new TBinaryProtocolFactory());

            /* Alternate server types
               ## For single-threaded server
               TSimpleServer server(processor, serverTransport, transportFactory, protocolFactory);
               ## For multi-threaded with thread pool and reusing threads
               TThreadPoolServer server(processor, serverTransport, transportFactory, protocolFactory, threadManager);
            */
            // For multi-threaded server
            TThreadedServer server(processor, serverTransport, transportFactory, protocolFactory);
            g_CatSleipnir().info("SeekRPC: Listen on port: %d", args.port);
            server.serve();
        }
        printf("Exiting Server ...\n");
    } catch(exception &err) {
        g_CatSleipnir().error("SeekRPC: Fatal Server Exception");
        cerr << print_exception_stack(err) << endl;
        return -1;
    } catch(...) {
        cout << "SeekRPC: Uncaught Exception\n";
    }
}
