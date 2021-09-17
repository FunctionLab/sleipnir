#include "stdafx.h"
#include "cmdline.h"
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>
#include <list>
#include <cassert>
#include "seekerror.h"
#include "seekhelper.h"
#include "seekcentral.h"

#include "PclQuery.h"

#define DEFAULT_CACHE_SIZE 100
#define NUM_THREADS 16
char THREAD_OCCUPIED[NUM_THREADS];
pthread_mutex_t mutexGet;

#define BACKLOG 20   // how many pending connections queue will hold
char const *PORT = "9000";

void sigchld_handler(int s) {
    while (waitpid(-1, NULL, WNOHANG) > 0);
}

// get sockaddr, IPv4 or IPv6:
void *get_in_addr(struct sockaddr *sa) {
    if (sa->sa_family == AF_INET) {
        return &(((struct sockaddr_in *) sa)->sin_addr);
    }
    return &(((struct sockaddr_in6 *) sa)->sin6_addr);
}


int main(int iArgs, char **aszArgs) {
    static const size_t c_iBuffer = 1024;
#ifdef WIN32
    pthread_win32_process_attach_np( );
#endif // WIN32

    gengetopt_args_info sArgs;

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }

    signal(SIGPIPE, SIG_IGN);
    size_t i, j;
    for (i = 0; i < NUM_THREADS; i++) {
        THREAD_OCCUPIED[i] = 0;
    }

    if (sArgs.port_arg != NULL) {
        PORT = sArgs.port_arg;
    }

    SeekSettings settings;
    settings.species = "human";
    settings.port = atol(PORT);
    settings.numThreads = NUM_THREADS;
    settings.numBufferedDBs = 1;

    CSeekDBSetting *dbSetting = new CSeekDBSetting(
            "NA", // gvar arg, argument not needed for PCLServer
            sArgs.sinfo_arg, sArgs.platform_arg, sArgs.prep_arg,
            "NA",  // DB arg, argument not needed for PCLServer
            sArgs.gene_arg, "NA", sArgs.quant_arg, sArgs.dset_arg, 
            "NA", // dataset size file, not needed for PCLServer
            99999 // num_db arg, argument not needed for PCLServer
    );
    dbSetting->setPclDir(sArgs.input_arg);
    settings.dbs.push_back(dbSetting);

    string add_db = sArgs.additional_db_arg;
    if (add_db != "NA") {
        // leagacyReadDBConfig will append additional dbs to the cc CSeekDBSetting vector
        bool bres = legacyReadDBConfigFile(add_db, settings.dbs, 
                                           CSeekDataset::CORRELATION, false);
        if (bres == false) {
            fprintf(stderr, "Error in legacyReadDBConfigFile %s\n", add_db.c_str());
            return 1;
        }
    }

    CSeekCentral seekCentral;
    seekCentral.InitializeFromSeekConfig(settings);
    cout << "PCL Dir: " << seekCentral.m_vecDBSetting[0]->pclDir << endl;

    // Make a cache from string to shared_ptr CPCL - that way
    //  a returned shared pointer can be used without worrying
    //  if it will be evicted from the cache during use.
    LRUCache <string, PclPtrS> pclCache(DEFAULT_CACHE_SIZE);

    //==================================================


    int sockfd, new_fd;
    struct addrinfo hints, *servinfo, *p;
    struct sockaddr_storage their_addr;
    socklen_t sin_size;
    struct sigaction sa;
    char s[INET6_ADDRSTRLEN];
    char buf[10];
    int rv;
    int yes = 1;

    memset(&hints, 0, sizeof(hints));
    hints.ai_family = AF_UNSPEC;
    hints.ai_socktype = SOCK_STREAM;
    hints.ai_flags = AI_PASSIVE;

    if ((rv = getaddrinfo(NULL, PORT, &hints, &servinfo)) != 0) {
        fprintf(stderr, "getaddrinfo: %s\n", gai_strerror(rv));
        return 1;
    }

    // loop through all the results and bind to the first we can
    for (p = servinfo; p != NULL; p = p->ai_next) {
        if ((sockfd = socket(p->ai_family, p->ai_socktype,
                             p->ai_protocol)) == -1) {
            perror("server: socket");
            continue;
        }
        if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &yes,
                       sizeof(int)) == -1) {
            perror("setsockopt");
            exit(1);
        }
        if (::bind(sockfd, p->ai_addr, p->ai_addrlen) == -1) {
            close(sockfd);
            perror("server: bind");
            continue;
        }
        break;
    }

    if (p == NULL) {
        fprintf(stderr, "server: failed to bind\n");
        return 2;
    }

    freeaddrinfo(servinfo); // all done with this structure

    if (listen(sockfd, BACKLOG) == -1) {
        perror("listen");
        exit(1);
    }

    sa.sa_handler = sigchld_handler; // reap all dead processes
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESTART;
    if (sigaction(SIGCHLD, &sa, NULL) == -1) {
        perror("sigaction");
        exit(1);
    }

    struct thread_data thread_arg[NUM_THREADS];
    pthread_t th[NUM_THREADS];
    pthread_attr_t attr[NUM_THREADS];
    int d = 0;

    pthread_mutex_init(&mutexGet, NULL);

    fprintf(stderr, "Finished initializations.\n");
    printf("server: waiting for connections...\n");

    while (1) {
        sin_size = sizeof their_addr;
        new_fd = accept(sockfd, (struct sockaddr *) &their_addr, &sin_size);
        if (new_fd == -1) {
            perror("accept");
            continue;
        }
        inet_ntop(their_addr.ss_family, get_in_addr(
                (struct sockaddr *) &their_addr), s, sizeof s);
        printf("server, got connection from %s\n", s);

        pthread_mutex_lock(&mutexGet);
        for (d = 0; d < NUM_THREADS; d++) {
            if (THREAD_OCCUPIED[d] == 0 ||
                thread_arg[d].isComplete == true) {
                    break;
                }
        }
        if (d == NUM_THREADS) {
            close(new_fd);
            continue;
        }
        THREAD_OCCUPIED[d] = 1;
        thread_arg[d].isComplete = false;

        pthread_mutex_unlock(&mutexGet);

        pthread_attr_init(&attr[d]);
        pthread_attr_setdetachstate(&attr[d], PTHREAD_CREATE_DETACHED);

        string mode;

        if (CSeekNetwork::Receive(new_fd, mode) != 0) {
            fprintf(stderr, "Error: receiving message\n");
            close(new_fd);
            continue;
        }

        //mode, 5 digits
        //1 - output coexpression calculation (with gene hubbiness removed),
        //require setting two groups of genes, one for query (1 or many),
        //and one other gene (to calculate coexpression on)
        //2 - output gene-normalized expression instead of original expression
        //3 - output gene expression (can be normalized or unnormalized, depends on 2)
        //4 - output query expression
        //5 - output query coexpression (compared to 1, which is gene-to-query
        //    coexpression. This is query-to-query coexpression)

        bool outputCoexpression = false;
        bool outputNormalized = false;
        bool outputExpression = false;
        bool outputQueryExpression = false;
        bool outputQueryCoexpression = false;
        if (mode[0] == '1')
            outputCoexpression = true;
        if (mode[1] == '1')
            outputNormalized = true;
        if (mode[2] == '1')
            outputExpression = true;
        if (mode[3] == '1')
            outputQueryExpression = true;
        if (mode[4] == '1')
            outputQueryCoexpression = true;

        vector <string> dsetName, geneName, queryName;
        string qname, gname, dname;
        string strrbp; //-1: disabled, 0-0.9999
        float rbp_p = -1;

        if (outputCoexpression || outputQueryCoexpression) {
            if (CSeekNetwork::Receive(new_fd, qname) != 0 ||
                CSeekNetwork::Receive(new_fd, gname) != 0 ||
                CSeekNetwork::Receive(new_fd, dname) != 0 ||
                CSeekNetwork::Receive(new_fd, strrbp) != 0) {
                fprintf(stderr, "Error: receiving message\n");
                close(new_fd);
                continue;
            }
            CMeta::Tokenize(dname.c_str(), dsetName);
            CMeta::Tokenize(gname.c_str(), geneName);
            CMeta::Tokenize(qname.c_str(), queryName);
            rbp_p = atof(strrbp.c_str());
        } else {
            if (CSeekNetwork::Receive(new_fd, gname) != 0 ||
                CSeekNetwork::Receive(new_fd, dname) != 0) {
                fprintf(stderr, "Error: receiving message\n");
                close(new_fd);
                continue;
            }
            CMeta::Tokenize(dname.c_str(), dsetName);
            CMeta::Tokenize(gname.c_str(), geneName);
        }

        thread_arg[d].new_fd = new_fd;
        thread_arg[d].geneNames = geneName;
        thread_arg[d].datasetNames = dsetName;
        thread_arg[d].queryGeneNames = queryName;
        thread_arg[d].outputNormalized = outputNormalized;
        thread_arg[d].outputCoexpression = outputCoexpression;
        thread_arg[d].outputExpression = outputExpression;
        thread_arg[d].outputQueryExpression = outputQueryExpression;
        thread_arg[d].outputQueryCoexpression = outputQueryCoexpression;
        thread_arg[d].rbp_p = rbp_p;
        thread_arg[d].seekCentral = &seekCentral;
        thread_arg[d].pclCache = &pclCache;
        thread_arg[d].isComplete = false;

        fprintf(stderr, "Arguments: %d %d %s %s\n", d, new_fd, dname.c_str(), gname.c_str());
        if (outputCoexpression) {
            fprintf(stderr, "Arguments: %s\n", qname.c_str());
        }

        int ret;
        pthread_create(&th[d], &attr[d], do_query, (void *) &thread_arg[d]);
        pthread_detach(th[d]);
        pthread_attr_destroy(&attr[d]);
    }


#ifdef WIN32
    pthread_win32_process_detach_np( );
#endif // WIN32
    return 0;
}