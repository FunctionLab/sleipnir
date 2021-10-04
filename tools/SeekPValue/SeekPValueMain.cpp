/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "stdafx.h"
#include "cmdline.h"
#include "SeekPValue.h"

#define BACKLOG 10   // how many pending connections queue will hold
char *PORT;

char THREAD_OCCUPIED[NUM_THREADS];
pthread_mutex_t mutexGet;


void sigchld_handler(int s) {
    while (waitpid(-1, NULL, WNOHANG) > 0);
}

// get sockaddr, IPv4 or IPv6:
void *get_in_addr(struct sockaddr *sa) {
    if (sa->sa_family == AF_INET)
        return &(((struct sockaddr_in *) sa)->sin_addr);
    return &(((struct sockaddr_in6 *) sa)->sin6_addr);
}


int main(int iArgs, char **aszArgs) {
    static const size_t c_iBuffer = 1024;
#ifdef WIN32
    pthread_win32_process_attach_np( );
#endif // WIN32
    gengetopt_args_info sArgs{};

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return -1;
    }

    signal(SIGPIPE, SIG_IGN);
    size_t i;
    for (i = 0; i < NUM_THREADS; i++) {
        THREAD_OCCUPIED[i] = 0;
    }

    PORT = sArgs.port_arg;
    string strMode = sArgs.mode_arg;
    float nan = sArgs.nan_arg; //only used for strMode=="genes"

    SeekSettings settings;
    settings.species = "human";
    settings.port = atol(PORT);
    settings.numThreads = NUM_THREADS;
    settings.numBufferedDBs = 1;

    // sArgs.input_arg is gene map file
    CSeekDBSetting *dbSetting = new CSeekDBSetting(
            "NA", // gvar dir arg, argument not needed for PCLServer
            "NA", // sinfo dir arg, argument not needed for PCLServer
            "NA", // platform dir arg, argument not needed for PCLServer
            "NA", // prep dir arg, argument not needed for PCLServer
            "NA",  // DB arg, argument not needed for PCLServer
            sArgs.input_arg, 
            "NA", // gene symbol file arg, argument not needed for PCLServer
            "NA", // quant file arg, argument not needed for PCLServer
            sArgs.dset_platform_arg, 
            "NA", // dataset size file, not needed for PCLServer
            99999 // num_db arg, argument not needed for PCLServer
    );
    dbSetting->setRandomDir(sArgs.random_dir_arg);
    settings.dbs.push_back(dbSetting);

    CSeekCentral seekCentral;
    seekCentral.InitializeFromSeekConfig(settings);
    initializePvalue(seekCentral);


    //find a free port and attempt binding to the port
    int sockfd, new_fd;
    struct addrinfo hints{}, *servinfo, *p;
    struct sockaddr_storage their_addr{};
    socklen_t sin_size;
    struct sigaction sa{};
    char s[INET6_ADDRSTRLEN];
    char buf[10];
    int rv;
    int yes = 1;

    memset(&hints, 0, sizeof(hints));
    hints.ai_family = AF_UNSPEC;
    hints.ai_socktype = SOCK_STREAM;
    hints.ai_flags = AI_PASSIVE;

    if ((rv = getaddrinfo(nullptr, PORT, &hints, &servinfo)) != 0) {
        fprintf(stderr, "getaddrinfo: %s\n", gai_strerror(rv));
        return -1;
    }

    // loop through all the results and bind to the first we can
    for (p = servinfo; p != nullptr; p = p->ai_next) {
        if ((sockfd = socket(p->ai_family, p->ai_socktype, p->ai_protocol)) == -1) {
            perror("server: socket");
            continue;
        }
        if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(int)) == -1) {
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

    if (p == nullptr) {
        fprintf(stderr, "server: failed to bind\n");
        return -2;
    }

    freeaddrinfo(servinfo); // all done with this structure

    if (listen(sockfd, BACKLOG) == -1) {
        perror("listen");
        exit(1);
    }

    sa.sa_handler = sigchld_handler; // reap all dead processes
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESTART;
    if (sigaction(SIGCHLD, &sa, nullptr) == -1) {
        perror("sigaction");
        exit(1);
    }

    printf("server: waiting for connections...\n");
    struct thread_data thread_arg[NUM_THREADS];
    pthread_t th[NUM_THREADS];

    pthread_mutex_init(&mutexGet, nullptr);

    while (1) {
        sin_size = sizeof their_addr;
        new_fd = accept(sockfd, (struct sockaddr *) &their_addr, &sin_size);
        if (new_fd == -1) {
            perror("accept");
            continue;
        }
        inet_ntop(their_addr.ss_family, get_in_addr((struct sockaddr *) &their_addr), s, sizeof s);
        printf("server, got connection from %s\n", s);

        int d = 0;
        pthread_mutex_lock(&mutexGet);
        for (d = 0; d < NUM_THREADS; d++) {
            if (THREAD_OCCUPIED[d] == 0 ||
                thread_arg[d].isComplete == true) {
                    break;
            }
        }

        if (d == NUM_THREADS) {
            close(new_fd);
            pthread_mutex_unlock(&mutexGet);
            continue;
        }

        THREAD_OCCUPIED[d] = 1;
        thread_arg[d].isComplete = false;
        pthread_mutex_unlock(&mutexGet);

        thread_arg[d].threadid = d;
        thread_arg[d].new_fd = new_fd;

        if (strMode == "genes") {
            string strQuery;
            vector<float> vf;
            vector <string> query;
            string sMode;
            int mode;

            if (CSeekNetwork::Receive(new_fd, sMode) == -1) {
                fprintf(stderr, "Error receiving from client\n");
            }

            if (sMode == "rank")
                mode = 0;
            else if (sMode == "score")
                mode = 1;

            if (CSeekNetwork::Receive(new_fd, strQuery) == -1) {
                fprintf(stderr, "Error receiving from client!\n");
            }

            if (CSeekNetwork::Receive(new_fd, vf) == -1) {
                fprintf(stderr, "Error receiving from client!\n");
            }

            CMeta::Tokenize(strQuery.c_str(), query, " ");
            //=========================================================
            thread_arg[d].section = 0; //genes section
            thread_arg[d].query = query;
            thread_arg[d].gene_score = vf;
            thread_arg[d].nan = nan;
            thread_arg[d].mode = mode;
        } else if (strMode == "datasets") {
            string strDataset;
            vector <string> dataset;
            vector<float> qsize;
            vector<float> vf;
            if (CSeekNetwork::Receive(new_fd, strDataset) == -1) {
                fprintf(stderr, "Error receiving from client!\n");
            }
            if (CSeekNetwork::Receive(new_fd, vf) == -1) {
                fprintf(stderr, "Error receiving from client!\n");
            }
            if (CSeekNetwork::Receive(new_fd, qsize) == -1) {
                fprintf(stderr, "Error receiving from client!\n");
            }
            vector<int> vi;
            vi.resize(qsize.size());
            for (int ki = 0; ki < qsize.size(); ki++)
                vi[ki] = (int) qsize[ki];

            CMeta::Tokenize(strDataset.c_str(), dataset, " ");
            //========================================================
            thread_arg[d].dset = dataset;
            thread_arg[d].dset_score = vf;
            thread_arg[d].dset_qsize = vi;
            thread_arg[d].section = 1;
        }
        thread_arg[d].seekCentral = &seekCentral;
        thread_arg[d].isComplete = false;

        int ret;
        pthread_create(&th[d], NULL, do_pvalue_query, (void *) &thread_arg[d]);
    }

#ifdef WIN32
    pthread_win32_process_detach_np( );
#endif // WIN32
    return 0;

}
