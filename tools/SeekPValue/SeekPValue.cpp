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

#define BACKLOG 10   // how many pending connections queue will hold
char *PORT;

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

#define NUM_THREADS 8
char THREAD_OCCUPIED[NUM_THREADS];

//mode is genes============================================
map<string, int> mapstrintGene;
vector <string> vecstrGenes;
vector <string> vecstrGeneID;
vector <vector<int>> randomRank;
vector <vector<float>> randomSc;
vector<int> querySize;
int numGenes;
//=========================================================
struct parameter {
    int size;
    double scale;
    double shape;
    double threshold;
    double portion;
    vector<double> quantile;
};
map<string, int> mapstrintDataset;
vector <string> vecstrDataset;
vector <vector<struct parameter>> dsetScore; //co-expression score for dataset i, query-size j
//=========================================================

struct thread_data {
    int section; //0 - genes, 1 - datasets (which mode to turn on)
    //Section "genes"
    vector <string> query;
    vector<float> gene_score;
    int mode; //0 - p-value on rank, 1 - p-value on score
    float nan;
    //Section "datasets"
    vector <string> dset;
    vector<float> dset_score; //scores to test
    vector<int> dset_qsize; //number of genes for which coexpression score is calculated, for all dset
    //============================================
    int threadid;
    int new_fd;
};

void *do_query(void *th_arg) {
    struct thread_data *my = (struct thread_data *) th_arg;
    int new_fd = my->new_fd;
    int threadid = my->threadid;
    float nan = my->nan;
    int mode = my->mode;
    int section = my->section;

    vector <string> dset = my->dset;
    vector<float> dset_score = my->dset_score;
    vector<int> dset_qsize = my->dset_qsize;

    vector <string> queryGenes = my->query;
    vector<float> geneScores = my->gene_score;

    size_t i, j, jj, k;

    if (section == 1) { //dataset
        vector<float> pval;
        for (i = 0; i < dset.size(); i++) {
            if (mapstrintDataset.find(dset[i]) == mapstrintDataset.end()) {
                fprintf(stderr, "Error: cannot find dataset %s\n", dset[i].c_str());
                pval.push_back(-1);
                continue;
            }
            if (dset_qsize[i] < 2) {
                pval.push_back(0.99);
                continue;
            }
            if (dset_score[i] <= 0) {
                pval.push_back(0.99);
                continue;
            }
            int pi = mapstrintDataset[dset[i]];
            float sc = log(dset_score[i]);
            int qsize = dset_qsize[i];
            if (dsetScore[pi].size() == 0) {
                //This dataset is a dummy, null dataset (because it is bad or low quality
                pval.push_back(0.99);
                continue;
            }
            if (qsize - 2 >= dsetScore[pi].size()) {
                //Don't have modeling for this query size, most likely because query size > 100
                pval.push_back(0.99);
                continue;
            }
            struct parameter &par = dsetScore[pi][qsize - 2];
            vector<double> &quantile = par.quantile;
            //fprintf(stderr, "threshold %.2e, scale %.2e, shape %.2e, qsize %d\n", par.threshold, par.scale, par.shape, qsize);

            if ((double) sc <= par.threshold) {
                double min_diff = 9999;
                int min_diff_index = -1;
                /*float p_val = 0.05;
                for(k=quantile.size()-1; k>=0; k--){
                    if( (double) sc > quantile[k] ){
                        break;
                    }else{
                        p_val+=0.05;
                    }
                }
                if(p_val >= 0.99){
                    p_val = 0.99;
                }*/
                int cut_index = 0;
                for (k = 0; k <= quantile.size(); k++) {
                    if (quantile[k] >= par.threshold) {
                        cut_index = k;
                        break;
                    }
                }

                for (k = 0; k <= cut_index; k++) {
                    double diff = -1;
                    if (k == cut_index) {
                        diff = fabs(par.threshold - (double) sc);
                        if (diff < min_diff) {
                            min_diff = diff;
                            min_diff_index = k;
                        }
                        continue;
                    }
                    diff = fabs(quantile[k] - (double) sc);
                    if (diff < min_diff) {
                        min_diff = diff;
                        min_diff_index = k;
                    }
                    //fprintf(stderr, "Diff %d %.2f %.2f %d\n", k, quantile[k], diff, qsize);
                }
                if (min_diff_index == -1) {
                    //fprintf(stderr, "Error, negative one!\n");
                }
                float pv = 1.0 - (0.05 + (float) min_diff_index * 0.01);
                pval.push_back(pv);
                //pval.push_back(p_val);
                //fprintf(stderr, "original %.2e pval %.2e\n", sc, pv);
                continue;
            }

            double limit = -1.0 * par.scale / par.shape;
            limit = floor(limit * 100.0) / 100.0;

            double diff = (double) sc - par.threshold;
            double cdf = 0;

            if (par.shape == 0) {
                cdf = 1.0 - exp(-1.0 * diff / par.scale);
            } else {
                if (par.shape < 0 && diff >= limit) { //check that difference is within range
                    diff = limit; //if difference is bigger than it can accept, assume limit
                }
                cdf = 1.0 - pow(1.0 + par.shape * diff / par.scale, -1.0 / par.shape);
            }
            double pv = par.portion * (1.0 - cdf);
            //fprintf(stderr, "before original %.2e pval %.2e\n", sc, pv);
            if (isnan(pv)) {
                pv = 2.0e-5;
            }
            pval.push_back((float) pv);
            //fprintf(stderr, "original %.2e pval %.2e\n", sc, pv);
        }

        if (CSeekNetwork::Send(new_fd, pval) == -1) {
            fprintf(stderr, "Error sending message to client!\n");
        }


    } else if (section == 0) { //genes
        vector <AResultFloat> sortedGenes;
        sortedGenes.resize(geneScores.size());
        for (i = 0; i < sortedGenes.size(); i++) {
            sortedGenes[i].i = i;
            sortedGenes[i].f = geneScores[i];
        }

        vector<int> queryGeneID;
        for (i = 0; i < queryGenes.size(); i++)
            queryGeneID.push_back(mapstrintGene[queryGenes[i]]);
        //Query genes themselves have lowest score, to prevent
        //them from being counted in PR
        for (i = 0; i < queryGeneID.size(); i++)
            sortedGenes[queryGeneID[i]].f = nan;

        sort(sortedGenes.begin(), sortedGenes.end());

        //comparison
        vector<int> geneRank;
        geneRank.resize(numGenes);
        for (jj = 0; jj < numGenes; jj++) {
            geneRank[sortedGenes[jj].i] = jj;
        }

        vector<float> pval;
        CSeekTools::InitVector(pval, geneScores.size(), (float) nan);

        for (jj = 0; jj < geneScores.size(); jj++) {
            int gene = sortedGenes[jj].i;
            int gene_rank = jj;
            float gene_score = sortedGenes[jj].f;
            if (gene_score == nan) break;
            //if(gene_score<0)
            //	continue;
            vector<int> &rR = randomRank[gene];
            vector<float> &rF = randomSc[gene];
            int kk = 0;
            if (mode == 1) {
                if (gene_score >= 0) {
                    for (kk = 0; kk < rF.size(); kk++) {
                        if (gene_score >= rF[kk] || kk == rF.size() - 1)
                            pval[gene] = (float) kk / (float) rF.size();
                        //fprintf(stderr, "%s\t%d\t%d\t%.5e\t%.5e\n", vecstrGenes[gene].c_str(),
                        //gene_rank, kk, gene_score, randomSc[gene][kk]);
                        if (gene_score >= rF[kk])
                            break;
                    }
                } else {
                    for (kk = rF.size() - 1; kk >= 0; kk--) {
                        if (gene_score <= rF[kk] || kk == 0)
                            pval[gene] = (float) (rF.size() - 1 - kk) / (float) rF.size();
                        //fprintf(stderr, "%s\t%d\t%d\t%.5e\t%.5e\n", vecstrGenes[gene].c_str(),
                        //gene_rank, rF.size()-1-kk, gene_score, randomSc[gene][kk]);
                        if (gene_score <= rF[kk])
                            break;
                    }
                }
            } else if (mode == 0) {
                if (gene_rank < 17600 / 2) {
                    for (kk = 0; kk < rR.size(); kk++) {
                        if (gene_rank <= rR[kk] || kk == rR.size() - 1)
                            pval[gene] = (float) kk / (float) rR.size();
                        //fprintf(stderr, "%s\t%d\t%d\t%.5e\t%.5e\n", vecstrGenes[gene].c_str(),
                        //gene_rank, kk, gene_score, randomSc[gene][kk]);
                        if (gene_rank <= rR[kk])
                            break;
                    }
                } else {
                    for (kk = rR.size() - 1; kk >= 0; kk--) {
                        if (gene_rank >= rR[kk] || kk == 0)
                            pval[gene] = (float) (rR.size() - 1 - kk) / (float) rF.size();
                        //fprintf(stderr, "%s\t%d\t%d\t%.5e\t%.5e\n", vecstrGenes[gene].c_str(),
                        //gene_rank, rR.size()-1-kk, gene_score, randomSc[gene][kk]);
                        if (gene_rank >= rR[kk])
                            break;
                    }
                }
            }
        }

        if (CSeekNetwork::Send(new_fd, pval) == -1) {
            fprintf(stderr, "Error sending message to client!\n");
        }
    }


    pthread_mutex_lock(&mutexGet);
    close(new_fd);
    THREAD_OCCUPIED[threadid] = 0;
    pthread_mutex_unlock(&mutexGet);
    int ret = 0;
    pthread_exit((void *) ret);
}
//mode ends=========================================================

bool ReadParameter(const string& param_file, vector<struct parameter> &v) {
    ifstream ifsm;
    ifsm.open(param_file.c_str());
    if (!ifsm.is_open()) {
        fprintf(stderr, "Error opening file %s\n", param_file.c_str());
        return false;
    }
    const int lineSize = 6084;
    char acBuffer[lineSize];
    utype c_iBuffer = lineSize;
    v.clear();

    while (!ifsm.eof()) {
        ifsm.getline(acBuffer, c_iBuffer - 1);
        if (acBuffer[0] == 0) break;
        acBuffer[c_iBuffer - 1] = 0;
        vector <string> tok;
        CMeta::Tokenize(acBuffer, tok);
        struct parameter par;
        par.size = strtol(tok[0].c_str(), nullptr, 10);
        par.portion = strtod(tok[2].c_str(), nullptr);
        par.threshold = strtod(tok[3].c_str(), nullptr);
        par.scale = strtod(tok[5].c_str(), nullptr);
        par.shape = strtod(tok[6].c_str(), nullptr);
        par.quantile = vector<double>();
        //for(int k=8; k<tok.size(); k++){
        for (int k = 7; k < tok.size(); k++) {
            par.quantile.push_back(strtod(tok[k].c_str(), nullptr));
            //fprintf(stderr, "This value is %d %.2f\n", k, atof(tok[k].c_str()));
        }
        v.push_back(par);
    }
    ifsm.close();
    return true;
}


int main(int iArgs, char **aszArgs) {
    static const size_t c_iBuffer = 1024;
#ifdef WIN32
    pthread_win32_process_attach_np( );
#endif // WIN32
    gengetopt_args_info sArgs{};

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }

    signal(SIGPIPE, SIG_IGN);
    size_t i;
    for (i = 0; i < NUM_THREADS; i++) {
        THREAD_OCCUPIED[i] = 0;
    }

    PORT = sArgs.port_arg;
    string strMode = sArgs.mode_arg;
    float nan = sArgs.nan_arg; //only used for strMode=="genes"

    //preparation=====================================================
    if (strMode == "genes") {
        if (!CSeekTools::ReadListTwoColumns(sArgs.input_arg, vecstrGeneID, vecstrGenes))
            return false;
        for (i = 0; i < vecstrGenes.size(); i++)
            mapstrintGene[vecstrGenes[i]] = (int) i;

        numGenes = vecstrGenes.size();

        string random_directory = sArgs.random_dir_arg;
        int num_random = sArgs.random_num_arg;
        int ii, jj;
        char ac[256];

        randomRank.resize(numGenes);
        randomSc.resize(numGenes);
        for (ii = 0; ii < numGenes; ii++) {
            randomRank[ii].resize(num_random);
            randomSc[ii].resize(num_random);
        }

        for (ii = 0; ii < num_random; ii++) {
            vector<float> randomScores;
            sprintf(ac, "%s/%d.gscore", random_directory.c_str(), ii);
            CSeekTools::ReadArray(ac, randomScores);
            /*vector<string> queryGenes;
            sprintf(ac, "%s/%d.query", random_directory.c_str(), ii);
            CSeekTools::ReadMultiGeneOneLine(ac, queryGenes);
            querySize.push_back(queryGenes.size());
            */
            vector <AResultFloat> sortedRandom;
            sortedRandom.resize(randomScores.size());
            for (jj = 0; jj < randomScores.size(); jj++) {
                sortedRandom[jj].i = jj;
                sortedRandom[jj].f = randomScores[jj];
            }
            sort(sortedRandom.begin(), sortedRandom.end());
            for (jj = 0; jj < randomScores.size(); jj++) {
                randomRank[sortedRandom[jj].i][ii] = jj;
                randomSc[sortedRandom[jj].i][ii] = sortedRandom[jj].f;
            }
        }

        for (jj = 0; jj < numGenes; jj++) {
            sort(randomRank[jj].begin(), randomRank[jj].end());
            sort(randomSc[jj].begin(), randomSc[jj].end(), std::greater<float>());
        }
    } else if (strMode == "datasets") {
        vector <string> vD, vDP;
        if (!CSeekTools::ReadListTwoColumns(sArgs.dset_platform_arg, vD, vDP))
            return false;
        for (i = 0; i < vD.size(); i++) {
            vecstrDataset.push_back(vD[i]);
        }
        for (i = 0; i < vecstrDataset.size(); i++) {
            mapstrintDataset[vecstrDataset[i]] = i;
        }
        string param_dir = sArgs.param_dir_arg;
        dsetScore.resize(vecstrDataset.size());
        for (i = 0; i < vecstrDataset.size(); i++) {
            string param_file = param_dir + "/" + vecstrDataset[i] + ".param";
            dsetScore[i] = vector<struct parameter>();
            if (!ReadParameter(param_file, dsetScore[i])) {
                fprintf(stderr, "Making this dataset null... (will always return insignificant)\n");
            }
        }
    }

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
        return 1;
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
            if (THREAD_OCCUPIED[d] == 0) break;
        }

        if (d == NUM_THREADS) {
            close(new_fd);
            pthread_mutex_unlock(&mutexGet);
            continue;
        }

        THREAD_OCCUPIED[d] = 1;
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

        int ret;
        pthread_create(&th[d], NULL, do_query, (void *) &thread_arg[d]);
    }

#ifdef WIN32
    pthread_win32_process_detach_np( );
#endif // WIN32
    return 0;

}
