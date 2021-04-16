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
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>
#include <list>

#define BACKLOG 20   // how many pending connections queue will hold
char *PORT;
int NUM_DSET_MEMORY = 100;
string pcl_input_dir;

//CPCL **pcl;
list<int> available;
char *loaded;
map<string, int> DNAME_MAP;
map<int, int> return_val;
map<int, string> DNAME_RMAP;

pthread_mutex_t mutexGet;

vector<CSeekDBSetting *> cc;

/*string strPrepInputDirectory;
string strSinfoInputDirectory;
string strDatasetInputDirectory;
string strPlatformInputDirectory;*/
map<string, int> mapstrintGene;

vector <string> vecstrGeneID;
vector <string> vecstrGenes;
vector <string> vecstrDatasets;
vector <string> vecstrDP;

map<string, int> mapstrintDatasetDB;
map <string, utype> mapstriPlatform;
map <string, string> mapstrstrDatasetPlatform;
map <string, utype> mapstrintDataset;
vector <string> vecstrPlatforms;
vector <CSeekPlatform> vp;

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

#define NUM_THREADS 16
char THREAD_OCCUPIED[NUM_THREADS];

int send_msg(int new_fd, char *c, int size) {
    int r = send(new_fd, c, size, 0);
    if (r == -1) {
        printf("client exists");
    }
    return r;
}

void cl(char *b, int size) {
    int i = 0;
    for (i = 0; i < size; i++) {
        b[i] = '\0';
    }
}

struct thread_data {
    vector <string> datasetName;
    vector <string> geneName;
    vector <string> queryName;
    int threadid;
    int new_fd;
    float rbp_p;
    bool outputNormalized;
    bool outputCoexpression;
    bool outputQueryCoexpression;
    bool outputExpression;
    bool outputQueryExpression;
    CPCL **pcl;
    CPCL **vc;
};

int cpy(char *d, char *s, int beg, int num) {
    int i;
    for (i = 0; i < num; i++) {
        d[beg + i] = s[i];
    }
    return beg + num;
}

int GetOpenSlot() {
    int i = -1;
    int size = 0;
    while (size < available.size()) {
        i = available.front();
        available.pop_front();
        available.push_back(i);
        if (loaded[i] == 0) break;
        size++;
    }
    if (size == available.size()) {
        return -1;
    }
    return i;
}

void *do_query(void *th_arg) {
    struct thread_data *my = (struct thread_data *) th_arg;
    vector <string> datasetName = my->datasetName;
    vector <string> geneName = my->geneName;
    vector <string> queryName = my->queryName;
    bool outputNormalized = my->outputNormalized;
    bool outputCoexpression = my->outputCoexpression;
    bool outputQueryCoexpression = my->outputQueryCoexpression;
    bool outputExpression = my->outputExpression;
    bool outputQueryExpression = my->outputQueryExpression;
    int new_fd = my->new_fd;
    int tid = my->threadid;
    float RBP_P = my->rbp_p;
    CPCL **pcl = my->pcl;
    CPCL **vc = my->vc;

    /*
    cl(buf, 5000);
    sprintf(buf, "Begin...\n");
    if(send_msg(new_fd, buf, 5000)==-1){
        THREAD_OCCUPIED[tid] = 0;
        pthread_exit(0);
    }*/

    //pthread_mutex_lock(&mutexGet);

    size_t i;

    vector<string>::const_iterator iterS = datasetName.begin();
    vector<int> occupied;

    fprintf(stderr, "before allocation...\n");
    vc = new CPCL *[datasetName.size()];
    for (i = 0; i < datasetName.size(); i++) {
        vc[i] = new CPCL();
        fprintf(stderr, "%p\n", vc[i]);
    }

    fprintf(stderr, "start processing...\n");

    pthread_mutex_lock(&mutexGet); //QZ disabled 1/31/2015

    for (i = 0; i < NUM_DSET_MEMORY; i++) {
        loaded[i] = 0;
    }
    int ind = -1;
    for (ind = 0; iterS != datasetName.end(); iterS++, ind++) {
        int n = -1;
        map<string, int>::const_iterator iterM = DNAME_MAP.find(*iterS);
        if (iterM != DNAME_MAP.end()) {
            n = iterM->second;
            fprintf(stderr, "found %d for dataset %s...\n", n, iterS->c_str());
            vc[ind]->Open(*pcl[n]);
            loaded[n] = 1;
            continue;
        }

        n = GetOpenSlot();
        map<int, string>::const_iterator iterRM = DNAME_RMAP.find(n);
        if (iterRM != DNAME_RMAP.end()) {
            DNAME_MAP.erase(iterRM->second);
            DNAME_RMAP.erase(n);
        }
        DNAME_MAP[*iterS] = n;
        DNAME_RMAP[n] = *iterS;
        loaded[n] = 1;

        fprintf(stderr, "acquired %d for dataset %s...\n", n, iterS->c_str());
        //pcl[n]->Reset();

        fprintf(stderr, "dataset reset \n");
        //string pcl_path = pcl_input_dir + "/" + *iterS; //for human-SEEK
        string pcl_path = pcl_input_dir + "/" + *iterS + ".bin"; //for model-organism-SEEK

        if (pcl[n] != NULL) {
            delete pcl[n];
        }
        pcl[n] = new CPCL();
        pcl[n]->Open(pcl_path.c_str());

        fprintf(stderr, "Allocating in created is called\n");
        fprintf(stderr, "%ld\n", CMeta::GetMemoryUsage());
        vc[ind]->Open(*pcl[n]);
        fprintf(stderr, "dataset opened\n");
    }
    for (i = 0; i < NUM_DSET_MEMORY; i++) {
        loaded[i] = 0;
    }


    pthread_mutex_unlock(&mutexGet); //QZ disabled 1/31

    int genes = geneName.size();
    int queries = queryName.size(); //(EXTRA)
    int datasets = datasetName.size();

    vector<float> vecG, vecQ, vecCoexpression, vecqCoexpression;
    vector<float> sizeD;

    sizeD.resize(datasets);

    vector <vector<float>> d_vecG, d_vecQ, d_vecCoexpression, d_vecqCoexpression;
    d_vecG.resize(datasets);
    d_vecQ.resize(datasets);
    d_vecCoexpression.resize(datasets);
    d_vecqCoexpression.resize(datasets);

    //Qian added
    //vector<CSeekDataset *> vd; //(EXTRA)
    //vd.resize(vc.size()); //(EXTRA)
    //CFullMatrix<float> *vCoexpression = new CFullMatrix<float>(); //(EXTRA)
    //vCoexpression->Initialize(genes, datasets); //(EXTRA)

    fprintf(stderr, "Reading data...\n");

    float NaN = -9999;
    //for each dataset

    vector<float> quant;

    //QUANT file must be consistent
    CSeekTools::ReadQuantFile(cc[0]->quantFile, quant);
    //CSeekTools::ReadQuantFile("/home/qzhu/Seek/quant2", quant);

#pragma omp parallel for \
    private(i) \
    firstprivate(genes, queries, datasets, outputCoexpression, outputQueryCoexpression, outputNormalized, outputExpression, \
    outputQueryExpression, NaN) \
    shared(pcl, vc, sizeD, datasetName, queryName, geneName, d_vecG, d_vecQ, d_vecCoexpression, d_vecqCoexpression, \
    mapstrintGene, mapstrstrDatasetPlatform, mapstriPlatform, vp) \
    schedule(dynamic)
    for (i = 0; i < datasets; i++) {
        CPCL *pp = vc[i];
        size_t j, k;
        int ps = pp->GetExperiments();
        int gs = pp->GetExperiments();

        CFullMatrix<float> *fq = NULL;

        CFullMatrix<float> *fq_RBP = NULL;

        CFullMatrix<float> *ff = NULL;
        vector<float> vCoexpression;
        vCoexpression.resize(genes);
        vector<float> vqCoexpression;
        vqCoexpression.resize(queryName.size());

        CSeekDataset *vd = NULL;
        CSeekPlatform *pl = NULL;
        sizeD[i] = (float) ps;
        d_vecG[i] = vector<float>();
        d_vecQ[i] = vector<float>();
        d_vecCoexpression[i] = vector<float>();
        d_vecqCoexpression[i] = vector<float>();

        set <string> absent;

        if (outputCoexpression || outputQueryCoexpression) {
            vd = new CSeekDataset();
            //string strFileStem = datasetName[i].substr(0, datasetName[i].find(".bin")); //for human-SEEK
            string strFileStem = datasetName[i]; //for model-organism-SEEK

            int dbID = mapstrintDatasetDB[strFileStem];

            string strAvgPath =
                    cc[dbID]->prepDir + "/" + strFileStem + ".gavg"; //avg and prep path share same directory
            string strPresencePath = cc[dbID]->prepDir + "/" + strFileStem + ".gpres";
            string strSinfoPath = cc[dbID]->sinfoDir + "/" + strFileStem + ".sinfo";


            vd->ReadGeneAverage(strAvgPath);
            vd->ReadGenePresence(strPresencePath);
            vd->ReadDatasetAverageStdev(strSinfoPath);
            vd->InitializeGeneMap();

            string strPlatform = mapstrstrDatasetPlatform.find(strFileStem)->second;
            utype platform_id = mapstriPlatform.find(strPlatform)->second;
            vd->SetPlatform(vp[platform_id]);
            pl = &vd->GetPlatform();

            fq = new CFullMatrix<float>();
            fq->Initialize(queryName.size(), ps);

            //NEW 3/5/2013====================================
            if (outputQueryCoexpression) {
                const vector <string> gNames = pp->GetGeneNames();
                fq_RBP = new CFullMatrix<float>();
                fq_RBP->Initialize(gNames.size(), ps);
                //fprintf(stderr, "%d %d X1\n", i, gNames.size());
                //#pragma omp parallel for \
				shared(pp, fq_RBP) private(k, j) firstprivate(gs) \
				schedule(dynamic)
                /*for(k=0; k<gNames.size(); k++){
                    int g = pp->GetGene(gNames[k]);
                    float *vv = pp->Get(g);
                    for(j=0; j<gs; j++)
                        fq_RBP->Set(g, j, vv[j]);
                    float mean = 0;
                    for(j=0; j<gs; j++)
                        mean+=fq_RBP->Get(g, j);
                    mean/=(float) (gs);
                    float stdev = 0;
                    for(j=0; j<gs; j++)
                        stdev+=(fq_RBP->Get(g, j) - mean) * (fq_RBP->Get(g, j) - mean);
                    stdev/=(float) (gs);
                    if(stdev<=0){
                        absent.insert(gNames[k]);
                    }
                    stdev = sqrt(stdev);
                    if(isnan(stdev) || isinf(stdev)){
                        fprintf(stderr, "Error:, standard deviation is zero\n");
                    }
                    for(j=0; j<gs; j++){
                        float t1 = fq_RBP->Get(g, j);
                        fq_RBP->Set(g, j, (t1 - mean) / stdev);

                    }
                }*/
                //fprintf(stderr, "%d X2\n", i);
            }
            //====================================================


            if (outputQueryExpression) {
                for (k = 0; k < queryName.size(); k++) {
                    int g = pp->GetGene(queryName[k]);
                    if (g == -1) { //does not contain the gene in the dataset
                        for (j = 0; j < gs; j++) {
                            fq->Set(k, j, NaN);
                            d_vecQ[i].push_back(fq->Get(k, j));
                        }
                        continue;
                    }
                    float *vv = pp->Get(g);
                    for (j = 0; j < gs; j++)
                        fq->Set(k, j, vv[j]);
                    if (!outputNormalized) {
                        for (j = 0; j < gs; j++)
                            d_vecQ[i].push_back(fq->Get(k, j));
                    }

                    //normalize
                    float mean = 0;
                    for (j = 0; j < gs; j++)
                        mean += fq->Get(k, j);
                    mean /= (float) (gs);
                    float stdev = 0;
                    for (j = 0; j < gs; j++)
                        stdev += (fq->Get(k, j) - mean) * (fq->Get(k, j) - mean);
                    stdev /= (float) (gs);
                    stdev = sqrt(stdev);
                    for (j = 0; j < gs; j++) {
                        float t1 = fq->Get(k, j);
                        fq->Set(k, j, (t1 - mean) / stdev);
                    }

                    if (outputNormalized) {
                        for (j = 0; j < gs; j++)
                            d_vecQ[i].push_back(fq->Get(k, j));
                    }
                }
            }
        }

        fprintf(stderr, "allocating space %lu %d...\n", geneName.size(),
                ps);
        ff = new CFullMatrix<float>();
        ff->Initialize(genes, ps);
        fprintf(stderr, "done allocating space.\n");

        if (outputExpression) {
            for (k = 0; k < geneName.size(); k++) {
                int g = pp->GetGene(geneName[k]);
                if (g == -1) {
                    for (j = 0; j < gs; j++) {
                        ff->Set(k, j, NaN);
                        d_vecG[i].push_back(ff->Get(k, j));
                    }
                    continue;
                }
                float *vv = pp->Get(g);
                for (j = 0; j < gs; j++)
                    ff->Set(k, j, vv[j]);

                if (!outputNormalized) {
                    for (j = 0; j < gs; j++) {
                        d_vecG[i].push_back(ff->Get(k, j));
                    }
                }

                //normalize
                float mean = 0;
                for (j = 0; j < gs; j++)
                    mean += ff->Get(k, j);
                mean /= (float) (gs);
                float stdev = 0;
                for (j = 0; j < gs; j++)
                    stdev += (ff->Get(k, j) - mean) * (ff->Get(k, j) - mean);
                stdev /= (float) (gs);
                stdev = sqrt(stdev);
                for (j = 0; j < gs; j++) {
                    float t1 = ff->Get(k, j);
                    ff->Set(k, j, (t1 - mean) / stdev);
                }

                if (outputNormalized) {
                    for (j = 0; j < gs; j++) {
                        d_vecG[i].push_back(ff->Get(k, j));
                    }
                }
            }
        }

        if (outputCoexpression) {
            int kk;
            CMeasurePearNorm pn;
            for (k = 0; k < geneName.size(); k++) {
                int g = pp->GetGene(geneName[k]);
                if (g == -1) {
                    vCoexpression[k] = NaN;
                    d_vecCoexpression[i].push_back(vCoexpression[k]);
                    continue;
                }
                float avgP = 0;
                int totalQueryPresent = queryName.size();
                for (kk = 0; kk < queryName.size(); kk++) {
                    int gg = pp->GetGene(queryName[kk]);
                    if (gg == -1) {
                        totalQueryPresent--;
                        continue;
                    }

                    float *x1 = NULL;
                    float *x2 = NULL;
                    if (g < gg) {
                        x1 = pp->Get(g);
                        x2 = pp->Get(gg);
                    } else {
                        x1 = pp->Get(gg);
                        x2 = pp->Get(g);
                    }

                    float p = (float) pn.Measure(x1, ps, x2, ps,
                                                 IMeasure::EMapCenter, NULL, NULL);

                    p = (p - vd->GetDatasetAverage()) / vd->GetDatasetStdev();

                    int gID = mapstrintGene[geneName[k]];
                    int qID = mapstrintGene[queryName[kk]];

                    int qb = CMeta::Quantize(p, quant);
                    p = quant[qb];

                    //fprintf(stderr, "Correlation %d %d %.5f %.5f %.5f %.5f\n", qID, gID, p, vd->GetGeneAverage(gID),
                    //	pl->GetPlatformAvg(qID), pl->GetPlatformStdev(qID));

                    //subtract hubbiness
                    p = (p - vd->GetGeneAverage(gID) - pl->GetPlatformAvg(qID)) / pl->GetPlatformStdev(qID);
                    //p = p - vd->GetGeneAverage(gID);
                    p = max((float) min(p, (float) 3.2), (float) -3.2);

                    avgP += p;
                }
                if (totalQueryPresent == 0)
                    avgP = NaN;
                else
                    avgP /= (float) (totalQueryPresent);

                vCoexpression[k] = avgP;
                d_vecCoexpression[i].push_back(avgP);
            }
        }

        if (outputQueryCoexpression) {
            int kk = 0;

            const vector <string> gNames = pp->GetGeneNames();
            vector<char> qMap;
            CSeekTools::InitVector(qMap, vecstrGenes.size(), (char) 0);

            int totQuery = 0;
            for (k = 0; k < queryName.size(); k++) {
                int g = pp->GetGene(queryName[k]);
                if (g == -1) continue;
                int mg = mapstrintGene[queryName[k]];
                qMap[mg] = 1;
                totQuery++;
            }

            CMeasurePearNorm pn;

            for (k = 0; k < queryName.size(); k++) {
                int g = pp->GetGene(queryName[k]);
                if (g == -1) {
                    vqCoexpression[k] = NaN;
                    d_vecqCoexpression[i].push_back(vqCoexpression[k]);
                    continue;
                }
                vector <AResult> ar;
                ar.resize(gNames.size());

                for (kk = 0; kk < gNames.size(); kk++) {
                    int gg = pp->GetGene(gNames[kk]);
                    ar[kk].i = (utype) mapstrintGene[gNames[kk]];
                    if (g == gg) {
                        ar[kk].f = 0;
                        continue;
                    }

                    float *x1 = NULL;
                    float *x2 = NULL;
                    if (g < gg) {
                        x1 = pp->Get(g);
                        x2 = pp->Get(gg);
                    } else {
                        x1 = pp->Get(gg);
                        x2 = pp->Get(g);
                    }

                    float p = (float) pn.Measure(x1, ps, x2, ps,
                                                 IMeasure::EMapCenter, NULL, NULL);

                    float px = p;
                    if (!(p < 5.0 && p > -5.0)) {
                        ar[kk].f = 0;
                        continue;
                    }

                    //get z-score (dataset wide)
                    p = (p - vd->GetDatasetAverage()) / vd->GetDatasetStdev();
                    int gID = mapstrintGene[gNames[kk]];
                    int qID = mapstrintGene[queryName[k]];

                    int qb = CMeta::Quantize(p, quant);
                    p = quant[qb];

                    //fprintf(stderr, "Correlation %d %d %.5f %.5f %.5f %.5f\n", qID, gID, p, vd->GetGeneAverage(gID),
                    //	pl->GetPlatformAvg(qID), pl->GetPlatformStdev(qID));

                    //subtract hubbiness
                    p = (p - vd->GetGeneAverage(gID) - pl->GetPlatformAvg(qID)) / pl->GetPlatformStdev(qID);
                    //p = p - vd->GetGeneAverage(gID);
                    p = max((float) min(p, (float) 3.2), (float) -3.2);

                    //fprintf(stderr, "%d error, infinite or nan! %.3f, %.3f, %.3f\n", i,
                    //	vd->GetDatasetAverage(), vd->GetDatasetStdev(),
                    //	vd->GetGeneAverage(mapstrintGene[gNames[kk]]));
                    //fprintf(stderr, "Correlation %d %d %.5f\n", qID, gID, p);
                    ar[kk].f = (utype)(p * 100.0 + 320);

                }

                int TOP = 1000;
                if (ar.size() < TOP) {
                    TOP = ar.size();
                }
                //fprintf(stderr, "%d H2\n", i);
                nth_element(ar.begin(), ar.begin() + TOP, ar.end());
                sort(ar.begin(), ar.begin() + TOP);
                //fprintf(stderr, "%d H3\n", i);

                float rbp = 0;
                //utype jj = 0;
                for (kk = 0; kk < TOP; kk++) {
                    if (qMap[ar[kk].i] == 0) continue;
                    if (ar[kk].f == 0) break;
                    rbp += pow(RBP_P, kk);
                    fprintf(stderr, "Sorted %d %d %d %.5f\n", i, kk, ar[kk].i, (ar[kk].f - 320) / 100.0f);
                    //jj++;
                    //fprintf(stderr, "Good %d %d\n", i, kk);
                    //jj++;
                }
                rbp *= (1.0 - RBP_P);

                rbp = rbp / totQuery * 1000;
                fprintf(stderr, "%zu %.3e\n", i, rbp);
                vqCoexpression[k] = rbp;
                d_vecqCoexpression[i].push_back(rbp);
            }


            /*for(k=0; k<queryName.size(); k++){
                int g = pp->GetGene(queryName[k]);
                if(g==-1){
                    vqCoexpression[k] = NaN;
                    vecqCoexpression.push_back(vqCoexpression[k]);
                    continue;
                }
                float avgP = 0;
                int totalQueryPresent = queryName.size() - 1;

                for(kk=0; kk<queryName.size(); kk++){
                    if(kk==k) continue;
                    int gg = pp->GetGene(queryName[kk]);
                    if(gg==-1){
                        totalQueryPresent--;
                        continue;
                    }

                    float *x1 = NULL;
                    float *x2 = NULL;
                    if(g<gg){
                        x1 = pp->Get(g);
                        x2 = pp->Get(gg);
                    }else{
                        x1 = pp->Get(gg);
                        x2 = pp->Get(g);
                    }

                    float p = (float) pn.Measure(x1, ps, x2, ps,
                        IMeasure::EMapCenter, NULL, NULL);
                    float px = p;

                    if(!(p<5.0 && p>-5.0)){
                        continue;
                    }

                    //get z-score (dataset wide)
                    p = (p - vd->GetDatasetAverage()) / vd->GetDatasetStdev();
                    int gID = mapstrintGene[queryName[kk]];
                    int qID = mapstrintGene[queryName[k]];

                    int qb = CMeta::Quantize(p, quant);
                    p = quant[qb];

                    //fprintf(stderr, "Correlation %d %d %.5f %.5f %.5f %.5f\n", qID, gID, p, vd->GetGeneAverage(gID),
                    //	pl->GetPlatformAvg(qID), pl->GetPlatformStdev(qID));

                    //subtract hubbiness
                    p = (p - vd->GetGeneAverage(gID) - pl->GetPlatformAvg(qID)) / pl->GetPlatformStdev(qID) ;
                    //p = p - vd->GetGeneAverage(gID);
                    p = max((float) min(p, (float) 3.2), (float)-3.2);

                    //float p = 0;
                    //for(j=2; j<gs; j++)
                    //	p+= fq->Get(k, j-2)*fq->Get(kk, j-2);
                    //p/=(float)(gs-2);
                    //p = 0.5 * log((1.0+p)/(1.0-p));
                    //p = max((float) min(p, (float) 3.2), (float)-3.2);
                    //get z-score (dataset wide)
                    //p = (p - vd->GetDatasetAverage()) / vd->GetDatasetStdev();
                    //subtract hubbiness
                    //p = p - vd->GetGeneAverage(mapstrintGene[queryName[kk]]);


                    avgP+=p;
                }
                if(totalQueryPresent==0)
                    avgP = NaN;
                else
                    avgP/=(float)(totalQueryPresent);

                vqCoexpression[k] = avgP;
                d_vecqCoexpression[i].push_back(avgP);
            }*/



        }

        if (outputCoexpression || outputQueryCoexpression) {
            delete vd;
            delete fq;
            if (outputQueryCoexpression)
                delete fq_RBP;
        }

        delete ff;
    }

    for (i = 0; i < datasets; i++) {
        vecG.insert(vecG.end(), d_vecG[i].begin(), d_vecG[i].end());
        vecQ.insert(vecQ.end(), d_vecQ[i].begin(), d_vecQ[i].end());
        vecCoexpression.insert(vecCoexpression.end(),
                               d_vecCoexpression[i].begin(), d_vecCoexpression[i].end());
        vecqCoexpression.insert(vecqCoexpression.end(),
                                d_vecqCoexpression[i].begin(), d_vecqCoexpression[i].end());
    }

    fprintf(stderr, "before deletion...\n");
    for (i = 0; i < datasets; i++) {
        fprintf(stderr, "%p\n", vc[i]);
    }

    for (i = 0; i < datasets; i++) {
        fprintf(stderr, "Deleting in vc[i] is called\n");
        fprintf(stderr, "%ld\n", CMeta::GetMemoryUsage());
        //getchar();
        delete vc[i];
        fprintf(stderr, "Deleting in vc[i] is finished\n");
        fprintf(stderr, "%ld\n", CMeta::GetMemoryUsage());
        //getchar();
    }

    delete[] vc;

    if (CSeekNetwork::Send(new_fd, sizeD) != 0) {
        fprintf(stderr, "Error sending messages\n");
    }

    if (outputExpression) {
        if (CSeekNetwork::Send(new_fd, vecG) != 0) {
            fprintf(stderr, "Error sending messages\n");
        }
    }

    if (outputQueryExpression) {
        if (CSeekNetwork::Send(new_fd, vecQ) != 0) {
            fprintf(stderr, "Error sending messages\n");
        }
    }

    if (outputCoexpression) {
        if (CSeekNetwork::Send(new_fd, vecCoexpression) != 0) {
            fprintf(stderr, "Error sending messages\n");
        }
    }

    if (outputQueryCoexpression) {
        if (CSeekNetwork::Send(new_fd, vecqCoexpression) != 0) {
            fprintf(stderr, "Error sending messages\n");
        }
    }

    /*for(i=0; i<NUM_DSET_MEMORY; i++){
        loaded[i] = 0;
    }*/

    pthread_mutex_lock(&mutexGet);
    THREAD_OCCUPIED[tid] = 0;
    pthread_mutex_unlock(&mutexGet);

    int ret = 0;
    close(new_fd);
    pthread_exit((void *) ret);
    //return void;
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

    PORT = sArgs.port_arg;

    CSeekDBSetting *dbSetting = new CSeekDBSetting(
            "NA", //default gvar arg, argument not needed for PCLServer
            sArgs.sinfo_arg, sArgs.platform_arg, sArgs.prep_arg,
            ".",  //default DB arg, argument not needed for PCLServer
            sArgs.gene_arg, "NA", sArgs.quant_arg, sArgs.dset_arg, "NA",
            21702 //default num_db arg, argument not needed for PCLServer
    );

    //vector<CSeekDBSetting*> cc;
    cc.push_back(dbSetting);

    pcl_input_dir = sArgs.input_arg;

    string add_db = sArgs.additional_db_arg;
    if (add_db != "NA") {
        ifstream ifsm;
        ifsm.open(add_db.c_str());
        const int lineSize = 1024;
        if (!ifsm.is_open()) {
            fprintf(stderr, "Error opening file %s\n", add_db.c_str());
            return false;
        }
        char acBuffer[lineSize];
        utype c_iBuffer = lineSize;
        vector <map<string, string>> parameters; //an array of CDatabase's
        while (!ifsm.eof()) {
            ifsm.getline(acBuffer, c_iBuffer - 1);
            if (acBuffer[0] == 0) break;
            acBuffer[c_iBuffer - 1] = 0;
            string strB = acBuffer;
            if (strB == "START") {
                map <string, string> p;
                while (!ifsm.eof()) {
                    ifsm.getline(acBuffer, c_iBuffer - 1);
                    if (acBuffer[0] == 0) {
                        fprintf(stderr, "Invalid line (empty)\n");
                        return 1;
                    }
                    strB = acBuffer;
                    if (strB == "END") break;
                    vector <string> tok;
                    CMeta::Tokenize(acBuffer, tok); //separator is tab
                    p[tok[0]] = tok[1];
                }
                parameters.push_back(p);
            }
        }
        ifsm.close();
        if (parameters.size() == 0) {
            fprintf(stderr, "Error, extra_db setting file must begin with START and end with END lines\n");
            return 1;
        }

        for (i = 0; i < parameters.size(); i++) {
            string sinfo_dir = "NA";
            string gvar_dir = "NA";
            string platform_dir = "NA";
            string prep_dir = "NA";
            string db_dir = "NA";
            string dset_map_file = "NA";
            string gene_map_file = "NA";
            string quant_file = "NA";
            int num_db = -1;

            if (parameters[i].find("SINFO_DIR")->second == "NA") {
                fprintf(stderr, "Please specify an sinfo directory for the extra db\n");
                return false;
            }
            sinfo_dir = parameters[i].find("SINFO_DIR")->second;
            if (parameters[i].find("GVAR_DIR") != parameters[i].end())
                gvar_dir = parameters[i].find("GVAR_DIR")->second;
            if (parameters[i].find("PREP_DIR") == parameters[i].end() ||
                parameters[i].find("PLATFORM_DIR") == parameters[i].end() ||
                parameters[i].find("DB_DIR") == parameters[i].end() ||
                parameters[i].find("DSET_MAP_FILE") == parameters[i].end() ||
                parameters[i].find("GENE_MAP_FILE") == parameters[i].end() ||
                parameters[i].find("QUANT_FILE") == parameters[i].end() ||
                parameters[i].find("NUMBER_OF_DB") == parameters[i].end()) {
                fprintf(stderr, "Some arguments are missing. Please make sure the following are provided:\n");
                fprintf(stderr, "PREP_DIR, DB_DIR, DSET_MAP_FILE, GENE_MAP_FILE, QUANT_FILE, NUMBER_OF_DB\n");
                return false;
            }

            platform_dir = parameters[i].find("PLATFORM_DIR")->second;
            db_dir = parameters[i].find("DB_DIR")->second;
            prep_dir = parameters[i].find("PREP_DIR")->second;
            dset_map_file = parameters[i].find("DSET_MAP_FILE")->second;
            gene_map_file = parameters[i].find("GENE_MAP_FILE")->second;
            quant_file = parameters[i].find("QUANT_FILE")->second;
            num_db = atoi(parameters[i].find("NUMBER_OF_DB")->second.c_str());

            CSeekDBSetting *dbSetting2 = new CSeekDBSetting(gvar_dir, sinfo_dir,
                                                            platform_dir, prep_dir, db_dir, 
                                                            gene_map_file, "NA", quant_file,
                                                            dset_map_file, "NA",
                                                            num_db);
            cc.push_back(dbSetting2);
        }
    }

    if (!CSeekTools::ReadListTwoColumns(sArgs.gene_arg, vecstrGeneID, vecstrGenes))
        return false;
    for (i = 0; i < vecstrGenes.size(); i++)
        mapstrintGene[vecstrGenes[i]] = (int) i;

    for (i = 0; i < cc.size(); i++) {
        vector <string> vD, vDP;
        if (!CSeekTools::ReadListTwoColumns(cc[i]->datasetFile, vD, vDP))
            return false;
        for (j = 0; j < vD.size(); j++) {
            vecstrDatasets.push_back(vD[j]);
            vecstrDP.push_back(vDP[j]);
            mapstrintDatasetDB[vD[j]] = (int) i;
        }

        SeekPlatforms platforms;
        platforms.loadPlatformDataFromFiles(cc[i]->platDir);
        vector<CSeekPlatform> &vpx = platforms.getCSeekPlatforms();
        map<string, utype> &mP = platforms.getPlatformMap();

        int cur = vp.size();
        for (map<string, utype>::iterator it = mP.begin();
             it != mP.end(); it++) {
            mapstriPlatform[it->first] = it->second + cur;
        }
        vp.resize(cur + vpx.size());
        for (j = 0; j < vpx.size(); j++)
            vp[cur + j].Copy(vpx[j]);
    }

    for (i = 0; i < vecstrDatasets.size(); i++) {
        mapstrstrDatasetPlatform[vecstrDatasets[i]] = vecstrDP[i];
        mapstrintDataset[vecstrDatasets[i]] = i;
    }

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

    printf("server: waiting for connections...\n");
    struct thread_data thread_arg[NUM_THREADS];
    pthread_t th[NUM_THREADS];
    pthread_attr_t attr[NUM_THREADS];
    int d = 0;

    pthread_mutex_init(&mutexGet, NULL);

    fprintf(stderr, "Start init.\n");

    CPCL **pcl = new CPCL *[NUM_DSET_MEMORY];

    CPCL ***vcx = new CPCL **[NUM_THREADS];

    loaded = new char[NUM_DSET_MEMORY];

    available.clear();
    for (i = 0; i < NUM_DSET_MEMORY; i++) {
        available.push_back(i);
        loaded[i] = 0;
        pcl[i] = new CPCL();
    }

    fprintf(stderr, "Finished initializations.\n");
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
            if (THREAD_OCCUPIED[d] == 0) break;
        }
        if (d == NUM_THREADS) {
            close(new_fd);
            continue;
        }
        THREAD_OCCUPIED[d] = 1;
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

        thread_arg[d].threadid = d;
        thread_arg[d].new_fd = new_fd;
        thread_arg[d].geneName = geneName;
        thread_arg[d].datasetName = dsetName;
        thread_arg[d].queryName = queryName;
        thread_arg[d].outputNormalized = outputNormalized;
        thread_arg[d].outputCoexpression = outputCoexpression;
        thread_arg[d].outputExpression = outputExpression;
        thread_arg[d].outputQueryExpression = outputQueryExpression;
        thread_arg[d].outputQueryCoexpression = outputQueryCoexpression;
        thread_arg[d].rbp_p = rbp_p;
        thread_arg[d].pcl = pcl;
        thread_arg[d].vc = vcx[d];

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
