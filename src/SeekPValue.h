#ifndef SEEKPVALUE_H
#define SEEKPVALUE_H

#include <string>
#include <vector>
#include <map>
#include <pthread.h>
#include "seekcentral.h"

using namespace std;

#define NUM_THREADS 8

struct parameter {
    int size;
    double scale;
    double shape;
    double threshold;
    double portion;
    vector<double> quantile;
};

typedef struct PValueData_s {
    vector <vector<int>> randomRank;
    vector <vector<float>> randomSc;
    vector <vector<struct parameter>> dsetScore;
} PValueData;


struct pvalue_thread_data {
    int queryType; //0 - genes, 1 - datasets (which mode to turn on)
    //Section "genes"
    bool rankBased; //true - p-value on rank, false - p-value on score
    vector <string> query;
    vector<string> gene_entrezIds;
    vector<double> gene_scores;
    vector<int> gene_ranks;
    float nan = -320;
    //Section "datasets"
    vector <string> dset;
    vector<float> dset_score; //scores to test
    vector<int> dset_qsize; //number of genes for which coexpression score is calculated, for all dset
    //============================================
    int new_fd;
    bool isComplete = false;
    vector<double> *resPvalues; // resulting pvalues
    PValueData *pvalueData;
    CSeekCentral *seekCentral;
};

// Unique new variables
extern vector <vector<int>> randomRank;
extern vector <vector<float>> randomSc;
extern vector <vector<struct parameter>> dsetScore; //co-expression score for dataset i, query-size j

// Functions Definitions
void *do_pvalue_query(void *th_arg);
bool ReadParameter(const string& param_file, vector<struct parameter> &v, PValueData &pvalueData);
bool initializePvalue(CSeekCentral &seekCentral, int numRandQueries, PValueData &pvalueData);
bool loadPvalueArrays(string dirname, PValueData &pvalueData);

#endif  // SEEKPVALUE_H