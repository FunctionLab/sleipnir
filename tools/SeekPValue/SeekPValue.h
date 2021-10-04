#ifndef SEEKPVALUE_H
#define SEEKPVALUE_H

#include <string>
#include <vector>
#include <map>
#include <pthread.h>

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
    bool isComplete = false;
    CSeekCentral *seekCentral;
};

// Unique new variables
extern vector <vector<int>> randomRank;
extern vector <vector<float>> randomSc;
extern vector <vector<struct parameter>> dsetScore; //co-expression score for dataset i, query-size j

// Functions Definitions
void *do_pvalue_query(void *th_arg);
bool ReadParameter(const string& param_file, vector<struct parameter> &v);
bool initializePvalue(CSeekCentral &seekCentral);

#endif  // SEEKPVALUE_H