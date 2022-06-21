#ifndef PCLSERVER_H
#define PCLSERVER_H

#include <string>
#include <list>
#include "pcl.h"
#include "seekcentral.h"

using namespace std;
using namespace Sleipnir;

using PclPtrS = shared_ptr<CPCL>;

// These all became part of the CSeekCentral datastructure that is initialized
// vector<CSeekDBSetting *> cc;
// vector <string> vecstrGenes;
// vector <string> vecstrDatasets;
// vector <string> vecstrDP;
// map <string, string> mapstrstrDatasetPlatform;
// map <string, utype> mapstrintDataset;
// map<string, int> mapstrintGene;
// vector <string> vecstrPlatforms; // Part of SeekPlatforms
// vector <CSeekPlatform> vp;  // Part of SeekPlatforms
// map <string, utype> mapstriPlatform;  // Part of SeekPlatforms
// map<string, int> mapstrintDatasetDB;  // Perhaps add this to CSeekCentral

struct pcl_thread_data {
    vector <string> datasetNames;
    vector <string> geneNames;
    vector <string> queryGeneNames;
    bool isComplete = false;
    bool outputNormalized;
    bool outputCoexpression;
    bool outputQueryCoexpression;
    bool outputExpression;
    bool outputQueryExpression;
    float rbp_p;
    int new_fd;
    CSeekCentral *seekCentral;
    LRUCache <string, PclPtrS> *pclCache;
    vector <int> *resDatasetSizes;
    vector <double> *resGeneExpression;
    vector <double> *resGeneCoexpression;
    vector <double> *resQueryExpression;
    vector <double> *resQueryCoexpression;
    vector <string> *resExperimentNames;
};

void *do_pcl_query(void *th_arg);

#endif  // PCLSERVER_H
