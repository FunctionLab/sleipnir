#ifndef SEEKINTERFACE_H
#define SEEKINTERFACE_H

#include "seekcentral.h"
#include "seekhelper.h"
#include "gen-cpp/SeekRPC.h"


using namespace std;
using namespace SeekRPC;

class SeekInterface {
  public:
    SeekInterface(vector<string> &configFiles);
    void seek_query(const SeekQuery &query, QueryResult &result);
    int32_t seek_query_async(const SeekQuery &query);
    void seek_get_result(int32_t task_id, QueryResult &result);
    string get_progress_message(int32_t task_id);
    int32_t ping();
    int32_t pvalue_genes();
    int32_t pvalue_datasets();
    int32_t pcl_data();
  private:
    void SeekQueryCommon(const SeekQuery &query, QueryResult &result);

    // Class variables
    map<string, SeekSettings> speciesConfigs;
    map<string, CSeekCentral> speciesSeekCentrals;
};

#endif  // SEEKINTERFACE_H