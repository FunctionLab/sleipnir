#ifndef SEEKINTERFACE_H
#define SEEKINTERFACE_H

#include <map>
#include <thread>
#include <mutex>
#include <queue>
#include <atomic>
#include <shared_mutex>
#include "seekcentral.h"
#include "seekhelper.h"
#include "PclQuery.h"
#include "gen-cpp/SeekRPC.h"


using namespace std;
using namespace SeekRPC;

enum QueryType {
  Seek,
  Pcl,
  Pvalue
};

class TaskInfo {
public:
    SeekQueryArgs seekQuery;
    SeekResult seekResult;
    PclQueryArgs pclQuery;
    PclResult pclResult;
    QueryType queryType;
    int64_t taskId;
    bool isComplete = false;
    time_t timestamp;
    mutex taskMutex;
    queue<string> messageLog;
    unique_ptr<thread> _thread;
};

using TaskInfoPtrS = shared_ptr<TaskInfo>;

class SeekInterface {
  public:
    SeekInterface(vector<string> &configFiles, uint32_t maxConcurreny, uint32_t taskTimeoutSec);
    void seekQuery(const SeekQueryArgs &query, SeekResult &result);
    int64_t seekQueryAsync(const SeekQueryArgs &query);
    bool isQueryComplete(int64_t task_id);
    void getSeekResult(int64_t task_id, bool block, SeekResult &result);
    string getProgressMessage(int64_t task_id);
    int32_t getRpcVersion();
    int32_t ping();
    int32_t pvalueGenes();
    int32_t pvalueDatasets();
    void pclQuery(const PclQueryArgs &query, PclResult &result);
    int64_t pclQueryAsync(const PclQueryArgs &query);
    void getPclResult(int64_t task_id, bool block, PclResult &result);
  private:
    void seekQueryCommon(const SeekQueryArgs &query, SeekResult &result, queue<string> &log);
    void pclQueryCommon(const PclQueryArgs &query, PclResult &result);
    int64_t commonAsync(TaskInfoPtrS task);
    void runQueryThread(TaskInfoPtrS task);
    void runCleanTasksThread(uint32_t intervalSec);
    bool cleanStaleTask(int64_t task_id);
    TaskInfoPtrS getTask(int64_t taskId);
    void joinTask(TaskInfo &task);
    int removeMappedTask(int64_t taskId);

    // Class variables
    map<string, SeekSettings> speciesConfigs;  // speciesName --> Config
    map<string, CSeekCentral> speciesSeekCentrals; // speciesName --> SeekCentralStruct
    map<string, LRUCache <string, PclPtrS>> speciesPclCache; // speciesName --> PclCache
    map<int64_t, TaskInfoPtrS> taskMap;  // task_id --> TaskInfo
    shared_mutex taskMapMutex;
    Semaphore querySemaphore;
    double maxTaskTimeSec;
    thread _cleanerThread;
    atomic_int64_t next_task_id = 1;  // first task_id
 };

#endif  // SEEKINTERFACE_H
