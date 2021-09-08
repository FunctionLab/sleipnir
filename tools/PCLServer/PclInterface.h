#ifndef PCLINTERFACE_H
#define PCLINTERFACE_H

#include <map>
#include <thread>
#include <mutex>
#include <queue>
#include <atomic>
#include <shared_mutex>
#include "PCLServer.h"
#include "seekcentral.h"
#include "seekhelper.h"
#include "gen-cpp/PclRPC.h"

using namespace std;
using namespace PclRPC;


class TaskInfo {
public:
    PclQueryArgs pclQueryArgs;
    PclResult pclResult;
    bool isComplete = false;
    time_t timestamp;
    mutex taskMutex;
    queue<string> messageLog;
    unique_ptr<thread> _thread;
};

using TaskInfoPtrS = shared_ptr<TaskInfo>;

class PclInterface {
  public:
    PclInterface(vector<string> &configFiles, uint32_t maxConcurreny, uint32_t taskTimeoutSec);
    void pclQuery(const PclQueryArgs &query, PclResult &result);
    int64_t pclQueryAsync(const PclQueryArgs &query);
    bool isQueryComplete(int64_t task_id);
    void getQueryResult(int64_t task_id, bool block, PclResult &result);
    string getProgressMessage(int64_t task_id);
    int32_t getRpcVersion();
    int32_t ping();
  private:
    void pclQueryCommon(const PclQueryArgs &query, PclResult &result);
    void runPclQueryThread(TaskInfoPtrS task);
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

#endif  // PCLINTERFACE_H