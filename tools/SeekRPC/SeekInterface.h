#ifndef SEEKINTERFACE_H
#define SEEKINTERFACE_H

#include <map>
#include <thread>
#include <mutex>
#include <atomic>
#include <shared_mutex>
#include "seekcentral.h"
#include "seekhelper.h"
#include "gen-cpp/SeekRPC.h"


using namespace std;
using namespace SeekRPC;


class TaskInfo {
public:
    SeekQuery seekQuery;
    QueryResult seekResult;
    bool isComplete = false;
    time_t timestamp;
    mutex taskMutex;
    ThreadSafeQueue<string> messageLog;
    unique_ptr<thread> _thread;
};

using TaskInfoPtrS = shared_ptr<TaskInfo>;

class SeekInterface {
  public:
    SeekInterface(vector<string> &configFiles, uint32_t maxConcurreny, uint32_t taskTimeoutSec);
    void seekQuery(const SeekQuery &query, QueryResult &result);
    int64_t seekQueryAsync(const SeekQuery &query);
    bool isQueryComplete(int64_t task_id);
    void getQueryResult(int64_t task_id, bool block, QueryResult &result);
    string getProgressMessage(int64_t task_id);
    int32_t getRpcVersion();
    int32_t ping();
    int32_t pvalueGenes();
    int32_t pvalueDatasets();
    int32_t pclData();
  private:
    void SeekQueryCommon(const SeekQuery &query, QueryResult &result, ThreadSafeQueue<string> &log);
    void runSeekQueryThread(TaskInfoPtrS task);
    void runCleanTasksThread(uint32_t intervalSec);
    bool cleanStaleTask(int64_t task_id);
    TaskInfoPtrS getTask(int64_t taskId);
    void joinTask(TaskInfo &task);
    int removeMappedTask(int64_t taskId);

    // Class variables
    map<string, SeekSettings> speciesConfigs;  // speciesName --> Config
    map<string, CSeekCentral> speciesSeekCentrals; // speciesName --> SeekCentralStruct
    map<int64_t, TaskInfoPtrS> taskMap;  // task_id --> TaskInfo
    shared_mutex taskMapMutex;
    Semaphore querySemaphore;
    double maxTaskTimeSec;
    thread _cleanerThread;
    atomic_int64_t next_task_id = 1;  // first task_id
 };

#endif  // SEEKINTERFACE_H
