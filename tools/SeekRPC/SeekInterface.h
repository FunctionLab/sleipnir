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
    unique_ptr<thread> _thread;
};

using TaskInfoPtrS = shared_ptr<TaskInfo>;

class SeekInterface {
  public:
    SeekInterface(vector<string> &configFiles, uint32_t maxConcurreny, uint32_t taskTimeoutSec);
    void seek_query(const SeekQuery &query, QueryResult &result);
    int64_t seek_query_async(const SeekQuery &query);
    void seek_get_result(int64_t task_id, QueryResult &result);
    string get_progress_message(int64_t task_id);
    int32_t ping();
    int32_t pvalue_genes();
    int32_t pvalue_datasets();
    int32_t pcl_data();
  private:
    void SeekQueryCommon(const SeekQuery &query, QueryResult &result);
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
