#include <iostream>
#include <cassert>
#include <filesystem>
#include <boost/algorithm/string/join.hpp>
#include "SeekInterface.h"
#include "seekcentral.h"
#include "seekerror.h"
#include "PclQuery.h"
#include "gen-cpp/seek_rpc_constants.h"

using namespace std;
using namespace SeekRPC;

string getLogMessages(ThreadSafeQueue<string> &messageLog);

SeekInterface::SeekInterface(vector<string> &configFiles,
                             uint32_t maxConcurreny,
                             uint32_t taskTimeoutSec) :
                                querySemaphore(maxConcurreny, maxConcurreny),
                                maxTaskTimeSec(taskTimeoutSec)
{
    getConfigs(configFiles, this->speciesConfigs);

    for( auto const& [speciesName, config] : this->speciesConfigs ) {
        // Initialize a seekCentral instance for each species.
        //  C++ automatically default-inits the mapped CSeekCentral() values on first reference
        try {
            cout << "Initialize " << speciesName << endl;
            this->speciesSeekCentrals[speciesName].InitializeFromSeekConfig(config);
            this->speciesPclCache.emplace(speciesName, config.pclCacheSize);
        } catch(exception &err) {
            throw_with_nested(config_error(FILELINE + "Error initializing CSeekCentral for species " + speciesName));
        }
    }

    // start the task cleaner thread
    this->_cleanerThread = thread(&SeekInterface::runCleanTasksThread, this, taskTimeoutSec);
}

void SeekInterface::seekQuery(const SeekQueryArgs &query, SeekResult &result)
{
    // spin off a thread to run this query but wait immediately for it
    int64_t taskId = this->seekQueryAsync(query);
    this->getSeekResult(taskId, true, result);
    return;
}

int64_t SeekInterface::seekQueryAsync(const SeekQueryArgs &query)
{
    TaskInfoPtrS task = make_shared<TaskInfo>();
    task->queryType = QueryType::Seek;
    task->seekQuery = query;
    return commonAsync(task);
}




void SeekInterface::pclQuery(const PclQueryArgs &query, PclResult &result)
{
    // spin off a thread to run this query but wait immediately for it
    int64_t taskId = this->pclQueryAsync(query);
    this->getPclResult(taskId, true, result);
    return;
}

int64_t SeekInterface::pclQueryAsync(const PclQueryArgs &query)
{
    TaskInfoPtrS task = make_shared<TaskInfo>();
    task->queryType = QueryType::Pcl;
    task->pclQuery = query;
    return commonAsync(task);
}

bool SeekInterface::isQueryComplete(int64_t task_id) {
    /* lookup the task */
    TaskInfoPtrS task = this->getTask(task_id);
    if (task == nullptr || task->isComplete) {
        return true;
    }
    return false;
  }

int32_t SeekInterface::getRpcVersion()
{
    return g_seek_rpc_constants.RPCVersion;
}

string SeekInterface::getProgressMessage(int64_t task_id)
{
    /* lookup the task */
    TaskInfoPtrS task = this->getTask(task_id);
    if (task == nullptr) {
        return "";
    }
    return getLogMessages(task->messageLog);
}

int32_t SeekInterface::ping()
{
    printf("ping\n");
    return 1;
}

int32_t SeekInterface::pvalueGenes()
{
    printf("pvalueGenes\n");
    return 0;
}

int32_t SeekInterface::pvalueDatasets()
{
    printf("pvalueDatasets\n");
    return 0;
}


int64_t SeekInterface::commonAsync(TaskInfoPtrS task)
{
    /* create a task and thread */
    /* Use lock_guard to lock the mutex and guarantee
        *   that it unlocks when the scope exits
        * Note: lock_guard or unique_lock could be used here.
        *   Since this is a shared_mutex I'll use unique_lock
        * and shared_lock operations.
        */
    unique_lock mlock(this->taskMapMutex);

    /* get next task_id using atomic increment counter */
    int64_t task_id = this->next_task_id++;
    task->taskId = task_id;
    if (this->taskMap.count(task_id) != 0) {
        /* taskMap for this id already exists
            * this should never happen, program logic error
            */
        throw state_error("taskMap already contains task_id " + to_string(task_id));
    } else {
        /* add a new task to the map and start the thread */
        task->timestamp = time(0);
        task->_thread = make_unique<thread>(&SeekInterface::runQueryThread, this, task);
        this->taskMap.emplace(task_id, move(task));
    }
    return task_id;
    /* mutex automatically unlocks when scope exits */
}

void SeekInterface::runQueryThread(TaskInfoPtrS task) {
    // block with completionFlag lock_guard
    {
        /* A lock_guard will automatically set the completionFlag
         *  to true when the scope exits. This will allow other threads
         *  such as the cleanerThread to know the execution has completed.
         */
        BoolFlag completionFlag(task->isComplete);
        lock_guard<BoolFlag> flag_lock(completionFlag);

        try {
            /* Wait on semaphore if max queries already running,
            *   the lock_guard will automatically free the semaphore
            *   resource when the scope exits.
            */
            lock_guard<Semaphore> sem_lock(this->querySemaphore);
            /* Run the query */
            if (task->queryType == QueryType::Seek) {
                this->seekQueryCommon(task->seekQuery, task->seekResult, task->messageLog);
            } else if (task->queryType == QueryType::Pcl) {
                this->pclQueryCommon(task->pclQuery, task->pclResult);
            }
        } catch (named_error &err) {
            string trace = print_exception_stack(err);
            task->seekResult.success = false;
            task->seekResult.status = QueryStatus::Error;
            task->seekResult.statusMsg = trace;
            task->seekResult.__isset.status = true;
            task->seekResult.__isset.statusMsg = true;
            /* (for future) can use exception_ptr to return the exception to main thread */
        }
    }
    assert(task->isComplete = true); // set by flag lock_guard
    return;
}

void SeekInterface::seekQueryCommon(const SeekQueryArgs &query, SeekResult &result, ThreadSafeQueue<string>  &log) {
    const SeekQueryParams &params = query.parameters;

    if (this->speciesSeekCentrals.find(query.species) == this->speciesSeekCentrals.end()) {
        // no matching species initialized
        throw query_error(FILELINE + "Invalid species name: " + query.species);
    }
    CSeekCentral &speciesSC = this->speciesSeekCentrals[query.species];

    bool bSubtractGeneAvg = false;
    bool bNormPlatform = false;
    enum CSeekDataset::DistanceMeasure eDM;
    if (params.distanceMeasure == DistanceMeasure::Correlation) {
        eDM = CSeekDataset::CORRELATION;
    } else if (params.distanceMeasure == DistanceMeasure::ZScore) {
        eDM = CSeekDataset::Z_SCORE;
    } else if (params.distanceMeasure == DistanceMeasure::ZScoreHubbinessCorrected) {
        eDM = CSeekDataset::Z_SCORE;
        bSubtractGeneAvg = true;
        bNormPlatform = true;
    } else {
        throw argument_error(FILELINE + "Unknown distance measure: " + to_string(params.distanceMeasure));
    }

    // Create the output directory if needed
    filesystem::create_directory(query.outputDir);

    vector<string> queryGenes(query.genes);
    if (params.useGeneSymbols == true) {
        // convert query genes from sybmol to entrez
        try {
            speciesSC.convertGenesSymbolToEntrez(query.genes, queryGenes);
        } catch(exception &err) {
            throw_with_nested(query_error(FILELINE + "Convert query genes from symbol to entrez id"));
        }
    }

    // seekcentral InitializeQuery expects a string with genes delimited by " " and multiple queries separated by "|"
    std::string joinedGenes = boost::algorithm::join(queryGenes, " ");
    std::string joinedDatasets;
    if (query.datasets.size() == 0) {
        // if query.datasets is empty then use all datasets, specified by "NA" string
        joinedDatasets = "NA";
    } else {
        joinedDatasets = boost::algorithm::join(query.datasets, " ");
    }

    int networkConnection = 0;  // no direct network connection to client

    CSeekCentral querySC;
    querySC.setUsingRPC(true, log);
    bool res = querySC.InitializeQuery(query.outputDir,
                                       joinedGenes,
                                       joinedDatasets,
                                       &speciesSC,
                                       networkConnection,
                                       params.minQueryGenesFraction,
                                       params.minGenomeFraction,
                                       eDM,
                                       bSubtractGeneAvg,
                                       bNormPlatform,
                                       params.useNegativeCorrelation,
                                       params.checkDatasetSize);
    if (res == false) {
        result.success = false;
        result.status = QueryStatus::Error;
        result.statusMsg = "Initialize query failed, check database settings";
        result.__isset.status = true;
        result.__isset.statusMsg = true;
        querySC.Destruct();
        return;
    }

    querySC.setSimulateWeightFlag(params.simulateWeights);

    if (params.searchMethod == SearchMethod::EqualWeighting) {
        querySC.EqualWeightSearch();
    } else if (params.searchMethod == SearchMethod::OrderStatistics) {
        querySC.OrderStatistics();
    } else {
        const gsl_rng_type *T;
        gsl_rng *rnd;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        rnd = gsl_rng_alloc(T);
        // gsl_rng_set(rnd, 100);
        utype FOLD = 5;
        //enum PartitionMode PART_M = CUSTOM_PARTITION;
        enum CSeekQuery::PartitionMode PART_M = CSeekQuery::LEAVE_ONE_IN;
        if (params.searchMethod == SearchMethod::CVCustom){
            if (query.guideGenes.size() == 0) {
                throw request_error(FILELINE + "No guide gene set specified for CVCustom search");
            }
            vector<vector<string> > vecGuideGeneSet;
            vecGuideGeneSet.push_back(query.guideGenes);
            querySC.CVCustomSearch(vecGuideGeneSet, rnd, PART_M, FOLD, params.rbpParam);
        } else { //"RBP"
            querySC.CVSearch(rnd, PART_M, FOLD, params.rbpParam);
        }
        gsl_rng_free(rnd);
    }


    if (querySC.numQueries() > 1) {
        // Error - should only be running one query per RPC request
        // TODO - throw state error?
        assert(querySC.numQueries() == 1);
    }

    // get the genes and scores and add them to the rpc result reply
    vector<StrDoublePair> geneResults = querySC.getGeneResult(0);
    int numGenes = geneResults.size();
    for (int i=0; i<numGenes; i++) {
        SeekRPC::StringDoublePair pair;
        if (params.useGeneSymbols == true) {
            string entrez = speciesSC.entrezToSymbol(geneResults[i].key);
            pair.__set_name(entrez);
        } else {
            pair.__set_name(geneResults[i].key);
        }
        pair.__set_value(geneResults[i].val);
        result.geneScores.push_back(pair);
    }

    // get the datasets and weigts and add them to the rpc result reply
    vector<StrDoublePair> datasetResults = querySC.getDatasetResult(0);
    int numDsets = datasetResults.size();
    for (int i=0; i<numDsets; i++) {
        SeekRPC::StringDoublePair pair;
        pair.__set_name(datasetResults[i].key);
        pair.__set_value(datasetResults[i].val);
        result.datasetWeights.push_back(pair);
    }
    result.__isset.datasetWeights = true;
    result.status = QueryStatus::Complete;
    result.__isset.status = true;
    result.success = true;
    result.statusMsg = getLogMessages(log);
    result.__isset.statusMsg = true;

    querySC.Destruct();
}


void SeekInterface::getSeekResult(int64_t task_id, bool block, SeekResult &result)
{
    /* lookup the task */
    TaskInfoPtrS task = this->getTask(task_id);
    if (task == nullptr) {
        result.success = false;
        result.status = QueryStatus::Error;
        result.statusMsg = "Task not found or timed out: task_id: " + to_string(task_id);
        result.__isset.status = true;
        result.__isset.statusMsg = true;
        return;
    }

    if (!block && !task->isComplete) {
        result.success = false;
        result.status = QueryStatus::Incomplete;
        result.__isset.status = true;
        return;
    }

    /* join the task thread */
    this->joinTask(*task);

    /* populate the return results */
    {
        lock_guard tlock(task->taskMutex);
        result = task->seekResult;
    }

    /* remove the task from the taskMap */
    this->removeMappedTask(task_id);

    return;
}

// TODO - not obvious how to combine getSeekResult (above) and getPclResult
void SeekInterface::getPclResult(int64_t task_id, bool block, PclResult &result)
{
    /* lookup the task */
    TaskInfoPtrS task = this->getTask(task_id);
    if (task == nullptr) {
        result.success = false;
        result.status = QueryStatus::Error;
        result.statusMsg = "Task not found or timed out: task_id: " + to_string(task_id);
        result.__isset.status = true;
        result.__isset.statusMsg = true;
        return;
    }

    if (!block && !task->isComplete) {
        result.success = false;
        result.status = QueryStatus::Incomplete;
        result.__isset.status = true;
        return;
    }

    /* join the task thread */
    this->joinTask(*task);

    /* populate the return results */
    {
        lock_guard tlock(task->taskMutex);
        result = task->pclResult;
    }

    /* remove the task from the taskMap */
    this->removeMappedTask(task_id);

    return;
}


void SeekInterface::pclQueryCommon(const PclQueryArgs &query, PclResult &result) {
    const PclSettings &settings = query.settings;

    if (this->speciesSeekCentrals.count(query.species) == 0) {
        // no matching species initialized
        throw query_error(FILELINE + "Invalid species name: " + query.species);
    }

    if (this->speciesPclCache.count(query.species) == 0) {
        // no matching species initialized
        throw query_error(FILELINE + "No cache found for species: " + query.species);
    }

    thread_data thread_arg;
    thread_arg.new_fd = -1; // used by original main() server, -1 indicates don't send
    thread_arg.isComplete = false;
    thread_arg.geneNames = query.genes;
    thread_arg.datasetNames = query.datasets;
    thread_arg.queryGeneNames = query.queryGenes;
    thread_arg.outputNormalized = settings.outputNormalized;
    thread_arg.outputExpression = settings.outputGeneExpression;
    thread_arg.outputCoexpression = settings.outputGeneCoexpression;
    thread_arg.outputQueryExpression = settings.outputQueryExpression;
    thread_arg.outputQueryCoexpression = settings.outputQueryCoexpression;
    thread_arg.rbp_p = settings.rbp;
    thread_arg.seekCentral = &this->speciesSeekCentrals[query.species];
    thread_arg.pclCache = &this->speciesPclCache.at(query.species);
    // structs to hold results if new_fd == -1 (so results not sent back within do_query)
    thread_arg.resDatasetSizes = &result.datasetSizes;
    thread_arg.resGeneExpression = &result.geneExpressions;
    thread_arg.resGeneCoexpression = &result.geneCoexpressions;
    thread_arg.resQueryExpression = &result.queryExpressions;
    thread_arg.resQueryCoexpression = &result.queryCoexpressions;

    try {
        do_query(&thread_arg);
    } catch (named_error &err) {
        string trace = print_exception_stack(err);
        result.success = false;
        result.status = QueryStatus::Error;
        result.statusMsg = trace;
        result.__isset.status = true;
        result.__isset.statusMsg = true;
        /* (for future) can use exception_ptr to return the exception to main thread */
        return;
    }

    result.success = true;
    result.status = QueryStatus::Complete;
    result.__isset.status = true;
    result.__isset.geneExpressions = settings.outputGeneExpression;
    result.__isset.geneCoexpressions = settings.outputGeneCoexpression;
    result.__isset.queryExpressions = settings.outputQueryExpression;
    result.__isset.queryCoexpressions = settings.outputQueryCoexpression;
}

TaskInfoPtrS SeekInterface::getTask(int64_t taskId) {
    shared_lock mlock(this->taskMapMutex);
    if (this->taskMap.count(taskId) > 0) {
        return this->taskMap[taskId];
    } else {
        return nullptr;
    }
}

void SeekInterface::joinTask(TaskInfo &task) {
    lock_guard tlock(task.taskMutex);
    if (task._thread->joinable()) {
        task._thread->join();
        assert(task.isComplete == true);
    }
}

int SeekInterface::removeMappedTask(int64_t taskId) {
    // remove the task from the taskMap
    unique_lock mlock(this->taskMapMutex);
    uint32_t numRemoved = this->taskMap.erase(taskId);
    assert(numRemoved == 0 || numRemoved == 1);
    return numRemoved;
}

/* This cleaner thread handles cases where the client never calls get_result.
 *   It will clean up abandoned threads and remove task from the taskMap
 */
void SeekInterface::runCleanTasksThread(uint32_t intervalSec) {
    while (true) {
        sleep(intervalSec);
        // cout << "### cleaner running ###" << endl;
        int64_t beginTaskId = 0;
        {
            /* find the oldest (lowest) taskId */
            shared_lock mlock(this->taskMapMutex);
            if (this->taskMap.size() == 0) { continue; } // empty
            auto iter = this->taskMap.lower_bound(beginTaskId);
            if (iter == this->taskMap.end()) { continue; } // empty
            beginTaskId = iter->first;
        }

        /* loop starting from oldest taskid checking if task is stale */
        for (int64_t idx = beginTaskId; idx < this->next_task_id; idx++) {
            bool isStale = this->cleanStaleTask(idx);
            if (isStale == false) {
                /* remaining threads are newer so within time window */
                break;
            }
        }
    }
}

bool SeekInterface::cleanStaleTask(int64_t task_id) {
    bool isTimedOut = false;

    TaskInfoPtrS task = this->getTask(task_id);
    if (task == nullptr) {
        // thread completed, return true to continue checking other threads
        return true;
    }
    // check elapsed time
    double elapsed_time = difftime(time(0), task->timestamp);
    // cout << "### Thread elapsed time: " << to_string(elapsed_time) << endl;
    if (elapsed_time > this->maxTaskTimeSec) {
        isTimedOut = true;
        if (task->isComplete == true) {
            // join the thread
            this->joinTask(*task);

            // remove the task from the taskMap
            uint32_t numRemoved = removeMappedTask(task_id);
            if (numRemoved > 0) {
                // cout << "### Cleaned Task: " << to_string(task_id) << endl;
            }
        }
    }
    return isTimedOut;
}

string getLogMessages(ThreadSafeQueue<string> &messageLog) {
    string joinedMessages;
    try {
        while (!messageLog.empty()) {
            string msg = messageLog.dequeue();
            joinedMessages += msg + "\n";
        }
    } catch (state_error e) {
        // catch if dequeue called on empty queue
    }
    return joinedMessages;
}

