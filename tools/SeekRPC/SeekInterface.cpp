#include <iostream>
#include <cassert>
#include <filesystem>
#include <boost/algorithm/string/join.hpp>
#include "SeekInterface.h"
#include "seekcentral.h"
#include "seekerror.h"
#include "PclQuery.h"
#include "SeekPValue.h"
#include "stdafx.h"
#include "gen-cpp/seek_rpc_constants.h"
// #include "/usr/local/Cellar/gperftools/2.9.1_1/include/gperftools/heap-profiler.h"

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
            g_CatSleipnir().info("SeekRPC: Initialize %s", speciesName.c_str());
            // Initialize the SeekCentral struct
            this->speciesSeekCentrals[speciesName].InitializeFromSeekConfig(config);
            // Initialize the PCL cache
            this->speciesPclCache.emplace(speciesName, config.pclCacheSize);
            // Initialize the Pvalue struct
            bool res;
            CSeekCentral &speciesSC = this->speciesSeekCentrals[speciesName];
            string pvalueDir = speciesSC.roAttr->m_vecDBSetting[0]->pvalueDir;
            // Try loading the pvalue metadata arrays
            res = loadPvalueArrays(pvalueDir, this->speciesPvalueData[speciesName]);
            if (res == false) {
                // Try creating the metadata from the raw random score outputs
                res = initializePvalue(speciesSC, -1, this->speciesPvalueData[speciesName]);
                if (res == false) {
                    // Disable pvalue queries for this species
                    cerr << "WARNING: PValue queries disabled for (" << speciesName;
                    cerr << "): unable to initialize data structures" << endl;
                    this->speciesPvalueData.erase(speciesName);
                }
            }
        } catch(exception &err) {
            throw_with_nested(config_error(FILELINE + "Error initializing CSeekCentral for species " + speciesName));
        }
    }

    // start the task cleaner thread
    this->_cleanerThread = thread(&SeekInterface::runCleanTasksThread, this, taskTimeoutSec);
    g_CatSleipnir().info("SeekRPC: Ready");
}

void SeekInterface::seekQuery(const SeekQueryArgs &query, SeekResult &result)
{
    // HeapProfilerStart("/tmp/heap");
    // spin off a thread to run this query but wait immediately for it
    int64_t taskId = this->seekQueryAsync(query);
    this->getSeekResult(taskId, true, result);
    // HeapProfilerDump("seekQuery");
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

void SeekInterface::pvalueGenes(const PValueGeneArgs& query, PValueResult& result) {

    if (this->speciesPvalueData.count(query.species) > 0) {
        int64_t taskId = this->pvalueGenesAsync(query);
        this->getPvalueResult(taskId, true, result);
    } else {
        result.success = false;
        result.status = QueryStatus::Error;
        result.statusMsg = "Pvalue data not initialized properly for species (" +
            query.species + "), check pvalue directory";
        result.__isset.status = true;
        result.__isset.statusMsg = true;
    }
    return;
}

int64_t SeekInterface::pvalueGenesAsync(const PValueGeneArgs& query) {
    TaskInfoPtrS task = make_shared<TaskInfo>();
    task->queryType = QueryType::Pvalue;
    task->pvalueGeneQuery = query;
    return commonAsync(task);
}

void SeekInterface::pvalueDatasets(const PValueDatasetArgs& query, PValueResult& result) {
    string errMsg = "pvalueDatasets query not implemented in RPC server yet";
    result.success = false;
    result.status = QueryStatus::Error;
    result.statusMsg = errMsg;
    result.__isset.status = true;
    result.__isset.statusMsg = true;
}

double SeekInterface::getRpcVersion()
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
        AtomicBoolFlag completionFlag(task->isComplete);
        lock_guard<AtomicBoolFlag> flag_lock(completionFlag);

        // For Logging
        string species = "NotSet";
        string firstQueryGene = "NotSet";
        int numQueryGenes = 0;
        ResultUnion_t resUnion;
        switch (task->queryType) {
            case QueryType::Seek:
                species = task->seekQuery.species;
                numQueryGenes = task->seekQuery.genes.size();
                if (numQueryGenes > 0) {
                    firstQueryGene = task->seekQuery.genes[0];
                }
                resUnion.seekResult = &task->seekResult;
                break;
            case QueryType::Pcl:
                species = task->pclQuery.species;
                numQueryGenes = task->pclQuery.genes.size();
                if (numQueryGenes > 0) {
                    firstQueryGene = task->pclQuery.genes[0];
                }
                resUnion.pclResult = &task->pclResult;
                break;
            case QueryType::Pvalue:
                species = task->pvalueGeneQuery.species;
                numQueryGenes = task->pvalueGeneQuery.genes.size();
                if (numQueryGenes > 0) {
                    firstQueryGene = task->pvalueGeneQuery.genes[0];
                }
                resUnion.pvalueResult = &task->pvalueResult;
                break;
        }
        g_CatSleipnir().info("%s(%ld): species %s: num_genes %d, query gene[0]: %s",
                            queryTypeName(task->queryType), task->taskId,
                            species.c_str(), numQueryGenes, firstQueryGene.c_str());

        try {
            /* Wait on semaphore if max queries already running,
            *   the lock_guard will automatically free the semaphore
            *   resource when the scope exits.
            */
            lock_guard<Semaphore> sem_lock(this->querySemaphore);

            /* update the start timestamp */
            task->timestamp = time(0);

            /* Run the query */
            if (task->queryType == QueryType::Seek) {
                this->seekQueryCommon(task->seekQuery, task->seekResult, task->messageLog);
            } else if (task->queryType == QueryType::Pcl) {
                this->pclQueryCommon(task->pclQuery, task->pclResult);
            } else if (task->queryType == QueryType::Pvalue) {
                this->pvalueGenesCommon(task->pvalueGeneQuery, task->pvalueResult);
            }
        } catch (named_error &err) {
            string traceMsg = print_exception_stack(err);
            setResultError(task->queryType, QueryStatus::Error, traceMsg, resUnion);
            /* (for future) can use exception_ptr to return the exception to main thread */
        }
        g_CatSleipnir().info("%s(%ld): Completed", queryTypeName(task->queryType), task->taskId);
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
    if (query.outputDir.length() > 0) {
        if (!filesystem::exists(query.outputDir)) {
            filesystem::create_directory(query.outputDir);
        }
    }

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
    float percentDatasetCoverage = 0.5;
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
                                       params.checkDatasetSize,
                                       percentDatasetCoverage);
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
        gsl_rng_set(rnd, 0);
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

void setResultError(QueryType qtype, QueryStatus::type status, string errMsg, ResultUnion_t res) {
    switch (qtype) {
        case QueryType::Seek:
            res.seekResult->success = false;
            res.seekResult->status = status;
            res.seekResult->statusMsg = errMsg;
            res.seekResult->__isset.status = true;
            res.seekResult->__isset.statusMsg = true;
            break;
        case QueryType::Pcl:
            res.pclResult->success = false;
            res.pclResult->status = status;
            res.pclResult->statusMsg = errMsg;
            res.pclResult->__isset.status = true;
            res.pclResult->__isset.statusMsg = true;
            break;
        case QueryType::Pvalue:
            res.pvalueResult->success = false;
            res.pvalueResult->status = status;
            res.pvalueResult->statusMsg = errMsg;
            res.pvalueResult->__isset.status = true;
            res.pvalueResult->__isset.statusMsg = true;
            break;
    }
}

void SeekInterface::getResultCommon(int64_t task_id, QueryType qtype, bool block, ResultUnion_t res)
{
    /* lookup the task */
    TaskInfoPtrS task = this->getTask(task_id);
    if (task == nullptr) {
        string errMsg = "Task not found or timed out: task_id: " + to_string(task_id);
        setResultError(qtype, QueryStatus::Error, errMsg, res);
        return;
    }

    if (task->queryType != qtype) {
        string errMsg = "Wrong getResult() function called for this TaskID's query type";
        setResultError(qtype, QueryStatus::Error, errMsg, res);
        return;
    }

    if (!block && !task->isComplete) {
        setResultError(qtype, QueryStatus::Incomplete, "", res);
        return;
    }

    /* join the task thread */
    this->joinTask(*task);

    /* populate the return results */
    {
        lock_guard tlock(task->taskMutex);
        switch (qtype) {
            case QueryType::Seek:
                *res.seekResult = task->seekResult;
                break;
            case QueryType::Pcl:
                *res.pclResult = task->pclResult;
                break;
            case QueryType::Pvalue:
                *res.pvalueResult = task->pvalueResult;
                break;
        }
    }

    /* remove the task from the taskMap */
    this->removeMappedTask(task->taskId);

    return;
}

void SeekInterface::getSeekResult(int64_t task_id, bool block, SeekResult &result)
{
    ResultUnion_t resUnion;
    resUnion.seekResult = &result;
    this->getResultCommon(task_id, QueryType::Seek, block, resUnion);

    return;
}

void SeekInterface::getPclResult(int64_t task_id, bool block, PclResult &result)
{
    ResultUnion_t resUnion;
    resUnion.pclResult = &result;
    this->getResultCommon(task_id, QueryType::Pcl, block, resUnion);

    return;
}

void SeekInterface::getPvalueResult(int64_t task_id, bool block, PValueResult &result)
{
    ResultUnion_t resUnion;
    resUnion.pvalueResult = &result;
    this->getResultCommon(task_id, QueryType::Pvalue, block,  resUnion);

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

    pcl_thread_data thread_arg;
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
        do_pcl_query(&thread_arg);
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


void SeekInterface::pvalueGenesCommon(const PValueGeneArgs &query, PValueResult &result) {
    if (this->speciesSeekCentrals.count(query.species) == 0) {
        // no matching species initialized
        throw query_error(FILELINE + "Invalid species name: " + query.species);
    }

    if (this->speciesPvalueData.count(query.species) == 0) {
        // no matching species initialized
        throw query_error(FILELINE + "No cache found for species: " + query.species);
    }

    int numReqGenes = query.genes.size();
    if (query.useRank == true) {
        if (numReqGenes > 0 && numReqGenes != query.geneRanks.size()) {
            throw request_error("Error: PValue rank-based query should have "
                                "equal num genes and num ranks provided. Perhaps "
                                "the useRanks flag was not set.");
        }
     } else {
        if (numReqGenes > 0 && numReqGenes != query.geneScores.size()) {
            throw request_error("Error: PValue score-based query should have "
                                "equal num genes and num scores provided");
        }
     }

    pvalue_thread_data thread_arg;
    thread_arg.new_fd = -1; // used by original main() server, -1 indicates don't send
    thread_arg.isComplete = false;
    thread_arg.useGeneMapOrder = (numReqGenes == 0); // true if no gene supplied
    thread_arg.queryType = 0; // genes
    thread_arg.gene_entrezIds = query.genes;
    thread_arg.gene_scores = query.geneScores;
    thread_arg.gene_ranks = query.geneRanks;
    thread_arg.rankBased = query.useRank;
    thread_arg.seekCentral = &this->speciesSeekCentrals[query.species];
    thread_arg.pvalueData = &this->speciesPvalueData.at(query.species);
    // structs to hold results if new_fd == -1 (so results not sent back within do_query)
    thread_arg.resPvalues = &result.pvalues;

    try {
        do_pvalue_query(&thread_arg);
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
    result.__isset.pvalues = true;
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

int SeekInterface::numTasksOutstanding() {
    shared_lock mlock(this->taskMapMutex);
    int numTasks = this->taskMap.size();
    return numTasks;
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
