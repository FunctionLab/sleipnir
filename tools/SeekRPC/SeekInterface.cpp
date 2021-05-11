#include <iostream>
#include <cassert>
#include <filesystem>
#include <boost/algorithm/string/join.hpp>
#include "SeekInterface.h"
#include "seekcentral.h"
#include "seekerror.h"

using namespace std;
using namespace SeekRPC;


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
        } catch(exception &err) {
            throw_with_nested(config_error(FILELINE + "Error initializing CSeekCentral for species " + speciesName));
        }
    }

    // start the task cleaner thread
    this->_cleanerThread = thread(&SeekInterface::runCleanTasksThread, this, taskTimeoutSec);
}

void SeekInterface::seek_query(const SeekQuery &query, QueryResult &result)
{
    // spin off a thread to run this query but wait immediately for it
    int64_t taskId = this->seek_query_async(query);
    this->seek_get_result(taskId, result);
    return;
}

int64_t SeekInterface::seek_query_async(const SeekQuery &query)
{
    /* get next task_id using atomic increment counter */
    int64_t task_id = this->next_task_id++;
    /* if counter ever loops then throw exception, allows 4 billion queries */
    assert(task_id > 0);

    /* create a task and thread */
    {
        /* Use lock_guard to lock the mutex and guarantee 
         *   that it unlocks when the scope exits
         */
        unique_lock mlock(this->taskMapMutex);
        if (this->taskMap.count(task_id) != 0) {
            /* taskMap for this id already exists
             * this should never happen, program logic error
             */
            throw state_error("taskMap already contains task_id " + to_string(task_id));
        } else {
            /* add a new task to the map and start the thread */
            TaskInfoPtrS task = make_shared<TaskInfo>();
            task->seekQuery = query;
            task->timestamp = time(0);
            task->_thread = make_unique<thread>(&SeekInterface::runSeekQueryThread, this, task);
            this->taskMap.emplace(task_id, move(task));
        }
        /* mutex automatically unlocks when scope exits */
    }

    return task_id;
}


void SeekInterface::seek_get_result(int64_t task_id, QueryResult &result)
{
    /* lookup the task */
    TaskInfoPtrS task = this->getTask(task_id);
    if (task == nullptr) {
        result.success = false;
        result.statusMsg = "Task not found or timed out: task_id: " + to_string(task_id);
        result.__isset.statusMsg = true;
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

string SeekInterface::get_progress_message(int64_t task_id)
{
    printf("get_progress_message\n");
    return "Status Message";
}

int32_t SeekInterface::ping()
{
    printf("ping\n");
    return 1;
}

int32_t SeekInterface::pvalue_genes()
{
    printf("pvalue_genes\n");
    return 0;
}

int32_t SeekInterface::pvalue_datasets()
{
    printf("pvalue_datasets\n");
    return 0;
}

int32_t SeekInterface::pcl_data()
{
    printf("pcl_data\n");
    return 0;
}

void SeekInterface::runSeekQueryThread(TaskInfoPtrS task) {
    BoolFlag completionFlag(task->isComplete);
    try {
         /* Wait on semaphore if max queries already running,
          *   the lock_guard will automatically free the semaphore
          *   resource when the scope exits.
          */
        lock_guard<Semaphore> sem_lock(this->querySemaphore);
        /* A lock_guard will automatically set the completionFlag
         *  to true when the scope exits. This will allow other threads
         *  such as the cleanerThread to know the execution has completed.
         */
        lock_guard<BoolFlag> flag_lock(completionFlag);
        /* Run the query */
        this->SeekQueryCommon(task->seekQuery, task->seekResult);
    } catch (named_error &err) {
        string trace = print_exception_stack(err);
        task->seekResult.success = false;
        task->seekResult.statusMsg = trace;
        task->seekResult.__isset.statusMsg = true;
        /* (for future) can use exception_ptr to return the exception to main thread */
    }
    assert(task->isComplete = true); // set by flag lock_guard
    return;
}

void SeekInterface::SeekQueryCommon(const SeekQuery &query, QueryResult &result) {
    const QueryParams &params = query.parameters;

    if (this->speciesSeekCentrals.find(query.species) == this->speciesSeekCentrals.end()) {
        // no matching species initialized
        throw query_error(FILELINE + "Invalid species name: " + query.species);
    }
    CSeekCentral &speciesSC = this->speciesSeekCentrals[query.species];

    bool bSubtractGeneAvg = false;
    bool bNormPlatform = false;
    enum CSeekDataset::DistanceMeasure eDM;
    if (params.distance_measure == "Correlation") {
        eDM = CSeekDataset::CORRELATION;
    } else if (params.distance_measure == "Zscore") {
        eDM = CSeekDataset::Z_SCORE;
    } else if (params.distance_measure == "ZscoreHubbinessCorrected") {
        eDM = CSeekDataset::Z_SCORE;
        bSubtractGeneAvg = true;
        bNormPlatform = true;
    } else {
        throw argument_error(FILELINE + "Unknown distance measure: " + params.distance_measure);
    }

    // Create the output directory if needed
    filesystem::create_directory(query.outputDir);

    vector<string> queryGenes(query.genes);
    if (params.use_gene_symbols == true) {
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
    querySC.setUsingRPC(true);
    bool res = querySC.InitializeQuery(query.outputDir,
                                       joinedGenes,
                                       joinedDatasets,
                                       &speciesSC,
                                       networkConnection,
                                       params.min_query_genes_fraction,
                                       params.min_genome_fraction,
                                       eDM,
                                       bSubtractGeneAvg,
                                       bNormPlatform,
                                       params.useNegativeCorrelation,
                                       params.check_dataset_size);
    if (res == false) {
        result.success = false;
        result.statusMsg = "Initialize query failed, check database settings";
        result.__isset.statusMsg = true;
        querySC.Destruct();
        return;
    }

    if (params.search_method == "EqualWeighting") {
        querySC.EqualWeightSearch();
    } else if (params.search_method == "OrderStatistics") {
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
        if(params.search_method == "CVCUSTOM"){
            if (query.guideGenes.size() == 0) {
                throw request_error(FILELINE + "No guide gene set specified for CVCustom search");
            }
            vector<vector<string> > vecGuideGeneSet;
            vecGuideGeneSet.push_back(query.guideGenes);
            querySC.CVCustomSearch(vecGuideGeneSet, rnd, PART_M, FOLD, params.rbp_param);
        } else { //"RBP"
            querySC.CVSearch(rnd, PART_M, FOLD, params.rbp_param);
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
        if (params.use_gene_symbols == true) {
            string entrez = speciesSC.entrezToSymbol(geneResults[i].key);
            pair.__set_name(entrez);
        } else {
            pair.__set_name(geneResults[i].key);
        }
        pair.__set_value(geneResults[i].val);
        result.gene_scores.push_back(pair);
    }

    // get the datasets and weigts and add them to the rpc result reply
    vector<StrDoublePair> datasetResults = querySC.getDatasetResult(0);
    int numDsets = datasetResults.size();
    for (int i=0; i<numDsets; i++) {
        SeekRPC::StringDoublePair pair;
        pair.__set_name(datasetResults[i].key);
        pair.__set_value(datasetResults[i].val);
        result.dataset_weights.push_back(pair);
    }
    result.__isset.dataset_weights = true;
    result.success = true;

    querySC.Destruct();
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
