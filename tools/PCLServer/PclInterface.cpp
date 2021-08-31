#include <iostream>
#include <cassert>
#include <filesystem>
#include "PclInterface.h"
#include "seekerror.h"
#include "gen-cpp/pcl_rpc_constants.h"

using namespace std;
using namespace PclRPC;

string getLogMessages(queue<string> &messageLog);

PclInterface::PclInterface(vector<string> &configFiles, 
                           uint32_t maxConcurreny,
                           uint32_t taskTimeoutSec) :
                              querySemaphore(maxConcurreny, maxConcurreny),
                              maxTaskTimeSec(taskTimeoutSec)
{
    getConfigs(configFiles, this->speciesConfigs);

    for( auto const& [speciesName, config] : this->speciesConfigs ) {
        // Initialize a seekCentral instance for each species.
        // C++ automatically default-inits the mapped CSeekCentral() values on first reference
        try {
            cout << "Initialize " << speciesName << endl;
            this->speciesSeekCentrals[speciesName].InitializeFromSeekConfig(config);
            // TODO - implement different cache size for different species
            this->speciesPclCache.emplace(speciesName, CACHE_SIZE);
        } catch(exception &err) {
            throw_with_nested(config_error(FILELINE + "Error initializing CSeekCentral for species " + speciesName));
        }
    }

    // start the task cleaner thread
    // this->_cleanerThread = thread(&PclInterface::runCleanTasksThread, this, taskTimeoutSec);
}

void PclInterface::pclQuery(const PclQueryArgs &query, PclResult &result)
{
    // TODO - replace with call to pclQueryAsync
    try {
        pclQueryCommon(query, result);
    } catch (named_error &err) {
        string trace = print_exception_stack(err);
        result.success = false;
        result.status = PclStatus::Error;
        result.statusMsg = trace;
        result.__isset.status = true;
        result.__isset.statusMsg = true;
        /* (for future) can use exception_ptr to return the exception to main thread */
    }
    return;
}

int64_t PclInterface::pclQueryAsync(const PclQueryArgs &query)
{
    return -1;
}

bool PclInterface::isQueryComplete(int64_t task_id)
{
    return false;
}

void PclInterface::getQueryResult(int64_t task_id, bool block, PclResult &result)
{
    result.success = false;
    result.status = PclStatus::Error;
    result.statusMsg = "Not Implemented";
    result.__isset.status = true;
    result.__isset.statusMsg = true;
}

int32_t PclInterface::getRpcVersion()
{
    return g_pcl_rpc_constants.PclRPCVersion;
}

int32_t PclInterface::ping()
{
    printf("ping\n");
    return 1;
}

void PclInterface::pclQueryCommon(const PclQueryArgs &query, PclResult &result) {
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
        result.status = PclStatus::Error;
        result.statusMsg = trace;
        result.__isset.status = true;
        result.__isset.statusMsg = true;
        /* (for future) can use exception_ptr to return the exception to main thread */
        return;
    }

    result.success = true;
    result.status = PclStatus::Complete;
    result.__isset.status = true;
    result.__isset.geneExpressions = settings.outputGeneExpression;
    result.__isset.geneCoexpressions = settings.outputGeneCoexpression;
    result.__isset.queryExpressions = settings.outputQueryExpression;
    result.__isset.queryCoexpressions = settings.outputQueryCoexpression;
}