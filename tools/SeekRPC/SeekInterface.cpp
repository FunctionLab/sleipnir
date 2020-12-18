#include <iostream>
#include <cassert>
#include <boost/algorithm/string/join.hpp>
#include "SeekInterface.h"
#include "seekcentral.h"
#include "seekerror.h"

using namespace std;
using namespace SeekRPC;

SeekInterface::SeekInterface(vector<string> &configFiles) {
    getConfigs(configFiles, this->speciesConfigs);

    for( auto const& [speciesName, config] : this->speciesConfigs ) {
        // Initialize a seekCentral instance for each species.
        //  C++ automatically default inits the mapped CSeekCentral() values on first reference
        try {
            this->speciesSeekCentrals[speciesName].InitializeFromSeekConfig(config);
        } catch(exception &err) {
            throw_with_nested(config_error(FILELINE + "Error initializing CSeekCentral for species " + speciesName));
        }
    }
}

void SeekInterface::seek_query(const SeekQuery &query, QueryResult &result)
{
    // TODO - spin off a thread to run this query
    try {
        this->SeekQueryCommon(query, result);
    } catch (named_error &err) {
        print_exception_stack(err);
        result.success = false;
        result.statusMsg = err.what();
        // TODO - print stack to string to return in result
        // TODO - add a result status to QueryResult and a statusString or errString
    }

    // printf("seek_query\n");
    // for (auto gene: query.genes)
    //   cout << gene << endl;
    // const QueryParams &params = query.parameters;
    // if (params.rbp_param) {
    //     cout << "RPB Set: " << params.rbp_param << endl;
    // }
    // if (params.__isset.rbp_param) {
    //     cout << "is set" << endl;
    // }
    // vector<string> rgenes;
    // rgenes.push_back("a_separate_result");
    // result.genes = rgenes;
    // // result.__set_genes(rgenes);
    // // result.genes.push_back("best_fit_gene");
    return;
}

int32_t SeekInterface::seek_query_async(const SeekQuery &query)
{ 
    printf("seek_query_async\n");
    return 0;
}

void SeekInterface::seek_get_result(int32_t task_id, QueryResult &result)
{
    printf("seek_get_result\n");
    return;
}

string SeekInterface::get_progress_message(int32_t task_id)
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

    // seekcentral InitializeQuery expects a string with genes delimited by " " and multiple queries separated by "|"
    std::string joinedGenes = boost::algorithm::join(query.genes, " ");
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
    bool r = querySC.InitializeQuery(query.outputDir,
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

    if (r) {
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
            gsl_rng_set(rnd, 100);
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
    }

    if (querySC.numQueries() > 1) {
        // Error - should only be running one query per RPC request
        // TODO - throw state error?
        assert(querySC.numQueries() == 1);
    }

    vector<double> gene_scores;
    result.__set_gene_scores(gene_scores);
    vector<string> datasets;
    result.__set_datasets(datasets);
    vector<double> weights;
    result.__set_dataset_weights(weights);

    // get the genes and scores and add them to the rpc result reply
    vector<PairedResult<string, float>> geneResults = querySC.getGeneResult(0);
    int numGenes = geneResults.size();
    for (int i=0; i<numGenes; i++) {
        result.genes.push_back(geneResults[i].key);
        result.gene_scores.push_back(geneResults[i].val);
    }

    // get the datasets and weigts and add them to the rpc result reply
    vector<PairedResult<string, float>> datasetResults = querySC.getDatasetResult(0);
    int numDsets = datasetResults.size();
    for (int i=0; i<numDsets; i++) {
        result.datasets.push_back(datasetResults[i].key);
        result.dataset_weights.push_back(datasetResults[i].val);
    }

    // result.__set_genes();
    // result.__set_gene_scores();
    result.success = true;

    querySC.Destruct();

}