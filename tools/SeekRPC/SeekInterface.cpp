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
            cout << "Initialize " << speciesName << endl;
            this->speciesSeekCentrals[speciesName].InitializeFromSeekConfig(config);
        } catch(exception &err) {
            throw_with_nested(config_error(FILELINE + "Error initializing CSeekCentral for species " + speciesName));
        }
    }
}

void SeekInterface::seek_query(const SeekQuery &query, QueryResult &result)
{
    // TODO - spin off a thread to run this query
    // call seek_query_async()
    // then seek_get_result()

    try {
        this->SeekQueryCommon(query, result);
    } catch (named_error &err) {
        string trace = print_exception_stack(err);
        result.success = false;
        result.statusMsg = trace;
        result.__isset.statusMsg = true;
    }
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

