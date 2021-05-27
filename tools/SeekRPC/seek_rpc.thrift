namespace cpp SeekRPC

// Define the possible values for search method
enum SearchMethod {
    CV = 1,
    CVCustom,
    EqualWeighting,
    OrderStatistics,
}

// Define the possible values for distance measure
enum DistanceMeasure {
    ZScore = 1,
    ZScoreHubbinessCorrected,
    Correlation,
}

struct QueryParams {
    1: optional SearchMethod search_method = SearchMethod.CV;
    2: optional DistanceMeasure distance_measure = DistanceMeasure.ZScoreHubbinessCorrected;
    3: optional double min_query_genes_fraction = 0.0;
    4: optional double min_genome_fraction = 0.0;
    5: optional double rbp_param = 0.99;
    6: optional bool useNegativeCorrelation = false;
    7: optional bool check_dataset_size = false;
    8: optional bool use_gene_symbols = false;
}

struct SeekQuery {
    1: required string species = "Unknown";
    2: required list<string> genes;
    3: optional list<string> datasets;
    4: optional QueryParams parameters;
    5: optional list<string> guideGenes;
    6: optional string outputDir = "/tmp/seek";
}

struct StringDoublePair {
    1: required string name;
    2: required double value;
}

enum QueryStatus {
    Complete = 1,
    Incomplete,
    Error,
}

struct QueryResult {
    1: required bool success;
    2: required list<StringDoublePair> gene_scores;
    3: optional list<StringDoublePair> dataset_weights;
    4: optional QueryStatus status;
    5: optional string statusMsg;
    // string dataset_availability;  /* these will come through status channel instead */
    // string query_availability;    /* these will come through status channel instead */
}

// Version to track compatibility across changes to the RPC interface
const i32 RPC_Version = 1;

/* Make seek_query an asynchronous call and have a get_status() call that returns
 *  status strings while the query is processing and then returns "Done Search"
 *  when it's complete.
 *  dataset_availability and query_availability strings will come through status channel.
 */

service SeekRPC {
    QueryResult seek_query(1: SeekQuery query);
    i64 seek_query_async(1: SeekQuery query);  // returns a task id
    QueryResult seek_get_result(1: i64 task_id, 2: bool block=true);  // returns result from an async task
    bool is_query_complete(1: i64 task_id),  // returns true if get_result won't block
    string get_progress_message(1: i64 task_id);  // to retrieve status info for async task given by id
    i32 get_rpc_version();
    i32 ping();  // returns monotonic increasing int
    i32 pvalue_genes();  // input and return types to be determined
    i32 pvalue_datasets();  // input and return types to be determined
    i32 pcl_data();  // input and return types to be determined
}
