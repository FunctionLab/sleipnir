namespace cpp SeekRPC

struct QueryParams {
    1: required string search_method;
    2: optional string distance_measure;
    3: optional double min_query_genes_fraction;
    4: optional double min_genome_fraction;
    5: optional string correlation_sign;
    6: optional double rpb_param;
}

struct SeekQuery {
    1: required string species;
    2: required list<string> genes;
    3: required list<string> datasets;
    4: required string outputDir;
    5: optional QueryParams parameters;
    6: optional list<string> guideGenes;
}

struct QueryResult {
    1: required list<string> genes;
    2: required list<double> gene_scores;
    3: required list<string> datasets;
    4: required list<double> dataset_weights;
    // string dataset_availability;  /* these will come through status channel instead */
    // string query_availability;    /* these will come through status channel instead */
}


/* Make seek_query an asynchronous call and have a get_status() call that returns
 *  status strings while the query is processing and then returns "Done Search"
 *  when it's complete.
 *  dataset_availability and query_availability strings will come through status channel.
 */

service SeekRPC {
    QueryResult seek_query(1: SeekQuery query);
    i32 seek_query_async(1: SeekQuery query);  // returns a task id
    QueryResult seek_get_result(1: i32 task_id);  // returns result from an async task
    string get_progress_message(1: i32 task_id);  // to retrieve status info for async task given by id
    i32 ping();  // returns monotonic increasing int
    i32 pvalue_genes();  // input and return types to be determined
    i32 pvalue_datasets();  // input and return types to be determined
    i32 pcl_data();  // input and return types to be determined
}
