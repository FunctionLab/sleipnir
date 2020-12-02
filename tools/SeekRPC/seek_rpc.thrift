namespace SeekRPC  // maybe need to use {} around items to put in the namespace

struct QParamsType {
    1: string search_method,
    2: double rpb_param,
    3: double min_query_genes_fraction,
    4: double min_genome_fraction,
    5: string distance_measure,
    6: string correlation_sign,
}

struct QueryType {
    1: list<string> genes,
    2: list<string> datasets,
    3: list<string> guideGenes,
    4: string outputDir,
    5: QParamsType parameters,
}

struct QResultType {
    1: list<string> genes,
    2: list<double> gene_scores,
    3: list<string> datasets,
    4: list<double> dataset_weights,
    // string dataset_availability,  /* these will come through status channel instead */
    // string query_availability,    /* these will come through status channel instead */
}


/* Make seek_query an asynchronous call and have a get_status() call that returns
 *  status strings while the query is processing and then returns "Done Search"
 *  when it's complete.
 *  dataset_availability and query_availability strings will come through status channel.
 */

service SeekRPC {
    QResultType seek_query(1: QueryType query),
    i32 seek_query_async(1: QueryType query),  // returns a task id
    QResultType seek_get_result(1: i32 task_id),  // returns result from an async task
    string get_progress_message(1: i32 task_id),  // to retrieve status info for async task given by id
    i32 ping(),  // returns monotonic increasing int
    i32 pvalue_genes(),  // input and return types to be determined
    i32 pvalue_datasets(),  // input and return types to be determined
    i32 pcl_data(),  // input and return types to be determined
}