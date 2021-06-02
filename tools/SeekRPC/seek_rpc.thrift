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
    1: optional SearchMethod searchMethod = SearchMethod.CV;
    2: optional DistanceMeasure distanceMeasure = DistanceMeasure.ZScoreHubbinessCorrected;
    3: optional double minQueryGenesFraction = 0.0;
    4: optional double minGenomeFraction = 0.0;
    5: optional double rbpParam = 0.99;
    6: optional bool useNegativeCorrelation = false;
    7: optional bool checkDatasetSize = false;
    8: optional bool useGeneSymbols = false;
    9: optional bool simulateWeights = false;
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
    2: required list<StringDoublePair> geneScores;
    3: optional list<StringDoublePair> datasetWeights;
    4: optional QueryStatus status;
    5: optional string statusMsg;
    // string datasetAvailability;  /* these will come through status channel instead */
    // string queryAvailability;    /* these will come through status channel instead */
}

// Version to track compatibility across changes to the RPC interface
const i32 RPCVersion = 1;

/* Make seekQuery an asynchronous call and have a getStatus() call that returns
 *  status strings while the query is processing and then returns "Done Search"
 *  when it's complete.
 *  datasetAvailability and queryAvailability strings will come through status channel.
 */

service SeekRPC {
    QueryResult seekQuery(1: SeekQuery query);
    i64 seekQueryAsync(1: SeekQuery query);  // returns a task id
    QueryResult getQueryResult(1: i64 taskId, 2: bool block=true);  // returns result from an async task
    bool isQueryComplete(1: i64 taskId),  // returns true if getResult won't block
    string getProgressMessage(1: i64 taskId);  // to retrieve status info for async task given by id
    i32 getRpcVersion();
    i32 ping();  // returns monotonic increasing int
    i32 pvalueGenes();  // input and return types to be determined
    i32 pvalueDatasets();  // input and return types to be determined
    i32 pclData();  // input and return types to be determined
}
