namespace cpp SeekRPC

struct StringDoublePair {
    1: required string name;
    2: required double value;
}

enum QueryStatus {
    Complete = 1,
    Incomplete,
    Error,
}

// Version to track compatibility across changes to the RPC interface
const double RPCVersion = 1.3;

// ### Seek Query Data Structures ###
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

// was QueryParams
struct SeekQueryParams {
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

// was SeekQuery
struct SeekQueryArgs {
    1: required string species = "Unknown";
    2: required list<string> genes;
    3: optional list<string> datasets;
    4: optional SeekQueryParams parameters;
    5: optional list<string> guideGenes;
    6: optional string outputDir = "";
}

// was QueryResult
struct SeekResult {
    1: required bool success;
    2: required list<StringDoublePair> geneScores;
    3: optional list<StringDoublePair> datasetWeights;
    4: optional QueryStatus status;
    5: optional string statusMsg;
    // string datasetAvailability;  /* these will come through status channel instead */
    // string queryAvailability;    /* these will come through status channel instead */
}

// ### PCL Data Structures ###
struct PclSettings {
    1: optional bool outputNormalized = false;
    2: optional bool outputGeneExpression = false;
    3: optional bool outputQueryExpression = false;
    4: optional bool outputGeneCoexpression = false;
    5: optional bool outputQueryCoexpression = false;
    6: optional double rbp = 0.99;
}

struct PclQueryArgs {
    1: required string species = "Unknown";
    2: required list<string> datasets;
    3: optional list<string> genes;
    4: optional list<string> queryGenes;
    5: optional PclSettings settings;
    6: optional string outputDir = "/tmp/seek";
}

struct PclResult {
    1: required bool success;
    2: required list<i32> datasetSizes;
    3: optional list<string> experimentNames;
    4: optional list<double> geneExpressions;
    5: optional list<double> geneCoexpressions;
    6: optional list<double> queryExpressions;
    7: optional list<double> queryCoexpressions;
    8: optional QueryStatus status;
    9: optional string statusMsg;
}

// Note that for backward compatibility if the geneScores for all
//  genes are passed in, then the gene entrezIDs don't need to be
//  provided and it will use the gene_map gene order to match
//  scores to genes. In addition, if all genes are passed in and
//  useRank=true the server will sort and rank the gene scores
//  and return the pvalues based on the ranks (in that case no
//  geneRanks need to be passed in).
struct PValueGeneArgs {
    1: required string species = "Unknown";
    2: optional list<string> genes;   // EntrezIDs
    3: optional list<double> geneScores;  // Seek query gene scores
    4: optional list<i32> geneRanks;  // Seek query gene ranks
    5: optional bool useRank = false;
}

struct PValueDatasetArgs {
    1: required string species = "Unknown";
    2: optional list<string> datasets;   // dataset names
    3: optional list<double> datasetWeights;  // Seek query result weights
}

struct PValueResult {
    1: required bool success;
    2: optional list<double> pvalues;
    3: optional QueryStatus status;
    4: optional string statusMsg;
}


/* Make seekQuery an asynchronous call and have a getStatus() call that returns
 *  status strings while the query is processing and then returns "Done Search"
 *  when it's complete.
 *  datasetAvailability and queryAvailability strings will come through status channel.
 */

service SeekRPC {
    SeekResult seekQuery(1: SeekQueryArgs query);
    i64 seekQueryAsync(1: SeekQueryArgs query);  // returns a task id
    SeekResult getSeekResult(1: i64 taskId, 2: bool block=true);  // returns result from an async task
    bool isQueryComplete(1: i64 taskId),  // returns true if getResult won't block
    string getProgressMessage(1: i64 taskId);  // to retrieve status info for async task given by id
    PclResult pclQuery(1: PclQueryArgs query);
    i64 pclQueryAsync(1: PclQueryArgs query);  // returns a task id
    PclResult getPclResult(1: i64 taskId, 2: bool block=true);  // returns result from an async task
    PValueResult pvalueGenes(1: PValueGeneArgs query);
    PValueResult pvalueDatasets(1: PValueDatasetArgs query);
    double getRpcVersion();
    i32 numTasksOutstanding();
    i32 ping();  // returns monotonic increasing int
}
