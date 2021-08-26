namespace cpp PclRPC

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

enum PclStatus {
    Complete = 1,
    Incomplete,
    Error,
}

struct PclResult {
    1: required bool success;
    2: required list<i32> datasetSizes;
    3: optional list<double> geneExpressions;
    4: optional list<double> geneCoexpressions;
    5: optional list<double> queryExpressions;
    6: optional list<double> queryCoexpressions;
    7: optional PclStatus status;
    8: optional string statusMsg;
}

// Version to track compatibility across changes to the RPC interface
const i32 PclRPCVersion = 1;

service PclRPC {
    PclResult pclQuery(1: PclQueryArgs query);
    i64 pclQueryAsync(1: PclQueryArgs query);  // returns a task id
    PclResult getQueryResult(1: i64 taskId, 2: bool block=true);  // returns result from an async task
    bool isQueryComplete(1: i64 taskId),  // returns true if getResult won't block
    i32 getRpcVersion();
    i32 ping();  // returns monotonic increasing int
}
