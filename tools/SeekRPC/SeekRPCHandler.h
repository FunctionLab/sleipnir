// This is copied and slightly modified from SeekRPC_server.skeleton.cpp to 
// serve as a link between the SeekRPCIf (interrface) and the SeekInterface.h class.
// This is basically a pass through that calls the corresponding functions in SeekInterface
// when an RPC request arrives. 

#include "gen-cpp/SeekRPC.h"
#include "SeekInterface.h"


using namespace std;
using namespace SeekRPC;


class SeekRPCHandler : virtual public SeekRPCIf {
 private:
  SeekInterface &seekInterface;

 public:
  SeekRPCHandler(SeekInterface &_seekInterface) : seekInterface(_seekInterface) {
  }

  int32_t getRpcVersion() {
    return seekInterface.getRpcVersion();
  }

  void seekQuery(QueryResult& _return, const SeekQuery& query) {
    return seekInterface.seekQuery(query, _return);
  }

  int64_t seekQueryAsync(const SeekQuery& query) {
    return seekInterface.seekQueryAsync(query);
  }

  bool isQueryComplete(const int64_t task_id) {
    return seekInterface.isQueryComplete(task_id);
  }

  void getQueryResult(QueryResult& _return, const int64_t task_id, const bool block) {
    return seekInterface.getQueryResult(task_id, block, _return);
  }

  void getProgressMessage(std::string& _return, const int64_t task_id) {
    _return = seekInterface.getProgressMessage(task_id);
    return;
  }

  int32_t ping() {
    return seekInterface.ping();
  }

  int32_t pvalueGenes() {
    return seekInterface.pvalueGenes();
  }

  int32_t pvalueDatasets() {
    return seekInterface.pvalueDatasets();
  }

  int32_t pclData() {
    // Your implementation goes here
    return seekInterface.pclData();
  }

};