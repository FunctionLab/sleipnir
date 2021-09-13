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

  void seekQuery(SeekResult& _return, const SeekQueryArgs& query) {
    return seekInterface.seekQuery(query, _return);
  }

  int64_t seekQueryAsync(const SeekQueryArgs& query) {
    return seekInterface.seekQueryAsync(query);
  }

  bool isQueryComplete(const int64_t taskId) {
    return seekInterface.isQueryComplete(taskId);
  }

  void getSeekResult(SeekResult& _return, const int64_t taskId, const bool block) {
    return seekInterface.getSeekResult(taskId, block, _return);
  }

  void getProgressMessage(std::string& _return, const int64_t taskId) {
    _return = seekInterface.getProgressMessage(taskId);
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

  void pclQuery(PclResult& _return, const PclQueryArgs& query) {
    return seekInterface.pclQuery(query, _return);
  }

  int64_t pclQueryAsync(const PclQueryArgs& query) {
    return seekInterface.pclQueryAsync(query);
  }

  void getPclResult(PclResult& _return, const int64_t taskId, const bool block) {
    return seekInterface.getPclResult(taskId, block, _return);
  }
};