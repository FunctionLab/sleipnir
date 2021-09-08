// This is copied and slightly modified from PclRPC_server.skeleton.cpp to 
// serve as a link between the PclRPCIf (interrface) and the PclInterface.h class.
// This is basically a pass through that calls the corresponding functions in PclInterface
// when an RPC request arrives. 

#include "PclInterface.h"
#include "gen-cpp/PclRPC.h"

using namespace std;
using namespace PclRPC;

class PclRPCHandler : virtual public PclRPCIf {
 private:
  PclInterface &pclInterface;

 public:
  PclRPCHandler(PclInterface &_pclInterface) : pclInterface(_pclInterface) {
  }

  void pclQuery(PclResult& _return, const PclQueryArgs& query) {
    return pclInterface.pclQuery(query, _return);
  }

  int64_t pclQueryAsync(const PclQueryArgs& query) {
    return pclInterface.pclQueryAsync(query);
  }

  void getQueryResult(PclResult& _return, const int64_t taskId, const bool block) {
    return pclInterface.getQueryResult(taskId, block, _return);
  }

  bool isQueryComplete(const int64_t taskId) {
    return pclInterface.isQueryComplete(taskId);
  }

  int32_t getRpcVersion() {
    return pclInterface.getRpcVersion();
  }

  int32_t ping() {
    return pclInterface.ping();
  }

};