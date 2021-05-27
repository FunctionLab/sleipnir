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

  int32_t get_rpc_version() {
    return seekInterface.get_rpc_version();
  }

  void seek_query(QueryResult& _return, const SeekQuery& query) {
    return seekInterface.seek_query(query, _return);
  }

  int64_t seek_query_async(const SeekQuery& query) {
    return seekInterface.seek_query_async(query);
  }

  bool is_query_complete(const int64_t task_id) {
    return seekInterface.is_query_complete(task_id);
  }

  void seek_get_result(QueryResult& _return, const int64_t task_id, const bool block) {
    return seekInterface.seek_get_result(task_id, block, _return);
  }

  void get_progress_message(std::string& _return, const int64_t task_id) {
    _return = seekInterface.get_progress_message(task_id);
    return;
  }

  int32_t ping() {
    return seekInterface.ping();
  }

  int32_t pvalue_genes() {
    return seekInterface.pvalue_genes();
  }

  int32_t pvalue_datasets() {
    return seekInterface.pvalue_datasets();
  }

  int32_t pcl_data() {
    // Your implementation goes here
    return seekInterface.pcl_data();
  }

};