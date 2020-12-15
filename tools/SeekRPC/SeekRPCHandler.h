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

  void seek_query(QueryResult& _return, const SeekQuery& query) {
    return seekInterface.seek_query(query, _return);
  }

  int32_t seek_query_async(const SeekQuery& query) {
    return seekInterface.seek_query_async(query);
  }

  void seek_get_result(QueryResult& _return, const int32_t task_id) {
    return seekInterface.seek_get_result(task_id, _return);
  }

  void get_progress_message(std::string& _return, const int32_t task_id) {
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