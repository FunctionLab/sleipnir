#include <iostream>
#include "SeekInterface.h"

using namespace std;
using namespace SeekRPC;

void SeekInterface::seek_query(const SeekQuery &query, QueryResult &result)
{
    printf("seek_query\n");
    for (auto gene: query.genes)
      cout << gene << endl;
    vector<string> rgenes;
    rgenes.push_back("a_separate_result");
    result.genes = rgenes;
    // result.__set_genes(rgenes);
    // result.genes.push_back("best_fit_gene");
    return;
}

int32_t SeekInterface::seek_query_async(const SeekQuery &query)
{ 
    printf("seek_query_async\n");
    return 0;
}

void SeekInterface::seek_get_result(int32_t task_id, QueryResult &result)
{
    printf("seek_get_result\n");
    return;
}

string SeekInterface::get_progress_message(int32_t task_id)
{
    printf("get_progress_message\n");
    return "Status Message";
}

int32_t SeekInterface::ping()
{
    printf("ping\n");
    return 1;
}

int32_t SeekInterface::pvalue_genes()
{
    printf("pvalue_genes\n");
    return 0;
}

int32_t SeekInterface::pvalue_datasets()
{
    printf("pvalue_datasets\n");
    return 0;
}

int32_t SeekInterface::pcl_data()
{
    printf("pcl_data\n");
    return 0;
}