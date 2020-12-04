#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <iostream>
#include <cassert>

using namespace std;
#include "databasei.h"

typedef struct {
  char *dblet1;
  char *dblet2;
  char *dbletCombined;
} Args;

bool parseArgs(int argc, char **argv, Args &args)
{
  string usage = "verifyMergedDBFiles -1 <dblet1> -2 <dblet2> -c <dbletCombined>\n"
    "Verfies that the combination of dblet1 with dblet2 would be dbletCombined\n";

  static struct option long_options[] = {
      {"dblet1",        required_argument, 0,  '1' },
      {"dblet2",        required_argument, 0,  '2' },
      {"dbletCombined", required_argument, 0,  'c' },
      {0,               0,                 0,  0   }
  };
  int opt = 0;
  int long_index = 0;

  while ((opt = getopt_long(argc, argv, "c:1:2:",
                            long_options, &long_index)) != -1)
  {
    switch (opt)
    {
    case '1':
      args.dblet1 = optarg;
      break;
    case '2':
      args.dblet2 = optarg;
      break;
    case 'c':
      args.dbletCombined = optarg;
      break;
    default:
      cout << "Error: unrecognized options" << endl;
      cout << usage << endl;
      return false;
    }
  }
  if (args.dblet1 == nullptr || args.dblet2 == nullptr || args.dbletCombined == nullptr) {
    cout << "Error: missing input args for dblets to compare" << endl;
    cout << usage << endl;
    return false;
  }
  return true;
}

int main(int argc, char** argv) 
{

  Args args;
  bool res = parseArgs(argc, argv, args);
  if (res == false) {
    exit(-1);
  }
  // input is 3 dblet names, orig1, orig2 and combined
  cout << args.dblet1 << " " << args.dblet2 << " " << args.dbletCombined << endl;

  // open dblets
  // database.cpp line 729
  Sleipnir::CDatabaselet dblet1(false);
  Sleipnir::CDatabaselet dblet2(false);
  Sleipnir::CDatabaselet dbletCombined(false);
  dblet1.Open(args.dblet1);
  dblet2.Open(args.dblet2);
  dbletCombined.Open(args.dbletCombined);

  // Verify number of datasets in combined is sum of parts
  size_t numDatasets1 = dblet1.GetDatasets();
  size_t numDatasets2 = dblet2.GetDatasets();
  size_t numDatasetsCombined = dbletCombined.GetDatasets();
  printf("Check num datasets: %zu + %zu = %zu\n", numDatasets1, numDatasets2, numDatasetsCombined);

  if(numDatasets1 + numDatasets2 != numDatasetsCombined) {
    printf("Num datasets don't match: %zu + %zu = %zu\n", numDatasets1, numDatasets2, numDatasetsCombined);
    exit(-1);
  }

  // Veirfy image size matches
  size_t imgSize1 = dblet1.GetImageSize();
  size_t imgSize2 = dblet2.GetImageSize();
  size_t imgSizeCombined = dbletCombined.GetImageSize();
  printf("Check gene image sizes: %zu + %zu = %zu\n", imgSize1, imgSize2, imgSizeCombined);
  if (imgSize1 + imgSize2 != imgSizeCombined) {
    printf("Gene image sizes don't match: %zu + %zu = %zu\n", imgSize1, imgSize2, imgSizeCombined);
    exit(-1);
  }

  // Verify gene names are the same
  size_t numGenes1 = dblet1.GetGenes();
  size_t numGenes2 = dblet2.GetGenes();
  size_t numGenesCombined = dbletCombined.GetGenes();
  assert(numGenes1 == numGenes2);
  assert(numGenes2 == numGenesCombined);
  for (int i = 0; i < numGenes1; i++) {
    assert(dblet1.GetGene(i) == dblet2.GetGene(i));
    assert(dblet2.GetGene(i) == dbletCombined.GetGene(i));
  }

  // Loop through gene data and compare them, make sure the data isn't all 0xFF
  vector<unsigned char> data1, data2, dataCombined;
  size_t iGenes = dblet1.GetTotalNumGenesInCollection();
  size_t numValidData = 0;
  for (int i = 0; i < numGenes1; i++) {
    printf("Compare data for gene %d:  ", i);
    for (int j = 0; j < iGenes; j++) {
      data1.clear();
      data2.clear();
      dataCombined.clear();

      dblet1.Get(i, j, data1);
      dblet2.Get(i, j, data2);
      dbletCombined.Get(i, j, dataCombined);
      // loop through the array comparing values for each dataset
      // dblet1 will be the first set of values in dbletCombined
      for(int x = 0; x < numDatasets1; x++) {
        if (data1[x] != dataCombined[x]) {
          printf("Gene data mismatch(1): gene: %d, co-gene: %d, dataset: %d", i, j, x);
          assert(false);
        }
        if (data1[x] != 0xFF) {
          numValidData++;
        }
      }
      /// dblet2 will be appened to dbletCombined after dblet1 data
      for(int y = 0; y < numDatasets2; y++) {
        if (data2[y] != dataCombined[numDatasets1 + y]) {
          printf("Gene data mismatch(2): gene: %d, co-gene: %d, dataset: %d", i, j, y);
          assert(false);
        }
        if (data2[y] != 0xFF) {
          numValidData++;
        }
      }
    }
    printf("numValidData: %zu\n", numValidData);
  }
  assert(numValidData != 0);
}

// Notes on CDatabaselet calls
// // database.cpp line 729
// CDatabaselet::Open(const string &strFile)

// // database.cpp line 685
// bool CDatabaselet::Get(size_t iOne, size_t iTwo, vector<unsigned char> &vecbData)

// // database.cpp line 699, set fReplace to false to overwrite vecbData rather than append to it
// bool CDatabaselet::Get(size_t iGene, vector<unsigned char> &vecbData, bool fReplace) const {

// Open a CDatabase consisting of CDatabaselets
// See database.cpp line 1173
//For SeekMiner
// bool CDatabase::Open(const char *db_dir, const vector <string> &vecstrGenes, const size_t &iDatasets, const size_t &iNumDBs) {
