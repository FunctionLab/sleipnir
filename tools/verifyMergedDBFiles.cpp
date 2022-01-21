#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <iostream>
#include <filesystem>
#include <stdexcept>
#include "seekhelper.h"
#include "databasei.h"

using namespace std;
namespace fs = std::filesystem;


typedef struct {
  char *inputDirList;
  char *dbFile;
  char *combinedDir;
  bool verbose = false;
} Args;

bool parseArgs(int argc, char **argv, Args &args)
{
  string usage = "verifyMergedDBFiles -1 <dblet1> -2 <dblet2> -c <dbletCombined>\n"
    "Verfies that the combination of dblet1 with dblet2 would be dbletCombined\n";

  static struct option long_options[] = {
      {"inputDirList",  required_argument, 0,  'i' },
      {"combinedDir",   required_argument, 0,  'c' },
      {"dbFile",        required_argument, 0,  'f' },
      {"verbose",       no_argument,       0,  'v'},
      {0,               0,                 0,  0   }
  };
  int opt = 0;
  int long_index = 0;

  while ((opt = getopt_long(argc, argv, "i:c:f:",
                            long_options, &long_index)) != -1)
  {
    switch (opt)
    {
    case 'i':
      args.inputDirList = optarg;
      break;
    case 'c':
      args.combinedDir = optarg;
      break;
    case 'f':
      args.dbFile = optarg;
      break;
    case 'v':
      args.verbose = true;
      break;
    default:
      cout << "Error: unrecognized options" << endl;
      cout << usage << endl;
      return false;
    }
  }
  if (args.inputDirList == nullptr || args.dbFile == nullptr || args.combinedDir == nullptr) {
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
  cout << "VerifyMergedDBFiles: " << args.inputDirList << " " << args.dbFile << " " << args.combinedDir << endl;

  // read in list of input dirs
  vector<string> inputDirs;
  loadOneColumnTextFile(args.inputDirList, inputDirs);

  // Open the individual input files
  // vector<Sleipnir::CDatabaselet> dblets(inputDirs.size(), CDatabaselet(false));
  // vector<Sleipnir::CDatabaselet> dblets;
  vector<unique_ptr<Sleipnir::CDatabaselet>> dblets;
  for (int i=0; i<inputDirs.size(); i++) {
    fs::path filePath(inputDirs[i]);
    filePath /= args.dbFile;
    // open dblets
    // database.cpp line 729
    auto dblet = make_unique<Sleipnir::CDatabaselet>(false);
    dblets.push_back(move(dblet));
    dblets[i]->Open(filePath);
  }

  // Open the combined output file
  Sleipnir::CDatabaselet dbletCombined(false);
  dbletCombined.Open(fs::path(args.combinedDir) / args.dbFile);

  // Verify number of datasets in combined is sum of parts
  size_t numDatasets = 0;
  size_t imgSize = 0;
  for (auto & dblet: dblets) {
    numDatasets += dblet->GetDatasets();
    imgSize += dblet->GetImageSize();
  }

  // Check num datasets match
  size_t numDatasetsCombined = dbletCombined.GetDatasets();
  if(numDatasets == numDatasetsCombined) {
    printf("Num datasets match: %zu == %zu\n", numDatasets, numDatasetsCombined);
  } else {
    printf("Num datasets don't match: %zu == %zu\n", numDatasets, numDatasetsCombined);
    exit(-1);
  }

  // Check size of datasets match
  size_t imgSizeCombined = dbletCombined.GetImageSize();
  if (imgSize == imgSizeCombined) {
    printf("Gene image sizes match: %zu == %zu\n", imgSize, imgSizeCombined);
  } else {
    printf("Gene image sizes don't match: %zu == %zu\n", imgSize, imgSizeCombined);
    exit(-1);
  }

  // Verify gene names are the same
  size_t numDbletGenes = dbletCombined.GetGenes();
  for (auto & dblet: dblets) {
    if (dblet->GetGenes() != numDbletGenes) {
      throw runtime_error("Num dbFile genes differ: " + string(args.dbFile));
    }
    for (int i = 0; i < numDbletGenes; i++) {
      if (dblet->GetGene(i) != dbletCombined.GetGene(i)) {
        throw runtime_error("Gene names differ: " + string(args.dbFile));
      }
    }
  }
  printf("Num genes in dbFile %s: %zu\n", args.dbFile, numDbletGenes);

  // Loop through gene data and compare them, make sure the data isn't all 0xFF
    vector<unsigned char> data1, dataCombined;
  size_t iGenes = dbletCombined.GetTotalNumGenesInCollection();
  size_t numValidData = 0;
  for (int i = 0; i < numDbletGenes; i++) {
    // For each gene stored in the db file
    if (args.verbose) {
      printf("Compare data for gene %d:  \n", i);
    }
    for (int j = 0; j < iGenes; j++) {
      // for each of the dataset wide genes
      dataCombined.clear();
      dbletCombined.Get(i, j, dataCombined);
      size_t dset_offset = 0;
      for (int k=0; k<dblets.size(); k++) {
        // for each of the partial db files to be combined
        // printf("Compare dblet %d:  ", k);
        data1.clear();
        dblets[k]->Get(i, j, data1);
        // loop through the array comparing values for each dataset
        // dblet[0] will be the first set of values in dbletCombined
        // followed by dblet[1] etc
        size_t numDatasets = dblets[k]->GetDatasets();
        /// dblet2 will be appened to dbletCombined after dblet1 data
        for(int x = 0; x < numDatasets; x++) {
          if (data1[x] != dataCombined[dset_offset + x]) {
            printf("Gene data mismatch(2): gene: %d, dblet %d, co-gene: %d, dataset: %d\n", i, k, j, x);
            throw runtime_error("Gene data mismatch: " + string(args.dbFile));
          }
          if (data1[x] != 0xFF) {
            numValidData++;
          }
        }
        dset_offset += numDatasets;
      }
    }
    if (args.verbose) {
      printf("numValidData: %zu\n", numValidData);
    }
  }
  if (numValidData == 0) {
    throw runtime_error("All data is NA: " + string(args.dbFile));
  }
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
