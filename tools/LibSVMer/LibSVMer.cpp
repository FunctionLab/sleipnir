#include <fstream>

#include <vector>
#include <queue>

/*****************************************************************************
 * This file is provided under the Creative Commons Attribution 3.0 license.
 *
 * You are free to share, copy, distribute, transmit, or adapt this work
 * PROVIDED THAT you attribute the work to the authors listed below.
 * For more information, please see the following web page:
 * http://creativecommons.org/licenses/by/3.0/
 *
 * This file is a component of the Sleipnir library for functional genomics,
 * authored by:
 * Curtis Huttenhower (chuttenh@princeton.edu)
 * Mark Schroeder
 * Maria D. Chikina
 * Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
 *
 * If you use this library, the included executable tools, or any related
 * code in your work, please cite the following publication:
 * Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
 * Olga G. Troyanskaya.
 * "The Sleipnir library for computational functional genomics"
 *****************************************************************************/
#include "stdafx.h"
#include "cmdline.h"
#include "statistics.h"

using namespace LIBSVM;

vector<LIBSVM::SVMLabel> ReadLabels(ifstream & ifsm) {

  static const size_t c_iBuffer = 1024;
  char acBuffer[c_iBuffer];
  vector<string> vecstrTokens;
  vector<LIBSVM::SVMLabel> vecLabels;
  size_t numPositives, numNegatives;
  numPositives = numNegatives = 0;
  while (!ifsm.eof()) {
    ifsm.getline(acBuffer, c_iBuffer - 1);
    acBuffer[c_iBuffer - 1] = 0;
    vecstrTokens.clear();
    CMeta::Tokenize(acBuffer, vecstrTokens);
    if (vecstrTokens.empty())
      continue;
    if (vecstrTokens.size() != 2) {
      cerr << "Illegal label line (" << vecstrTokens.size() << "): "
        << acBuffer << endl;
      continue;
    }
    vecLabels.push_back(LIBSVM::SVMLabel(vecstrTokens[0], atof(
            vecstrTokens[1].c_str())));
    if (vecLabels.back().Target > 0)
      numPositives++;
    else
      numNegatives++;
  }
  return vecLabels;
}


struct SortResults {

  bool operator()(const LIBSVM::Result& rOne, const LIBSVM::Result & rTwo) const {
    return (rOne.Value > rTwo.Value);
  }
};


size_t PrintResults(vector<LIBSVM::Result> vecResults, ofstream & ofsm) {
  sort(vecResults.begin(), vecResults.end(), SortResults());
  int LabelVal;
  for (size_t i = 0; i < vecResults.size(); i++) {
    ofsm << vecResults[i].GeneName << '\t' << vecResults[i].Target << '\t'
      << vecResults[i].Value << endl;
  }
};

struct ParamStruct {
  vector<float> vecK, vecTradeoff;
  vector<size_t> vecLoss;
  vector<char*> vecNames;
};

int main(int iArgs, char** aszArgs) {

  gengetopt_args_info sArgs;

  CPCL PCL;//data
  LIBSVM::CLIBSVM SVM;//model

  size_t i, j, iGene, jGene;
  ifstream ifsm;

  bool added;
  added = false;

  if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
    cmdline_parser_print_help();
    return 1;
  }


  //Set model parameters

  if (sArgs.cross_validation_arg < 1){
    cerr << "cross_valid is <1. Must be set at least 1" << endl;
    return 1;
  }
  else if(sArgs.cross_validation_arg < 2){
    cerr << "cross_valid is set to 1. No cross validation holdouts will be run." << endl;
  }

  if (sArgs.num_cv_runs_arg < 1){
    cerr << "number of cv runs is < 1. Must be set at least 1" << endl;
    return 1;
  }

  SVM.SetTradeoff(sArgs.tradeoff_arg);
  SVM.SetNu(sArgs.nu_arg);
  SVM.SetSVMType(sArgs.svm_type_arg);
  SVM.SetBalance(sArgs.balance_flag);

  if (!SVM.parms_check()) {
    cerr << "Sanity check failed, see above errors" << endl;
    return 1;
  }

  //TODO: allow multiple PCL files
  //size_t iFile; //TODO
  // vector<string> PCLs; //TODO

  //check data file
  if (sArgs.input_given) {
    if (!PCL.Open(sArgs.input_arg, sArgs.skip_arg, sArgs.mmap_flag)) {
      cerr << "Could not open input PCL" << endl;
      return 1;
    }
  }

  //read label files
  vector<LIBSVM::SVMLabel> vecLabels;
  set<string> setLabeledGenes;
  if (sArgs.labels_given) {
    ifsm.clear();
    ifsm.open(sArgs.labels_arg);
    if (ifsm.is_open())
      vecLabels = ReadLabels(ifsm);
    else {
      cerr << "Could not read label file" << endl;
      return 1;
    }
    for (i = 0; i < vecLabels.size(); i++)
      setLabeledGenes.insert(vecLabels[i].GeneName);
  }

  LIBSVM::SAMPLE* pTrainSample;

  size_t numSample;
  numSample = sArgs.cross_validation_arg * sArgs.num_cv_runs_arg;
  vector<LIBSVM::SVMLabel> pTrainVector[numSample];
  vector<LIBSVM::SVMLabel> pTestVector[numSample];
  vector<LIBSVM::Result> AllResults;
  vector<LIBSVM::Result> tmpAllResults;

  if (sArgs.model_given && sArgs.labels_given) { //learn once and write to file
    //TODO
    /*
    pTrainSample = CLIBSVM::CreateSample(PCL, vecLabels);
    SVM.Learn(*pTrainSample);
    SVM.WriteModel(sArgs.model_arg);
    */

  } else if (sArgs.model_given && sArgs.output_given) { //read model and classify all
    //TODO: test
    /*
    vector<SVMLabel> vecAllLabels;
    for (size_t i = 0; i < PCL.GetGenes(); i++)
      vecAllLabels.push_back(SVMLabel(PCL.GetGene(i), 0));

    SVM.ReadModel(sArgs.model_arg);
    AllResults = SVM.Classify(PCL, vecAllLabels);
    ofstream ofsm;
    ofsm.open(sArgs.output_arg);
    if (ofsm.is_open())
      PrintResults(AllResults, ofsm);
    else {
      cerr << "Could not open output file" << endl;
    }
    */

  } else if (sArgs.output_given && sArgs.labels_given) {
    size_t ii, index;
    //do learning and classifying with cross validation
    if( sArgs.cross_validation_arg > 1 && sArgs.num_cv_runs_arg >= 1 ){
      for (ii = 0; ii < sArgs.num_cv_runs_arg; ii++) {                    
        std::random_shuffle(vecLabels.begin(), vecLabels.end());

        for (i = 0; i < sArgs.cross_validation_arg; i++) {                  
          index = sArgs.cross_validation_arg * ii + i;
          pTestVector[index].reserve((size_t) vecLabels.size()
              / sArgs.cross_validation_arg + sArgs.cross_validation_arg);
          pTrainVector[index].reserve((size_t) vecLabels.size()
              / (sArgs.cross_validation_arg)
              * (sArgs.cross_validation_arg - 1)
              + sArgs.cross_validation_arg);
          for (j = 0; j < vecLabels.size(); j++) {
            if (j % sArgs.cross_validation_arg == i) {
              pTestVector[index].push_back(vecLabels[j]);
            } else {
              pTrainVector[index].push_back(vecLabels[j]);
            }
          }
        }

      }
    }  
    else{ // if you have less than 2 fold cross, no cross validation is done, all train genes are used and predicted

      // no holdout so train is the same as test gene set
      pTestVector[0].reserve((size_t) vecLabels.size() + sArgs.cross_validation_arg);
      pTrainVector[0].reserve((size_t) vecLabels.size() + sArgs.cross_validation_arg);

      for (j = 0; j < vecLabels.size(); j++) {
        pTestVector[0].push_back(vecLabels[j]);		      
        pTrainVector[0].push_back(vecLabels[j]);		    
      }
    }


    vector<SVMLabel> vec_allUnlabeledLabels;
    vector<Result> vec_allUnlabeledResults;
    vector<Result> vec_tmpUnlabeledResults;
    if (sArgs.all_flag) {
      vec_allUnlabeledLabels.reserve(PCL.GetGenes());
      vec_allUnlabeledResults.reserve(PCL.GetGenes());
      for (i = 0; i < PCL.GetGenes(); i++) {
        if (setLabeledGenes.find(PCL.GetGene(i))
            == setLabeledGenes.end()) {
          vec_allUnlabeledLabels.push_back(
              SVMLabel(PCL.GetGene(i), 0));
          vec_allUnlabeledResults.push_back(Result(PCL.GetGene(i)));
        }
      }
    }


    for (i = 0; i < sArgs.cross_validation_arg * sArgs.num_cv_runs_arg; i++) {
      pTrainSample = LIBSVM::CLIBSVM::CreateSample(PCL, //TODO: make more efficient
          pTrainVector[i]);

      cerr << "Cross Validation Trial " << i << endl;

      SVM.Learn(*pTrainSample);
      cerr << "Learned" << endl;


      tmpAllResults = SVM.Classify(PCL, pTestVector[i]);
      cerr << "Classified " << tmpAllResults.size() << " examples" << endl;
      for(std::vector<LIBSVM::Result>::iterator it = tmpAllResults.begin() ; it != tmpAllResults.end() ; it ++){
        added = false;
        for(std::vector<LIBSVM::Result>::iterator ita = AllResults.begin() ; ita != AllResults.end() ; ita ++){
          if ( (*it).GeneName.compare((*ita).GeneName) == 0 ){
            (*ita).Value += (*it).Value;
            added = true;
            break;
          }

        }

        if(!added)
          AllResults.push_back((*it));

      }
      tmpAllResults.resize(0);
      if (sArgs.all_flag) {
        vec_tmpUnlabeledResults = SVM.Classify(
            PCL, vec_allUnlabeledLabels);
        for (j = 0; j < vec_tmpUnlabeledResults.size(); j++)
          vec_allUnlabeledResults[j].Value
            += vec_tmpUnlabeledResults[j].Value;

      }
      LIBSVM::CLIBSVM::PrintSample(*pTrainSample);

      size_t mem = CMeta::GetMemoryUsage();
      cerr << "before free: " << mem << endl;

      if (i > 0) {
        //LIBSVM::CLIBSVM::FreeSample(*pTrainSample);
        free(pTrainSample);
      }

      mem = CMeta::GetMemoryUsage();
      cerr << "after free: " << mem << endl;
      cerr << "end of a cv run" << endl;
    }

    for(std::vector<LIBSVM::Result>::iterator it = AllResults.begin();
        it != AllResults.end(); ++ it){
      (*it).Value /= sArgs.num_cv_runs_arg;

    }



    if (sArgs.all_flag) { //add the unlabeled results
      for (j = 0; j < vec_allUnlabeledResults.size(); j++)
        vec_allUnlabeledResults[j].Value
          /= (sArgs.cross_validation_arg * sArgs.num_cv_runs_arg);
      AllResults.insert(AllResults.end(),
          vec_allUnlabeledResults.begin(),
          vec_allUnlabeledResults.end());
    }

    ofstream ofsm;
    ofsm.clear();
    ofsm.open(sArgs.output_arg);
    PrintResults(AllResults, ofsm);
    return 0;

  } else {
    cerr << "More options are needed" << endl;
  }

}

