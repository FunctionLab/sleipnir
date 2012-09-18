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

using namespace SVMLight;

struct ParamStruct {
	vector<float> vecK, vecTradeoff;
	vector<size_t> vecLoss;
	vector<char*> vecNames;
};

ParamStruct ReadParamsFromFile(ifstream& ifsm, string outFile) {
	static const size_t c_iBuffer = 1024;
	char acBuffer[c_iBuffer];
	char* nameBuffer;
	vector<string> vecstrTokens;
	size_t extPlace;
	string Ext, FileName;
	if ((extPlace = outFile.find_first_of(".")) != string::npos) {
		FileName = outFile.substr(0, extPlace);
		Ext = outFile.substr(extPlace, outFile.size());
	} else {
		FileName = outFile;
		Ext = "";
	}
	ParamStruct PStruct;
	size_t index = 0;
	while (!ifsm.eof()) {
		ifsm.getline(acBuffer, c_iBuffer - 1);
		acBuffer[c_iBuffer - 1] = 0;
		vecstrTokens.clear();
		CMeta::Tokenize(acBuffer, vecstrTokens);
		if (vecstrTokens.empty())
			continue;
		if (vecstrTokens.size() != 3) {
			cerr << "Illegal params line (" << vecstrTokens.size() << "): "
					<< acBuffer << endl;
			continue;
		}
		if (acBuffer[0] == '#') {
			cerr << "skipping " << acBuffer << endl;
		} else {
			PStruct.vecLoss.push_back(atoi(vecstrTokens[0].c_str()));
			PStruct.vecTradeoff.push_back(atof(vecstrTokens[1].c_str()));
			PStruct.vecK.push_back(atof(vecstrTokens[2].c_str()));
			PStruct.vecNames.push_back(new char[c_iBuffer]);
			if (PStruct.vecLoss[index] == 4 || PStruct.vecLoss[index] == 5)
				sprintf(PStruct.vecNames[index], "%s_l%d_c%4.6f_k%4.3f%s",
						FileName.c_str(), PStruct.vecLoss[index],
						PStruct.vecTradeoff[index], PStruct.vecK[index],
						Ext.c_str());
			else
				sprintf(PStruct.vecNames[index], "%s_l%d_c%4.6f%s",
						FileName.c_str(), PStruct.vecLoss[index],
						PStruct.vecTradeoff[index], Ext.c_str());
			index++;
		}

	}
	return PStruct;
}

int main(int iArgs, char** aszArgs) {
	gengetopt_args_info sArgs;	
	SVMLight::CSVMPERF SVM;
	
	size_t i, j, iGene, jGene;
	ifstream ifsm;
	float d;
	int   iRet;
	map<string, size_t>	mapstriZeros, mapstriDatasets;
	vector<string>  vecstrDatasets;
	DIR* dp;
	struct dirent* ep;
	
	if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
		cmdline_parser_print_help();
		return 1;
	}
	SVM.SetVerbosity(sArgs.verbosity_arg);
	SVM.SetLossFunction(sArgs.error_function_arg);
	if (sArgs.k_value_arg > 1) {
		cerr << "k_value is >1. Setting default 0.5" << endl;
		SVM.SetPrecisionFraction(0.5);
	} else if (sArgs.k_value_arg <= 0) {
		cerr << "k_value is <=0. Setting default 0.5" << endl;
		SVM.SetPrecisionFraction(0.5);
	} else {
		SVM.SetPrecisionFraction(sArgs.k_value_arg);
	}

	
	if (sArgs.cross_validation_arg < 1){
	  cerr << "cross_valid is <1. Must be set at least 1" << endl;
	  return 1;
	}
	else if(sArgs.cross_validation_arg < 2){
	  cerr << "cross_valid is set to 1. No cross validation holdouts will be run." << endl;
	}
	
	SVM.SetTradeoff(sArgs.tradeoff_arg);
	if (sArgs.slack_flag)
		SVM.UseSlackRescaling();
	else
		SVM.UseMarginRescaling();

	
	if (!SVM.parms_check()) {
		cerr << "Sanity check failed, see above errors" << endl;
		return 1;
	}

	//  cout << "there are " << vecLabels.size() << " labels processed" << endl;
	/*
	if (sArgs.input_given) {
		if (!PCL.Open(sArgs.input_arg, sArgs.skip_arg, sArgs.mmap_flag)) {
		cerr << "Could not open input PCL" << endl;
		return 1;
		}
	}
	*/
	
	// read in the list of datasets
	if(sArgs.directory_arg ) {
	  dp = opendir (sArgs.directory_arg);
	  if (dp != NULL){
	    while (ep = readdir (dp)){
	      // skip . .. files and temp files with ~
	      if (ep->d_name[0] == '.' || ep->d_name[strlen(ep->d_name)-1] == '~') 
		continue;
	      
	      // currently opens all files. Add filter here if want pick file extensions
	      vecstrDatasets.push_back((string)sArgs.directory_arg + "/" + ep->d_name);	      
	    }
	    (void) closedir (dp);	    
	    
	    cerr << "Input Dir contrains # datasets: " << vecstrDatasets.size() << '\n';
	  }
	  else{
	    cerr << "Couldn't open the directory: " << sArgs.directory_arg << '\n';
	    return 1;
	  }	  
	}
	
	///######################
	// Chris added
	vector<SVMLight::SVMLabelPair*> vecLabels;
	/*********
	 ** PCL gold standard later
	 
	CPCL Labels;
	if (sArgs.labels_given) {	  
	  if (!Labels.Open(sArgs.labels_arg, sArgs.skip_arg, sArgs.mmap_flag)) {
	    cerr << "Could not open input PCL" << endl;
	    return 1;
	  }
	  
	  for(size_t i = 0; i < PCL.GetGenes(); i++)
	    for(size_t j = 0; j < PCL.GetExperiments(); j++)
	      if (!CMeta::IsNaN(d = PCL.Get(i, j)))  
		vecLabels.push_back(SVMLight::SVMLabelPairPair(d, i, j));
	  
	}
	****/
	CDat Labels;
	CDat Results;
	if (sArgs.labels_given) {	  
	  if (!Labels.Open(sArgs.labels_arg, sArgs.mmap_flag)) {
	    cerr << "Could not open input labels Dat" << endl;
	    return 1;
	  }
	  
	  // set all NaN values to negatives
	  if( sArgs.nan2neg_given ){
	    for(i = 0; i < Labels.GetGenes(); i++)
	      for(j = (i+1); j < Labels.GetGenes(); j++)
		if (CMeta::IsNaN(d = Labels.Get(i, j)))  
		  Labels.Set(i, j, -1);
	  }
	  
	  // 
	  //if( sArgs.geneq_arg ) {
	    
	  //}
	  	  
	  for(i = 0; i < Labels.GetGenes(); i++)
	    for(j = (i+1); j < Labels.GetGenes(); j++)
	      if (!CMeta::IsNaN(d = Labels.Get(i, j)))  
		vecLabels.push_back(new SVMLight::SVMLabelPair(d, i, j));
	  
	}
	
	SVMLight::SAMPLE* pTrainSample;
	SVMLight::SAMPLE* pAllSample;
	vector<SVMLight::SVMLabelPair*> pTrainVector[sArgs.cross_validation_arg];
	vector<SVMLight::SVMLabelPair*> pTestVector[sArgs.cross_validation_arg];
	//vector<SVMLight::Result> AllResults;
	//vector<SVMLight::Result> tmpAllResults;
	
	if (sArgs.model_given && sArgs.labels_given) { //learn once and write to file
	  //DEBUG implement
	  //pTrainSample = CSVMPERF::CreateSample(vecstrDatasets, vecLabels, Labels.GetGeneNames());
	  //SVM.Learn(*pTrainSample);
	  //SVM.WriteModel(sArgs.model_arg);
	} else if (sArgs.model_given && sArgs.output_given) { //read model and classify all
	  
	  // DEBUG need to figure out which/where to select all genes to predict
	  // Also, need to paralize the full predictions
	  
	  SVM.ReadModel(sArgs.model_arg);
	  //AllResults = SVM.Classify(PCL, vecAllLabels);
	  
	  // DEBUG output predictions
	  // sArgs.output_arg
	  
	} else if (sArgs.output_given && sArgs.labels_given) {
		//do learning and classifying with cross validation
	        if( sArgs.cross_validation_arg > 1){
		  for (i = 0; i < sArgs.cross_validation_arg; i++) {
		    pTestVector[i].reserve((size_t) vecLabels.size()
					   / sArgs.cross_validation_arg + sArgs.cross_validation_arg);
		    pTrainVector[i].reserve((size_t) vecLabels.size()
					    / (sArgs.cross_validation_arg)
					    * (sArgs.cross_validation_arg - 1)
					    + sArgs.cross_validation_arg);
		    for (j = 0; j < vecLabels.size(); j++) {
		      if (j % sArgs.cross_validation_arg == i) {
			pTestVector[i].push_back(vecLabels[j]);
		      } else {
			pTrainVector[i].push_back((vecLabels[j]));
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
				
		// initalize the results
		Results.Open(Labels.GetGeneNames(), true);
		
		//if (sArgs.all_flag) {
		  // DEBUG, probably don't allow all flag
		  // orig for classify all genes
		//}
		if (sArgs.params_given) { //reading paramters from file
			ifsm.close();
			ifsm.clear();
			ifsm.open(sArgs.params_arg);
			if (!ifsm.is_open()) {
				cerr << "Could not open: " << sArgs.params_arg << endl;
				return 1;
			}
			ParamStruct PStruct;
			string outFile(sArgs.output_arg);
			PStruct = ReadParamsFromFile(ifsm, outFile);

			size_t iParams;
			ofstream ofsm;
			SVMLight::SAMPLE * ppTrainSample[sArgs.cross_validation_arg];
			
			//build all the samples since they are being reused
			/*
			for (i = 0; i < sArgs.cross_validation_arg; i++)
			  ppTrainSample[i] = SVMLight::CSVMPERF::CreateSample(vecstrDatasets,
									      pTrainVector[i]);
			
			for (iParams = 0; iParams < PStruct.vecTradeoff.size(); iParams++) {
			  // DEBUG need to implement			  
			}
			*/
		} else { //run once
		  
		  // create sample for all labels
		  // DEBUG need to implement
		  cerr << "CreateDocs!"<< endl;
		  SVMLight::CSVMPERF::CreateDoc(vecstrDatasets,
						vecLabels,
						Labels.GetGeneNames());
		  
		  for (i = 0; i < sArgs.cross_validation_arg; i++) {
		    cerr << "Cross validation: " << i << endl;
		    
		    // DEBUG need to implement
		    pTrainSample = SVMLight::CSVMPERF::CreateSample(pTrainVector[i]);
		    cerr << "Cross Validation Trial " << i << endl;		    
		    
		    SVM.Learn(*pTrainSample);
		    cerr << "Learned" << endl;
		    SVM.Classify(Results,
				 pTestVector[i]);
		    //cerr << "Classified " << tmpAllResults.size() << " examples"
		    //<< endl;
		    //AllResults.insert(AllResults.end(), tmpAllResults.begin(),
		    //				      tmpAllResults.end());
		    //tmpAllResults.resize(0);
		    ///if (sArgs.all_flag) {
		      // DEBUG, probably don't allow all flag
		      // orig for classify all genes				  
		    //}
		    
		    // MAIRA, Also, why did you not start at i==0???
		    
		    // DEBUG
		    SVMLight::CSVMPERF::FreeSample_leave_Doc(*pTrainSample);
		    free(pTrainSample);
		  }
		  
		  //if (sArgs.all_flag) { //add the unlabeled results
		    // DEBUG, probably don't allow all flag
		    // orig for classify all genes
		  //}
		  
		  Results.Save(sArgs.output_arg);
		  return 0;
		}
	} else {
		cerr << "More options are needed" << endl;
	}

}

