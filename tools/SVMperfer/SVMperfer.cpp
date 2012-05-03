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
//#include "../../extlib/svm_light/svm_light/kernel.h"

vector<SVMLight::SVMLabel> ReadLabels(ifstream & ifsm) {

	static const size_t c_iBuffer = 1024;
	char acBuffer[c_iBuffer];
	vector<string> vecstrTokens;
	vector<SVMLight::SVMLabel> vecLabels;
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
		vecLabels.push_back(SVMLight::SVMLabel(vecstrTokens[0], atof(
				vecstrTokens[1].c_str())));
		if (vecLabels.back().Target > 0)
			numPositives++;
		else
			numNegatives++;
	}
	return vecLabels;
}

struct SortResults {

	bool operator()(const SVMLight::Result& rOne, const SVMLight::Result & rTwo) const {
		return (rOne.Value > rTwo.Value);
	}
};

size_t PrintResults(vector<SVMLight::Result> vecResults, ofstream & ofsm) {
	sort(vecResults.begin(), vecResults.end(), SortResults());
	int LabelVal;
	for (size_t i = 0; i < vecResults.size(); i++) {
		ofsm << vecResults[i].GeneName << '\t' << vecResults[i].Target << '\t'
				<< vecResults[i].Value << endl;
	}
}
;

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

	CPCL PCL;
	SVMLight::CSVMPERF SVM;

	size_t i, j, iGene, jGene;
	ifstream ifsm;
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
	size_t iFile;
	vector<string> PCLs;
	if (sArgs.input_given) {
		if (!PCL.Open(sArgs.input_arg, sArgs.skip_arg, sArgs.mmap_flag)) {
			cerr << "Could not open input PCL" << endl;
			return 1;
		}
	}

	vector<SVMLight::SVMLabel> vecLabels;
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

	SVMLight::SAMPLE* pTrainSample;
	vector<SVMLight::SVMLabel> pTrainVector[sArgs.cross_validation_arg];
	vector<SVMLight::SVMLabel> pTestVector[sArgs.cross_validation_arg];
	vector<SVMLight::Result> AllResults;
	vector<SVMLight::Result> tmpAllResults;

	if (sArgs.model_given && sArgs.labels_given) { //learn once and write to file
		pTrainSample = CSVMPERF::CreateSample(PCL, vecLabels);
		SVM.Learn(*pTrainSample);
		SVM.WriteModel(sArgs.model_arg);
	} else if (sArgs.model_given && sArgs.output_given) { //read model and classify all
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
			for (i = 0; i < sArgs.cross_validation_arg; i++)
				ppTrainSample[i] = SVMLight::CSVMPERF::CreateSample(PCL,
						pTrainVector[i]);
			
			for (iParams = 0; iParams < PStruct.vecTradeoff.size(); iParams++) {
				SVM.SetLossFunction(PStruct.vecLoss[iParams]);
				SVM.SetTradeoff(PStruct.vecTradeoff[iParams]);
				SVM.SetPrecisionFraction(PStruct.vecK[iParams]);
				for (j = 0; j < vec_allUnlabeledResults.size(); j++)
					vec_allUnlabeledResults[j].Value = 0;
				for (i = 0; i < sArgs.cross_validation_arg; i++) {
					cerr << "Cross Validation Trial " << i << endl;
					SVM.Learn(*ppTrainSample[i]);
					
					cerr << "Learned" << endl;					
					
					tmpAllResults = SVM.Classify(PCL, pTestVector[i]);
					cerr << "Classified " << tmpAllResults.size()
							<< " examples" << endl;
					AllResults.insert(AllResults.end(), tmpAllResults.begin(),
							tmpAllResults.end());
					tmpAllResults.resize(0);
					if (sArgs.all_flag && vec_allUnlabeledLabels.size() > 0) {
						vec_tmpUnlabeledResults = SVM.Classify(PCL,
								vec_allUnlabeledLabels);
						for (j = 0; j < vec_tmpUnlabeledResults.size(); j++)
							vec_allUnlabeledResults[j].Value
									+= vec_tmpUnlabeledResults[j].Value;
					}

				}


				ofsm.open(PStruct.vecNames[iParams]);
				if (sArgs.all_flag) { //add the unlabeled results
					for (j = 0; j < vec_tmpUnlabeledResults.size(); j++)
						vec_allUnlabeledResults[j].Value
								/= sArgs.cross_validation_arg;
					AllResults.insert(AllResults.end(),
							vec_allUnlabeledResults.begin(),
							vec_allUnlabeledResults.end());
				}

				PrintResults(AllResults, ofsm);
				ofsm.close();
				ofsm.clear();
				if (i > 0 || iParams > 0)
					SVM.FreeModel();
				AllResults.resize(0);
			}
		} else { //run once
			for (i = 0; i < sArgs.cross_validation_arg; i++) {
				pTrainSample = SVMLight::CSVMPERF::CreateSample(PCL,
						pTrainVector[i]);

				cerr << "Cross Validation Trial " << i << endl;

				SVM.Learn(*pTrainSample);
				cerr << "Learned" << endl;
				tmpAllResults = SVM.Classify(PCL,
						pTestVector[i]);
				cerr << "Classified " << tmpAllResults.size() << " examples"
						<< endl;
				AllResults.insert(AllResults.end(), tmpAllResults.begin(),
						tmpAllResults.end());
				tmpAllResults.resize(0);
				if (sArgs.all_flag) {
					vec_tmpUnlabeledResults = SVM.Classify(
							PCL, vec_allUnlabeledLabels);
					for (j = 0; j < vec_tmpUnlabeledResults.size(); j++)
						vec_allUnlabeledResults[j].Value
								+= vec_tmpUnlabeledResults[j].Value;

				}
				if (i > 0) {
					SVMLight::CSVMPERF::FreeSample(*pTrainSample);
				}
			}

			if (sArgs.all_flag) { //add the unlabeled results
				for (j = 0; j < vec_allUnlabeledResults.size(); j++)
					vec_allUnlabeledResults[j].Value
							/= sArgs.cross_validation_arg;
				AllResults.insert(AllResults.end(),
						vec_allUnlabeledResults.begin(),
						vec_allUnlabeledResults.end());
			}

			ofstream ofsm;
			ofsm.clear();
			ofsm.open(sArgs.output_arg);
			PrintResults(AllResults, ofsm);
			return 0;
		}
	} else {
		cerr << "More options are needed" << endl;
	}

}

