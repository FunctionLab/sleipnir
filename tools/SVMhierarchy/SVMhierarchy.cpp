#include <fstream>
#include <iostream>
#include <iterator>
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

using namespace SVMArc;
//#include "../../extlib/svm_light/svm_light/kernel.h"







int main(int iArgs, char** aszArgs) {
	gengetopt_args_info sArgs;

	CPCL PCL;
	SVMArc::CSVMSTRUCTTREE SVM;

	size_t i, j, k , iGene, jGene;
	double bestscore;
	;
	ifstream ifsm;
	if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
		cmdline_parser_print_help();
		return 1;
	}

	//Set Parameters
	SVM.SetLearningAlgorithm(sArgs.learning_algorithm_arg);
	SVM.SetVerbosity(sArgs.verbosity_arg);
	SVM.SetLossFunction(sArgs.loss_function_arg);
	SVM.SetNThreads(sArgs.threads_arg);
	cerr << "SetLossFunction: " <<sArgs.loss_function_arg<< endl;

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

	SVM.ReadOntology(sArgs.ontoparam_arg); // Read Ontology File

	if (!SVM.parms_check()) {
		cerr << "Parameter check not passed, see above errors" << endl;
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



	//Read labels from file
	vector<SVMArc::SVMLabel> vecLabels;
	set<string> setLabeledGenes;
	if (sArgs.labels_given) {
		ifsm.clear();
		ifsm.open(sArgs.labels_arg);
		if (ifsm.is_open())
			vecLabels = SVM.ReadLabels(ifsm);
		else {
			cerr << "Could not read label file" << endl;
			exit(1);
		}
		for (i = 0; i < vecLabels.size(); i++)
			setLabeledGenes.insert(vecLabels[i].GeneName);
	}
	cerr << "Read labels from file" << endl;

	SVM.InitializeLikAfterReadLabels();

	//Training
	SAMPLE* pTrainSample;
	vector<SVMArc::SVMLabel> pTrainVector[sArgs.cross_validation_arg];
	vector<SVMArc::SVMLabel> pTestVector[sArgs.cross_validation_arg];
	vector<SVMArc::Result> AllResults;
	vector<SVMArc::Result> tmpAllResults;

	if (sArgs.model_given && sArgs.output_given && (!sArgs.labels_given) ) { //read model and classify all
		vector<SVMLabel> vecAllLabels;

		for (size_t i = 0; i < PCL.GetGenes(); i++)
			vecAllLabels.push_back(SVMLabel(PCL.GetGene(i), 0));

		SVM.ReadModel(sArgs.model_arg);
		AllResults = SVM.Classify(PCL, vecAllLabels);
		ofstream ofsm;
		ofsm.open(sArgs.output_arg);
		if (ofsm.is_open())
			SVM.PrintResults(AllResults, ofsm);
		else {
			cerr << "Could not open output file" << endl;
		}
	} else if (sArgs.output_given && sArgs.labels_given) {
		//do learning and classifying with cross validation
		//set up training data
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
		//set up training data done

		//set up validation data
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
		//run once
		for (i = 0; i < sArgs.cross_validation_arg; i++) {
			pTrainSample = SVM.CreateSample(PCL,
				pTrainVector[i]);

			cerr << "Cross Validation Trial " << i << endl;
			SVM.Learn(*pTrainSample);
			cerr << "Learned" << endl;
			if (i > -1) {
				SVMArc::CSVMSTRUCTTREE::FreeSample(*pTrainSample);
			}
			tmpAllResults = SVM.Classify(PCL,	pTestVector[i]);
			cerr << "Classified " << tmpAllResults.size() << " examples"<< endl;



			AllResults.insert(AllResults.end(), tmpAllResults.begin(), tmpAllResults.end());
			tmpAllResults.resize(0);

			if(i == (sArgs.cross_validation_arg-1)){
				if (sArgs.all_flag || sArgs.model_given ) {
					if(sArgs.cross_validation_arg!=1){
						pTrainSample = SVM.CreateSample(PCL, vecLabels);
						cerr << "Train with All Labeled Data " <<  endl;
						SVM.Learn(*pTrainSample);
						cerr << "Learned" << endl;
						if (i > -1) {
							SVMArc::CSVMSTRUCTTREE::FreeSample(*pTrainSample);
						}
					}
					if (sArgs.model_given ){  //learn once and write to file
						SVM.WriteModel(sArgs.model_arg);
						cerr <<" Model Writen to file "<<sArgs.model_arg<<endl;
					}
					if(sArgs.all_flag){
						vec_allUnlabeledResults = SVM.Classify(PCL, vec_allUnlabeledLabels);
						cerr << "Classified " << vec_allUnlabeledResults.size() << " examples"<< endl;
					}

				}
			}






		}

		if (sArgs.all_flag) { //add the unlabeled results

			AllResults.insert(AllResults.end(),
				vec_allUnlabeledResults.begin(),
				vec_allUnlabeledResults.end());
		}

		ofstream ofsm;
		ofsm.clear();
		ofsm.open(sArgs.output_arg);
		SVM.PrintResults(AllResults, ofsm);
		return 0;

	} else {
		cerr << "More options are needed" << endl;
	}

}

