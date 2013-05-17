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

#ifndef NO_SVM_STRUCT
#ifndef SVMSTRUCTI_H
#define SVMSTRUCTI_H
#include "pclset.h"
#include "meta.h"
#include "dat.h"
#ifndef NO_SVM_STRUCT
#define SVMSTRUCT_H
extern "C" {

#define class Class

#include <svm_multiclass/svm_light/svm_common.h>
#include <svm_multiclass/svm_light/svm_learn.h>
#include <svm_multiclass/svm_struct_api_types.h>
#include <svm_multiclass/svm_struct/svm_struct_common.h>
#include <svm_multiclass/svm_struct_api.h>
#include <svm_multiclass/svm_struct/svm_struct_learn.h>
#undef class
	//#include "svm_struct_api.h"

}
#endif

#include <stdio.h>
using namespace Sleipnir;
using namespace std;

/* removed to support cygwin */
//#include <execinfo.h>

namespace SVMArc {
	class SVMLabel {
	public:
		string GeneName;
		size_t Target; //Save single integer label; used for single label classification (0-1, or multiclass)
		vector<char> TargetM; //Save multiple labels; used for hierarchical multi-label classification;

		size_t index;
		bool hasIndex;
		SVMLabel(std::string name, size_t target) {
			GeneName = name;
			Target = target;
			hasIndex = false;
			index = -1;
		}

		SVMLabel(std::string name, vector<char> cl) {
			GeneName = name;
			TargetM = cl;
			hasIndex = false;
			index = -1;
		}
		SVMLabel() {
			GeneName = "";
			Target = 0;
		}
		void SetIndex(size_t i) {
			index = i;
			hasIndex = true;
		}
	};

	class Result {
	public:
		std::string GeneName;
		int Target; //for single label prediction
		int Value; //for single label prediction
		vector<char> TargetM;//for multi label prediction
		vector<char> ValueM; //for multi label prediction
		vector<double> Scores;
		int num_class;
		int CVround;
		int Rank;
		Result() {
			GeneName = "";
			Target = 0;
			Value = -1;
		}

		Result(std::string name, int cv = -1) {
			GeneName = name;
			Target = 0;
			Value = 0;
			CVround = cv;
			Rank = -1;
			num_class = 0;

		}
		string toString() {
			stringstream ss;
			ss << GeneName << '\t' << Target << '\t' << Value << '\t' << "CV"
				<< CVround;
			if (Rank != -1) {
				ss << '\t' << Rank;
			}
			return ss.str();
		}
		string toStringMC() {
			stringstream ss;
			ss << GeneName << '\t' << Target << '\t' << Value << '\t';
			for(size_t j=1;j<=num_class;j++)
				ss << Scores[j]<<'\t';
			return ss.str();
		}
		string toStringTREE(map<int, string>* ponto_map_rev, int returnindex) {
			stringstream ss;
			int mark=1;
			ss << GeneName << '\t';
			for(size_t j=0;j<num_class;j++){
				if(TargetM[j])
					if(mark){
						if(returnindex)
							ss<<j;
						else
							ss <<(*ponto_map_rev)[j];
						mark = 0;
					}
					else
						ss <<','<<(*ponto_map_rev)[j];
			}
			if(mark==1)
				ss<<"??"<<'\t';
			else
				ss<<'\t';

			mark=1;
			for(size_t j=0;j<num_class;j++){
				if(ValueM[j])
					if(mark){
						if(returnindex)
							ss<<j;
						else
							ss <<(*ponto_map_rev)[j];
						mark = 0;
					}
					else
						ss <<','<<(*ponto_map_rev)[j];
			}
			if(mark)
				ss<<"??";
			ss <<'\t';
			for(size_t j=0;j<num_class;j++)
				ss << Scores[j]<<'\t';
			return ss.str();
		}
	};

	enum EFilter {
		EFilterInclude = 0, EFilterExclude = EFilterInclude + 1,
	};

	class CSVMSTRUCT{
		/* This base class is solely intended to serve as a common template for different SVM Struct implementations
		A few required functions are not defined here because their parameter type or return type has to differ 
		among different implementations, but I listed them in comments. */
	public:
		virtual vector<Result> Classify(Sleipnir::CPCL& PCL, vector<SVMLabel> SVMLabels) = 0;
		virtual void SetTradeoff(double tradeoff)=0;
		virtual void SetLossFunction(size_t loss_f)=0;
		virtual void SetLearningAlgorithm(int alg)=0;
		virtual void UseSlackRescaling()=0;
		virtual void UseMarginRescaling()=0;
		virtual void ReadModel(char* model_file)=0;
		virtual void WriteModel(char* model_file)=0;
		virtual vector<SVMLabel> ReadLabels(ifstream & ifsm)=0;
		virtual void SetVerbosity(size_t V)=0;
		virtual bool parms_check() = 0;
		virtual bool initialize() = 0;

		/*The following functions should also be implemented
		SAMPLE* CreateSample(Sleipnir::CPCL &PCL, vector<SVMLabel> SVMLabels);
		static void FreeSample(sample s)
		void Learn(SAMPLE &sample)
		*/
	};




	//this class encapsulates the model and parameters and has no associated data


	//class for SVMStruct
	class CSVMSTRUCTMC : CSVMSTRUCT{

	public:
		LEARN_PARM learn_parm;
		KERNEL_PARM kernel_parm;
		STRUCT_LEARN_PARM struct_parm;
		STRUCTMODEL structmodel;
		int Alg;
		CSVMSTRUCTMC() {
			initialize();
			//set_struct_verbosity(5);
		}

		void SetLossFunction(size_t loss_f) {
			struct_parm.loss_function = loss_f;
		}

		void SetTradeoff(double tradeoff) {
			struct_parm.C = tradeoff;
		}
		void SetLearningAlgorithm(int alg) {
			Alg = alg;
		}
		void SetKernel(int K) {
			kernel_parm.kernel_type = K;
		}
		void SetPolyD(int D) {
			kernel_parm.poly_degree = D;
		}

		//void UseCPSP() {
		//	Alg = 9;
		//	struct_parm.preimage_method = 2;
		//	struct_parm.sparse_kernel_size = 500;
		//	struct_parm.bias = 0;
		//}

		//void SetRBFGamma(double g) {
		//	kernel_parm.rbf_gamma = g;
		//	UseCPSP();
		//}

		void UseSlackRescaling() {
			struct_parm.loss_type = SLACK_RESCALING;
		}

		void UseMarginRescaling() {
			struct_parm.loss_type = MARGIN_RESCALING;
		}



		void ReadModel(char* model_file) {
			FreeModel();
			structmodel = read_struct_model(model_file, &struct_parm);
		}

		void WriteModel(char* model_file) {
			if (kernel_parm.kernel_type == LINEAR) {
				ofstream ofsm;
				ofsm.open(model_file);
				for (size_t i = 0; i < structmodel.sizePsi; i++) {
					ofsm << structmodel.w[i+1] << endl;
				}
			} else {
				write_struct_model(model_file, &structmodel, &struct_parm);
			}
		}

		void WriteWeights(ostream& osm) {
			osm << structmodel.w[0];
			for (size_t i = 1; i < structmodel.sizePsi + 1; i++)
				osm << '\t' << structmodel.w[i];
			osm << endl;
		}

		static void FreePattern(pattern x) {
			free_pattern(x);
		}

		static void FreeLabel(label y) {
			free_label(y);
		}

		void FreeModel() {
			free_struct_model(structmodel);
		}

		static void FreeSample(sample s) {
			free_struct_sample(s);
		}

		static void FreeDoc(DOC* pDoc) {
			free_example(pDoc, true);
		}
		void SetVerbosity(size_t V);

		//static members process data
		//single gene predictions


		//creates a Doc for a given gene index in a single microarray
		static DOC* CreateDoc(Sleipnir::CPCL &PCL, size_t iGene, size_t iDoc);


		//read in labels
		vector<SVMLabel> ReadLabels(ifstream & ifsm);

		//Creates a sample using a single PCL and SVMlabels Looks up genes by name.
		SAMPLE
			* CreateSample(Sleipnir::CPCL &PCL, vector<SVMLabel> SVMLabels);

		//Classify single genes
		vector<Result> Classify(Sleipnir::CPCL& PCL, vector<SVMLabel> SVMLabels);

		//MEMBER functions wraps learning
		void Learn(SAMPLE &sample) {
			cerr << "SLACK NORM =" << struct_parm.slack_norm << endl;
			/*  if (kernel_parm.kernel_type==CUSTOM)
			svm_learn_struct_joint_custom(sample, &struct_parm, &learn_parm, &kernel_parm, &structmodel);
			else*/


			cerr << "ALG=" << Alg << endl;

			if(Alg == 0)
				svm_learn_struct(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,NSLACK_ALG);
			else if(Alg == 1)
				svm_learn_struct(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,NSLACK_SHRINK_ALG);
			else if(Alg == 2)
				svm_learn_struct_joint(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,ONESLACK_PRIMAL_ALG);
			else if(Alg == 3)
				svm_learn_struct_joint(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,ONESLACK_DUAL_ALG);
			else if(Alg == 4)
				svm_learn_struct_joint(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,ONESLACK_DUAL_CACHE_ALG);
			else if(Alg == 9)
				svm_learn_struct_joint_custom(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel);
			else
				exit(1);
			//
		}

		struct SortResults {

			bool operator()(const Result& rOne, const Result & rTwo) const {
				return (rOne.Value < rTwo.Value);
			}
		};

		size_t PrintResults(vector<Result> vecResults, ofstream & ofsm) {
			sort(vecResults.begin(), vecResults.end(), SortResults());
			int LabelVal;
			for (size_t i = 0; i < vecResults.size(); i++) {
				ofsm << vecResults[i].GeneName << '\t' << vecResults[i].Target << '\t'
					<< vecResults[i].Value<<'\t';
				for(size_t j=1;j<=vecResults[i].num_class;j++)
					ofsm << vecResults[i].Scores[j]<<'\t';
				ofsm<< endl;

			}
		};

		bool parms_check();
		bool initialize();



		// free the sample but don't free the Docs
		static void FreeSample_leave_Doc(SAMPLE s);



		STRUCTMODEL read_struct_model_w_linear(char *file, STRUCT_LEARN_PARM *sparm);
	};


};


#endif // NO_SVM_SVMSTRUCT
#endif // SVMSTRUCT_H
