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

#ifndef NO_LIBSVM
#ifndef LIBSVMI_H
#define LIBSVMI_H
#include "pclset.h"
#include "meta.h"
#include "dat.h"

#include <stdio.h>

/* removed to support cygwin */
//#include <execinfo.h>

//#include <svm.h>

namespace LIBSVM {


extern "C" {
#define class Class2
#include <libsvm/svm.h>
#undef class

}

typedef struct sample { /* a sample is a set of examples */
   size_t     n;            /* n is the total number of examples */
   size_t  numFeatures; 
   struct svm_problem *problems;
//   struct svm_node *x_space;
   sample() {
     n = 0;
     numFeatures = 0;
     problems = NULL;
//     x_space = NULL;
   }
   
   ~sample(){
     //no destructor for problem struct
     free(problems->y);
     free(problems->x);
//     free(x_space);
     problems = NULL;
   }
} SAMPLE;
 

class SVMLabel {
public:
	string GeneName;
	double Target;
	size_t index;
	bool hasIndex;
        
	SVMLabel(std::string name, double target) {
		GeneName = name;
		Target = target;
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
	double Target;
	double Value;
	int CVround;
	int Rank;
	Result() {
		GeneName = "";
		Target = 0;
		Value = Sleipnir::CMeta::GetNaN();
	}

	Result(std::string name, int cv = -1) {
		GeneName = name;
		Target = 0;
		Value = 0;
		CVround = cv;
		Rank = -1;
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

};

enum EFilter {
	EFilterInclude = 0, EFilterExclude = EFilterInclude + 1,
};

//this class encapsulates the model and parameters and has no associated data

class CLIBSVM {
public:
  //struct svm_parameter parm;
  struct svm_model* model;
  struct svm_parameter parm;
  int balance;

  static struct svm_node *x_space;

  CLIBSVM() {
    initialize();
  }

  ~CLIBSVM() {
    svm_free_and_destroy_model( &model );
    model = NULL;
  }

  void SetBalance(int bal){
    balance = bal;
  }

  void SetSVMType(int type) {
    parm.svm_type = type;
  }

  void SetTradeoff(double tradeoff) {
    parm.C = tradeoff; //TODO: only applicable for vanilla svm
  }

  void SetKernel(int K) {
    parm.kernel_type = K;
  }

  void SetPolyD(int D) {
    parm.degree = D;
  }

  void SetRBFGamma(double g) {
    parm.gamma = g;
  }

  void SetNu(double nu) {
    parm.nu = nu;
  }

  void ReadModel(char* model_file) {
    FreeModel();
    model = svm_load_model(model_file); 
  }

  void FreeModel() {
    svm_free_and_destroy_model(&model);
  }

  void WriteModel(char* model_file) {
    svm_save_model(model_file, model);
  }
  

  //static members process data
  //
  
  static void SetXSpace(Sleipnir::CPCL& PCL);

  //
  //single gene predictions

  //TODO: add functions to handle PCL files
  //creates a svm_problem for a given gene index in a microarray set
  //static svm_problem* CreateProblem(Sleipnir::CPCLSet &PCLSet, size_t iGene, size_t iProblem);

  //creates a svm_problem for a given gene in a Dat file using all other genes as features
  //static svm_problem* CreateProblem(Sleipnir::CDat& Dat, size_t iGene, size_t iProblem);

  //Creates a sample using a PCLset and SVMlabels Looks up genes by name.
  //static SAMPLE* CreateSample(Sleipnir::CPCLSet &PCLSet,
  //			vector<SVMLabel> SVMLabels);
  
  //Creates a sample of svm_problems using a single PCL and SVMlabels Looks up genes by name.
  static SAMPLE* CreateSample(Sleipnir::CPCL &PCL, vector<SVMLabel> SVMLabels);

  //Same as above except creates bootstrap samples and does not duplicate data
  //static SAMPLE** CreateSampleBootStrap(Sleipnir::CPCL &PCL,
  //	vector<SVMLabel>& SVMLabels, 
  //      vector<vector<size_t> > vecvecIndex);

  //Creates a sample using a Dat and SVMlabels. Looks up genes by name
  static SAMPLE* CreateSample(Sleipnir::CDat& CDat,
			vector<SVMLabel> SMVLabels);

  //Classify single genes
  vector<Result> Classify(Sleipnir::CPCL& PCL, vector<SVMLabel> SVMLabels);
  
  //vector<Result> Classify(Sleipnir::CPCLSet& PCLSet,
  //			vector<SVMLabel> SVMLabels);
  //vector<Result> Classify(Sleipnir::CDat& Dat, vector<SVMLabel> SVMLabels);

  //MEMBER functions wraps learning
  void Learn(SAMPLE &sample) {
    //only L2 for LibSVM
    //cerr << "SLACK NORM =" << struct_parm.slack_norm << endl;
    //slack_norm = type of regularization

    //Take care of the labels here
    size_t i;
    size_t numn, nump;

    struct svm_problem* prob = sample.problems;

    numn = nump = 0;

    for(i = 0; i < sample.n; i++){
      if (((*prob).y)[i] > 0){
        nump ++;
      }else{
        numn ++;
      }
    }

    if (balance) {
cerr << "balancing the weights between postivies and negatives. " << endl;
      parm.nr_weight = 2;
      parm.weight_label = (int *) realloc(parm.weight_label, sizeof(int)*parm.nr_weight);
      parm.weight = (double *) realloc(parm.weight, sizeof(double)*parm.nr_weight);
      parm.weight_label[0] = 1;
      parm.weight[0] = numn;
      parm.weight_label[1] = -1;
      parm.weight[1] = nump;
    }

    if(parms_check()){
      model = svm_train(prob,&parm);
    }else{
    }
    prob = NULL;

  }
/*
  static void FreeSample(SAMPLE s){
    FreeProblem(s.problems);
    free(s.x_space);
  }

  static void FreeProblem(svm_problem *prob){
    int i ;
    free(prob->y);
    free(prob->x);
    return;
  }
*/
  static void PrintSample(SAMPLE s){
    PrintProblem(s.problems);
  }

  static void PrintProblem(svm_problem *prob){
    size_t i, j ;
    i = j = 0;

    for(i = 0 ; i < 3 ; i++){
      for(j = 0 ; j < 2 ; j ++){
        PrintNode((prob->x)[i][j]);
      }
    }

    return;
  }

  static void PrintNode(svm_node node){
    cerr << "index: " << node.index << endl;
    cerr << "value: " << node.value << endl;
  }


  //no pairwise learning for libSVM wrapper

  bool parms_check();
  bool initialize();
	
  // functions to convert probablity
  //void sigmoid_train(Sleipnir::CDat& Results, vector<SVMLabelPair*>& SVMLabels, float& A, float& B);
  //void sigmoid_predict(Sleipnir::CDat& Results, vector<SVMLabelPair*>& SVMLabels, float A, float B);
        
        //not sure exactly what this does in svmperf compare to just ReadModel
	// read in a SVM model file that's only has the w vector written out for linear kernel
/*	void ReadModelLinear(char* model_file) {
	  FreeModel();
	  structmodel = read_struct_model_w_linear(model_file, &struct_parm);
	}*/
	
	//STRUCTMODEL read_struct_model_w_linear(char *file, STRUCT_LEARN_PARM *sparm);
};
}

#endif // NO_SVM_LIBSVM
#endif // LIBSVM_H
