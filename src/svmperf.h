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

#ifndef NO_SVM_PERF
#ifndef SVMPERFI_H
#define SVMPERFI_H
#include "pclset.h"
#include "meta.h"
#include "dat.h"

#include <stdio.h>

/* removed to support cygwin */
//#include <execinfo.h>

namespace SVMLight {
extern "C" {

#define class Class

#include <svm_light/svm_common.h>
#include <svm_light/svm_learn.h>
#include <svm_struct_api_types.h>
#include <svm_struct/svm_struct_common.h>
#include <svm_struct_api.h>
#include <svm_struct/svm_struct_learn.h>
#undef class
//#include "svm_struct_api.h"

}

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

class SVMLabelPair {
public:
	double Target;
	size_t iidx;
	size_t jidx;
	bool hasIndex;
	DOC* pDoc;
	
	SVMLabelPair(double target, size_t i, size_t j) {
		Target = target;
		hasIndex = true;
		iidx = i;
		jidx = j;
	}

	SVMLabelPair() {
		Target = 0;
		hasIndex = false;
	}
	
	void SetIndex(size_t i, size_t j) {
		iidx = i;
		jidx = j;
		hasIndex = true;
	}
	
	void SetDoc(DOC* inDoc) {
	  pDoc = inDoc;
	}
	DOC* GetDoc() {
	  return pDoc;
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

class CSVMPERF {
public:
	LEARN_PARM learn_parm;
	KERNEL_PARM kernel_parm;
	STRUCT_LEARN_PARM struct_parm;
	STRUCTMODEL structmodel;
	int Alg;
	CSVMPERF(int a = 3) {
		Alg = a;
		initialize();
		//set_struct_verbosity(5);
	}

	void SetLossFunction(size_t loss_f) {
		struct_parm.loss_function = loss_f;
	}

	void SetTradeoff(double tradeoff) {
		struct_parm.C = tradeoff;

	}
	void SetKernel(int K) {
		kernel_parm.kernel_type = K;
	}
	void SetPolyD(int D) {
		kernel_parm.poly_degree = D;
	}

	void UseCPSP() {
		Alg = 9;
		struct_parm.preimage_method = 2;
		struct_parm.sparse_kernel_size = 500;
		struct_parm.bias = 0;
	}

	void SetRBFGamma(double g) {
		kernel_parm.rbf_gamma = g;
		UseCPSP();
	}

	void UseSlackRescaling() {
		struct_parm.loss_type = SLACK_RESCALING;
	}

	void UseMarginRescaling() {
		struct_parm.loss_type = MARGIN_RESCALING;
	}

	void SetPrecisionFraction(double frac) {
		struct_parm.prec_rec_k_frac = frac;
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

	//creates a Doc for a given gene index in a microarray set
	static DOC* CreateDoc(Sleipnir::CPCLSet &PCLSet, size_t iGene, size_t iDoc);

	//creates a Doc for a given gene index in a single microarray
	static DOC* CreateDoc(Sleipnir::CPCL &PCL, size_t iGene, size_t iDoc);

	//creates a Doc for a given gene in a Dat file using all other genes as features
	static DOC* CreateDoc(Sleipnir::CDat& Dat, size_t iGene, size_t iDoc);

	//Creates a sample using a PCLset and SVMlabels Looks up genes by name.
	static SAMPLE* CreateSample(Sleipnir::CPCLSet &PCLSet,
			vector<SVMLabel> SVMLabels);

	//Creates a sample using a single PCL and SVMlabels Looks up genes by name.
	static SAMPLE
	* CreateSample(Sleipnir::CPCL &PCL, vector<SVMLabel> SVMLabels);

	//Same as above except creates bootstrap samples and does not duplicate data
	static SAMPLE** CreateSampleBootStrap(Sleipnir::CPCL &PCL,
			vector<SVMLabel>& SVMLabels, vector<vector<size_t> > vecvecIndex);

	//Creates a sample using a Dat and SVMlabels. Looks up genes by name
	static SAMPLE* CreateSample(Sleipnir::CDat& CDat,
			vector<SVMLabel> SMVLabels);

	//Classify single genes
	vector<Result> Classify(Sleipnir::CPCL& PCL, vector<SVMLabel> SVMLabels);
	vector<Result> Classify(Sleipnir::CPCLSet& PCLSet,
			vector<SVMLabel> SVMLabels);
	vector<Result> Classify(Sleipnir::CDat& Dat, vector<SVMLabel> SVMLabels);

	//MEMBER functions wraps learning
	void Learn(SAMPLE &sample) {
		cerr << "SLACK NORM =" << struct_parm.slack_norm << endl;
		/*  if (kernel_parm.kernel_type==CUSTOM)
		 svm_learn_struct_joint_custom(sample, &struct_parm, &learn_parm, &kernel_parm, &structmodel);
		 else*/
		//Take care of the labels here
		size_t numn, nump;
		numn = nump = 0;
		size_t i;
		for (i = 0; i < sample.examples[0].x.totdoc; i++) {
			if (sample.examples[0].y.Class[i] > 0)
				nump++;
			else
				numn++;
		}
		//make scaling appropriate for the loss function being used
		if ((struct_parm.loss_function == ZEROONE)
				|| (struct_parm.loss_function == FONE)
				|| (struct_parm.loss_function == ERRORRATE)
				|| (struct_parm.loss_function == PRBEP)
				|| (struct_parm.loss_function == PREC_K)
				|| (struct_parm.loss_function == REC_K)) {
			for (i = 0; i < sample.examples[0].x.totdoc; i++) {
				if (sample.examples[0].y.Class[i] > 0)
					sample.examples[0].y.Class[i] = 0.5 * 100.0 / (numn + nump);
				else
					sample.examples[0].y.Class[i] = -0.5 * 100.0
							/ (numn + nump);
			}
		}
		/* change label value for easy computation of rankmetrics (i.e. ROC-area) */
		if (struct_parm.loss_function == SWAPPEDPAIRS) {
			for (i = 0; i < sample.examples[0].x.totdoc; i++) {
				/*      if(sample.examples[0].y.class[i]>0)
				 sample.examples[0].y.class[i]=numn*0.5;
				 else
				 sample.examples[0].y.class[i]=-nump*0.5; */
				if (sample.examples[0].y.Class[i] > 0)
					sample.examples[0].y.Class[i] = 0.5 * 100.0 / nump;
				else
					sample.examples[0].y.Class[i] = -0.5 * 100.0 / numn;
				/*				cerr << sample.examples[0].x.doc[i]->fvec->words[0].weight
				 << '\t' << sample.examples[0].y.Class[i] << endl;*/
			}
		}
		if (struct_parm.loss_function == AVGPREC) {
			for (i = 0; i < sample.examples[0].x.totdoc; i++) {
				if (sample.examples[0].y.Class[i] > 0)
					sample.examples[0].y.Class[i] = numn;
				else
					sample.examples[0].y.Class[i] = -nump;
			}
		}
		cerr << "ALG=" << Alg << endl;
		svm_learn_struct_joint(sample, &struct_parm, &learn_parm, &kernel_parm,
				       &structmodel, Alg);
		//
	}

	// Pair learning

	//creates a Doc for a pair of genes by feature-wise multiplication
	static DOC* CreateDoc(Sleipnir::CPCL &PCL, size_t iGene, size_t jGene,
			size_t iDoc);

	//creates the complete Sample given  Dat of answers and a list of CV Genes
	static SAMPLE* CreateSample(Sleipnir::CPCL &PCL, Sleipnir::CDat& Answers,
			const vector<string>& CVGenes);

	//Classify  gene pairs from the given CV List
	void Classify(Sleipnir::CPCL& PCL, Sleipnir::CDat& Answers,
			Sleipnir::CDat& Values, Sleipnir::CDat& Counts,
			const vector<string>& CVGenes);

	//Will classify all pairs EXCEPT ones that touch the CV list
	void ClassifyAll(Sleipnir::CPCL& PCL, Sleipnir::CDat& Values,
			Sleipnir::CDat& Counts, const vector<string>& CVGenes);
	//Same as above but won't keep track of the counts, saves on memory
	void ClassifyAll(Sleipnir::CPCL& PCL, Sleipnir::CDat& Values, const vector<
			string>& CVGenes);
	bool parms_check();
	bool initialize();
	
	
	//Pair & Multiple dabs learning
	static bool CreateDoc(vector<string>& vecstrDatasets,
			      vector<SVMLabelPair*>& vecLabels,
			      const vector<string>& LabelsGene,
			      Sleipnir::CDat::ENormalize eNormalize = Sleipnir::CDat::ENormalizeNone);
	
	static SAMPLE* CreateSample(vector<SVMLabelPair*>& SVMLabels);	
	void Classify(Sleipnir::CDat &Results,
		      vector<SVMLabelPair*>& SVMLabels);
	
	// free the sample but don't free the Docs
	static void FreeSample_leave_Doc(SAMPLE s);
};
}

#endif // NO_SVM_SVMPERF
#endif // SVMPERF_H
