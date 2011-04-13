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
#include "svmperf.h"
#include "pclset.h"
#include "dataset.h"
#include "meta.h"
#include "genome.h"
#include "compactmatrix.h"
#include <vector>
#include <set>

#define  SLACK_RESCALING    1
#define  MARGIN_RESCALING   2

namespace SVMLight {
#include "svmperf.h"
extern "C" {
//    void free_struct_model(STRUCTMODEL sm);
void free_struct_sample(SAMPLE s);
//    void svm_learn_struct_joint_custom(SAMPLE sample,
//            STRUCT_LEARN_PARM *sparm,
//            LEARN_PARM *lparm, KERNEL_PARM *kparm,
//            STRUCTMODEL *sm);
//    SAMPLE read_struct_examples_sleipnir(DOC **all_docs, double*all_labels, int example_size, int total_features, STRUCT_LEARN_PARM *sparm);
//    void free_struct_model(STRUCTMODEL sm);
//    void free_struct_sample(SAMPLE s);
//    void set_struct_verbosity(long verb);
//    double estimate_r_delta_average(DOC **, long, KERNEL_PARM *);
//    MODEL *read_model(char *);
LABEL classify_struct_example(PATTERN x, STRUCTMODEL *sm,
		STRUCT_LEARN_PARM *sparm);
DOC* create_example(long, long, long, double, SVECTOR *);
SVECTOR * create_svector(WORD *, char *, double);
void set_struct_verbosity(long verb);

}

void CSVMPERF::SetVerbosity(size_t V) {
	struct_verbosity = (long) V;
}

bool CSVMPERF::initialize() {

	//set directionality


	/* set default */

	//Learn_parms
	struct_parm.C = 0.01;
	struct_parm.slack_norm = 1;
	struct_parm.epsilon = DEFAULT_EPS;
	struct_parm.custom_argc = 0;
	struct_parm.loss_function = 4;
	struct_parm.loss_type = DEFAULT_RESCALING;
	struct_parm.newconstretrain = 100;
	struct_parm.ccache_size = 5;
	struct_parm.batch_size = 100;

	//Learn_parms
	//strcpy (learn_parm.predfile, "trans_predictions");
	strcpy(learn_parm.alphafile, "");
	//verbosity=0;/*verbosity for svm_light*/
	//struct_verbosity = 1; /*verbosity for struct learning portion*/
	learn_parm.biased_hyperplane = 1;
	learn_parm.remove_inconsistent = 0;
	learn_parm.skip_final_opt_check = 0;
	learn_parm.svm_maxqpsize = 10;
	learn_parm.svm_newvarsinqp = 0;
	learn_parm.svm_iter_to_shrink = -9999;
	learn_parm.maxiter = 1000;
	learn_parm.kernel_cache_size = 40;
	learn_parm.svm_c = 99999999; /* overridden by struct_parm->C */
	learn_parm.eps = 0.001; /* overridden by struct_parm->epsilon */
	learn_parm.transduction_posratio = -1.0;
	learn_parm.svm_costratio = 1.0;
	learn_parm.svm_costratio_unlab = 1.0;
	learn_parm.svm_unlabbound = 1E-5;
	learn_parm.epsilon_crit = 0.001;
	learn_parm.epsilon_a = 1E-15; /* changed from 1e-15 */
	learn_parm.compute_loo = 0;
	learn_parm.rho = 1.0;
	learn_parm.xa_depth = 0;

	//kernel parms
	kernel_parm.kernel_type = 0;
	kernel_parm.poly_degree = 3;
	kernel_parm.rbf_gamma = 1.0;
	kernel_parm.coef_lin = 1;
	kernel_parm.coef_const = 1;
	strcpy(kernel_parm.custom, "empty");

	if (learn_parm.svm_iter_to_shrink == -9999) {
		learn_parm.svm_iter_to_shrink = 100;
	}

	if ((learn_parm.skip_final_opt_check)
			&& (kernel_parm.kernel_type == LINEAR)) {
		printf(
				"\nIt does not make sense to skip the final optimality check for linear kernels.\n\n");
		learn_parm.skip_final_opt_check = 0;
	}

	//struct parms
	struct_parm.bias = 0;
	struct_parm.prec_rec_k_frac = 0.5;
	struct_parm.sparse_kernel_type = LINEAR;
	struct_parm.sparse_kernel_size = 500;
	struct_parm.shrinking = 1;
	strcpy(struct_parm.sparse_kernel_file, "");
	/* set number of features to -1, indicating that it will be computed
	 in init_struct_model() */
	struct_parm.num_features = -1;

	return true;
}

bool CSVMPERF::parms_check() {
	if ((learn_parm.skip_final_opt_check) && (learn_parm.remove_inconsistent)) {
		fprintf(
				stderr,
				"\nIt is necessary to do the final optimality check when removing inconsistent \nexamples.\n");
		return false;
	}
	if ((learn_parm.svm_maxqpsize < 2)) {
		fprintf(
				stderr,
				"\nMaximum size of QP-subproblems not in valid range: %ld [2..]\n",
				learn_parm.svm_maxqpsize);
		return false;
	}
	if ((learn_parm.svm_maxqpsize < learn_parm.svm_newvarsinqp)) {
		fprintf(
				stderr,
				"\nMaximum size of QP-subproblems [%ld] must be larger than the number of\n",
				learn_parm.svm_maxqpsize);
		fprintf(
				stderr,
				"new variables [%ld] entering the working set in each iteration.\n",
				learn_parm.svm_newvarsinqp);
		return false;
	}
	if (learn_parm.svm_iter_to_shrink < 1) {
		fprintf(
				stderr,
				"\nMaximum number of iterations for shrinking not in valid range: %ld [1,..]\n",
				learn_parm.svm_iter_to_shrink);
		return false;
	}
	if (struct_parm.C < 0) {
		fprintf(
				stderr,
				"\nTrade-off between training error and margin is not set (C<0)!\nC value will be set to default value. Clight = Cpef * 100 / n \n");
	}
	if (learn_parm.transduction_posratio > 1) {
		fprintf(stderr,
				"\nThe fraction of unlabeled examples to classify as positives must\n");
		fprintf(stderr, "be less than 1.0 !!!\n\n");
		return false;
	}
	if (learn_parm.svm_costratio <= 0) {
		fprintf(stderr,
				"\nThe COSTRATIO parameter must be greater than zero!\n\n");
		return false;
	}
	if (struct_parm.epsilon <= 0) {
		fprintf(stderr,
				"\nThe epsilon parameter must be greater than zero!\n\n");
		return false;
	}
	if ((struct_parm.slack_norm < 1) || (struct_parm.slack_norm > 2)) {
		fprintf(stderr,
				"\nThe norm of the slacks must be either 1 (L1-norm) or 2 (L2-norm)!\n\n");
		return false;
	}

	if ((struct_parm.loss_type != SLACK_RESCALING) && (struct_parm.loss_type
			!= MARGIN_RESCALING)) {
		fprintf(
				stderr,
				"\nThe loss type must be either 1 (slack rescaling) or 2 (margin rescaling)!\n\n");
		return false;
	}

	if (learn_parm.rho < 0) {
		fprintf(stderr,
				"\nThe parameter rho for xi/alpha-estimates and leave-one-out pruning must\n");
		fprintf(stderr,
				"be greater than zero (typically 1.0 or 2.0, see T. Joachims, Estimating the\n");
		fprintf(stderr,
				"Generalization Performance of an SVM Efficiently, ICML, 2000.)!\n\n");
		return false;
	}
	if ((learn_parm.xa_depth < 0) || (learn_parm.xa_depth > 100)) {
		fprintf(stderr,
				"\nThe parameter depth for ext. xi/alpha-estimates must be in [0..100] (zero\n");
		fprintf(stderr,
				"for switching to the conventional xa/estimates described in T. Joachims,\n");
		fprintf(
				stderr,
				"Estimating the Generalization Performance of an SVM Efficiently, ICML, 2000.)\n");
		return false;
	}

	/* Note that the validity of the value for struct_parm->prec_rec_k_frac in
	 relation to #pos is checked in read_struct_examples() */
	if (struct_parm.prec_rec_k_frac < 0) {
		fprintf(stderr,
				"\nThe value of option --k must be greater then zero!\n\n");
		return false;
	}

	return true;
}

DOC* CSVMPERF::CreateDoc(Sleipnir::CPCLSet &PCLSet, size_t iGene, size_t iDoc) {
	WORD* aWords;
	size_t i, j, iWord, iWords, iPCL, iExp;
	float d;
	DOC* pRet;
	pRet->fvec->words[0].weight;
	//get number of features
	iWords = 0;
	for (i = 0; i < PCLSet.GetPCLs(); i++) {
		//	  cout<<"CD:PCLSET= "<<i<<endl;
		//  cout<<"CD:numExp= "<<PCLSet.Get(i).GetExperiments()<<endl;
		iWords += PCLSet.Get(i).GetExperiments();
	}
	//      cout << "CD:iwords=" << iWords << endl;
	aWords = new WORD[iWords + 1];
	//number the words
	for (i = 0; i < iWords; ++i) {
		//   cout<<i<<endl;
		aWords[i].wnum = i + 1;
		// asWords[ i ].wnum = 0;
	}
	aWords[i].wnum = 0;
	//get the values;
	iWord = 0;
	for (i = 0; i < PCLSet.GetPCLs(); i++) {
		iExp = PCLSet.Get(i).GetExperiments();
		for (j = 0; j < iExp; j++) {
			//     cout<<"CD:iWord="<<iWord<<endl;
			if (!Sleipnir::CMeta::IsNaN(d = PCLSet.Get(i, iGene, j))) {
				//   if (i==0 && j==0)
				//       cout<<"First value is "<<d<<endl;
				aWords[iWord].weight = d;
			} else
				aWords[iWord].weight = 0;
			iWord++;
		}
	}
	pRet = create_example(iDoc, 0, 0, 1, create_svector(aWords, "", 1));
	delete[] aWords;
	// cout<<"done creating DOC"<<endl;
	return pRet;
}
//For single genes usign single PCL

DOC* CSVMPERF::CreateDoc(Sleipnir::CPCL &PCL, size_t iGene, size_t iDoc) {
	WORD* aWords;
	size_t i, j, iWord, iWords, iPCL, iExp;
	float d;
	DOC* pRet;
	pRet->fvec->words[0].weight;
	//get number of features
	iWords = PCL.GetExperiments();
	//cerr<<"Newing WORDS "<<(iWords+1)*sizeof(WORD)<<endl;
	aWords = new WORD[iWords + 1];
	//set the words
	for (i = 0; i < iWords; ++i) {
		aWords[i].wnum = i + 1;
		if (!Sleipnir::CMeta::IsNaN(d = PCL.Get(iGene, i)))
			aWords[i].weight = d;
		else
			aWords[i].weight = 0;
	}
	aWords[i].wnum = 0;
	// cerr<<"START Create Example"<<endl;
	pRet = create_example(iDoc, 0, 0, 1, create_svector(aWords, "", 1));
	//cerr<<"END create example"<<endl;
	delete[] aWords;
	return pRet;
}
//Single Gene using a Dat for data

DOC* CSVMPERF::CreateDoc(Sleipnir::CDat& Dat, size_t iGene, size_t iDoc) {
	WORD* aWords;
	size_t i, j, iWord, iWords;
	float d;
	DOC* pRet;
	pRet->fvec->words[0].weight;
	//get number of features
	iWords = Dat.GetGenes();
	//      cout << "CD:iwords=" << iWords << endl;
	aWords = new WORD[iWords + 1];
	//number the words
	for (i = 0; i < iWords; ++i) {
		//   cout<<i<<endl;
		aWords[i].wnum = i + 1;
		// asWords[ i ].wnum = 0;
	}
	aWords[i].wnum = 0;
	//get the values;
	iWord = 0;
	for (i = 0; i < Dat.GetGenes(); i++) {
		if (!Sleipnir::CMeta::IsNaN(d = Dat.Get(iGene, i))) {
			//   if (i==0 && j==0)
			//       cout<<"First value is "<<d<<endl;
			aWords[iWord].weight = d;
		} else
			aWords[iWord].weight = 0;
		iWord++;
	}
	pRet = create_example(iDoc, 0, 0, 1, create_svector(aWords, "", 1));
	delete[] aWords;
	// cout<<"done creating DOC"<<endl;
	return pRet;
}

//Create Sample functions for single genes

SAMPLE* CSVMPERF::CreateSample(Sleipnir::CPCLSet &PCLSet,
		vector<SVMLabel> SVMLabels) {
	size_t i, j, iGene, iDoc;
	vector<DOC*> vec_pDoc;
	vector<double> vecClass;
	vector<size_t> veciGene;
	iDoc = 0;
	float numPositives, numNegatives;
	numPositives = numNegatives = 0;
	for (i = 0; i < SVMLabels.size(); i++) {
		//     cout<< "processing gene " << SVMLabels[i].GeneName << endl;
		if (!SVMLabels[i].hasIndex) {
			SVMLabels[i].SetIndex(PCLSet.GetGene(SVMLabels[i].GeneName));
		}
		iGene = SVMLabels[i].index;
		//   cout << SVMLabels[i].GeneName<<" gene at location "<<iGene << endl;
		if (iGene != -1) {
			//       cout << "creating doc" << endl;
			iDoc++;
			vec_pDoc.push_back(CreateDoc(PCLSet, iGene, iDoc - 1));
			vecClass.push_back(SVMLabels[i].Target);
		}
	}

	DOC** ppDoc;
	ppDoc = new DOC*[vec_pDoc.size()];
	copy(vec_pDoc.begin(), vec_pDoc.end(), ppDoc);
	vec_pDoc.clear();
	PATTERN* pPattern = new PATTERN;
	pPattern->doc = ppDoc;

	pPattern->totdoc = iDoc;
	//   cout << "number of document=" << pPattern->totdoc << endl;
	LABEL* pLabel = new LABEL;
	double* aClass;
	aClass = new double[vecClass.size()];
	copy(vecClass.begin(), vecClass.end(), aClass);
	vecClass.clear();
	pLabel->Class = aClass;
	pLabel->totdoc = iDoc;

	EXAMPLE* aExample;
	aExample = new EXAMPLE[1];
	//cout<<"aExample @"<<aExample<<endl;
	aExample[0].x = *pPattern;
	aExample[0].y = *pLabel;
	SAMPLE* pSample = new SAMPLE;
	pSample->n = 1;
	pSample->examples = aExample;
	/* cout << "examples @" << pSample->examples << endl;
	 cout<< "ppDoc="<<ppDoc<<endl;
	 cout << "docs @" << pSample->examples[0].x.doc << endl;
	 cout<<"done creating sample"<<endl;
	 cout<<"sample @ "<<pSample<<endl;*/
	return pSample;
}

SAMPLE* CSVMPERF::CreateSample(Sleipnir::CPCL &PCL, vector<SVMLabel> SVMLabels) {
	size_t i, j, iGene, iDoc;
	// cerr<<"CREATE pDoc vector"<<endl;
	vector<DOC*> vec_pDoc;
	//cerr << "RESERVE pDoc"<<endl;
	vec_pDoc.reserve(SVMLabels.size());
	//cerr<<"CREATE class"<<endl;
	vector<double> vecClass;
	//cerr<<"RESERVE class"<<endl;
	vecClass.reserve(SVMLabels.size());
	iDoc = 0;
	float numPositives, numNegatives;
	numPositives = numNegatives = 0;
	for (i = 0; i < SVMLabels.size(); i++) {
		//     cout<< "processing gene " << SVMLabels[i].GeneName << endl;
		if (!SVMLabels[i].hasIndex) {
			SVMLabels[i].SetIndex(PCL.GetGene(SVMLabels[i].GeneName));
		}
		iGene = SVMLabels[i].index;
		//   cout << SVMLabels[i].GeneName<<" gene at location "<<iGene << endl;
		if (iGene != -1) {
			//       cout << "creating doc" << endl;
			iDoc++;
			vec_pDoc.push_back(CreateDoc(PCL, iGene, iDoc - 1));
			vecClass.push_back(SVMLabels[i].Target);
		}
	}

	DOC** ppDoc;
	//cerr<<"NEW ppDoc"<<endl;
	ppDoc = new DOC*[vec_pDoc.size()];
	copy(vec_pDoc.begin(), vec_pDoc.end(), ppDoc);
	vec_pDoc.clear();
	//cerr<<"NEW Pattern"<<endl;
	PATTERN* pPattern = new PATTERN;
	pPattern->doc = ppDoc;

	pPattern->totdoc = iDoc;
	LABEL* pLabel = new LABEL;
	double* aClass;
	cerr << "NEW Class array" << endl;
	aClass = new double[vecClass.size()];
	copy(vecClass.begin(), vecClass.end(), aClass);
	vecClass.clear();
	pLabel->Class = aClass;
	pLabel->totdoc = iDoc;

	EXAMPLE* aExample;
	//cerr<<"NEW Example"<<endl;
	aExample = new EXAMPLE[1];
	//cout<<"aExample @"<<aExample<<endl;
	aExample[0].x = *pPattern;
	aExample[0].y = *pLabel;
	//cerr<<"NEW Sample"<<endl;
	SAMPLE* pSample = new SAMPLE;
	pSample->n = 1;
	pSample->examples = aExample;
	/* cout << "examples @" << pSample->examples << endl;
	 cout<< "ppDoc="<<ppDoc<<endl;
	 cout << "docs @" << pSample->examples[0].x.doc << endl;
	 cout<<"done creating sample"<<endl;
	 cout<<"sample @ "<<pSample<<endl;*/
	return pSample;
	//cerr<<"DONE CreateSample"<<endl;
}

SAMPLE** CSVMPERF::CreateSampleBootStrap(Sleipnir::CPCL &PCL,
		vector<SVMLabel>& SVMLabels, vector<vector<size_t> > vecvecIndex) {

	size_t i, j, iGene, iDoc;
	vector<DOC*> vec_pDoc;
	vec_pDoc.reserve(SVMLabels.size());
	vector<vector<double> > vecvecClass;
	vector<double> vecClass;
	vecClass.reserve(SVMLabels.size());
	iDoc = 0;
	float numPositives, numNegatives;
	numPositives = numNegatives = 0;
	//Creat all the docs once
	for (i = 0; i < SVMLabels.size(); i++) {
		//the labels will have an index
		iGene = SVMLabels[i].index;
		if (iGene != -1) {
			iDoc++;
			vec_pDoc.push_back(CreateDoc(PCL, iGene, iDoc - 1));
			vecClass.push_back(SVMLabels[i].Target);

		}
	}
	size_t numBootStraps = vecvecIndex.size();
	SAMPLE** ppSample = new SAMPLE*[numBootStraps];
	DOC** ppDoc;
	for (i = 0; i < numBootStraps; i++) {
		//get a new ppDoc
		ppDoc = new DOC*[vecvecIndex[i].size()];
		for (j = 0; j < vecvecIndex[i].size(); j++) {
			ppDoc[j] = vec_pDoc[vecvecIndex[i][j]]; //assign the pointer
		}
		//set up the pattern
		PATTERN* pPattern = new PATTERN;
		pPattern->doc = ppDoc;
		pPattern->totdoc = iDoc;

		//set up the labels
		LABEL* pLabel = new LABEL;
		double* aClass;
		aClass = new double[vecvecIndex[i].size()];
		for (j = 0; j < vecvecIndex[i].size(); j++) {
			aClass[j] = SVMLabels[vecvecIndex[i][j]].Target;
		}

		pLabel->Class = aClass;
		pLabel->totdoc = iDoc;

		//set up the Example
		EXAMPLE* aExample;
		aExample = new EXAMPLE[1];
		aExample[0].x = *pPattern;
		aExample[0].y = *pLabel;

		//set up the Sample
		ppSample[i] = new SAMPLE;
		ppSample[i]->n = 1;
		ppSample[i]->examples = aExample;
	}
	return ppSample;
}

SAMPLE * CSVMPERF::CreateSample(Sleipnir::CDat& Dat, vector<SVMLabel> SVMLabels) {
	size_t i, j, iGene, iDoc;
	vector<DOC*> vec_pDoc;
	vector<double> vecClass;
	vector<size_t> veciGene;
	iDoc = 0;
	float numPositives, numNegatives;
	numPositives = numNegatives = 0;
	for (i = 0; i < SVMLabels.size(); i++) {
		//     cout<< "processing gene " << SVMLabels[i].GeneName << endl;
		iGene = Dat.GetGene(SVMLabels[i].GeneName);
		//   cout << SVMLabels[i].GeneName<<" gene at location "<<iGene << endl;
		if (iGene != -1) {
			//       cout << "creating doc" << endl;
			iDoc++;
			vec_pDoc.push_back(CreateDoc(Dat, iGene, iDoc - 1));
			vecClass.push_back(SVMLabels[i].Target);
		}
	}

	DOC** ppDoc;
	ppDoc = new DOC*[vec_pDoc.size()];
	copy(vec_pDoc.begin(), vec_pDoc.end(), ppDoc);
	vec_pDoc.clear();
	PATTERN* pPattern = new PATTERN;
	pPattern->doc = ppDoc;

	pPattern->totdoc = iDoc;
	//   cout << "number of document=" << pPattern->totdoc << endl;
	LABEL* pLabel = new LABEL;
	double* aClass;
	aClass = new double[vecClass.size()];
	copy(vecClass.begin(), vecClass.end(), aClass);
	vecClass.clear();
	pLabel->Class = aClass;
	pLabel->totdoc = iDoc;

	EXAMPLE* aExample;
	aExample = new EXAMPLE[1];
	//cout<<"aExample @"<<aExample<<endl;
	aExample[0].x = *pPattern;
	aExample[0].y = *pLabel;
	SAMPLE* pSample = new SAMPLE;
	pSample->n = 1;
	pSample->examples = aExample;
	/* cout << "examples @" << pSample->examples << endl;
	 cout<< "ppDoc="<<ppDoc<<endl;
	 cout << "docs @" << pSample->examples[0].x.doc << endl;
	 cout<<"done creating sample"<<endl;
	 cout<<"sample @ "<<pSample<<endl;*/
	return pSample;
}

//Single gene classification

vector<Result> CSVMPERF::Classify(Sleipnir::CPCL &PCL,
		vector<SVMLabel> SVMLabels) {
	size_t i, j, iGene, iDoc;
	vector<DOC*> vec_pDoc;
	vector<double> vecClass;
	vector<Result> vecResult;
	iDoc = 0;
	DOC** ppDoc;
	ppDoc = new DOC*[1];
	PATTERN pattern;
	pattern.doc = ppDoc;
	pattern.totdoc = 1;
	//cerr << "CLASSIFY classifying " << endl;
	LABEL label;
	for (i = 0; i < SVMLabels.size(); i++) {
		if (!SVMLabels[i].hasIndex) {
			SVMLabels[i].SetIndex(PCL.GetGene(SVMLabels[i].GeneName));
		}
		iGene = SVMLabels[i].index;
		//   cout << "CLASS gene=" << iGene << endl;
		if (iGene != -1) {
			iDoc++;

			//  cout << "CLASS iDOC=" << iDoc << endl;
			ppDoc[0] = CreateDoc(PCL, iGene, iDoc);
			label
					= classify_struct_example(pattern, &structmodel,
							&struct_parm);
			vecClass.push_back(SVMLabels[i].Target);
			vecResult.resize(iDoc);
			vecResult[iDoc - 1].GeneName = SVMLabels[i].GeneName;
			vecResult[iDoc - 1].Target = SVMLabels[i].Target;
			vecResult[iDoc - 1].Value = label.Class[0];
			//cerr<<"CLASSIFY Called FreeDoc"<<endl;
			FreeDoc(ppDoc[0]);
			//cerr<<"CLASSIFY End FreeDoc"<<endl;
		}
	}
	//    cerr << "copying" << endl;
	for (i = 0; i < label.totdoc; i++) {
		vecResult[i].Value = label.Class[i];
		//     cout << "CLASS: i=" << i << " value=" << vecResult[i].Value << endl;
	}
	//cerr << "CLASSIFY:done copying" << endl;
	//  FreePattern(pattern);
	delete ppDoc;
	return vecResult;
}

vector<Result> CSVMPERF::Classify(Sleipnir::CPCLSet &PCLSet,
		vector<SVMLabel> SVMLabels) {
	size_t i, j, iGene, iDoc;
	vector<DOC*> vec_pDoc;
	vector<double> vecClass;
	vector<Result> vecResult;
	iDoc = 0;
	for (i = 0; i < SVMLabels.size(); i++) {
		iGene = PCLSet.GetGene(SVMLabels[i].GeneName);
		//   cout << "CLASS gene=" << iGene << endl;
		if (iGene != -1) {
			iDoc++;
			//  cout << "CLASS iDOC=" << iDoc << endl;
			vec_pDoc.push_back(CreateDoc(PCLSet, iGene, iDoc));
			vecClass.push_back(SVMLabels[i].Target);
			vecResult.resize(iDoc);
			vecResult[iDoc - 1].GeneName = SVMLabels[i].GeneName;
			vecResult[iDoc - 1].Target = SVMLabels[i].Target;
		}
	}
	DOC** ppDoc;
	ppDoc = new DOC*[vec_pDoc.size()];
	copy(vec_pDoc.begin(), vec_pDoc.end(), ppDoc);
	vec_pDoc.clear();
	PATTERN pattern;
	pattern.doc = ppDoc;
	pattern.totdoc = iDoc;

	LABEL label = classify_struct_example(pattern, &structmodel, &struct_parm);
	//   cout << "label totdoc=" << label.totdoc << endl;
	for (i = 0; i < label.totdoc; i++) {
		vecResult[i].Value = label.Class[i];
		//     cout << "CLASS: i=" << i << " value=" << vecResult[i].Value << endl;
	}
	FreePattern(pattern);
	return vecResult;
}

vector<Result> CSVMPERF::Classify(Sleipnir::CDat &Dat,
		vector<SVMLabel> SVMLabels) {
	size_t i, j, iGene, iDoc;
	vector<DOC*> vec_pDoc;
	vector<double> vecClass;
	vector<Result> vecResult;
	iDoc = 0;
	for (i = 0; i < SVMLabels.size(); i++) {
		iGene = Dat.GetGene(SVMLabels[i].GeneName);
		//   cout << "CLASS gene=" << iGene << endl;
		if (iGene != -1) {
			iDoc++;
			//  cout << "CLASS iDOC=" << iDoc << endl;
			vec_pDoc.push_back(CreateDoc(Dat, iGene, iDoc));
			vecClass.push_back(SVMLabels[i].Target);
			vecResult.resize(iDoc);
			vecResult[iDoc - 1].GeneName = SVMLabels[i].GeneName;
			vecResult[iDoc - 1].Target = SVMLabels[i].Target;
		}
	}
	DOC** ppDoc;
	ppDoc = new DOC*[vec_pDoc.size()];
	copy(vec_pDoc.begin(), vec_pDoc.end(), ppDoc);
	vec_pDoc.clear();
	PATTERN pattern;
	pattern.doc = ppDoc;
	pattern.totdoc = iDoc;

	LABEL label = classify_struct_example(pattern, &structmodel, &struct_parm);
	//   cout << "label totdoc=" << label.totdoc << endl;
	for (i = 0; i < label.totdoc; i++) {
		vecResult[i].Value = label.Class[i];
		//     cout << "CLASS: i=" << i << " value=" << vecResult[i].Value << endl;
	}
	FreePattern(pattern);
	return vecResult;
}

//Create DOC for pais of genes
//For pairs of genes using PCLSet

DOC* CSVMPERF::CreateDoc(Sleipnir::CPCL &PCL, size_t iGene, size_t jGene,
		size_t iDoc) {
	WORD* aWords;
	size_t i, j, iWord, iWords, iPCL, iExp;
	float d, e;
	DOC* pRet;
	pRet->fvec->words[0].weight;
	//get number of features

	iWords = PCL.GetExperiments();

	//      cout << "CD:iwords=" << iWords << endl;
	aWords = new WORD[iWords + 1];
	//number the words
	for (i = 0; i < iWords; ++i) {
		//   cout<<i<<endl;
		aWords[i].wnum = i + 1;
		// asWords[ i ].wnum = 0;
	}
	aWords[i].wnum = 0;
	//get the values;
	iWord = 0;

	for (j = 0; j < PCL.GetExperiments(); j++) {
		//     cout<<"CD:iWord="<<iWord<<endl;
		if (!Sleipnir::CMeta::IsNaN(d = PCL.Get(iGene, j))
				&& !Sleipnir::CMeta::IsNaN(e = PCL.Get(jGene, j))) {
			//   if (i==0 && j==0)
			//       cout<<"First value is "<<d<<endl;
			aWords[iWord].weight = d * e;
		} else
			aWords[iWord].weight = 0;
		iWord++;
	}

	pRet = create_example(iDoc, 0, 0, 1, create_svector(aWords, "", 1));
	delete[] aWords;
	// cout<<"done creating DOC"<<endl;
	return pRet;
}

//Create Sample for pairs of genes, wraps the above 2 functions

SAMPLE* CSVMPERF::CreateSample(Sleipnir::CPCL &PCL, Sleipnir::CDat& Answers,
		const vector<string>& CVGenes) {
	size_t i, j, iGene, jGene, iDoc;
	vector<DOC*> vec_pDoc;
	vector<double> vecClass;
	vector<size_t> veciGene;
	Sleipnir::CPCL* pP = &PCL;
	iDoc = 0;
	float numPositives, numNegatives;
	numPositives = numNegatives = 0;
	set<string> setGenes;
	set<string>::iterator iterSet;
	float w;
	for (i = 0; i < CVGenes.size(); i++) {
			setGenes.insert(CVGenes[i]);
	}
	for (i = 0; i < Answers.GetGenes() - 1; i++) {
		if ((setGenes.find(Answers.GetGene(i)) != setGenes.end()) && ((iGene
				= PCL.GetGene(Answers.GetGene(i))) != -1)) {
			for (j = i + 1; j < Answers.GetGenes(); j++) {
				if ((setGenes.find(Answers.GetGene(j)) != setGenes.end())
						&& ((jGene = PCL.GetGene(Answers.GetGene(j))) != -1)) {
					if (!Sleipnir::CMeta::IsNaN(w = Answers.Get(i, j))) {
						iDoc++;
						vec_pDoc.push_back(CreateDoc(PCL, iGene, jGene, iDoc
								- 1));

					}
				}
			}
		}
	}

	DOC** ppDoc;
	ppDoc = new DOC*[vec_pDoc.size()];
	copy(vec_pDoc.begin(), vec_pDoc.end(), ppDoc);
	vec_pDoc.clear();
	PATTERN* pPattern = new PATTERN;
	pPattern->doc = ppDoc;

	pPattern->totdoc = iDoc;
	LABEL* pLabel = new LABEL;
	double* aClass;
	aClass = new double[vecClass.size()];
	copy(vecClass.begin(), vecClass.end(), aClass);
	vecClass.clear();
	pLabel->Class = aClass;
	pLabel->totdoc = iDoc;

	EXAMPLE* aExample;
	aExample = new EXAMPLE[1];
	aExample[0].x = *pPattern;
	aExample[0].y = *pLabel;
	SAMPLE* pSample = new SAMPLE;
	pSample->n = 1;
	pSample->examples = aExample;
	return pSample;
}

void CSVMPERF::Classify(Sleipnir::CPCL& PCL, Sleipnir::CDat& Answers,
		Sleipnir::CDat& Values, Sleipnir::CDat& Counts,
		const vector<string>& CVGenes) {
	size_t i, j, iGene, jGene, iDoc;
	set<string> setGenes;
	set<string>::iterator iterSet;
	vector<DOC*> vec_pDoc;
	vector<pair<size_t, size_t> > vecPairIndex;
	iDoc = 0;
	for (i = 0; i < CVGenes.size(); i++) {
		setGenes.insert(CVGenes[i]);
	}

	cout << "the number of genes  to be classified is " << setGenes.size()
			<< endl;
	for (i = 0; i < Answers.GetGenes() - 1; i++) {
		if ((setGenes.find(Answers.GetGene(i)) != setGenes.end()) && ((iGene
				= PCL.GetGene(Answers.GetGene(i))) != -1)) {
			for (j = i + 1; j < Answers.GetGenes(); j++) {
				if ((setGenes.find(Answers.GetGene(j)) != setGenes.end())
						&& ((jGene = PCL.GetGene(Answers.GetGene(j))) != -1)) {
					if (!Sleipnir::CMeta::IsNaN(Answers.Get(i, j))) {
						iDoc++;
						vec_pDoc.push_back(CreateDoc(PCL, iGene, jGene, iDoc
								- 1));
						vecPairIndex.push_back(pair<size_t, size_t> (i, j));
					}
				}
			}
		}
	}
	DOC** ppDoc;
	ppDoc = new DOC*[vec_pDoc.size()];
	copy(vec_pDoc.begin(), vec_pDoc.end(), ppDoc);
	vec_pDoc.clear();
	PATTERN pattern;
	pattern.doc = ppDoc;
	pattern.totdoc = iDoc;

	LABEL label = classify_struct_example(pattern, &structmodel, &struct_parm);
	pair<size_t, size_t> tmpPair;
	float d;
	for (i = 0; i < vecPairIndex.size(); i++) {
		tmpPair = vecPairIndex[i];
		if (!Sleipnir::CMeta::IsNaN(d = Values.Get(tmpPair.first,
				tmpPair.second)))
			Values.Set(tmpPair.first, tmpPair.second, d + label.Class[i]);
		else
			Values.Set(tmpPair.first, tmpPair.second, label.Class[i]);
		if (!Sleipnir::CMeta::IsNaN(d = Counts.Get(tmpPair.first,
				tmpPair.second)))
			Counts.Set(tmpPair.first, tmpPair.second, d + 1);
		else
			Counts.Set(tmpPair.first, tmpPair.second, 1);
	}
		FreePattern(pattern);

}

void CSVMPERF::ClassifyAll(Sleipnir::CPCL& PCL, Sleipnir::CDat& Values,
		Sleipnir::CDat& Counts, const vector<string>& CVGenes) {
	size_t iGene, jGene, i, j, k;
	string strGeneOne, strGeneTwo;
	set<string> setGenes;
	set<string>::iterator iterSet;
	vector<DOC*> vec_pDoc;
	vector<pair<size_t, size_t> > vecPairIndex;
	size_t iDoc;
	for (i = 0; i < CVGenes.size(); i++) {
		setGenes.insert(CVGenes[i]);
	}

	for (i = 0; i < PCL.GetGenes() - 1; i++) {
		if (setGenes.find(PCL.GetGene(i)) == setGenes.end()) {
			iDoc = 0;
			for (j = i + 1; j < PCL.GetGenes() - 1; j++) {
				if (setGenes.find(PCL.GetGene(j)) == setGenes.end()) {
					iDoc++;
					vec_pDoc.push_back(CreateDoc(PCL, i, j, iDoc - 1));
					vecPairIndex.push_back(pair<size_t, size_t> (i, j));
				}
			}
			DOC** ppDoc;
			ppDoc = new DOC*[vec_pDoc.size()];
			copy(vec_pDoc.begin(), vec_pDoc.end(), ppDoc);
			vec_pDoc.clear();
			PATTERN pattern;
			pattern.doc = ppDoc;
			pattern.totdoc = iDoc;

			LABEL label = classify_struct_example(pattern, &structmodel,
					&struct_parm);
			pair<size_t, size_t> tmpPair;
			float d;
			for (k = 0; k < vecPairIndex.size(); k++) {
				tmpPair = vecPairIndex[k];
				if (!Sleipnir::CMeta::IsNaN(d = Values.Get(tmpPair.first,
						tmpPair.second)))
					Values.Set(tmpPair.first, tmpPair.second, d
							+ label.Class[k]);
				else
					Values.Set(tmpPair.first, tmpPair.second, label.Class[k]);
				if (!Sleipnir::CMeta::IsNaN(d = Counts.Get(tmpPair.first,
						tmpPair.second)))
					Counts.Set(tmpPair.first, tmpPair.second, d + 1);
				else
					Counts.Set(tmpPair.first, tmpPair.second, 1);
			}
			vecPairIndex.resize(0);
			FreePattern(pattern);
		}
	}
}

void CSVMPERF::ClassifyAll(Sleipnir::CPCL& PCL, Sleipnir::CDat& Values,
		const vector<string>& CVGenes) {
	size_t iGene, jGene, i, j, k;
	string strGeneOne, strGeneTwo;
	set<string> setGenes;
	set<string>::iterator iterSet;
	vector<DOC*> vec_pDoc;
	vector<pair<size_t, size_t> > vecPairIndex;
	size_t iDoc;
	for (i = 0; i < CVGenes.size(); i++) {
		setGenes.insert(CVGenes[i]);
	}

	for (i = 0; i < PCL.GetGenes() - 1; i++) {
		if (setGenes.find(PCL.GetGene(i)) == setGenes.end()) {
			iDoc = 0;
			for (j = i + 1; j < PCL.GetGenes() - 1; j++) {
				if (setGenes.find(PCL.GetGene(j)) == setGenes.end()) {
					iDoc++;
					vec_pDoc.push_back(CreateDoc(PCL, i, j, iDoc - 1));
					vecPairIndex.push_back(pair<size_t, size_t> (i, j));
				}
			}
			DOC** ppDoc;
			ppDoc = new DOC*[vec_pDoc.size()];
			copy(vec_pDoc.begin(), vec_pDoc.end(), ppDoc);
			vec_pDoc.clear();
			PATTERN pattern;
			pattern.doc = ppDoc;
			pattern.totdoc = iDoc;

			LABEL label = classify_struct_example(pattern, &structmodel,
					&struct_parm);
			pair<size_t, size_t> tmpPair;
			float d;
			for (k = 0; k < vecPairIndex.size(); k++) {
				tmpPair = vecPairIndex[k];
				Values.Set(tmpPair.first, tmpPair.second, d + label.Class[k]);
			}
			vecPairIndex.resize(0);
			FreePattern(pattern);
		}
	}
}

}

