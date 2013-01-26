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

			//cout << "CLASS iDOC=" << iDoc << endl;
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


// populate docs for each label gene pair from a vector of dat file names
bool CSVMPERF::CreateDoc(vector<string>& vecstrDatasets,
			 vector<SVMLabelPair*>& vecLabels,
			 const vector<string>& LabelsGene,
			 Sleipnir::CDat::ENormalize eNormalize){
  
  size_t i, j, k, iGene, jGene, iDoc, iWord, iWords;
  float d;
  vector<size_t> veciGene;  
  vector<WORD*> vec_pWord;
  vector<size_t> labelg2dat;  
  
  WORD* aWords;
  DOC* pRet;	
  //Sleipnir::CDat Dat;
  
  iWords = vecstrDatasets.size();
  
  vec_pWord.reserve(vecLabels.size());
  
  // initialize all the WORDs
  for(i=0; i < vecLabels.size(); i++){
    aWords = new WORD[iWords + 1];
    
    // set wnum values
    for (k = 0; k < iWords; ++k) {
      //   cout<<i<<endl;
      aWords[k].wnum = k + 1;
      // asWords[ i ].wnum = 0;
    }
    aWords[k].wnum = 0;    
    vec_pWord[i] = aWords;
  }
  
  // initialize the gene mappings 
  labelg2dat.resize(LabelsGene.size());
  
  // now open up all datasets
  for(i=0; i < vecstrDatasets.size(); i++){
    Sleipnir::CDat Dat;
    if(!Dat.Open(vecstrDatasets[i].c_str())) {
      cerr << vecstrDatasets[i].c_str() << endl;
      cerr << "Could not open: " << vecstrDatasets[i] << endl;
      return false;
    }
    
    // normalize dat file
    if( eNormalize != Sleipnir::CDat::ENormalizeNone ){
      cerr << "Normalize input data" << endl;      
      Dat.Normalize( eNormalize );
    }
    
    cerr << "Open: " << vecstrDatasets[i] << endl;
    
    // construct gene name mapping
    for(k=0; k < LabelsGene.size(); k++){
      labelg2dat[k] = Dat.GetGene(LabelsGene[k]);
    }    
    /////
    
    for (j = 0; j < vecLabels.size(); j++) {
      aWords = vec_pWord[j];
            
      if( ((iGene = labelg2dat[vecLabels[j]->iidx]) == -1 ) || 
	  ((jGene = labelg2dat[vecLabels[j]->jidx]) == -1 )){
	aWords[i].weight = 0;
	continue;
      }      
      
      if (!Sleipnir::CMeta::IsNaN(d = Dat.Get(iGene, jGene))) {
	aWords[i].weight = d;
      } else
	aWords[i].weight = 0;            
    }
  }
  
  // now create a Doc per label
  for (j = 0; j < vecLabels.size(); j++) {
    //pRet->fvec->words[0].weight;
    aWords = vec_pWord[j];
    pRet = create_example(j, 0, 0, 1, create_svector(aWords, "", 1));
    vecLabels[j]->pDoc = pRet;
    delete[] aWords;
  }
  
  return true;
}

SAMPLE* CSVMPERF::CreateSample(vector<SVMLabelPair*>& SVMLabels) {
	size_t i, j, iGene, iDoc;
	vector<DOC*> vec_pDoc;
	vector<double> vecClass;
	iDoc = 0;
	for (i = 0; i < SVMLabels.size(); i++) {
	  iDoc++;
	  vec_pDoc.push_back(SVMLabels[i]->pDoc);
	  vecClass.push_back(SVMLabels[i]->Target);
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

void CSVMPERF::Classify(Sleipnir::CDat &Results,
				  vector<SVMLabelPair*>& SVMLabels) {
	size_t i, iGene, iDoc;
	iDoc = 0;
	DOC** ppDoc;
	ppDoc = new DOC*[1];
	PATTERN pattern;
	pattern.doc = ppDoc;
	pattern.totdoc = 1;
	//cerr << "CLASSIFY classifying " << endl;
	LABEL label;
	for (i = 0; i < SVMLabels.size(); i++) {
	  ppDoc[0] = SVMLabels[i]->pDoc;
	  label
	    = classify_struct_example(pattern, &structmodel,
				      &struct_parm);
	  
	  Results.Set(SVMLabels[i]->iidx, SVMLabels[i]->jidx, label.Class[0]);
	}
	
	
	delete ppDoc;
}

void CSVMPERF::FreeSample_leave_Doc(SAMPLE s){
  /* Frees the memory of sample s. */
  int i;
  for(i=0;i<s.n;i++) {
    free(s.examples[i].x.doc);
    free_label(s.examples[i].y);
  }
  free(s.examples);
}

// Platt's binary SVM Probablistic Output
// Assume dec_values and labels have same dimensions and genes
void CSVMPERF::sigmoid_train(Sleipnir::CDat& Results,
			     vector<SVMLabelPair*>& SVMLabels,
			     float& A, float& B){
	double prior1=0, prior0 = 0;
	size_t i, j, idx;
	float d, lab;
	
	int max_iter=100;	// Maximal number of iterations
	double min_step=1e-10;	// Minimal step taken in line search
	double sigma=1e-12;	// For numerically strict PD of Hessian
	double eps=1e-5;
	vector<double> t;
	double fApB,p,q,h11,h22,h21,g1,g2,det,dA,dB,gd,stepsize;
	double newA,newB,newf,d1,d2;
	int iter; 
	vector<float> dec_values;
	
	// Negatives are values less than 0
	for(i = 0; i < SVMLabels.size(); i++)
	  if( SVMLabels[i]->Target > 0 )
	    prior1 += 1;
	  else if(SVMLabels[i]->Target < 0 )
	    prior0 += 1;

	// initialize size
	t.resize(prior0+prior1);
	
	
	// classify training examples
	LABEL label;
	size_t iGene, iDoc;
	iDoc = 0;
	DOC** ppDoc;
	ppDoc = new DOC*[1];
	PATTERN pattern;
	pattern.doc = ppDoc;
	pattern.totdoc = 1;
	
	dec_values.resize(t.size());
	
	idx = 0;
	for (i = 0; i < SVMLabels.size(); i++) {
	  if( SVMLabels[i]->Target == 0 )
	    continue;
	  if( Sleipnir::CMeta::IsNaN(d = Results.Get(SVMLabels[i]->iidx, SVMLabels[i]->jidx)))
	    continue;
	  
	  dec_values[idx] = d;
	  idx++;
	}
	
	// Initial Point and Initial Fun Value
	A=0.0; B=log((prior0+1.0)/(prior1+1.0));
	double hiTarget=(prior1+1.0)/(prior1+2.0);
	double loTarget=1/(prior0+2.0);			
	double fval = 0.0;
	
	idx = 0;
	for(i = 0; i < SVMLabels.size(); i++){
	  if ( SVMLabels[i]->Target > 0) t[idx]=hiTarget;
	  else if( SVMLabels[i]->Target < 0 ) t[idx]=loTarget;
	  else 
	    continue;
	  
	  fApB = dec_values[idx]*A+B;
	  if (fApB>=0)
	    fval += t[idx]*fApB + log(1+exp(-fApB));
	  else
	    fval += (t[idx] - 1)*fApB +log(1+exp(fApB));
	  
	  idx++;
	}
	
	for (iter=0;iter<max_iter;iter++)
	{
		// Update Gradient and Hessian (use H' = H + sigma I)
		h11=sigma; // numerically ensures strict PD
		h22=sigma;
		h21=0.0;g1=0.0;g2=0.0;
		for (i=0; i < dec_values.size(); i++)
		{
			fApB = dec_values[i]*A+B;
			if (fApB >= 0)
			{
				p=exp(-fApB)/(1.0+exp(-fApB));
				q=1.0/(1.0+exp(-fApB));
			}
			else
			{
				p=1.0/(1.0+exp(fApB));
				q=exp(fApB)/(1.0+exp(fApB));
			}
			d2=p*q;
			h11+=dec_values[i]*dec_values[i]*d2;
			h22+=d2;
			h21+=dec_values[i]*d2;
			d1=t[i]-p;
			g1+=dec_values[i]*d1;
			g2+=d1;
		}

		// Stopping Criteria
		if (fabs(g1)<eps && fabs(g2)<eps)
			break;

		// Finding Newton direction: -inv(H') * g
		det=h11*h22-h21*h21;
		dA=-(h22*g1 - h21 * g2) / det;
		dB=-(-h21*g1+ h11 * g2) / det;
		gd=g1*dA+g2*dB;


		stepsize = 1;		// Line Search
		while (stepsize >= min_step)
		{
			newA = A + stepsize * dA;
			newB = B + stepsize * dB;

			// New function value
			newf = 0.0;
			for (i=0; i < dec_values.size(); i++)
			{
				fApB = dec_values[i]*newA+newB;
				if (fApB >= 0)
					newf += t[i]*fApB + log(1+exp(-fApB));
				else
					newf += (t[i] - 1)*fApB +log(1+exp(fApB));
			}
			// Check sufficient decrease
			if (newf<fval+0.0001*stepsize*gd)
			{
				A=newA;B=newB;fval=newf;
				break;
			}
			else
				stepsize = stepsize / 2.0;
		}

		if (stepsize < min_step)
		{
		  cerr << "Line search fails in two-class probability estimates\n";
		  break;
		}
	}
	
	if (iter>=max_iter)
	  cerr << "Reaching maximal iterations in two-class probability estimates\n";
}

void CSVMPERF::sigmoid_predict(Sleipnir::CDat& Results, vector<SVMLabelPair*>& SVMLabels, float A, float B){
  size_t i;
  float d, fApB;
  
  for(i = 0; i < SVMLabels.size(); i++){
    if (!Sleipnir::CMeta::IsNaN(d = Results.Get(SVMLabels[i]->iidx, SVMLabels[i]->jidx))){
	fApB = d*A+B;
	// 1-p used later; avoid catastrophic cancellation
	if (fApB >= 0)
	  Results.Set(SVMLabels[i]->iidx, SVMLabels[i]->jidx, exp(-fApB)/(1.0+exp(-fApB)));
	else
	  Results.Set(SVMLabels[i]->iidx, SVMLabels[i]->jidx, 1.0/(1+exp(fApB)));      
    }    
  }
}
  
}

