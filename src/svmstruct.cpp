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
#include "svmstruct.h"
#include "pclset.h"
#include "dataset.h"
#include "meta.h"
#include "genome.h"
#include "compactmatrix.h"
#include <vector>
#include <set>

#define  SLACK_RESCALING    1
#define  MARGIN_RESCALING   2

namespace SVMArc {
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

	void CSVMSTRUCTMC::SetVerbosity(size_t V) {
		struct_verbosity = (long) V;
	}

	bool CSVMSTRUCTMC::initialize() {

		//set directionality


		/* set default */
		Alg = DEFAULT_ALG_TYPE;
		//Learn_parms
		struct_parm.C=0.01;
		struct_parm.slack_norm=1;
		struct_parm.epsilon=DEFAULT_EPS;
		struct_parm.custom_argc=0;
		struct_parm.loss_function=DEFAULT_LOSS_FCT;
		struct_parm.loss_type=DEFAULT_RESCALING;
		struct_parm.newconstretrain=100;
		struct_parm.ccache_size=5;
		struct_parm.batch_size=100;
		//Learn_parms
		//strcpy (learn_parm.predfile, "trans_predictions");
		strcpy(learn_parm.alphafile, "");
		//verbosity=0;/*verbosity for svm_light*/
		//struct_verbosity = 1; /*verbosity for struct learning portion*/
		learn_parm.biased_hyperplane=1;
		learn_parm.remove_inconsistent=0;
		learn_parm.skip_final_opt_check=0;
		learn_parm.svm_maxqpsize=10;
		learn_parm.svm_newvarsinqp=0;
		learn_parm.svm_iter_to_shrink=-9999;
		learn_parm.maxiter=100000;
		learn_parm.kernel_cache_size=40;
		learn_parm.svm_c=99999999;  /* overridden by struct_parm.C */
		learn_parm.eps=0.001;       /* overridden by struct_parm.epsilon */
		learn_parm.transduction_posratio=-1.0;
		learn_parm.svm_costratio=1.0;
		learn_parm.svm_costratio_unlab=1.0;
		learn_parm.svm_unlabbound=1E-5;
		learn_parm.epsilon_crit=0.001;
		learn_parm.epsilon_a=1E-10;  /* changed from 1e-15 */
		learn_parm.compute_loo=0;
		learn_parm.rho=1.0;
		learn_parm.xa_depth=0;
		kernel_parm.kernel_type=0;
		kernel_parm.poly_degree=3;
		kernel_parm.rbf_gamma=1.0;
		kernel_parm.coef_lin=1;
		kernel_parm.coef_const=1;
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

		/* set number of features to -1, indicating that it will be computed
		in init_struct_model() */
		struct_parm.num_features = -1;

		return true;
	}

	bool CSVMSTRUCTMC::parms_check() {
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



		return true;
	}



	DOC* CSVMSTRUCTMC::CreateDoc(Sleipnir::CPCL &PCL, size_t iGene, size_t iDoc) {
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


	vector<SVMLabel> CSVMSTRUCTMC::ReadLabels(ifstream & ifsm) {

		static const size_t c_iBuffer = 1024;
		char acBuffer[c_iBuffer];
		vector<string> vecstrTokens;
		vector<SVMLabel> vecLabels;
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
			vecLabels.push_back(SVMArc::SVMLabel(vecstrTokens[0], atoi(
				vecstrTokens[1].c_str())));
			if (vecLabels.back().Target > 0)
				numPositives++;
			else
				numNegatives++;
		}
		return vecLabels;
	}


	SAMPLE* CSVMSTRUCTMC::CreateSample(Sleipnir::CPCL &PCL, vector<SVMLabel> SVMLabels) {
		size_t i, j, iGene, iDoc;
		int     n;       /* number of examples */
		int *target;
		long num_classes=0;
		SAMPLE* pSample = new SAMPLE;
		EXAMPLE* examples;
		DOC** docs;
		vector<DOC*> vec_pDoc;
		vec_pDoc.reserve(SVMLabels.size());
		vector<int> vecClass;
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

		//copy patterns and labels to new vector
		docs = new DOC*[vec_pDoc.size()];
		n = vec_pDoc.size();
		//cout << "Read in " << n << "Standards"<<endl;
		copy(vec_pDoc.begin(), vec_pDoc.end(), docs);
		vec_pDoc.clear();

		//cerr << "NEW Class array" << endl;
		target = new int[vecClass.size()];
		copy(vecClass.begin(), vecClass.end(), target);
		vecClass.clear();



		examples=(EXAMPLE *)my_malloc(sizeof(EXAMPLE)*n);
		for(i=0;i<n;i++)     /* find highest class label */
			if(num_classes < target[i]) 
				num_classes=target[i];

		for(i=0;i<n;i++)     /* make sure all class labels are positive */
			if(target[i]<1) {
				printf("\nERROR: The class label '%d' of example number %ld is not greater than '1'!\n",target[i],i+1);
				exit(1);
			} 
			for(i=0;i<n;i++) {          /* copy docs over into new datastructure */
				examples[i].x.doc=docs[i];
				examples[i].y.Class=target[i];
				examples[i].y.scores=NULL;
				examples[i].y.num_classes=num_classes;
			}
			free(target);
			free(docs);
			pSample->n=n;
			pSample->examples=examples;

			if(struct_verbosity>=0)
				printf(" (%d examples) ",pSample->n);




			return pSample;
			//cerr<<"DONE CreateSample"<<endl;
	}

	//Single gene classification

	vector<Result> CSVMSTRUCTMC::Classify(Sleipnir::CPCL &PCL,
		vector<SVMLabel> SVMLabels) {
			size_t i, j,k, iGene, iDoc;
			vector<int> vecClass;
			vector<Result> vecResult;
			iDoc = 0;
			PATTERN pattern;
			pattern.totdoc = 1;
			cerr << "CLASSIFY classifying " << endl;
			LABEL label;
			for (i = 0; i < SVMLabels.size(); i++) {
				if (!SVMLabels[i].hasIndex) {
					SVMLabels[i].SetIndex(PCL.GetGene(SVMLabels[i].GeneName));
				}
				iGene = SVMLabels[i].index;
				//cout << "CLASS gene=" << iGene << endl;
				if (iGene != -1) {
					iDoc++;

					//cout << "CLASS iDOC=" << iDoc << endl;
					pattern.doc = CreateDoc(PCL, iGene, iDoc);
					//cerr<<"Doc Created"<<endl;
					label	= classify_struct_example(pattern, &structmodel,
						&struct_parm);
					//cerr<<"CLASSIED"<<endl;
					vecClass.push_back(SVMLabels[i].Target);
					vecResult.resize(iDoc);
					vecResult[iDoc - 1].GeneName = SVMLabels[i].GeneName;
					vecResult[iDoc - 1].Target = SVMLabels[i].Target;
					vecResult[iDoc - 1].Value = label.Class;
					vecResult[iDoc - 1].num_class=struct_parm.num_classes;
					//vecResult[iDoc - 1].Scores.reserve(label.num_classes);
					for (k = 1; k <= struct_parm.num_classes; k++)
						vecResult[iDoc - 1].Scores.push_back(label.scores[k]);
					//cerr<<"CLASSIFY Called FreeDoc"<<endl;
					FreeDoc(pattern.doc);
					//cerr<<"CLASSIFY End FreeDoc"<<endl;
				}
			}

			return vecResult;
	}


	void CSVMSTRUCTMC::FreeSample_leave_Doc(SAMPLE s){
		/* Frees the memory of sample s. */
		int i;
		for(i=0;i<s.n;i++) {
			free(s.examples[i].x.doc);
			free_label(s.examples[i].y);
		}
		free(s.examples);
	}



}

