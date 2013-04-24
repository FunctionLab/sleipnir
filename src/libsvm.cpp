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
#include "libsvm.h"
#include "pclset.h"
#include "dataset.h"
#include "meta.h"
#include "genome.h"
#include "compactmatrix.h"
#include <vector>
#include <set>

#define Malloc(type,n) (type *)malloc((n)*sizeof(type)) 

//#include <libsvm.h>

namespace LIBSVM {

#include "libsvm.h"

bool CLIBSVM::initialize() {

	/* set default */
    //cout << "sadfsajdfsalkdjfalskdjfkalsjdf" << endl;

        parm.cache_size = 100;
        parm.C = 0.01;
	parm.eps = 1e-3;
        parm.svm_type = C_SVC;
        parm.p = 0.1;
        parm.shrinking = 1;
        parm.nr_weight = 0;
        parm.weight_label = NULL;
        parm.weight = NULL;
        parm.probability = 0;
        parm.nu = 0.5;
        parm.coef0 = 0;
        parm.gamma = 0;
        parm.degree = 3;
        parm.kernel_type = LINEAR;

	return true;
}

bool CLIBSVM::parms_check() {
	if (parm.C < 0) {
		fprintf(
				stderr,
				"\nTrade-off between training error and margin is not set (C<0)!\nC value will be set to default value. Clight = Cpef * 100 / n \n");
		fprintf(stderr, "be less than 1.0 !!!\n\n");
		return false;
	}
	if (parm.eps <= 0) {
		fprintf(stderr,
				"\nThe epsilon parameter must be greater than zero!\n\n");
		return false;
	}

        //TODO: add more parameter checks 

	return true;
}

SAMPLE * CLIBSVM::CreateSample(Sleipnir::CPCL& PCL, vector<SVMLabel> SVMLabels) {
	size_t i, j, k, iGene, iProblem, numFeatures, max_index;
        float d;

        struct svm_problem* prob;
        struct svm_node* x_space;

        prob = Malloc(struct svm_problem,1);

        prob->l = 0;//number of labels in PCL
        numFeatures = PCL.GetExperiments();

        cout << "number of features: " << numFeatures << endl;
//        return NULL;
//	vector<double> vecClass;
//	vector<size_t> veciGene;
	
        iProblem = 0;

	for (i = 0; i < SVMLabels.size(); i++) {
//		     cout<< "processing gene " << SVMLabels[i].GeneName << endl;
                if (!SVMLabels[i].hasIndex){
//                  cout << "setindex: " << PCL.GetGene(SVMLabels[i].GeneName) << endl;
                  SVMLabels[i].SetIndex(PCL.GetGene(SVMLabels[i].GeneName));
                }
		iGene = SVMLabels[i].index;
//		   cout << SVMLabels[i].GeneName<<" gene at location "<<iGene << endl;
		if (iGene != -1) {
//cout << iGene << endl;
                  (prob->l)++;
		}
	}

 
cout << "number of labels in data: " << prob->l << endl;
cout << "number of elements: " << (1+numFeatures) * prob->l << endl; 

        prob->y = Malloc(double,prob->l);
        prob->x = Malloc(struct svm_node *, prob->l);
        x_space = Malloc(struct svm_node, (1+numFeatures) * prob->l);

        max_index = numFeatures;

        j = 0;//element index

        for (i = 0; i < SVMLabels.size(); i++) {
            iGene = SVMLabels[i].index;
            if (iGene != -1){
              (prob->x)[i] = &x_space[j];
              (prob->y)[i] = SVMLabels[i].Target;

//cout << SVMLabels[i].Target << endl;
//cout << prob.y[i] << endl;

              for(k = 0; k < numFeatures; k++){
                x_space[j].index = k;
                if (!Sleipnir::CMeta::IsNaN(d = PCL.Get(iGene, k))) {
                  x_space[j].value = d;
                }else{
                  x_space[j].value = 1; //TODO: make this a flag!!!
                  //if missing value??? SVMPerf imputes 0 ... what about gene i to gene i ? should impute 1?
                }
//cout << k << endl;
//cout << "index: " << x_space[j].index << endl;
//cout << "value: " << x_space[j].value << endl;
//
//double blah = x_space[j].index + x_space[j].value;
                j++;
              }         
              x_space[j].index = -1;
// cout << k << endl;
//cout << "index: " << x_space[j].index << endl;
//cout << "value: " << x_space[j].value << endl;
              
//              x_space[j].index = -1;       
              j++;
            }
        }
//cout << "here" << endl;
        SAMPLE* pSample;
        pSample = Malloc(SAMPLE,1);

        pSample->n = prob->l;
        pSample->problems = prob;
//sleep(15);
//

//        cout << ((pSample->problems)->x)
        cout << ((pSample->problems)->y)[0] << endl;
        cout << ((pSample->problems)->y)[1] << endl;


//cout << "an svm_node index: " << (((*prob).x)[0][0]).index << endl;
//cout << "an svm_node value: " << (((*prob).x)[0][0]).value << endl;
//sleep(15);
	return pSample;
}



SAMPLE * CLIBSVM::CreateSample(Sleipnir::CDat& Dat, vector<SVMLabel> SVMLabels) {
	size_t i, j, k, iGene, iProblem, numFeatures, max_index;
        float d;

        struct svm_problem* prob;
        struct svm_node *x_space;

        prob->l = 0;//number of labels in Dat
        numFeatures = Dat.GetGenes();

//	vector<double> vecClass;
//	vector<size_t> veciGene;
	
        iProblem = 0;

	for (i = 0; i < SVMLabels.size(); i++) {
		//     cout<< "processing gene " << SVMLabels[i].GeneName << endl;
		iGene = Dat.GetGene(SVMLabels[i].GeneName);
		//   cout << SVMLabels[i].GeneName<<" gene at location "<<iGene << endl;
		if (iGene != -1) {
                    
			//       cout << "creating doc" << endl;
//			iProblem++;
//			vec_pproblem.push_back(CreateDoc(Dat, iGene, iDoc - 1));
//			vecClass.push_back(SVMLabels[i].Target);

                  (prob->l)++;
		}
	}

 
        prob->y = Malloc(double,prob->l);
        prob->x = Malloc(struct svm_node *, prob->l);
        x_space = Malloc(struct svm_node, numFeatures * (prob->l));

        max_index = numFeatures;
        j = 0;

        for (i = 0; i < SVMLabels.size(); i++) {
            iGene = Dat.GetGene(SVMLabels[i].GeneName);
            if (iGene != -1){
              (prob->x)[i] = &x_space[j];
              (prob->y)[i] = SVMLabels[i].Target;
              for(k = 0; k < numFeatures; k++){
                x_space[j].index = k;
                if (!Sleipnir::CMeta::IsNaN(d = Dat.Get(iGene, k))) {
                  x_space[j].value = d;
                }else{
                  x_space[j].value = 1; //TODO: make this a flag!!!
                  //if missing value??? SVMPerf imputes 0 ... what about gene i to gene i ? should impute 1?
                }
              }
              j++;
            }
        }

        SAMPLE* pSample;
        pSample = Malloc(SAMPLE,1);
        
        pSample->n = prob->l;
        pSample->problems = prob;



	return pSample;
}


vector<Result> CLIBSVM::Classify(Sleipnir::CPCL &PCL,
        vector<SVMLabel> SVMLabels) {
    size_t i, j, iGene;
    double predict_label;
    double* dec_values;
    double dec_value;
    struct svm_node *x;

    SAMPLE* pSample;
    vector<Result> vecResult;

    pSample = CLIBSVM::CreateSample(PCL, SVMLabels);

    int svm_type = svm_get_svm_type(model);
    int nr_class = svm_get_nr_class(model);

cerr << nr_class << endl;
cerr << pSample->n << endl;;

    dec_values = Malloc(double, nr_class*(nr_class-1)/2);
    vecResult.resize(pSample->n);

cerr << "number of samples: " << vecResult.size() << endl;
cerr << "length of svm labels: " << SVMLabels.size() << endl;

    j= 0; //pSample index
    for(i = 0 ; i < SVMLabels.size() ; i++){
      if(!SVMLabels[i].hasIndex){//assume createSample sequentially added to pSample TODO: currently hacky
        SVMLabels[i].SetIndex(PCL.GetGene(SVMLabels[i].GeneName));
      }
      iGene = SVMLabels[i].index;

      if (iGene != -1 ) {    
        x = pSample->problems->x[j];
        predict_label = svm_predict_values(model,x, dec_values);
        dec_value = dec_values[0]; //assume that positive class is the first class TODO: currently hackly
//cerr << "Gene: " << iGene << endl;
//cerr << dec_value << endl;        


cerr << j << endl;

        vecResult[j].GeneName = SVMLabels[i].GeneName;
        vecResult[j].Target = SVMLabels[i].Target;
        vecResult[j].Value = dec_value;
        
        j++;

      }
    }

//    delete dec_values;

    return vecResult;
}

    
//TODO: classify for dab/dat files
//

}
