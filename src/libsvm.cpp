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
bool CLIBSVM::posFeatOnly = false;

bool CLIBSVM::initialize() {

	/* set default */

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

        if (parm.nu < 0 | parm.nu > 1) {
            fprintf(stderr, "nu parameter must be between 0 and 1");
            return false;
        }

        //TODO: add more parameter checks 

	return true;
}

SAMPLE * CLIBSVM::CreateSample(Sleipnir::CPCL& PCL, vector<SVMLabel> SVMLabels) {
	size_t i, j, k, s, iGene, iProblem, numFeatures, numLabels, max_index;
        float d;
        bool posFeatOnly;

        struct svm_problem* prob;
        struct svm_node* x_space;
        vector<size_t> iPosFeats;

        prob = Malloc(struct svm_problem,1);

        numFeatures = PCL.GetExperiments();
        numLabels = 0;
        
        posFeatOnly = CLIBSVM::posFeatOnly; 
// cerr << "in create sample: " << posFeatOnly << endl;       

	
        iProblem = 0;

	for (i = 0; i < SVMLabels.size(); i++) {
                if (!SVMLabels[i].hasIndex){
                  SVMLabels[i].SetIndex(PCL.GetGene(SVMLabels[i].GeneName));
                }
		iGene = SVMLabels[i].index;
		if (iGene != -1) {
                  numLabels++;
                  if(posFeatOnly){
                    if(SVMLabels[i].Target > 0){
                      iPosFeats.push_back(iGene);
                    }
                  }
		}
	}

        if(posFeatOnly){
          numFeatures = iPosFeats.size();
        }

cerr << "number of features used: " << numFeatures << endl;
cerr << "number of labels given: " << SVMLabels.size() << endl;
cerr << "number of labels in data: " << numLabels << endl;

        prob->l = numLabels;
        prob->y = Malloc(double,numLabels);
        prob->x = Malloc(struct svm_node *, numLabels);
        x_space = Malloc(struct svm_node, (1+numFeatures) * numLabels);

        max_index = numFeatures;

        j = 0;//element index
        s = 0;//sample index
        for (i = 0; i < SVMLabels.size(); i++) {
            iGene = SVMLabels[i].index;

            if (iGene != -1){
              (prob->x)[s] = &x_space[j];
              (prob->y)[s] = SVMLabels[i].Target;

              for(k = 0; k < numFeatures; k++){

                if(posFeatOnly){
                  // ignore non-positive features
                  if(find(iPosFeats.begin(),iPosFeats.end(),k) != iPosFeats.end()){
                    continue; 
                  }
                }

                x_space[j].index = k;
                if (!Sleipnir::CMeta::IsNaN(d = PCL.Get(iGene, k))) {
                  x_space[j].value = d;
                }else{
                  // impute 0 for missing values
                  x_space[j].value = 0;
                }
                j++;
              }         
              x_space[j].index = -1;
              j++;
              s++;
            }
        }

        SAMPLE* pSample = new SAMPLE;

        pSample->n = prob->l;//number of labels
        pSample->problems = prob;
        pSample->numFeatures = numFeatures;
        pSample->x_space = x_space; 
	return pSample;
}

//TODO: create sample for dab/dat files
//

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

    dec_values = Malloc(double, nr_class*(nr_class-1)/2);
    vecResult.resize(pSample->n);

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

        vecResult[j].GeneName = SVMLabels[i].GeneName;
        vecResult[j].Target = SVMLabels[i].Target;
        vecResult[j].Value = dec_value;
        
        j++;

      }
    }
    FreeSample( *pSample );
    //delete pSample ;
    free(dec_values);
    x = NULL;

    return vecResult;
}

    
//TODO: classify for dab/dat files
//

}
