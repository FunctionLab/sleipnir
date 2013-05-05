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
struct svm_node* CLIBSVM::x_space = NULL;

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
        
        model = NULL;

        balance = 0;//Off

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

void CLIBSVM::SetXSpace(Sleipnir::CPCL& PCL) {
	size_t j, k, iGene, numFeatures, numLabels;
        std::vector<std::string> vecGeneNames;
        string GeneName;
        float d;

        numFeatures = PCL.GetExperiments();
        numLabels = PCL.GetGenes();
        vecGeneNames = PCL.GetGeneNames();
        
cerr << "number of features used: " << numFeatures << endl;
cerr << "number of labels in data: " << numLabels << endl;
cerr << "number of gene names: " << vecGeneNames.size() << endl;

        x_space = Malloc(struct svm_node, (1+numFeatures) * numLabels);

        j = 0;//element index
        
        for ( std::vector<std::string>::iterator it = vecGeneNames.begin(); it != vecGeneNames.end(); ++it) {
          GeneName = *it;
          iGene = PCL.GetGene(GeneName); 

          for ( k = 0; k < numFeatures; k++){
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
        }


}

SAMPLE * CLIBSVM::CreateSample(Sleipnir::CPCL& PCL, vector<SVMLabel> SVMLabels) {
	size_t i, j, k, s, iGene, iProblem, numFeatures, numLabels, max_index;
        float d;

        struct svm_problem* prob;

        prob = Malloc(struct svm_problem,1);

        numFeatures = PCL.GetExperiments();
        numLabels = 0;
        
        iProblem = 0;

	for (i = 0; i < SVMLabels.size(); i++) {
                if (!SVMLabels[i].hasIndex){
                  SVMLabels[i].SetIndex(PCL.GetGene(SVMLabels[i].GeneName));
                }
		iGene = SVMLabels[i].index;
		if (iGene != -1) {
                  numLabels++;
		}
	}

cerr << "number of features used: " << numFeatures << endl;
cerr << "number of labels given: " << SVMLabels.size() << endl;
cerr << "number of labels in data: " << numLabels << endl;

        prob->l = numLabels;
        prob->y = Malloc(double,numLabels);
        prob->x = Malloc(struct svm_node *, numLabels);

        if(x_space == NULL) {
          SetXSpace(PCL);
        }

        max_index = numFeatures;

        s = 0;//sample index
        for (i = 0; i < SVMLabels.size(); i++) {
            iGene = SVMLabels[i].index;

            if (iGene != -1){
              (prob->x)[s] = &x_space[iGene*(1+numFeatures)];
              (prob->y)[s] = SVMLabels[i].Target;
              s++;
            }
        }

        SAMPLE* pSample = new SAMPLE;

        pSample->n = prob->l;//number of labels
        pSample->problems = prob;
        pSample->numFeatures = numFeatures;
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
        dec_value = dec_values[0]; //assume that positive class is the first class TODO: currently hacklyi

//cerr << dec_value << endl;

        vecResult[j].GeneName = SVMLabels[i].GeneName;
        vecResult[j].Target = SVMLabels[i].Target;
        vecResult[j].Value = dec_value;
        
        j++;

      }
    }

    free(pSample);
    //delete pSample ;
    free(dec_values);
//    x = NULL;
//    dec_values = NULL;
//    pSample = NULL;

    return vecResult;
}

    
//TODO: classify for dab/dat files
//

}
