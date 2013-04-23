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


SAMPLE * CLIBSVM::CreateSample(Sleipnir::CDat& Dat, vector<SVMLabel> SVMLabels) {
	size_t i, j, k, iGene, iProblem, numFeatures, max_index;
        float d;

        struct svm_problem prob;
        struct svm_node *x_space;

        prob.l = 0;//number of labels in Dat
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

                  prob.l++;
		}
	}

 
        prob.y = Malloc(double,prob.l);
        prob.x = Malloc(struct svm_node *, prob.l);
        x_space = Malloc(struct svm_node, numFeatures * prob.l);

        max_index = numFeatures;
        j = 0;

        for (i = 0; i < SVMLabels.size(); i++) {
            iGene = Dat.GetGene(SVMLabels[i].GeneName);
            if (iGene != -1){
              prob.x[i] = &x_space[j];
              prob.y[i] = SVMLabels[i].Target;
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
        SAMPLE* pSample = new SAMPLE;
        pSample->n = prob.l;
        pSample->problems = &prob;

	return pSample;
}

}

