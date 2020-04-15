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

#ifndef NO_SVM_STRUCTTREE
#ifndef SVMSTRUCTTREEI_H
#define SVMSTRUCTTREEI_H

#include "pclset.h"
#include "meta.h"
#include "dat.h"

#include <stdio.h>
#include <map>

#ifndef NO_SVM_STRUCT
#define SVMSTRUCT_H
extern "C" {

#define class Class

#include <svm_hierarchy/svm_light/svm_common.h>
#include <svm_hierarchy/svm_light/svm_learn.h>
#include <svm_hierarchy/svm_struct_api_types.h>
#include <svm_hierarchy/svm_struct/svm_struct_common.h>
#include <svm_hierarchy/svm_struct_api.h>
#include <svm_hierarchy/svm_struct/svm_struct_learn.h>
#undef class

}
#endif

#include "svmstruct.h"

/* removed to support cygwin */
//#include <execinfo.h>

namespace SVMArc {

    struct TEMPNODE { // temporary intermediate data stucture
        set<ONTONODE *> children;
    };

    //class for SVMStruct
    class CSVMSTRUCTTREE : CSVMSTRUCTBASE {

    public:
        LEARN_PARM learn_parm;
        KERNEL_PARM kernel_parm;
        STRUCT_LEARN_PARM struct_parm;
        STRUCTMODEL structmodel;
        map<string, int> onto_map;
        map<int, string> onto_map_rev;


        int Alg;

        CSVMSTRUCTTREE() {
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

        void SetEpsilon(float eps) {
            struct_parm.epsilon = eps;
        }

        void UseSlackRescaling() {
            struct_parm.loss_type = SLACK_RESCALING;
        }

        void UseMarginRescaling() {
            struct_parm.loss_type = MARGIN_RESCALING;
        }

        void SetNThreads(int n) {
            struct_parm.n_threads = n;
        }

        void ReadModel(char *model_file) {

            structmodel = read_struct_model(model_file, &struct_parm);
            if (structmodel.svm_model->kernel_parm.kernel_type == LINEAR) { /* linear kernel */
                /* compute weight vector */
                add_weight_vector_to_linear_model(structmodel.svm_model);
                structmodel.w = structmodel.svm_model->lin_weights;
            }
        }

        void WriteModel(char *model_file) {
            //if (kernel_parm.kernel_type == LINEAR) {
            //	ofstream ofsm;
            //	ofsm.open(model_file);
            //	for (size_t i = 0; i < structmodel.sizePsi; i++) {
            //		ofsm << structmodel.w[i+1] << endl;
            //	}
            //} else {
            write_struct_model(model_file, &structmodel, &struct_parm);
            /*}*/
        }

        void WriteWeights(ostream &osm) {
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

        static void FreeDoc(DOC *pDoc) {
            free_example(pDoc, true);
        }

        void SetVerbosity(size_t V);


        void ReadOntology(const char *treefile);

        //creates a Doc for a given gene index in a single microarray
        static DOC *CreateDoc(Sleipnir::CPCL &PCL, size_t iGene, size_t iDoc);

        //static DOC* CreateDoc(Sleipnir::CDat& Dat, size_t iGene, size_t iDoc);
        //read labels
        vector <SVMLabel> ReadLabels(ifstream &ifsm);

        void vecsetZero(ONTONODE *node, vector<char> *ybar0, char zero);

        void preprocessLabel(vector<char> *multilabels);

        void InitializeLikAfterReadLabels();

        //Creates a sample using a single PCL and SVMlabels Looks up genes by name.
        SAMPLE *CreateSample(Sleipnir::CPCL &PCL, vector <SVMLabel> SVMLabels);
        //SAMPLE* CreateSample(Sleipnir::CDat& Dat, vector<SVMLabel> SVMLabels);

        //Classify single genes
        vector <Result> Classify(Sleipnir::CPCL &PCL, vector <SVMLabel> SVMLabels);

        //MEMBER functions wraps learning
        void Learn(SAMPLE &sample) {
            //cerr << "SLACK NORM =" << struct_parm.slack_norm << endl;
            cerr << "Algorithm " << Alg << " selected." << endl;

            if (Alg == 0)
                svm_learn_struct(sample, &struct_parm, &learn_parm, &kernel_parm, &structmodel, NSLACK_ALG);
            else if (Alg == 1)
                svm_learn_struct(sample, &struct_parm, &learn_parm, &kernel_parm, &structmodel, NSLACK_SHRINK_ALG);
            else if (Alg == 2)
                svm_learn_struct_joint(sample, &struct_parm, &learn_parm, &kernel_parm, &structmodel,
                                       ONESLACK_PRIMAL_ALG);
            else if (Alg == 3)
                svm_learn_struct_joint(sample, &struct_parm, &learn_parm, &kernel_parm, &structmodel,
                                       ONESLACK_DUAL_ALG);
            else if (Alg == 4)
                svm_learn_struct_joint(sample, &struct_parm, &learn_parm, &kernel_parm, &structmodel,
                                       ONESLACK_DUAL_CACHE_ALG);
            else if (Alg == 9)
                svm_learn_struct_joint_custom(sample, &struct_parm, &learn_parm, &kernel_parm, &structmodel);
            else
                exit(1);
            //
        }


        struct SortResults {
            bool operator()(const Result &rOne, const Result &rTwo) const {
                return (rOne.GeneName < rTwo.GeneName); //sort results by name
            }
        };

        size_t PrintResults(vector <Result> vecResults, ofstream &ofsm) {
            sort(vecResults.begin(), vecResults.end(), SortResults());
            int LabelVal;
            for (size_t i = 0; i < vecResults.size(); i++) {
                ofsm << vecResults[i].toStringTREE(&onto_map_rev, 0) << endl;
            }
        };


        bool parms_check();

        bool initialize();


        // free the sample but don't free the Docs
        static void FreeSample_leave_Doc(SAMPLE s);


        STRUCTMODEL read_struct_model_w_linear(char *file, STRUCT_LEARN_PARM *sparm);
    };


};


#endif // NO_SVM_SVMSTRUCTTREE
#endif // SVMSTRUCTTREE_H
