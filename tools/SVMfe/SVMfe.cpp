#include <fstream>

#include <vector>

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
#include "cmdline.h"
#include "statistics.h"

using namespace SVMLight;
//#include "../../extlib/svm_light/svm_light/kernel.h"

vector <SVMLight::SVMLabel> ReadLabels(ifstream &ifsm) {

    static const size_t c_iBuffer = 1024;
    char acBuffer[c_iBuffer];
    vector <string> vecstrTokens;
    vector <SVMLight::SVMLabel> vecLabels;
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
        vecLabels.push_back(SVMLight::SVMLabel(vecstrTokens[0], atof(
                vecstrTokens[1].c_str())));
    }
    return vecLabels;
}

struct SortResults {

    bool operator()(const SVMLight::Result &rOne, const SVMLight::Result &rTwo) const {
        return (rOne.Value > rTwo.Value);
    }
};

void PrintResults(vector <SVMLight::Result> vecResults, ofstream &ofsm) {
    sort(vecResults.begin(), vecResults.end(), SortResults());
    int LabelVal;
    for (size_t i = 0; i < vecResults.size(); i++) {
        ofsm << vecResults[i].toString() << endl;
    }
}

void EliminateSvec(SVMLight::SVECTOR *a, bool *p_bUseFeature) {
    register SVMLight::WORD *ai;
    ai = a->words;
    register bool *pUse = p_bUseFeature;
    size_t zeroed = 0;
    pUse++; //the first one is junk since feature numbering starts at 1
    while (ai->wnum) {
        if (!*pUse) {
            ai->weight = 0.0;
            zeroed++;
        }
        ai++;
        pUse++;
    }
    //cerr<<"Zeroed features ="<<zeroed<<endl;
}

void PrintSVec(SVMLight::SVECTOR *a) {
    register SVMLight::WORD *ai;
    ai = a->words;
    while (ai->wnum) {

        cerr << ' ' << ai->weight;
        ai++;
    }
    cerr << endl;
}

void EliminateSample(SVMLight::SAMPLE &sample, bool *p_bUseFeature) {
    for (size_t i = 0; i < sample.examples[0].x.totdoc; i++) {
        EliminateSvec(sample.examples[0].x.doc[i]->fvec, p_bUseFeature);
    }
}

void ZeroWordVec(vector <SVMLight::WORD> &vecWords) {
    for (size_t i = 0; i < vecWords.size(); i++) {
        vecWords[i].weight = 0;
    }
}

void AddToWordVec(vector <SVMLight::WORD> &vecWords, double *a_fW) {
    for (size_t i = 0; i < vecWords.size(); i++) {
        //vecWords[i].weight += svec.words[i].weight;
        vecWords[i].weight += a_fW[i];
    }
}

void CombineWithWordVecMin(vector <SVMLight::WORD> &vecWords, double *a_fW) {
    for (size_t i = 0; i < vecWords.size(); i++) {
        if (vecWords[i].weight * vecWords[i].weight > a_fW[i] * a_fW[i]) {
            vecWords[i].weight = a_fW[i];
        }
    }
}

void CombineWithWordVecMax(vector <SVMLight::WORD> &vecWords, double *a_fW) {
    for (size_t i = 0; i < vecWords.size(); i++) {
        if (vecWords[i].weight * vecWords[i].weight < a_fW[i] * a_fW[i]) {
            vecWords[i].weight = a_fW[i];
        }
    }
}


vector<float> MeanDiff(SAMPLE &sample, size_t iNumFeat) {
    size_t i, j, npos, nneg;
    npos = nneg = 0;
    for (i = 0; i < sample.n; i++) {
        if (sample.examples->y.Class[i] > 0) {
            npos++;
        } else if (sample.examples->y.Class[i] < 0) {
            nneg++;
        }

    }
    vector<float> vecMeanDiff;
    vecMeanDiff.resize(iNumFeat);
    fill(vecMeanDiff.begin(), vecMeanDiff.end(), 0.0f);
    float nfrac, pfrac;
    nfrac = 1.0f / nneg;
    pfrac = 1.0f / npos;
    for (i = 0; i < sample.n; i++) {
        if (sample.examples->y.Class[i] > 0) {
            for (j = 0; j < iNumFeat; j++) {
                vecMeanDiff[j]
                        += sample.examples->x.doc[0]->fvec->words[j].weight;
            }
        } else if (sample.examples->y.Class[i] < 0) {
            for (j = 0; j < iNumFeat; j++) {
                vecMeanDiff[j] -= nfrac
                                  * sample.examples->x.doc[0]->fvec->words[j].weight;
            }
        }

    }
    return vecMeanDiff;
}

bool SortWords(const SVMLight::WORD &first, const SVMLight::WORD &second) {
    if (first.weight * first.weight < second.weight * second.weight)
        return 1;
    else
        return 0;
}

struct ParamStruct {
    vector<float> vecK, vecTradeoff;
    vector <size_t> vecLoss;
    vector<char *> vecNames;
};

ParamStruct ReadParamsFromFile(ifstream &ifsm, string outFile) {
    static const size_t c_iBuffer = 1024;
    char acBuffer[c_iBuffer];
    char *nameBuffer;
    vector <string> vecstrTokens;
    size_t extPlace;
    string Ext, FileName;
    if ((extPlace = outFile.find_first_of(".")) != string::npos) {
        FileName = outFile.substr(0, extPlace);
        Ext = outFile.substr(extPlace, outFile.size());
    } else {
        FileName = outFile;
        Ext = "";
    }
    ParamStruct PStruct;
    size_t index = 0;
    while (!ifsm.eof()) {
        ifsm.getline(acBuffer, c_iBuffer - 1);
        acBuffer[c_iBuffer - 1] = 0;
        vecstrTokens.clear();
        CMeta::Tokenize(acBuffer, vecstrTokens);
        if (vecstrTokens.empty())
            continue;
        if (vecstrTokens.size() != 3) {
            cerr << "Illegal params line (" << vecstrTokens.size() << "): "
                 << acBuffer << endl;
            continue;
        }
        if (acBuffer[0] == '#') {
            cout << "skipping " << acBuffer << endl;
        } else {
            PStruct.vecLoss.push_back(atoi(vecstrTokens[0].c_str()));
            PStruct.vecTradeoff.push_back(atof(vecstrTokens[1].c_str()));
            PStruct.vecK.push_back(atof(vecstrTokens[2].c_str()));
            PStruct.vecNames.push_back(new char[c_iBuffer]);
            if (PStruct.vecLoss[index] == 4 || PStruct.vecLoss[index] == 5)
                sprintf(PStruct.vecNames[index], "%s_l%d_c%4.6f_k%4.3f%s",
                        FileName.c_str(), PStruct.vecLoss[index],
                        PStruct.vecTradeoff[index], PStruct.vecK[index],
                        Ext.c_str());
            else
                sprintf(PStruct.vecNames[index], "%s_l%d_c%4.6f%s",
                        FileName.c_str(), PStruct.vecLoss[index],
                        PStruct.vecTradeoff[index], Ext.c_str());
            index++;
        }

    }
    cout << "DONE reading params file" << endl;
    return PStruct;
}

int main(int iArgs, char **aszArgs) {
    gengetopt_args_info sArgs;
    CPCL PCL;
    SVMLight::CSVMPERF SVM;
    size_t i, j, iGene, jGene;
    ifstream ifsm;
    size_t elimSoFar, toElim;
    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }
    CMeta Meta(sArgs.verbosity_arg, sArgs.random_arg);
    SVM.SetVerbosity(sArgs.svm_verbosity_arg);
    SVM.SetLossFunction(sArgs.error_function_arg);
    if (sArgs.k_value_arg > 1) {
        cerr << "k_value is >1. Setting default 0.5" << endl;
        SVM.SetPrecisionFraction(0.5);
    } else if (sArgs.k_value_arg <= 0) {
        cerr << "k_value is <=0. Setting default 0.5" << endl;
        SVM.SetPrecisionFraction(0.5);
    } else {
        SVM.SetPrecisionFraction(sArgs.k_value_arg);
    }
    SVM.SetTradeoff(sArgs.tradeoff_arg);
    if (sArgs.slack_flag)
        SVM.UseSlackRescaling();

    if (!SVM.parms_check()) {
        cerr << "Sanity check failed, see above errors" << endl;
        return 1;
    }

    cerr << "SVM params are:" << endl;
    cerr << "loss function\t" << SVM.struct_parm.loss_function << endl;
    cerr << "k value\t" << SVM.struct_parm.prec_rec_k_frac << endl;
    cerr << "tradeoff\t" << SVM.struct_parm.C << endl;

    //  cout << "there are " << vecLabels.size() << " labels processed" << endl;
    size_t iFile;
    vector <string> PCLs;
    if (sArgs.input_given) {
        if (!PCL.Open(sArgs.input_arg, sArgs.skip_arg, sArgs.mmap_flag)) {
            cerr << "Could not open input PCL" << endl;
            return 1;
        }
    }

    vector <SVMLight::SVMLabel> vecLabels;
    set <string> setLabeledGenes;
    if (sArgs.labels_given) {
        ifsm.clear();
        ifsm.open(sArgs.labels_arg);
        if (ifsm.is_open())
            vecLabels = ReadLabels(ifsm);
        else {
            cerr << "could not read label file" << endl;
            return 1;
        }
        for (i = 0; i < vecLabels.size(); i++)
            setLabeledGenes.insert(vecLabels[i].GeneName);
    }

    vector <SVMLabel> vecLabelsIndex;
    size_t iIndex;
    for (i = 0; i < vecLabels.size(); i++) {
        iIndex = PCL.GetGene(vecLabels[i].GeneName);
        if (iIndex != -1) {
            vecLabelsIndex.push_back(vecLabels[i]);
            vecLabelsIndex.back().index = iIndex;
        }
    }
    vecLabels = vecLabelsIndex;

    vector <SVMLight::SVMLabel> pTestVector[sArgs.cross_validation_arg];
    vector <SVMLight::SVMLabel> pTrainVector[sArgs.cross_validation_arg];

    size_t nCV = sArgs.cross_validation_arg;
    size_t nLabels = (nCV) * (nCV - 1) + nCV;
    size_t iC, iB;

    for (i = 0; i < nCV; i++) {
        pTestVector[i].reserve((size_t) vecLabels.size() / nCV + nCV);
        pTrainVector[i].reserve((size_t) vecLabels.size() / (nCV) * (nCV - 1)
                                + nCV);
        for (j = 0; j < vecLabels.size(); j++) {
            if (j % sArgs.cross_validation_arg == i) {
                pTestVector[i].push_back(vecLabels[j]);
            } else {
                pTrainVector[i].push_back((vecLabels[j]));
            }
        }
    }

    vector <vector<size_t>> *pvecvecIndex[nCV];
    for (iC = 0; iC < nCV; iC++) {
        pvecvecIndex[iC] = new vector <vector<size_t>>;
        pvecvecIndex[iC]->resize(sArgs.bootstrap_arg);
        cerr << "numer of vectors is " << pvecvecIndex[iC]->size() << endl;
    }

    //use for bootstrap
    if (sArgs.bootstrap_arg > 1) {
        for (iC = 0; iC < nCV; iC++) {
            nLabels = pTrainVector[iC].size();

            for (iB = 0; iB < sArgs.bootstrap_arg; iB++) {
                //	pvecvecIndex[iC][iB].reserve(nLabels);
                for (i = 0; i < nLabels; i++) {
                    iIndex = (size_t)(rand() / (float) RAND_MAX * nLabels);
                    (*pvecvecIndex[iC])[iB].push_back(iIndex);
                }
            }
        }
    } else {
        for (iC = 0; iC < nCV; iC++) {
            nLabels = pTrainVector[iC].size();
            for (i = 0; i < nLabels; i++) {
                (*pvecvecIndex[iC])[0].push_back(i);
            }
        }
    }
    //build all the samples since they are being reused



    SVMLight::SAMPLE **pppTrainSample[nCV];

    for (iC = 0; iC < nCV; iC++) {
        cerr << iC << '\t' << pTrainVector[iC].size() << endl;

        cerr << iC << '\t' << pvecvecIndex[iC]->size() << endl;
        pppTrainSample[iC] = SVMLight::CSVMPERF::CreateSampleBootStrap(PCL,
                                                                       pTrainVector[iC], (*pvecvecIndex[iC]));
//		for (iB = 0; iB < sArgs.bootstrap_arg; iB++) {
//			*pppTrainSample[iC][iB];
//		}
    }



    //feature elimination data
    size_t numFeat = PCL.GetExperiments() + 1;
    bool *p_bUseFeature = new bool[numFeat + 1];

    // file names
    stringstream ss;
    ss.setf(ios_base::fmtflags(), ios_base::floatfield);
    ss << sArgs.output_arg;
    ss << "_cv" << sArgs.cross_validation_arg;
    ss << "_b" << sArgs.bootstrap_arg;
    ss << "_l" << sArgs.error_function_arg;
    if (sArgs.k_value_given)
        ss << "_k" << sArgs.k_value_arg;
    ss << "_c" << sArgs.tradeoff_arg;
    string strBaseFile = ss.str();
    static const size_t c_iBuffer = 1024;
    char outFile[c_iBuffer];
    ss << "_elim";
    ofstream ofsmElim;
    ofsmElim.open(ss.str().c_str());
    vector < ofstream * > vecpOfsmIter;

    size_t numIter;
    size_t **ppFeatureInfo = new size_t *[nCV];
    for (i = 0; i < nCV; i++)
        ppFeatureInfo[i] = new size_t[numFeat];


    vector <SVMLight::WORD> vecWav;
    vector <SVMLight::WORD> vecWfe;
    vecWav.resize(numFeat);
    vecWfe.resize(numFeat);
    //vectors to hold the result file
    vector <SVMLight::Result> AllResults;

    for (iC = 0; iC < nCV; iC++) {

        cerr << "Cross validation round " << iC << endl;
        size_t currNumFeat = numFeat;
        elimSoFar = 0;
        size_t iRound = 0;
        size_t iOffset = 0;
        set <size_t> sElim;
        fill(p_bUseFeature, p_bUseFeature + numFeat, true);
        while (currNumFeat >= sArgs.min_left_arg && iRound
                                                    <= sArgs.iter_max_arg) { //feature elimination loop
            //maintenance
            //reset the wnum values from the last sort
            if (iC == 0) {
                sprintf(outFile, "%s_round%d", strBaseFile.c_str(), iRound);
                cerr << outFile << endl;
                vecpOfsmIter.resize(iRound + 1);
                vecpOfsmIter[iRound] = new ofstream(outFile);

            }
            cerr << "Elimination round " << iRound << endl;
            for (i = 0; i < numFeat; i++)
                vecWfe[i].wnum = (int32_t) i;
            AllResults.resize(0); //clear results
            ZeroWordVec(vecWfe);
            ZeroWordVec(vecWav); //clear W

            //do the learning and collecting

            for (iB = 0; iB < sArgs.bootstrap_arg; iB++) {
                cerr << "Bootstrap round " << iB << endl;
                SVM.Learn(*pppTrainSample[iC][iB]);
                if (sArgs.min_flag && iB) {
                    CombineWithWordVecMin(vecWfe, SVM.structmodel.w);
                } else if (sArgs.max_flag && iB) {
                    CombineWithWordVecMax(vecWfe, SVM.structmodel.w);
                } else {
                    AddToWordVec(vecWfe, SVM.structmodel.w);
                }
                AddToWordVec(vecWav, SVM.structmodel.w);
                if (sArgs.verbosity_arg > 2) {
                    cerr << "W:";
                    SVM.WriteWeights(cerr);
                }
            }

            SVM.ReplaceModel(vecWav);
            AllResults = SVM.Classify(PCL, pTestVector[iC]);
            cerr << "Classified " << AllResults.size() << endl;
            for (j = 0; j < AllResults.size(); j++)
                AllResults[j].CVround = iC;
            /*AllResults.insert(AllResults.end(), tmpAllResults.begin(),
                    tmpAllResults.end());*/
            PrintResults(AllResults, (*vecpOfsmIter[iRound]));

            //			if (iRound > 0 || i > 0) //free the model if it exists
            //			SVM.FreeModel();

            if (sArgs.verbosity_arg > 1) {
                cerr << "SUM:";
                for (i = 0; i < vecWfe.size(); i++)
                    cerr << ' ' << vecWfe[i].weight;
                cerr << endl;
            }
            //do the sorting
            sort(vecWfe.begin(), vecWfe.end(), SortWords);
            iRound++;
            //eliminate (half for now)
            if (sArgs.verbosity_arg > 4) {
                cerr << "Before:";
                for (i = 0; i < numFeat; i++)
                    cerr << (int) p_bUseFeature[i] << ' ';
                cerr << endl;
            }
            i = 0;
            while (!vecWfe[i].weight) {
                i++;
            }
            //i--;
            iOffset = i;
            cerr << "Offset by number of features " << iOffset << endl;
            cerr << "Number of Features " << currNumFeat << endl;
            toElim = currNumFeat * sArgs.elim_fraction_arg;
            cerr << "will eliminate " << toElim << " Features" << endl;

            for (i = iOffset; i < toElim + iOffset; i++) {
                elimSoFar++;
                if (p_bUseFeature[(size_t) vecWfe[i].wnum]) {
                    if (sArgs.verbosity_arg > 1)
                        cerr << "Elimintating: " << vecWfe[i].wnum
                             << " with weight " << vecWfe[i].weight << endl;
                    p_bUseFeature[(size_t) vecWfe[i].wnum] = false;
                    if (i)
                        ofsmElim << iC << '\t' << vecWfe[i].wnum << '\t'
                                 << PCL.GetExperiment(vecWfe[i].wnum - 1)
                                 << '\t' << iRound << endl;
                }
            }
            cerr << "So far eliminated " << elimSoFar << " features" << endl;
            if (sArgs.verbosity_arg > 4) {
                cerr << "After:";
                for (i = 0; i < numFeat; i++)
                    cerr << (int) p_bUseFeature[i] << ' ';
                cerr << endl;
            }
            currNumFeat = numFeat - elimSoFar;
            //= (size_t) (currNumFeat * (1 - sArgs.elim_fraction_arg));
            for (iB = 0; iB < sArgs.bootstrap_arg; iB++) {
                EliminateSample(*pppTrainSample[iC][iB], p_bUseFeature);
            }
            //EliminateSample(*ppTestSample[i], p_bUseFeature);
        }
        iRound++;
        for (i = 1; i < numFeat; i++) {
            if (p_bUseFeature[i])
                ofsmElim << iC << '\t' << i << '\t' << PCL.GetExperiment(i - 1)
                         << '\t' << iRound << endl;
        }
        for (i = 1; i < numFeat; i++) {
            if (p_bUseFeature[i])
                ofsmElim << iC << '\t' << i << '\t' << PCL.GetExperiment(i - 1)
                         << '\t' << "final" << endl;
        }

    }

    return 0;
}

