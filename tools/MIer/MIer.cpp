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

struct SSorter {
    const vector <size_t> &m_veciSizes;

    SSorter(const vector <size_t> &veciSizes) : m_veciSizes(veciSizes) {}

    bool operator()(size_t iOne, size_t iTwo) const {

        return (m_veciSizes[iOne] > m_veciSizes[iTwo]);
    }
};

float find_value(size_t iOne, size_t iTwo, const CDat *pDat) {
    return (((iOne == -1) || (iTwo == -1)) ? CMeta::GetNaN() : pDat->Get(iOne, iTwo));
}

size_t find_value(float dValue, const CDataPair *pDat, size_t iDefault, bool fRandom) {
    if (pDat->IsContinuous())
        return 0;

    return (CMeta::IsNaN(dValue) ? (fRandom ? (rand() % pDat->GetValues()) : iDefault) :
            pDat->Quantize(dValue));
}

size_t find_value(size_t iOne, size_t iTwo, const CDataPair *pDat, size_t iDefault, bool fRandom) {
    float fValue = find_value(iOne, iTwo, pDat);
    size_t iValue = find_value(fValue, pDat, iDefault, fRandom);
    return iValue;
}

int main(int iArgs, char **aszArgs) {
    gengetopt_args_info sArgs;
    vector < CDataPair * > vpDatsBigmem;
    vector<float *> vecpadVals;
    size_t iCountJoint, iJoint, iValueOne, iValueTwo;
    vector <size_t> veciDefaults, veciSizes;
    vector <string> vecstrInputs, vecstrGenes;
    float dSubsample;
    map <string, size_t> mapZeros;
    IMeasure *pMeasure;
    CMeasurePearson Pearson;
    CMeasureEuclidean Euclidean;
    CMeasureKendallsTau KendallsTau;
    CMeasureKolmogorovSmirnov KolmSmir;
    CMeasureHypergeometric Hypergeom;
    CMeasureQuickPearson PearQuick;
    CMeasureInnerProduct InnerProd;
    CMeasureBinaryInnerProduct BinInnerProd;

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }
    CMeta Meta(sArgs.verbosity_arg, sArgs.random_arg);

    CMeasureSigmoid EuclideanSig(&Euclidean, false, 1.0f / sArgs.inputs_num);
    IMeasure *apMeasures[] = {&Pearson, &EuclideanSig, &KendallsTau,
                              &KolmSmir, &Hypergeom, &PearQuick, &InnerProd,
                              &BinInnerProd, nullptr
    };


    vecstrInputs.resize(sArgs.inputs_num);
    copy(sArgs.inputs, sArgs.inputs + sArgs.inputs_num, vecstrInputs.begin());
    if (sArgs.only_arg == -1) {
        vector <size_t> veciIndices, veciSizes;

        veciSizes.resize(vecstrInputs.size());
        veciIndices.resize(vecstrInputs.size());
        for (size_t i = 0; i < vecstrInputs.size(); ++i) {
            CDat Dat;
            ifstream ifsm;

            veciIndices[i] = i;
            ifsm.open(vecstrInputs[i].c_str());
            Dat.OpenGenes(ifsm, true);
            veciSizes[i] = Dat.GetGenes();
        }
        sort(veciIndices.begin(), veciIndices.end(), SSorter(veciSizes));
        CMeta::Permute(vecstrInputs, veciIndices);
    }
    {
        CDataset Data;

        Data.OpenGenes(vecstrInputs);
        vecstrGenes.resize(Data.GetGenes());
        copy(Data.GetGeneNames().begin(), Data.GetGeneNames().end(), vecstrGenes.begin());
    }
    veciSizes.resize(vecstrInputs.size());
    for (size_t i = 0; i < veciSizes.size(); ++i) {
        CDataPair DatSize;
        DatSize.OpenQuants(vecstrInputs[i].c_str());
        veciSizes[i] = DatSize.GetValues();
    }

    pMeasure = nullptr;
    if (sArgs.distance_arg)
        for (size_t i = 0; apMeasures[i]; ++i)
            if (!strcmp(apMeasures[i]->GetName(), sArgs.distance_arg)) {
                pMeasure = apMeasures[i];
                break;
            }

    if (sArgs.zeros_arg) {
        ifstream ifsm;
        vector <string> vecstrZeros;
        char acLine[1024];

        ifsm.open(sArgs.zeros_arg);
        if (!ifsm.is_open()) {
            cerr << "Could not open: " << sArgs.zeros_arg << endl;
            return 1;
        }
        while (!ifsm.eof()) {
            ifsm.getline(acLine, ARRAYSIZE(acLine) - 1);
            acLine[ARRAYSIZE(acLine) - 1] = 0;
            vecstrZeros.clear();
            CMeta::Tokenize(acLine, vecstrZeros);
            if (vecstrZeros.empty())
                continue;
            mapZeros[vecstrZeros[0]] = strtol(vecstrZeros[1].c_str(), nullptr, 10);
        }
    }
    veciDefaults.resize(vecstrInputs.size());
    for (size_t i = 0; i < veciDefaults.size(); ++i) {
        map<string, size_t>::const_iterator iterZero;

        if ((iterZero = mapZeros.find(CMeta::Deextension(CMeta::Basename(
                vecstrInputs[i].c_str())))) != mapZeros.end())
            veciDefaults[i] = iterZero->second;
        else
            veciDefaults[i] = sArgs.randomize_flag ? -1 :
                              (sArgs.zero_flag ? 0 : veciSizes[i]++);
    }

    dSubsample = (sArgs.subsample_arg == -1) ? 1 : (2.0f * sArgs.subsample_arg / vecstrGenes.size() /
                                                    (vecstrGenes.size() - 1));
    vector <pair<string, string>> vecpairstrChosen;
    for (size_t i = 0; i < vecstrGenes.size(); ++i) {
        for (size_t j = (i + 1); j < vecstrGenes.size(); ++j) {
            if (((float) rand() / RAND_MAX) < dSubsample) {
                vecpairstrChosen.push_back(make_pair(vecstrGenes[i], vecstrGenes[j]));
            }
        }
    }

    if (sArgs.bigmem_flag) {
        vecpadVals.resize(vecstrInputs.size(), NULL);
        vpDatsBigmem.resize(vecstrInputs.size(), NULL);
#pragma omp parallel for num_threads(sArgs.threads_arg)
        for (size_t iDat = 0; iDat < vecstrInputs.size(); ++iDat) {
            CDataPair *pDat = new CDataPair;
#pragma omp critical
            {
                if (!(pDat->Open(vecstrInputs[iDat].c_str(), false, !!sArgs.memmap_flag) ||
                      pDat->Open(vecstrInputs[iDat].c_str(), true, !!sArgs.memmap_flag))) {
                    cerr << "Could not open: " << vecstrInputs[iDat] << endl;
                    exit(1);
                }
            }
            if (pMeasure) {
                float *adVals = new float[vecpairstrChosen.size()];
                for (size_t i = 0; i < vecpairstrChosen.size(); ++i) {
                    pair <string, string> pairChosen = vecpairstrChosen[i];
                    size_t iGeneOne = pDat->GetGene(pairChosen.first);
                    size_t jGeneOne = pDat->GetGene(pairChosen.second);
                    float dValueOne = find_value(iGeneOne, jGeneOne, pDat);
                    adVals[i] = dValueOne;
                }
                vecpadVals[iDat] = adVals;
                delete pDat;
            } else {
                vpDatsBigmem[iDat] = pDat;
            }
        }
    }


    //Calculate pairwise measures for dataset
    vector<float> vecdMeasures(vecstrInputs.size() * (vecstrInputs.size() - 1) / 2 + vecstrInputs.size(), 0);
    size_t iCompletedPairs = 0;
    size_t iNumDatasets = vecstrInputs.size();
    for (size_t iDatOne = 0; iDatOne < iNumDatasets; ++iDatOne) {
        cerr << "Processing " << iDatOne + 1 << '/' << iNumDatasets << endl;
        CDataPair *DatOne;
        if (sArgs.bigmem_flag) {
            if (!pMeasure) {
                DatOne = vpDatsBigmem[iDatOne];
            }
        } else {
            DatOne = new CDataPair;
            if (!(DatOne->Open(vecstrInputs[iDatOne].c_str(), false, !!sArgs.memmap_flag) ||
                  DatOne->Open(vecstrInputs[iDatOne].c_str(), true, !!sArgs.memmap_flag))) {
                cerr << "Could not open: " << vecstrInputs[iDatOne] << endl;
                return 1;
            }
        }
#pragma omp parallel for num_threads(sArgs.threads_arg)
        for (size_t iDatTwo = iDatOne; iDatTwo <
                                       iNumDatasets; ++iDatTwo) { //Change loop to also include self pair for all measures, simplifies code a great deal
            CDataPair *DatTwo;
            vector <vector<size_t>> vecveciJoint;
            vector <size_t> veciOne, veciTwo;
            size_t iCountOne, iCountTwo;
            float dMeasure = 0;
            float *adOne;
            float *adTwo;

            //Load second dataset
            if (sArgs.bigmem_flag) {
                if (!pMeasure) {
                    DatTwo = vpDatsBigmem[iDatTwo];
                }
            } else {
                DatTwo = new CDataPair;
                if (!(DatTwo->Open(vecstrInputs[iDatTwo].c_str(), false, !!sArgs.memmap_flag) ||
                      DatTwo->Open(vecstrInputs[iDatTwo].c_str(), true, !!sArgs.memmap_flag))) {
                    cerr << "Could not open: " << vecstrInputs[iDatTwo] << endl;
                    exit(1);
                }
            }

            if (pMeasure) {//Calculate with measure
                if (sArgs.bigmem_flag) {
                    adOne = vecpadVals[iDatOne];
                    adTwo = vecpadVals[iDatTwo];
                } else {
                    adOne = new float[vecpairstrChosen.size()];
                    adTwo = new float[vecpairstrChosen.size()];
                    for (size_t i = 0; i < vecpairstrChosen.size(); ++i) {
                        pair <string, string> pairChosen = vecpairstrChosen[i];
                        size_t iGeneOne = DatOne->GetGene(pairChosen.first);
                        size_t jGeneOne = DatOne->GetGene(pairChosen.second);
                        size_t iGeneTwo = DatTwo->GetGene(pairChosen.first);
                        size_t jGeneTwo = DatTwo->GetGene(pairChosen.second);
                        float dValueOne = find_value(iGeneOne, jGeneOne, DatOne);
                        float dValueTwo = find_value(iGeneTwo, jGeneTwo, DatTwo);
                        adOne[i] = dValueOne;
                        adTwo[i] = dValueTwo;
                    }
                }
                dMeasure = (float) pMeasure->Measure(adOne, vecpairstrChosen.size(), adTwo, vecpairstrChosen.size());
                if (!sArgs.bigmem_flag) {
                    delete[] adTwo;
                    delete[] adOne;
                }
            } else {//use MI calculation
                veciOne.resize(veciSizes[iDatOne]);
                fill(veciOne.begin(), veciOne.end(), 0);

                veciTwo.resize(veciSizes[iDatTwo]);
                fill(veciTwo.begin(), veciTwo.end(), 0);

                vecveciJoint.resize(veciOne.size());
                for (size_t i = 0; i < vecveciJoint.size(); ++i) {
                    vecveciJoint[i].resize(veciTwo.size());
                    fill(vecveciJoint[i].begin(), vecveciJoint[i].end(), 0);
                }
                iCountOne = 0;
                iCountTwo = 0;
                iCountJoint = 0;
                for (size_t i = 0; i < vecpairstrChosen.size(); ++i) {
                    pair <string, string> pairChosen = vecpairstrChosen[i];
                    size_t iGeneOne = DatOne->GetGene(pairChosen.first);
                    size_t jGeneOne = DatOne->GetGene(pairChosen.second);
                    size_t iGeneTwo = DatTwo->GetGene(pairChosen.first);
                    size_t jGeneTwo = DatTwo->GetGene(pairChosen.second);
                    float dValueOne = find_value(iGeneOne, jGeneOne, DatOne);
                    float dValueTwo = find_value(iGeneTwo, jGeneTwo, DatTwo);
                    //Get information for MI calculation
                    if ((iValueTwo = find_value(dValueTwo, DatTwo, veciDefaults[iDatTwo], !!sArgs.randomize_flag)) !=
                        -1) {
                        iCountTwo++;
                        veciTwo[iValueTwo]++;
                        if ((iValueOne = find_value(dValueOne, DatOne, veciDefaults[iDatOne],
                                                    !!sArgs.randomize_flag)) != -1) {
                            iCountOne++;
                            iCountJoint++;
                            vecveciJoint[iValueOne][iValueTwo]++;
                        }
                    }
                }

                dMeasure = 0;
                for (size_t i = 0; i < veciOne.size(); ++i) {
                    float dOne = (float) veciOne[i] / iCountOne;
                    for (size_t j = 0; j < veciTwo.size(); ++j)
                        if ((iJoint = vecveciJoint[i][j])) {
                            float dJoint = (float) iJoint / iCountJoint;
                            dMeasure += dJoint * log(dJoint * iCountTwo / dOne / veciTwo[j]);
                        }
                }
                // This corrects for bias introduced by empirical estimation:
                // http://robotics.stanford.edu/~gal/Research/Redundancy-Reduction/Neuron_suppl/node2.html
                dMeasure -= (veciOne.size() - 1) * (veciTwo.size() - 1) / (2.0f * (iCountOne +
                                                                                   iCountTwo));
                dMeasure = (dMeasure < 0) ? 0 : (dMeasure / log(2.0f));
            }

            //store for later output
            vecdMeasures[iCompletedPairs + iDatTwo - iDatOne] = dMeasure;
            if (!sArgs.bigmem_flag) { //Clean up after ourselves if we're not bigmem
                if (iDatOne != iDatTwo) { //don't delete DatOne
                    delete DatTwo;
                }
            }
        }
        iCompletedPairs += iNumDatasets - iDatOne;
        if (!sArgs.bigmem_flag) {
            delete DatOne; //Clean up after ourselves
        }
    }

    //Now on to output
    iCompletedPairs = 0;
    if (sArgs.table_flag) { // Write out table
        for (size_t iDat = 0; iDat < iNumDatasets; ++iDat) {
            cout << '\t' << vecstrInputs[iDat];
        }
        cout << endl;
        for (size_t iDatOne = 0; iDatOne < iNumDatasets; ++iDatOne) {
            cout << vecstrInputs[iDatOne];
            for (size_t iDatTwo = 0; iDatTwo < iNumDatasets; ++iDatTwo) {
                if (iDatTwo < iDatOne) {
                    cout << '\t';
                } else {
                    cout << '\t' << vecdMeasures[iCompletedPairs + iDatTwo - iDatOne];
                }
            }
            cout << endl;
            iCompletedPairs += iNumDatasets - iDatOne;
        }
    } else { // write as pairs
        for (size_t iDatOne = 0; iDatOne < iNumDatasets; ++iDatOne) {
            for (size_t iDatTwo = iDatOne; iDatTwo < iNumDatasets; ++iDatTwo) {
                cout << vecstrInputs[iDatOne] << '\t' << vecstrInputs[iDatTwo] << '\t'
                     << vecdMeasures[iCompletedPairs + iDatTwo - iDatOne] << endl;
            }
            iCompletedPairs += iNumDatasets - iDatOne;
        }
    }

    return 0;
}
