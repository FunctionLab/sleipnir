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
#include <random>

enum ETFPN {
    ETFPN_TP = 0,
    ETFPN_FP = ETFPN_TP + 1,
    ETFPN_TN = ETFPN_FP + 1,
    ETFPN_FN = ETFPN_TN + 1
};

struct SDatum {
    float m_dValue;
    size_t m_iOne;
    size_t m_iTwo;
    float m_dAnswer;

    SDatum(float dValue, size_t iOne, size_t iTwo, float dAnswer) : m_dValue(dValue), m_iOne(iOne), m_iTwo(iTwo),
                                                                    m_dAnswer(dAnswer) {}
};

struct SSorter {
    bool m_fInvert;

    explicit SSorter(bool fInvert) : m_fInvert(fInvert) {}

    bool operator()(const SDatum &sOne, const SDatum &sTwo) const {

        return (m_fInvert ? (sOne.m_dValue > sTwo.m_dValue) : (sOne.m_dValue < sTwo.m_dValue));
    }
};

double AUCMod(const CDat &, const CDat &, const vector<bool> &, const vector<bool> &, const gengetopt_args_info &, bool,
              float);

int main(int iArgs, char **aszArgs) {
    CDat Answers, Data, wDat;
    gengetopt_args_info sArgs{};
    size_t i, j, k, m, iOne, iTwo, iGenes, iPositives, iNegatives, iBins, iRand;
    vector <size_t> veciGenes, veciRec, veciRecTerm;
    CFullMatrix<bool> MatGenes;
    CFullMatrix <size_t> MatResults;
    ETFPN eTFPN;
    int iMax;
    float dAnswer, dValue;
    vector<bool> vecfHere, vecfUbik;
    vector<float> vecdScores, vecdSSE, vecdBinValue, vecGeneWeights;
    vector <size_t> veciPositives, veciNegatives, veciGenesTerm;
    CGenome Genome;
    CGenes wGenes(Genome);
    ofstream ofsm;
    ostream *postm;
    map<float, size_t> mapValues;
    bool fMapAnswers, isDatWeighted = false;
    CGenes GenesTm(Genome), GenesUbik(Genome);

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }
    CMeta Meta(sArgs.verbosity_arg, sArgs.random_arg);

    fMapAnswers = !!sArgs.memmap_flag && !(sArgs.genes_arg || sArgs.genet_arg || sArgs.genex_arg || sArgs.genee_arg);
    if (!Answers.Open(sArgs.answers_arg, fMapAnswers)) {
        cerr << "Couldn't open: " << sArgs.answers_arg << endl;
        return 1;
    }
    if (sArgs.genes_arg && !Answers.FilterGenes(sArgs.genes_arg, CDat::EFilterInclude)) {
        cerr << "Couldn't open: " << sArgs.genes_arg << endl;
        return 1;
    }
    if (sArgs.genee_arg && !Answers.FilterGenes(sArgs.genee_arg, CDat::EFilterEdge)) {
        cerr << "Couldn't open: " << sArgs.genee_arg << endl;
        return 1;
    }
    if (sArgs.genet_arg) {
        if (!(Answers.FilterGenes(sArgs.genet_arg, CDat::EFilterTerm) &&
              GenesTm.Open(sArgs.genet_arg))) {
            cerr << "Couldn't open: " << sArgs.genet_arg << endl;
            return 1;
        }
        veciGenesTerm.reserve(GenesTm.GetGenes());
        for (i = 0; i < GenesTm.GetGenes(); ++i)
            if ((j = Answers.GetGene(GenesTm.GetGene(i).GetName())) != -1)
                veciGenesTerm.push_back(j);
    }
    if (sArgs.genex_arg && !Answers.FilterGenes(sArgs.genex_arg, CDat::EFilterExclude)) {
        cerr << "Couldn't open: " << sArgs.genex_arg << endl;
        return 1;
    }

    if (sArgs.genec_given) {
        CGenome tGenome;
        CGenes Genes(tGenome);
        CGenome uGenome;
        CGenes Genesubiq(uGenome);
        vector<bool> vecfctxt;
        vector<bool> vecfubiq;
        float d;

        if (!Genes.Open(sArgs.genec_arg)) {
            cerr << "Couldn't open: " << sArgs.genec_arg << endl;
            return 1;
        }

        vecfctxt.resize(Answers.GetGenes());
        for (i = 0; i < vecfctxt.size(); ++i) {
            vecfctxt[i] = Genes.IsGene(Answers.GetGene(i));
        }

        if (sArgs.ubiqg_given) {
            if (!Genesubiq.Open(sArgs.ubiqg_arg)) {
                cerr << "Could not open: " << sArgs.ubiqg_arg << endl;
                return 1;
            }
            vecfubiq.resize(Answers.GetGenes());
            for (i = 0; i < vecfubiq.size(); ++i) {
                vecfubiq[i] = Genesubiq.IsGene(Answers.GetGene(i));
            }
        }


        for (i = 0; i < Answers.GetGenes(); ++i) {
            for (j = (i + 1); j < Answers.GetGenes(); ++j) {

                if (CMeta::IsNaN(d = Answers.Get(i, j)))
                    continue;

                // remove all original negatives
                if (d < 1) {
                    Answers.Set(i, j, CMeta::GetNaN());
                    continue;
                }

                if (!vecfctxt[i] && !vecfctxt[j])
                    Answers.Set(i, j, 0);
                else if (vecfctxt[i] && vecfctxt[j])
                    continue;
                else if (sArgs.ubiqg_given &&
                         ((vecfubiq[i] && !vecfctxt[j]) ||
                          (vecfubiq[j] && !vecfctxt[i])))
                    Answers.Set(i, j, 0);
                else
                    Answers.Set(i, j, CMeta::GetNaN());
            }
        }
    }

    if (!Data.Open(sArgs.input_arg, !!sArgs.memmap_flag)) {
        cerr << "Couldn't open: " << sArgs.input_arg << endl;
        return 1;
    }

    if (sArgs.singlegene_flag) {
        float d;
        vector <size_t> veciPos, veciNeg, veciIndex, veciDgenes;
        size_t pj, nj, x, y, iOne, iTwo;
        CDat sAnswers;

        veciDgenes.resize(Answers.GetGenes());
        for (i = 0; i < Answers.GetGenes(); ++i)
            veciDgenes[i] = Data.GetGene(Answers.GetGene(i));

        //copy over
        for (i = 0; i < Answers.GetGenes(); ++i) {
            iOne = veciDgenes[i];
            for (j = (i + 1); j < Answers.GetGenes(); ++j) {
                iTwo = veciDgenes[j];
                if (iOne == -1 || iTwo == -1 || CMeta::IsNaN(d = Data.Get(iOne, iTwo)))
                    Answers.Set(i, j, CMeta::GetNaN());
            }
        }

        sAnswers.Open(Answers);
        for (i = 0; i < Answers.GetGenes(); ++i) {
            for (j = (i + 1); j < Answers.GetGenes(); ++j) {
                Answers.Set(i, j, CMeta::GetNaN());
            }
        }

        veciIndex.resize(sAnswers.GetGenes());
        for (i = 0; i < sAnswers.GetGenes(); ++i) {
            veciIndex[i] = i;
        }

        random_device rand_dev;
        mt19937 rand_gen(rand_dev());
        std::shuffle(veciIndex.begin(), veciIndex.end(), rand_gen);

        for (x = 0; x < sAnswers.GetGenes(); ++x) {
            i = veciIndex[x];

            veciPos.clear();
            veciNeg.clear();

            for (j = 0; j < sAnswers.GetGenes(); ++j) {
                if (i == j || CMeta::IsNaN(d = sAnswers.Get(i, j)))
                    continue;

                if (d == 0) {
                    veciNeg.push_back(j);
                } else if (d == 1) {
                    veciPos.push_back(j);
                }
            }

            pj = -1;
            nj = -1;
            if (!veciPos.empty())
                pj = veciPos[rand() % veciPos.size()];
            if (!veciNeg.empty())
                nj = veciNeg[rand() % veciNeg.size()];

            if (pj != -1)
                Answers.Set(i, pj, 1);
            if (nj != -1)
                Answers.Set(i, nj, 0);

            // remove gene rows that have already been sampled
            for (j = 0; j < sAnswers.GetGenes(); ++j) {
                if (i != j && j != pj && j != nj)
                    sAnswers.Set(i, j, CMeta::GetNaN());
            }

            if (pj != -1) {
                for (j = 0; j < sAnswers.GetGenes(); ++j) {
                    if (j == pj || CMeta::IsNaN(d = sAnswers.Get(j, pj)))
                        continue;
                    if (d == 1)
                        sAnswers.Set(j, pj, CMeta::GetNaN());
                }
            }
            if (nj != -1) {
                for (j = 0; j < sAnswers.GetGenes(); ++j) {
                    if (j == nj || CMeta::IsNaN(d = sAnswers.Get(j, nj)))
                        continue;
                    if (d == 0)
                        sAnswers.Set(j, nj, CMeta::GetNaN());
                }
            }
        }
    }

    if (sArgs.normalize_flag)
        Data.Normalize(CDat::ENormalizeMinMax);


    if (sArgs.weights_arg) {
        ifstream ifsm;
        ifsm.open(sArgs.weights_arg);
        if (!wGenes.OpenWeighted(ifsm)) {
            if (!wDat.Open(sArgs.weights_arg, !!sArgs.memmap_flag)) {
                cerr << "Couldn't open: " << sArgs.inputs[i] << endl;
                return 1;
            } else {
                isDatWeighted = true;
            }
        } else {
            vecGeneWeights.resize(Answers.GetGenes());
            for (i = 0; i < vecGeneWeights.size(); ++i) {
                vecGeneWeights[i] = wGenes.GetGeneWeight(wGenes.GetGene(Answers.GetGene(i)));
            }
        }
    }

    if (sArgs.abs_arg) {
        float d;
        for (i = 0; i < Data.GetGenes(); ++i) {
            for (j = (i + 1); j < Data.GetGenes(); ++j) {
                if (!CMeta::IsNaN((d = Data.Get(i, j)))) {
                    Data.Set(i, j, fabs(d - sArgs.abs_arg));
                }
            }
        }
    }

    veciGenes.resize(Answers.GetGenes());
    for (i = 0; i < Answers.GetGenes(); ++i)
        veciGenes[i] = Data.GetGene(Answers.GetGene(i));
    for (iRand = 0; iRand <= (size_t) sArgs.randomize_arg; ++iRand) {
        if (iRand)
            Data.Randomize();
        if (sArgs.finite_flag) {
            vector<float> vecdValues;
            {
                set<float> setdValues;

                for (i = 0; i < Answers.GetGenes(); ++i) {
                    if ((iOne = veciGenes[i]) == -1)
                        continue;
                    for (j = (i + 1); j < Answers.GetGenes(); ++j) {
                        if (((iTwo = veciGenes[j]) == -1) ||
                            CMeta::IsNaN(dValue = Data.Get(iOne, iTwo)) ||
                            CMeta::IsNaN(Answers.Get(i, j)))
                            continue;
                        if (sArgs.invert_flag)
                            dValue = 1 - dValue;
                        setdValues.insert(dValue);
                    }
                }
                vecdValues.resize(setdValues.size());
                copy(setdValues.begin(), setdValues.end(), vecdValues.begin());
            }
            sort(vecdValues.begin(), vecdValues.end());
            for (i = 0; i < vecdValues.size(); ++i)
                mapValues[vecdValues[i]] = i;
            iBins = mapValues.size();
        } else
            iBins = sArgs.bins_arg;
        MatResults.Initialize(iBins ? (iBins + 1) :
                              (size_t)((sArgs.max_arg - sArgs.min_arg) / sArgs.delta_arg) + 1, 4);
        MatGenes.Initialize(veciGenes.size(), MatResults.GetRows());

        for (iGenes = 0; !sArgs.inputs_num || (iGenes < sArgs.inputs_num); ++iGenes) {
            MatResults.Clear();
            MatGenes.Clear();

            if (sArgs.genep_arg) {
                CGenes Genes(Genome);
                if (!Genes.Open(sArgs.genep_arg)) {
                    cerr << "Couldn't open: " << sArgs.genep_arg << endl;
                    return 1;
                }
                vecfHere.resize(Answers.GetGenes());
                for (i = 0; i < vecfHere.size(); ++i) {
                    vecfHere[i] = Genes.IsGene(Answers.GetGene(i));
                }
            }
            if (sArgs.ubiqg_arg) {
                if (!GenesUbik.Open(sArgs.ubiqg_arg)) {
                    cerr << "Could not open: " << sArgs.ubiqg_arg << endl;
                    return 1;
                }
                vecfUbik.resize(Answers.GetGenes());
                for (i = 0; i < vecfUbik.size(); ++i) {
                    vecfUbik[i] = GenesUbik.IsGene(Answers.GetGene(i));
                }
            }

            if (!mapValues.empty()) {
                for (i = 0; i < Answers.GetGenes(); ++i) {
                    if ((iOne = veciGenes[i]) == -1)
                        continue;
                    for (j = (i + 1); j < Answers.GetGenes(); ++j) {
                        if (((iTwo = veciGenes[j]) == -1) ||
                            CMeta::IsNaN(dValue = Data.Get(iOne, iTwo)) ||
                            CMeta::IsNaN(dAnswer = Answers.Get(i, j)))
                            continue;
                        if (CMeta::SkipEdge(!!dAnswer, i, j, vecfHere, vecfUbik, sArgs.ctxtpos_flag, sArgs.ctxtneg_flag,
                                            sArgs.bridgepos_flag, sArgs.bridgeneg_flag, sArgs.outpos_flag,
                                            sArgs.outneg_flag)) {
                            continue;
                        }
                        if (sArgs.invert_flag)
                            dValue = 1 - dValue;
                        for (k = 0; k <= mapValues[dValue]; ++k) {
                            MatGenes.Set(i, k, true);
                            MatGenes.Set(j, k, true);
                            MatResults.Get(k, dAnswer ? ETFPN_TP : ETFPN_FP)++;
                        }
                        for (; k < MatResults.GetRows(); ++k)
                            MatResults.Get(k, dAnswer ? ETFPN_FN : ETFPN_TN)++;
                    }
                }
            } else if (iBins) {
                vector <SDatum> vecsData;
                size_t iChunk;

                for (iPositives = iNegatives = i = 0; i < Answers.GetGenes(); ++i) {
                    if ((iOne = veciGenes[i]) == -1)
                        continue;
                    for (j = (i + 1); j < Answers.GetGenes(); ++j) {
                        if (((iTwo = veciGenes[j]) == -1) ||
                            CMeta::IsNaN(dAnswer = Answers.Get(i, j)) ||
                            CMeta::IsNaN(dValue = Data.Get(iOne, iTwo)))
                            continue;
                        if (CMeta::SkipEdge(!!dAnswer, i, j, vecfHere, vecfUbik, sArgs.ctxtpos_flag, sArgs.ctxtneg_flag,
                                            sArgs.bridgepos_flag, sArgs.bridgeneg_flag, sArgs.outpos_flag,
                                            sArgs.outneg_flag)) {
                            continue;
                        }

                        MatGenes.Set(i, 0, true);
                        MatGenes.Set(j, 0, true);
                        if (dAnswer)
                            iPositives++;
                        else
                            iNegatives++;
                        vecsData.emplace_back(dValue, i, j, dAnswer);
                    }
                }
                sort(vecsData.begin(), vecsData.end(), SSorter(!!sArgs.invert_flag));
                //instead of putting all of the uneveness in one bin, spread it out into each bin.
                //N.B. only the part without the sse_flag is fixed in this regard
                auto perChunk = (size_t)(vecsData.size() / MatResults.GetRows());
                auto chnkRem = (size_t)(vecsData.size() % MatResults.GetRows());
                iChunk = (size_t)(0.5 + ((float) vecsData.size() / (MatResults.GetRows())));
                if (sArgs.sse_flag) {
                    vecdSSE.resize(MatResults.GetRows());
                    veciPositives.resize(vecdSSE.size());
                    for (i = 1, j = 0; i < vecdSSE.size(); ++i, j += iChunk) {
                        veciPositives[veciPositives.size() - i - 1] = veciPositives[veciPositives.size() - i];
                        vecdSSE[vecdSSE.size() - i - 1] = vecdSSE[vecdSSE.size() - i];
                        for (k = 0; k < iChunk; ++k) {
                            if ((j + k) >= vecsData.size())
                                break;
                            const SDatum &sDatum = vecsData[vecsData.size() - (j + k) - 1];

                            for (m = 0; m < (vecdSSE.size() - i); ++m) {
                                MatGenes.Set(sDatum.m_iOne, m, true);
                                MatGenes.Set(sDatum.m_iTwo, m, true);
                            }
                            dValue = sDatum.m_dValue - sDatum.m_dAnswer;
                            veciPositives[veciPositives.size() - i - 1]++;
                            vecdSSE[vecdSSE.size() - i - 1] += dValue * dValue;
                        }
                    }
                } else {
                    veciPositives.resize(MatResults.GetRows() - 1);
                    veciNegatives.resize(veciPositives.size());
                    vecdBinValue.resize(veciPositives.size() + 1);
                    for (i = 0; i < veciNegatives.size(); ++i)
                        veciNegatives[i] = veciPositives[i] = 0;
                    for (i = j = 0; i < veciPositives.size() + 1; ++i, j += k) {
                        size_t thisChunk = (i < chnkRem) ? (perChunk + 1) : (perChunk);
                        for (k = 0; k < thisChunk; ++k) {
                            //if( ( j + k ) >= vecsData.size( ) )
                            //	break;
                            const SDatum &sDatum = vecsData[j + k];
                            vecdBinValue[i] = sDatum.m_dValue;

                            for (m = i; m > 0; --m) {
                                MatGenes.Set(sDatum.m_iOne, m, true);
                                MatGenes.Set(sDatum.m_iTwo, m, true);
                            }
                            if (i >= veciPositives.size())
                                continue;
                            if (Answers.Get(sDatum.m_iOne, sDatum.m_iTwo))
                                veciPositives[i]++;
                            else
                                veciNegatives[i]++;
                        }
                    }
                    MatResults.Set(0, ETFPN_TP, iPositives);
                    MatResults.Set(0, ETFPN_FP, iNegatives);
                    MatResults.Set(0, ETFPN_TN, 0);
                    MatResults.Set(0, ETFPN_FN, 0);
                    for (i = 1; i < MatResults.GetRows(); ++i) {
                        MatResults.Set(i, ETFPN_TP, MatResults.Get(i - 1, ETFPN_TP) - veciPositives[i - 1]);
                        MatResults.Set(i, ETFPN_FP, MatResults.Get(i - 1, ETFPN_FP) - veciNegatives[i - 1]);
                        MatResults.Set(i, ETFPN_TN, MatResults.Get(i - 1, ETFPN_TN) + veciNegatives[i - 1]);
                        MatResults.Set(i, ETFPN_FN, MatResults.Get(i - 1, ETFPN_FN) +
                                                    veciPositives[i - 1]);
                    }
                }
            } else
                for (i = 0; i < Answers.GetGenes(); ++i) {
                    if (!(i % 1000))
                        cerr << "Processing gene " << i << '/' << Answers.GetGenes() << endl;
                    if ((iOne = veciGenes[i]) == -1)
                        continue;
                    for (j = (i + 1); j < Answers.GetGenes(); ++j) {
                        if (((iTwo = veciGenes[j]) == -1) ||
                            CMeta::IsNaN(dAnswer = Answers.Get(i, j)) ||
                            CMeta::IsNaN(dValue = Data.Get(iOne, iTwo)))
                            continue;
                        if (CMeta::SkipEdge(!!dAnswer, i, j, vecfHere, vecfUbik, sArgs.ctxtpos_flag, sArgs.ctxtneg_flag,
                                            sArgs.bridgepos_flag, sArgs.bridgeneg_flag, sArgs.outpos_flag,
                                            sArgs.outneg_flag)) {
                            continue;
                        }

                        if (sArgs.invert_flag)
                            dValue = 1 - dValue;

                        iMax = (int) ceil((dValue - sArgs.min_arg) / sArgs.delta_arg);
                        if (iMax > (int) MatResults.GetRows())
                            iMax = (int) MatResults.GetRows();
                        eTFPN = (ETFPN) !dAnswer;
                        for (k = 0; (int) k < iMax; ++k) {
                            MatResults.Get(k, eTFPN)++;
                            MatGenes.Set(i, k, true);
                            MatGenes.Set(j, k, true);
                        }
                        eTFPN = (ETFPN) (2 + !eTFPN);
                        for (; k < (int) MatResults.GetRows(); ++k)
                            MatResults.Get(k, eTFPN)++;
                    }
                }
            for (iPositives = iNegatives = i = 0; i < Answers.GetGenes(); ++i)
                for (j = (i + 1); j < Answers.GetGenes(); ++j) {
                    if (CMeta::IsNaN(dAnswer = Answers.Get(i, j))) {
                        continue;
                    }

                    if (CMeta::SkipEdge(!!dAnswer, i, j, vecfHere, vecfUbik, sArgs.ctxtpos_flag, sArgs.ctxtneg_flag,
                                        sArgs.bridgepos_flag, sArgs.bridgeneg_flag, sArgs.outpos_flag,
                                        sArgs.outneg_flag)) {
                        continue;
                    }

                    if (dAnswer)
                        iPositives++;
                    else {
                        iNegatives++;
                    }
                }

            veciRec.resize(MatResults.GetRows());
            veciRecTerm.resize(MatResults.GetRows());
            for (i = 0; i < veciRec.size(); ++i) {
                veciRec[i] = veciRecTerm[i] = 0;
                for (j = 0; j < MatGenes.GetRows(); ++j)
                    if (MatGenes.Get(j, i)) {
                        veciRec[i]++;
                        if (vecfHere.size() && vecfHere[j])
                            veciRecTerm[i]++;
                    }
                for (j = 0; j < veciGenesTerm.size(); ++j)
                    if (MatGenes.Get(veciGenesTerm[j], i) &&
                        (vecfHere.empty() || !vecfHere[veciGenesTerm[j]]))
                        veciRecTerm[i]++;
            }

            if (sArgs.inputs_num) {
                ofsm.open(((string) sArgs.directory_arg + '/' +
                           CMeta::Basename(sArgs.inputs[iGenes]) + ".bins").c_str());
                postm = &ofsm;
            } else
                postm = &cout;

            if (!sArgs.sse_flag) { //Weighted context is currently only supported for AUC calculation
                *postm << "#	P	" << iPositives << endl;
                *postm << "#	N	" << iNegatives << endl;
            }

            if (!sArgs.weights_arg) {
                *postm << "Cut	Genes	"
                       << (sArgs.sse_flag ? "Pairs	SSE" : "TP	FP	TN	FN	RC	PR	VALUE") << endl;
                for (i = 0; i < MatResults.GetRows(); ++i) {
                    *postm << (iBins ? i : (sArgs.min_arg + (i * sArgs.delta_arg))) << '\t' <<
                           veciRec[i];
                    if (sArgs.sse_flag)
                        *postm << '\t' << veciPositives[i] << '\t' << vecdSSE[i];
                    else
                        for (j = 0; j < MatResults.GetColumns(); ++j)
                            *postm << '\t' << MatResults.Get(i, j);
                    if (veciGenesTerm.size() || vecfHere.size())
                        *postm << '\t' << veciRecTerm[i];

                    // print precision/recall
                    *postm << '\t' << (float) MatResults.Get(i, 0) / (MatResults.Get(0, 0));

                    if ((MatResults.Get(i, 1) + MatResults.Get(i, 0)) != 0)
                        *postm << '\t' << (float) MatResults.Get(i, 0) / (MatResults.Get(i, 1) + MatResults.Get(i, 0));
                    else
                        *postm << '\t' << 0.0;
                    *postm << '\t' << vecdBinValue[i];

                    *postm << endl;
                }
            } else {
                *postm << "AUC is calculated using weighted context." << endl;
                if (sArgs.flipneg_flag)
                    *postm << "Flipneg ON" << endl;
                else
                    *postm << "Flipneg OFF" << endl;
            }

            if (!sArgs.sse_flag) {
                if (sArgs.weights_arg) {
                    if (isDatWeighted)
                        *postm << "#	AUC	"
                               << CStatistics::WilcoxonRankSum(Data, Answers, wDat, sArgs.flipneg_flag) << endl;
                    else
                        *postm << "#	AUC	"
                               << CStatistics::WilcoxonRankSum(Data, Answers, vecGeneWeights, sArgs.flipneg_flag)
                               << endl;
                } else {
                    if (sArgs.auc_arg)
                        *postm << "#	AUC	"
                               << AUCMod(Data, Answers, vecfHere, vecfUbik, sArgs, !!sArgs.invert_flag, sArgs.auc_arg)
                               << endl;
                    else {
                        *postm << "#	AUC	"
                               << CStatistics::WilcoxonRankSum(Data, Answers, vecfHere, vecfUbik, !!sArgs.ctxtpos_flag,
                                                               !!sArgs.ctxtneg_flag, !!sArgs.bridgepos_flag,
                                                               !!sArgs.bridgeneg_flag, !!sArgs.outpos_flag,
                                                               !!sArgs.outneg_flag, !!sArgs.invert_flag) << endl;

                    }
                }
            }

            if (sArgs.inputs_num)
                ofsm.close();
            else
                cout.flush();

            if (!sArgs.inputs_num)
                break;
        }
    }

    return 0;
}

struct SSorterMod {
    const vector<float> &m_vecdValues;

    explicit SSorterMod(const vector<float> &vecdValues) : m_vecdValues(vecdValues) {}

    bool operator()(size_t iOne, size_t iTwo) {

        return (m_vecdValues[iTwo] < m_vecdValues[iOne]);
    }
};

double AUCMod(const CDat &DatData, const CDat &DatAnswers, const vector<bool> &vecfHere, const vector<bool> &vecfUbik,
              const gengetopt_args_info &sArgs, bool fInvert,
              float dAUC) {
    size_t i, j, iPos, iNeg, iPosCur, iNegCur, iOne, iTwo, iIndex, iAUC;
    vector<float> vecdValues, vecdAnswers;
    vector <size_t> veciGenes, veciIndices;
    bool fAnswer;
    double dRet;
    float d, dAnswer, dSens, dSpec, dSensPrev, dSpecPrev;

    veciGenes.resize(DatAnswers.GetGenes());
    for (i = 0; i < veciGenes.size(); ++i)
        veciGenes[i] = DatData.GetGene(DatAnswers.GetGene(i));

    for (iPos = iNeg = i = 0; i < DatAnswers.GetGenes(); ++i) {
        if ((iOne = veciGenes[i]) == -1)
            continue;
        for (j = (i + 1); j < DatAnswers.GetGenes(); ++j) {
            if (((iTwo = veciGenes[j]) == -1) ||
                CMeta::IsNaN(dAnswer = DatAnswers.Get(i, j)) ||
                CMeta::IsNaN(d = DatData.Get(iOne, iTwo)))
                continue;
            fAnswer = dAnswer > 0;

            if (CMeta::SkipEdge(fAnswer, i, j, vecfHere, vecfUbik, sArgs.ctxtpos_flag, sArgs.ctxtneg_flag,
                                sArgs.bridgepos_flag, sArgs.bridgeneg_flag, sArgs.outpos_flag, sArgs.outneg_flag)) {
                continue;
            }

            if (fAnswer)
                iPos++;
            else
                iNeg++;
            if (fInvert)
                d = 1 - d;
            vecdAnswers.push_back(dAnswer);
            vecdValues.push_back(d);
        }
    }

    veciIndices.resize(vecdValues.size());
    for (i = 0; i < veciIndices.size(); ++i)
        veciIndices[i] = i;
    sort(veciIndices.begin(), veciIndices.end(), SSorterMod(vecdValues));

    iAUC = (size_t)((dAUC < 1) ? (dAUC * iNeg) : dAUC);
    dRet = dSensPrev = dSpecPrev = 0;
    for (iPosCur = iNegCur = i = 0; (i < veciIndices.size()) && (iNegCur < iAUC); ++i) {
        iIndex = veciIndices[i];
        if (vecdAnswers[iIndex] > 0)
            iPosCur++;
        else
            iNegCur++;
        dSens = (float) iPosCur / iPos;
        dSpec = 1 - (float) (iNeg - iNegCur) / iNeg;
        if (dSpec > dSpecPrev) {
            dRet += (dSpec - dSpecPrev) * dSens;
            dSensPrev = dSens;
            dSpecPrev = dSpec;
        }
    }
    dRet *= max(1.0f, (float) iNeg / iAUC);

    return dRet;
}
