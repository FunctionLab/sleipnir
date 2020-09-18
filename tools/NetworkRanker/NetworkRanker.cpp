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
#include "statistics.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <algorithm>


template<class tType>
struct SCompareRank {
    const vector <tType> &m_vecData;

    SCompareRank(const vector <tType> &vecData) : m_vecData(vecData) {}

    bool operator()(size_t iOne, size_t iTwo) const {
        return (m_vecData[iOne] < m_vecData[iTwo]);
    }
};

template<class tType>
struct SCompare {
    const vector <tType> &m_vecData;

    SCompare(const vector <tType> &vecData) : m_vecData(vecData) {}

    bool operator()(size_t iOne, size_t iTwo) const {
        if (CMeta::IsNaN(m_vecData[iOne]))
            return 0;
        if (CMeta::IsNaN(m_vecData[iTwo]))
            return 1;
        return (m_vecData[iOne] < m_vecData[iTwo]);
    }
};


double WilcoxonRankSum(const vector<float> vecdValues, const vector<float> vecAnswers) {
    std::vector <size_t> veciIndices;
    std::vector<float> vecdRanks;
    size_t iIndex, iCount, iPos, iNeg, i, j;
    double dSum, d;

    veciIndices.resize(vecdValues.size());
    for (i = 0; i < vecdValues.size(); ++i)
        veciIndices[i] = i;
    // Sort ranks by values
    std::sort(veciIndices.begin(), veciIndices.end(), SCompareRank<float>(vecdValues));
    vecdRanks.resize(veciIndices.size());
    for (i = 0; i < vecdRanks.size(); ++i) {
        iIndex = veciIndices[i];
        // Handle ties
        if (!i || (vecdValues[iIndex] != vecdValues[veciIndices[i - 1]])) {
            for (iCount = 0, j = i; j < veciIndices.size(); ++j) {
                if (vecdValues[veciIndices[j]] != vecdValues[iIndex])
                    break;
                iCount++;
            }
            d = i + (iCount - 1) / 2.0f;
        }
        vecdRanks[iIndex] = d;
    }

    for (i = 0, dSum = 0, iPos = 0, iNeg = 0; i < vecAnswers.size(); i++) {
        if (vecAnswers[i] > 0) {
            iPos += 1;
            dSum += vecdRanks[i];
        } else {
            iNeg += 1;
        }
    }

    dSum -= (iPos * (iPos - 1)) / 2;
    return (dSum / iPos / iNeg);
}

float score_correction(float d, size_t iNet, size_t iGene1, size_t iGene2, size_t iNetTotal,
                       CFullMatrix<float> &GeneDeg, CFullMatrix<float> &GeneStd,
                       vector<float> &RefGeneDeg, CFullMatrix<float> &MatSum,
                       CDat &DatRef, CDat &DatSum,
                       CDat &DatStd, CDat &DatMean,
                       CPCL &GeneExpPCL, const size_t pclNet, const size_t pclGene1, const size_t pclGene2,
                       float cutoff, float expVal,
                       bool gnorm, bool nnorm, bool enorm, bool backg, bool refnet, bool zscore) {


    float dRaw = d;

    // Correct for background in network (i) 
    if (gnorm) {
        d /= GeneDeg.Get(iGene1, iNet);
    } else if (nnorm) {
        d /= GeneDeg.Get(iGene2, iNet);
    } else if (enorm) {
        d /= sqrt(GeneDeg.Get(iGene1, iNet) * GeneDeg.Get(iGene2, iNet));
    } else if (zscore) {
        //float cMean = ( GeneDeg.Get( iGene1, iNet ) + GeneDeg.Get( iGene2, iNet ) ) / ( 2 * ( GeneDeg.GetRows() -1 ) );
        //float cStd = sqrt( pow( GeneStd.Get( iGene1, iNet ), 2 ) + pow( GeneStd.Get( iGene2, iNet ), 2 ) );
        //d = ( d - cMean ) / cStd;
        float cMean1 = GeneDeg.Get(iGene1, iNet) / (GeneDeg.GetRows() - 1);
        float cMean2 = GeneDeg.Get(iGene2, iNet) / (GeneDeg.GetRows() - 1);
        d = ((d - cMean1) / GeneStd.Get(iGene1, iNet) + (d - cMean2) / GeneStd.Get(iGene2, iNet)) / 1.4142;
    }


    // Divide by average edge weight across all networks
    if (backg) {
        float fDiv;
        if (refnet) {
            float refd;
            refd = DatRef.Get(iGene1, iGene2);
            if (gnorm)
                fDiv = RefGeneDeg[iGene1];
            else if (enorm)
                fDiv = sqrt(RefGeneDeg[iGene1] * RefGeneDeg[iGene2]);
            else
                fDiv = RefGeneDeg[iGene2];

            refd /= fDiv;

            d -= fDiv;
            if (d < 0) d = 0;
        } else if (gnorm) {
            fDiv = (MatSum.Get(iGene1, iGene2) / iNetTotal);
            if (d) d /= fDiv;
        } else if (nnorm) {
            fDiv = (MatSum.Get(iGene2, iGene1) / iNetTotal);
            if (d) d /= fDiv;
        } else if (zscore) {
            float zMeanBack = (dRaw - DatMean.Get(iGene1, iGene2)) / DatStd.Get(iGene1, iGene2);
            d = (d + zMeanBack) / 1.4142;
        } else {
            fDiv = (DatSum.Get(iGene1, iGene2) / iNetTotal);
            if (d) d /= fDiv;
        }
    }


    //if ( GeneExpPCL.GetExperiments() && GeneExpPCL.GetGenes() ) {
    //size_t iNet = GeneExpPCL.GetGene( netName );
    //size_t iGene1 = GeneExpPCL.GetExperiment( strGene1 );
    //size_t iGene2 = GeneExpPCL.GetExperiment( strGene2 );
    if (pclNet != -1 && pclGene1 != -1 && pclGene2 != -1) {
        float exp1 = GeneExpPCL.Get(pclNet, pclGene1);
        float exp2 = GeneExpPCL.Get(pclNet, pclGene2);
        if (!CMeta::IsNaN(exp1) && !CMeta::IsNaN(exp2) && exp1 < cutoff && exp2 < cutoff) {
            if (expVal == -1)
                expVal = CMeta::GetNaN();
            d = expVal;
        }

        //cout << pclNet << ", " << pclGene1 << ", " << pclGene2 << endl;
    }
    //}


    return d;
}

int main(int iArgs, char **aszArgs) {

    gengetopt_args_info sArgs{};
    int iRet = 0;
    size_t i, j, k, l, iGene, iExp, iAnswer;
    float d, fScore, fSum, fAnswer, fDeg;
    DIR *dp;
    struct dirent *ep;
    CDat DatMean, DatStd, DatSum, DatCur, DatRef, Data, Answers;

    vector <size_t> veciGenesCur, veciAnnotIdx, veciGenesetIdx;
    vector<int> veciAnswers;
    vector<float> vecdValues, vecfAnswers;
    vector <vector<size_t>> veciIdxBins;
    vector < vector < vector < size_t > > > veciBinsIdx;

    CPCL PCL, AnnotPCL, GenesetsPCL, GeneExpPCL;
    std::vector <string> input_files;
    std::string dab_dir;
    CFullMatrix<float> GeneDeg = CFullMatrix<float>();
    CFullMatrix<float> GeneStd = CFullMatrix<float>();
    vector<float> RefGeneDeg = vector<float>();
    CFullMatrix<float> MatSum = CFullMatrix<float>();
    vector <string> features;
    vector<float> ExpPCLGenes = vector<float>();

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }
    CMeta Meta(sArgs.verbosity_arg);

    // check if directory valid
    if (sArgs.directory_arg) {
        dp = opendir(sArgs.directory_arg);
        if (dp != nullptr) {
            (void) closedir(dp);
            dab_dir = sArgs.directory_arg;
        } else {
            cerr << "Couldn't open the directory: " << sArgs.directory_arg << '\n';
            return 1;
        }
    }

    dp = opendir(dab_dir.c_str());
    if (dp != nullptr) {
        while ((ep = readdir(dp))) {
            // skip . .. files and temp files with ~
            if (ep->d_name[0] == '.' || ep->d_name[strlen(ep->d_name) - 1] == '~')
                continue;
            if (std::string(ep->d_name).substr(strlen(ep->d_name) - 4, 4) == string(".dab"))
                input_files.push_back((string) sArgs.directory_arg + "/" + ep->d_name);
        }
        (void) closedir(dp);
    }

    if (sArgs.genesets_arg) {
        GenesetsPCL.Open(sArgs.genesets_arg, 0);
    }

    if (sArgs.refnet_arg) {
        DatRef.Open(sArgs.refnet_arg);
        RefGeneDeg.resize(DatRef.GetGenes());
        // Calculate gene degrees
        for (j = 0; j < DatRef.GetGenes(); ++j) {
            fSum = 0;
            for (k = 0; k < DatRef.GetGenes(); ++k) {
                if (j == k) continue;
                fSum += DatRef.Get(j, k);
            }
            RefGeneDeg[j] = fSum / (DatRef.GetGenes() - 1);
        }
    }

    if (sArgs.barcode_arg) {
        GeneExpPCL.Open(sArgs.barcode_arg, 0);
    }

    if (sArgs.annot_arg) {
        AnnotPCL.Open(sArgs.annot_arg, 0);
    } else if (sArgs.enorm_flag || sArgs.gnorm_flag || sArgs.nnorm_flag || sArgs.pval_flag || sArgs.zscore_flag) {
        if (sArgs.pval_flag) {
            veciIdxBins.resize(input_files.size());
            veciBinsIdx.resize(input_files.size());
        }

        for (i = 0; i < input_files.size(); i++) {
            vector<float> vecfDegs;
            vector <size_t> veciIndices;

            cout << "calculating degree - opening: " << input_files[i] << endl;

            // open dat/dab network
            if (!DatCur.Open(input_files[i].c_str())) {
                cerr << "Couldn't open: " << input_files[i] << endl;
                continue;
            }

            vecfDegs.resize(DatCur.GetGenes());

            if (GeneDeg.GetRows() == 0) {
                GeneDeg.Initialize(DatCur.GetGenes(), input_files.size());
                GeneStd.Initialize(DatCur.GetGenes(), input_files.size());
            }

            // Calculate gene degrees
            for (j = 0; j < DatCur.GetGenes(); ++j) {
                float fSum = 0;
                vector<float> vecfNeigh = vector<float>();

                vecfNeigh.resize(DatCur.GetGenes() - 1);
                for (k = 0; k < DatCur.GetGenes(); ++k) {
                    if (j == k) continue;
                    fSum += DatCur.Get(j, k);
                    vecfNeigh[k] = DatCur.Get(j, k);
                }
                float fDeg = fSum / (DatCur.GetGenes() - 1);
                GeneDeg.Set(j, i, fDeg);
                GeneStd.Set(j, i, sqrt(CStatistics::Variance(vecfNeigh)));
                vecfDegs[j] = fDeg;
            }

            if (sArgs.pval_flag) {
                veciIndices.resize(vecfDegs.size());
                for (j = 0; j < vecfDegs.size(); ++j)
                    veciIndices[j] = j;
                // Sort ranks by values
                std::sort(veciIndices.begin(), veciIndices.end(), SCompare<float>(vecfDegs));

                veciIdxBins[i].resize(vecfDegs.size());
                veciBinsIdx[i].resize(ceil(vecfDegs.size() / 100.0));

                for (j = 0; j < vecfDegs.size(); j += 1) {
                    size_t bin = j / 100;
                    veciIdxBins[i][veciIndices[j]] = bin;
                    veciBinsIdx[i][bin].push_back(veciIndices[j]);
                }
            }
        }
    }


    if (!sArgs.annot_given) {
        // now iterate dat/dab networks
        for (i = 0; i < input_files.size(); ++i) {
            // open dat/dab network
            if (!DatCur.Open(input_files[i].c_str())) {
                cerr << "Couldn't open: " << input_files[i] << endl;
                return 1;
            }
            cerr << "opened: " << input_files[i] << endl;

            // if open first network, we will just add edge weights to this CDat
            if (i == 0) {
                DatSum.Open(DatCur.GetGeneNames());
                DatMean.Open(DatCur.GetGeneNames());
                DatStd.Open(DatCur.GetGeneNames());
                MatSum.Initialize(DatCur.GetGenes(), DatCur.GetGenes());
            }

            // now add edges to Dat
            for (j = 0; j < DatSum.GetGenes(); ++j) {
                for (k = (j + 1); k < DatSum.GetGenes(); ++k) {
                    if (CMeta::IsNaN(d = DatCur.Get(j, k))) {
                        cerr << "ERROR: missing values" << input_files[i] << endl;
                        return 1;
                    }

                    if (sArgs.gnorm_flag || sArgs.nnorm_flag) {
                        MatSum.Set(j, k, (i == 0 ? 0 : MatSum.Get(j, k)) + d / GeneDeg.Get(j, i));
                        MatSum.Set(k, j, (i == 0 ? 0 : MatSum.Get(k, j)) + d / GeneDeg.Get(k, i));
                        continue;
                    } else if (sArgs.enorm_flag) {
                        d /= sqrt(GeneDeg.Get(j, i) * GeneDeg.Get(k, i));
                    }

                    DatSum.Set(j, k, (i == 0 ? 0 : DatSum.Get(j, k)) + d);

                    float fMean = DatMean.Get(j, k);
                    float fVar = DatStd.Get(j, k);
                    if (CMeta::IsNaN(fMean)) fMean = 0;
                    if (CMeta::IsNaN(fVar)) fVar = 0;

                    float fDelta = d - fMean;
                    fMean = fMean + fDelta / (i + 1);
                    fVar = fVar + fDelta * (d - fMean);

                    DatMean.Set(j, k, fMean);
                    DatStd.Set(j, k, fVar);

                }
            }
        }

        for (j = 0; j < DatSum.GetGenes(); ++j) {
            for (k = (j + 1); k < DatSum.GetGenes(); ++k) {
                DatStd.Set(j, k, sqrt(DatStd.Get(j, k) / (input_files.size() - 1)));

            }
        }
    }

    for (i = 0; i < input_files.size(); ++i) {
        if (!DatCur.Open(input_files[i].c_str())) {
            cerr << "Couldn't open: " << input_files[i] << endl;
            return 1;
        }
        cerr << input_files[i] << endl;
        if (i == 0) {
            ExpPCLGenes.resize(DatCur.GetGenes());
            for (size_t x = 0; x < ExpPCLGenes.size(); x++) {
                ExpPCLGenes[x] = GeneExpPCL.GetExperiment(DatCur.GetGene(x));
            }

            if (sArgs.genesets_arg) {
                PCL.Open(GenesetsPCL.GetExperimentNames(), input_files, features);
            } else
                PCL.Open(DatCur.GetGeneNames(), input_files, features);
        }


        string netName = CMeta::Basename(input_files[i].c_str());
        netName = CMeta::Deextension(netName);

        cerr << netName << endl;

        size_t iNetExp = -1;
        if (GeneExpPCL.GetGenes()) {
            iNetExp = GeneExpPCL.GetGene(netName);
        }

        if (sArgs.annot_arg && !i) {
            veciAnnotIdx.resize(DatCur.GetGenes());
            for (j = 0; j < DatCur.GetGenes(); j++) {
                size_t idx = AnnotPCL.GetGene(DatCur.GetGene(j));
                veciAnnotIdx[j] = idx;
            }
        }
        if (sArgs.genesets_arg && !i) {
            veciGenesetIdx.resize(DatCur.GetGenes());
            for (j = 0; j < DatCur.GetGenes(); j++) {
                size_t idx = GenesetsPCL.GetGene(DatCur.GetGene(j));
                veciGenesetIdx[j] = idx;
            }
        }

        size_t iGenes, iGenesMax;
        iGenesMax = 1;
        if (sArgs.genesets_arg)
            iGenesMax = GenesetsPCL.GetExperiments();
        for (iGenes = sArgs.geneset_idx_arg; iGenes < std::max((size_t) sArgs.geneset_idx_arg, iGenesMax); iGenes++) {
            float fSetScore = 0;
            size_t iSet = 0;
            size_t iGenesetSize = 0;

            for (j = 0; j < DatCur.GetGenes(); ++j) {
                fScore = 0;
                float fGeneScore = 0;
                size_t iGeneScores = 0;

                if (sArgs.genesets_arg &&
                    (veciGenesetIdx[j] == -1 || GenesetsPCL.Get(veciGenesetIdx[j], iGenes) != 1))
                    continue;

                if (sArgs.annot_arg) {

                    iExp = AnnotPCL.GetExperiment(netName);
                    string strName = DatCur.GetGene(j);

                    vecdValues.clear();
                    vecfAnswers.clear();

                    for (k = 0; k < DatCur.GetGenes(); ++k) {
                        if (j == k) continue;
                        d = DatCur.Get(j, k);

                        if ((iGene = veciAnnotIdx[k]) == -1)
                            continue;

                        if ((fAnswer = AnnotPCL.Get(iGene, iExp)) == 0)
                            continue;

                        vecfAnswers.push_back(fAnswer);
                        vecdValues.push_back(d);
                    }

                    fScore = WilcoxonRankSum(vecdValues, vecfAnswers);
                    PCL.Set(j, i, fScore);

                    continue;
                }

                for (k = 0; k < DatCur.GetGenes(); ++k) {
                    if (j == k) continue;
                    if (sArgs.genesets_arg && k <= i) continue;
                    if (sArgs.genesets_arg &&
                        (veciGenesetIdx[k] == -1 || GenesetsPCL.Get(veciGenesetIdx[k], iGenes) != 1))
                        continue;

                    iGenesetSize++;

                    //cout << iNetExp << ", " << ExpPCLGenes[ j ] << ", " << ExpPCLGenes[ k ] << endl;

                    d = DatCur.Get(j, k);
                    d = score_correction(d, i, j, k, input_files.size(),
                                         GeneDeg, GeneStd, RefGeneDeg, MatSum, DatRef, DatSum,
                                         DatMean, DatStd,
                                         GeneExpPCL, iNetExp, ExpPCLGenes[j], ExpPCLGenes[k], sArgs.exp_cut_arg,
                                         sArgs.no_exp_arg,
                                         sArgs.gnorm_flag, sArgs.nnorm_flag, sArgs.enorm_flag,
                                         sArgs.backg_flag, sArgs.refnet_arg, sArgs.zscore_flag);

                    if (CMeta::IsNaN(d)) { continue; }

                    if (sArgs.genesets_arg) {
                        fSetScore += d;
                        iSet += 1;
                        fGeneScore += d;
                        iGeneScores += 1;
                    }

                        // Use connectivity score to weight edge, or use score itself
                    else if (sArgs.log_weight_flag) {
                        fScore += DatCur.Get(j, k) * log2(d);
                    } else if (sArgs.weight_flag) {
                        fScore += DatCur.Get(j, k) * d;
                    } else {
                        fScore += d;
                    }
                }

                if (!sArgs.genesets_arg) {
                    fScore /= (DatCur.GetGenes() - 1);
                    PCL.Set(j, i, fScore);
                }
            }

            fSetScore /= (float) iSet;

            if (sArgs.genesets_arg && sArgs.pval_flag && iGenesetSize) {
                float fRelRank = 0;
                vector<float> vecfRandVals;

                for (size_t p = 0; p < 100; p++) {

                    vector <size_t> ivecRand;
                    for (j = 0; j < DatCur.GetGenes(); ++j) {
                        if (veciGenesetIdx[j] != -1 &&
                            GenesetsPCL.Get(veciGenesetIdx[j], iGenes) == 1) {
                            size_t b = veciIdxBins[i][j];
                            size_t idx = veciBinsIdx[i][b][rand() % veciBinsIdx[i][b].size()];
                            ivecRand.push_back(idx);
                        }
                    }

                    float fSetScoreRand = 0;
                    size_t iSetRand = 0;

                    for (j = 0; j < ivecRand.size(); ++j) {
                        for (k = j + 1; k < ivecRand.size(); ++k) {

                            d = DatCur.Get(ivecRand[j], ivecRand[k]);
                            d = score_correction(d, i, ivecRand[j], ivecRand[k], input_files.size(),
                                                 GeneDeg, GeneStd, RefGeneDeg, MatSum, DatRef, DatSum,
                                                 DatMean, DatStd,
                                                 GeneExpPCL, iNetExp, ExpPCLGenes[j], ExpPCLGenes[k], sArgs.exp_cut_arg,
                                                 sArgs.no_exp_arg,
                                                 sArgs.gnorm_flag, sArgs.nnorm_flag, sArgs.enorm_flag,
                                                 sArgs.backg_flag, sArgs.refnet_arg, sArgs.zscore_flag);

                            if (CMeta::IsNaN(d)) continue;

                            fSetScoreRand += d;
                            iSetRand += 1;
                        }
                    }
                    fSetScoreRand /= (float) iSetRand;

                    if (fSetScoreRand < fSetScore)
                        fRelRank += .01;//1.0/100.0;

                    vecfRandVals.push_back(fSetScoreRand);
                }

                //if ( sArgs.zscore_flag ) {
                //    float std = sqrt( CStatistics::Variance( vecfRandVals ) );
                //    float mean = CStatistics::Average( vecfRandVals );
                //    fSetScore = (fSetScore - mean) / std;
                //}
                //else {
                fSetScore = fRelRank;
                //}
            }

            if (sArgs.genesets_arg)
                PCL.Set(iGenes, i, fSetScore);
        }
    }

    PCL.Save(sArgs.output_arg);

    return iRet;
}
