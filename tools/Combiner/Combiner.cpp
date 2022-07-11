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
#include <random>
static int MainDABs(const gengetopt_args_info &);

static int MainDATs(const gengetopt_args_info &);

static int MainRevDATs(const gengetopt_args_info &);

static int MainPCLs(const gengetopt_args_info &);

static int MainModules(const gengetopt_args_info &);

static const TPFnCombiner c_apfnCombiners[] = {MainPCLs, MainDATs, MainDABs, MainModules, MainRevDATs, NULL};
static const char *c_aszCombiners[] = {"pcl", "dat", "dad", "module", "revdat", NULL};
static const char c_acDab[] = ".dab";
static const char c_acDat[] = ".dat";
static const char c_acQDab[] = ".qdab";

// store all input file names
vector <string> input_files;

enum EMethod {
    EMethodBegin = 0,
    EMethodMean = EMethodBegin,
    EMethodSum = EMethodMean + 1,
    EMethodGMean = EMethodSum + 1,
    EMethodHMean = EMethodGMean + 1,
    EMethodMax = EMethodHMean + 1,
    EMethodMin = EMethodMax + 1,
    EMethodDiff = EMethodMin + 1,
    EMethodMeta = EMethodDiff + 1,
    EMethodQMeta = EMethodMeta + 1,
    EMethodProd = EMethodQMeta + 1,
    EMethodEnd = EMethodProd + 1
};

static const char *c_aszMethods[] = {
        "mean", "sum", "gmean", "hmean", "max", "min", "diff", "meta", "qmeta", "prod", nullptr
};

int main(int iArgs, char **aszArgs) {
    gengetopt_args_info sArgs{};
    int iRet = 0;
    size_t i;
    DIR *dp;
    struct dirent *ep;

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }
    CMeta Meta(sArgs.verbosity_arg);

    // now collect the data files from directory if given
    if (sArgs.directory_arg) {
        dp = opendir(sArgs.directory_arg);
        if (dp != nullptr) {
            while ((ep = readdir(dp))) {
                // skip . .. files and temp files with ~
                if (ep->d_name[0] == '.' || ep->d_name[strlen(ep->d_name) - 1] == '~')
                    continue;

                // currently opens all files. Add filter here if want pick file extensions
                input_files.push_back((string) sArgs.directory_arg + "/" + ep->d_name);
            }
            (void) closedir(dp);
        } else {
            cerr << "Couldn't open the directory: " << sArgs.directory_arg << '\n';
            return 1;
        }
    } else {
        input_files.resize(sArgs.inputs_num);
        copy(sArgs.inputs, sArgs.inputs + sArgs.inputs_num, input_files.begin());
    }

    for (i = 0; c_aszCombiners[i]; ++i)
        if (!strcmp(c_aszCombiners[i], sArgs.type_arg)) {
            iRet = c_apfnCombiners[i](sArgs);
            break;
        }

    return iRet;
}

int MainPCLs(const gengetopt_args_info &sArgs) {
    CPCL PCL, PCLNew;
    size_t i, j, iArg, iExp, iGene;
    vector <string> vecstrGenes, vecstrExps, vecstrFeatures;
    set <string> setstrGenes;
    set<string>::const_iterator iterGenes;
    ifstream ifsm;
    ofstream ofsm;
    float d;

    for (iArg = 0; iArg < input_files.size(); ++iArg) {
        if (!PCL.Open(input_files[iArg].c_str(), sArgs.skip_arg, !!sArgs.memmap_flag)) {
            cerr << "Could not open: " << input_files[iArg] << endl;
            return 1;
        }
        if (!iArg)
            for (i = 0; i < PCL.GetFeatures(); ++i)
                vecstrFeatures.push_back(PCL.GetFeature(i));
        for (i = 0; i < PCL.GetExperiments(); ++i)
            vecstrExps.push_back(PCL.GetExperiment(i));
        for (i = 0; i < PCL.GetGenes(); ++i)
            setstrGenes.insert(PCL.GetGene(i));
        ifsm.close();
    }
    vecstrGenes.resize(setstrGenes.size());
    copy(setstrGenes.begin(), setstrGenes.end(), vecstrGenes.begin());

    PCLNew.Open(vecstrGenes, vecstrExps, vecstrFeatures);
    iExp = 0;
    for (iArg = 0; iArg < input_files.size(); ++iArg) {
        cerr << "Processing " << input_files[iArg] << "..." << endl;
        ifsm.clear();
        ifsm.open(input_files[iArg].c_str());
        PCL.Open(input_files[iArg].c_str(), sArgs.skip_arg, !!sArgs.memmap_flag);
        if (sArgs.normalize_flag)
            PCL.Normalize(CPCL::ENormalizeRow);

        for (i = 0; i < PCLNew.GetGenes(); ++i)
            if ((iGene = PCL.GetGene(vecstrGenes[i])) != -1) {
                if (!iArg)
                    for (j = 1; j < PCLNew.GetFeatures(); ++j)
                        PCLNew.SetFeature(i, j, PCL.GetFeature(iGene, j));
                for (j = 0; j < PCL.GetExperiments(); ++j)
                    PCLNew.Set(i, iExp + j, PCL.Get(iGene, j));
            } else if (sArgs.zero_flag) {
                for (j = 0; j < PCL.GetExperiments(); ++j)
                    PCLNew.Set(i, iExp + j, 0.0f);
            }

        iExp += PCL.GetExperiments();
        ifsm.close();
    }

    if (sArgs.output_arg) {
        ofsm.open(sArgs.output_arg);
        PCLNew.Save(ofsm);
        ofsm.close();
    } else {
        PCLNew.Save(cout);
        cout.flush();
    }

    return 0;
}

struct SCallbackMeta {
    CDat m_DatWei;
    CDat m_DatYBare;
    CDat m_DatQe;
    CDat m_DatDeltaeD;
};

struct SCallbackVars {
    size_t m_iDataset;
    size_t m_iOne;
    size_t m_iTwo;
    float m_dValue;
    float m_dWeight;
    EMethod m_eMethod;
    size_t m_iDatasets;
    CDat *m_pDatOut;
    CDat *m_pDatCur;
    CHalfMatrix<float> *m_pMatCounts;
    const vector <set<size_t>> *m_pvecsetiGenes;
    const vector <vector<size_t>> *m_pvecveciGenes;
    SCallbackMeta m_sCallbackMeta;
};

struct SCallback {
    const gengetopt_args_info &m_sArgs;
    const CGenes &m_GenesIn;
    const CPCL &m_PCLWeights;
    const vector <string> &m_vecstrTerms;
    const vector<CGenes *> &m_vecpTerms;
    SCallbackVars m_sCallbackVars;

    SCallback(const gengetopt_args_info &sArgs, const CGenes &GenesIn, const CPCL &PCLWeights,
              const vector <string> &vecstrTerms, const vector<CGenes *> &vecpTerms) :
            m_sArgs(sArgs), m_GenesIn(GenesIn), m_PCLWeights(PCLWeights), m_vecstrTerms(vecstrTerms),
            m_vecpTerms(vecpTerms) {}
};

int iterate_inputs(SCallback &sCallback, void (*pfnCallback)(SCallbackVars &),
                   void (*pfnInitialize)(SCallbackVars &) = NULL) {
    SCallbackVars &sCallbackVars = sCallback.m_sCallbackVars;
    CDat &DatOut = *sCallbackVars.m_pDatOut;
    const gengetopt_args_info &sArgs = sCallback.m_sArgs;
    const CGenes &GenesIn = sCallback.m_GenesIn;
    const CPCL &PCLWeights = sCallback.m_PCLWeights;
    const vector <string> &vecstrTerms = sCallback.m_vecstrTerms;
    const vector<CGenes *> &vecpTerms = sCallback.m_vecpTerms;
    vector <set<size_t>> vecsetiGenes;
    vector <vector<size_t>> vecveciGenes;
    size_t i, j, k, iOne, iTwo, iA, iB;
    CDat DatCur;
    float d;

    sCallbackVars.m_iDatasets = input_files.size();
    sCallbackVars.m_pvecveciGenes = &vecveciGenes;
    sCallbackVars.m_pvecsetiGenes = &vecsetiGenes;
    vecsetiGenes.resize(DatOut.GetGenes());
    for (i = 0; i < input_files.size(); ++i) {
        if (!DatCur.Open(input_files[i].c_str(), !!sArgs.memmap_flag && !sArgs.normalize_flag && !GenesIn.GetGenes())) {
            cerr << "Couldn't open: " << input_files[i] << endl;
            return 1;
        }
        sCallbackVars.m_iDataset = i;
        sCallbackVars.m_pDatCur = &DatCur;
        if (PCLWeights.GetGenes()) {
            if ((j = PCLWeights.GetGene(CMeta::Deextension(CMeta::Basename(input_files[i].c_str())))) == -1) {
                cerr << "Ignoring unweighted graph: " << input_files[i] << endl;
                continue;
            }
            sCallbackVars.m_dWeight = PCLWeights.Get(j, 0);
        } else
            sCallbackVars.m_dWeight = 1;
        cerr << "Opened: " << input_files[i] << endl;
        if (sArgs.zero_flag) {
            vector <string> vecstrGenes;

            for (j = 0; j < DatOut.GetGenes(); ++j)
                if (DatCur.GetGene(DatOut.GetGene(j)) == -1)
                    vecstrGenes.push_back(DatOut.GetGene(j));
            DatCur.AddGenes(vecstrGenes);
            for (j = 0; j < DatCur.GetGenes(); ++j)
                for (k = (j + 1); k < DatCur.GetGenes(); ++k)
                    if (CMeta::IsNaN(DatCur.Get(j, k)))
                        DatCur.Set(j, k, 0);
        }
        if (sArgs.normalize_flag)
            DatCur.Normalize(CDat::ENormalizeZScore);
        if (sArgs.quantiles_arg > 0)
            DatCur.NormalizeQuantiles(sArgs.quantiles_arg);
        if (GenesIn.GetGenes())
            DatCur.FilterGenes(GenesIn, CDat::EFilterInclude);
        vecveciGenes.resize(DatCur.GetGenes());
        for (j = 0; j < vecveciGenes.size(); ++j)
            vecveciGenes[j].clear();
        for (j = 0; j < DatOut.GetGenes(); ++j) {
            vecsetiGenes[j].clear();
            if (vecstrTerms.empty()) {
                if ((iOne = DatCur.GetGene(DatOut.GetGene(j))) != -1)
                    vecveciGenes[iOne].push_back(j);
            } else
                for (k = 0; k < vecpTerms[j]->GetGenes(); ++k)
                    if ((iOne = DatCur.GetGene(vecpTerms[j]->GetGene(k).GetName())) != -1) {
                        vecveciGenes[iOne].push_back(j);
                        vecsetiGenes[j].insert(iOne);
                    }
        }
        if (pfnInitialize)
            pfnInitialize(sCallbackVars);
        for (j = 0; j < DatCur.GetGenes(); ++j)
            for (k = (j + 1); k < DatCur.GetGenes(); ++k) {
                if (CMeta::IsNaN(d = DatCur.Get(j, k)))
                    continue;
                sCallbackVars.m_dValue = d;
                for (iA = 0; iA < vecveciGenes[j].size(); ++iA) {
                    iOne = vecveciGenes[j][iA];
                    if (vecsetiGenes[iOne].find(k) != vecsetiGenes[iOne].end())
                        continue;
                    sCallbackVars.m_iOne = iOne;
                    for (iB = 0; iB < vecveciGenes[k].size(); ++iB) {
                        iTwo = vecveciGenes[k][iB];
                        if (vecsetiGenes[iTwo].find(j) != vecsetiGenes[iTwo].end())
                            continue;
                        sCallbackVars.m_iTwo = iTwo;
                        pfnCallback(sCallbackVars);
                    }
                }
            }
    }

    return 0;
}

float weight_weistar(const SCallbackVars &sCallback) {
    const SCallbackMeta &sCallbackMeta = sCallback.m_sCallbackMeta;
    float dDelta2;

    if ((dDelta2 = sCallbackMeta.m_DatDeltaeD.Get(sCallback.m_iOne, sCallback.m_iTwo))) {
        if ((dDelta2 = (sCallbackMeta.m_DatQe.Get(sCallback.m_iOne, sCallback.m_iTwo) -
                        sCallback.m_iDatasets + 1) / dDelta2) < 0)
            dDelta2 = 0;
    }
/*
if( ( (float)rand( ) / RAND_MAX ) < 0.000001 )
cerr << ( 1 / ( ( 1 / sCallbackMeta.m_DatWei.Get( sCallback.m_iOne, sCallback.m_iTwo ) ) + dDelta2 ) ) << '\t' <<
sCallbackMeta.m_DatWei.Get( sCallback.m_iOne, sCallback.m_iTwo ) << '\t' << dDelta2 << endl;
//*/
    return (1 / ((1 / sCallbackMeta.m_DatWei.Get(sCallback.m_iOne, sCallback.m_iTwo)) + dDelta2));
}

void callback_andersondarling(SCallbackVars &sCallback) {
    static const float c_dFrac = 0.001f;
    vector<float> vecdValues;
    size_t i, j, iGenes;
    float d;

    std::default_random_engine generator;
    std::uniform_real_distribution<float> dist(0, (float) RAND_MAX);

    iGenes = sCallback.m_pDatCur->GetGenes();
    vecdValues.reserve((size_t)(iGenes * iGenes * c_dFrac));
    for (i = 0; i < iGenes; ++i)
        for (j = (i + 1); j < iGenes; ++j)
            if (!CMeta::IsNaN(d = sCallback.m_pDatCur->Get(i, j)) &&
                (dist(generator) < c_dFrac))
                vecdValues.push_back(d);
    sCallback.m_dWeight = (float) (log(CStatistics::AndersonDarlingScore<float>(
            vecdValues.begin(), vecdValues.end())) / log(2.0));
}

void callback_combine(SCallbackVars &sCallback) {

    if (sCallback.m_eMethod == EMethodMeta)
        sCallback.m_dWeight = weight_weistar(sCallback);
    switch (sCallback.m_eMethod) {
        case EMethodGMean:
            sCallback.m_pDatOut->Get(sCallback.m_iOne, sCallback.m_iTwo) *= pow(sCallback.m_dValue,
                                                                                sCallback.m_dWeight);
            sCallback.m_pMatCounts->Get(sCallback.m_iOne, sCallback.m_iTwo) += sCallback.m_dWeight;
            break;

        case EMethodHMean:
            sCallback.m_pDatOut->Get(sCallback.m_iOne, sCallback.m_iTwo) += sCallback.m_dWeight / sCallback.m_dValue;
            sCallback.m_pMatCounts->Get(sCallback.m_iOne, sCallback.m_iTwo) += sCallback.m_dWeight;
            break;

        case EMethodMax:
            if (sCallback.m_dValue > sCallback.m_pDatOut->Get(sCallback.m_iOne, sCallback.m_iTwo))
                sCallback.m_pDatOut->Set(sCallback.m_iOne, sCallback.m_iTwo, sCallback.m_dValue);
            break;

        case EMethodMin:
            if (sCallback.m_dValue < sCallback.m_pDatOut->Get(sCallback.m_iOne, sCallback.m_iTwo))
                sCallback.m_pDatOut->Set(sCallback.m_iOne, sCallback.m_iTwo, sCallback.m_dValue);
            break;

        case EMethodProd:
            sCallback.m_pDatOut->Get(sCallback.m_iOne, sCallback.m_iTwo) *= pow(sCallback.m_dValue,
                                                                                sCallback.m_dWeight);
            sCallback.m_pMatCounts->Get(sCallback.m_iOne, sCallback.m_iTwo) += sCallback.m_dWeight;
            break;

        default:
            sCallback.m_pDatOut->Get(sCallback.m_iOne, sCallback.m_iTwo) += sCallback.m_dWeight * sCallback.m_dValue *
                                                                            ((sCallback.m_iDataset &&
                                                                              (sCallback.m_eMethod == EMethodDiff)) ? -1
                                                                                                                    : 1);
            sCallback.m_pMatCounts->Get(sCallback.m_iOne, sCallback.m_iTwo) += sCallback.m_dWeight;
    }
}

void callback_wei(SCallbackVars &sCallback) {
    static const float c_dFrac = 0.001f;
    SCallbackMeta &sCallbackMeta = sCallback.m_sCallbackMeta;
    const vector <vector<size_t>> &vecveciGenes = *sCallback.m_pvecveciGenes;
    const vector <set<size_t>> &vecsetiGenes = *sCallback.m_pvecsetiGenes;
    size_t i, j, iOne, iTwo, iA, iB, iGenes;
    float d, dSum, dSumSqs;
    vector<float> vecdSums, vecdSumSqs, vecdValues;

    iGenes = sCallback.m_pDatCur->GetGenes();
    vecdValues.reserve((size_t)(iGenes * (iGenes - 1) * c_dFrac));
    vecdSums.resize(sCallback.m_pDatOut->GetGenes());
    fill(vecdSums.begin(), vecdSums.end(), 0.0f);
    vecdSumSqs.resize(sCallback.m_pDatOut->GetGenes());
    fill(vecdSumSqs.begin(), vecdSumSqs.end(), 0.0f);
    dSum = dSumSqs = 0;
    for (i = 0; i < iGenes; ++i)
        for (j = (i + 1); j < iGenes; ++j) {
            if (CMeta::IsNaN(d = sCallback.m_pDatCur->Get(i, j)))
                continue;
            for (iA = 0; iA < vecveciGenes[i].size(); ++iA) {
                iOne = vecveciGenes[i][iA];
                if (vecsetiGenes[iOne].find(j) != vecsetiGenes[iOne].end())
                    continue;
                for (iB = 0; iB < vecveciGenes[j].size(); ++iB) {
                    iTwo = vecveciGenes[j][iB];
                    if (vecsetiGenes[iTwo].find(i) != vecsetiGenes[iTwo].end())
                        continue;
                    if (!(iA || iB) && (((float) rand() / (float) RAND_MAX) < c_dFrac))
                        vecdValues.push_back(d);
                    dSum += d;
                    vecdSums[iOne] += d;
                    vecdSums[iTwo] += d;
                    d *= d;
                    dSumSqs += d;
                    vecdSumSqs[iOne] += d;
                    vecdSumSqs[iTwo] += d;
                }
            }
        }
    i = iGenes * (iGenes - 1) / 2;
    dSum /= i;
    dSumSqs = (dSumSqs / (i - 1)) - (dSum * dSum);
    for (i = 0; i < vecdSums.size(); ++i) {
        d = (vecdSums[i] /= iGenes - 1);
        vecdSumSqs[i] = (vecdSumSqs[i] / (iGenes - 2)) - (d * d);
    }

    d = (float) (log(CStatistics::AndersonDarlingScore<float>(vecdValues.begin(), vecdValues.end())) /
                 log(2.0));
    sCallbackMeta.m_DatWei.Open(sCallback.m_pDatOut->GetGeneNames());
    sCallbackMeta.m_DatWei.Clear(0);
    for (i = 0; i < vecdSumSqs.size(); ++i)
        for (j = (i + 1); j < vecdSumSqs.size(); ++j)
            sCallbackMeta.m_DatWei.Set(i, j,
//				2 / ( vecdSumSqs[i] + vecdSumSqs[j] )			// pooled variance of adjacent genes
//				3 / ( vecdSumSqs[i] + vecdSumSqs[j] + dSumSqs )	// unweighted pool of adjanced genes + whole network
                                       3 * d / (vecdSumSqs[i] + vecdSumSqs[j] + dSumSqs)
            );
}

void callback_deltaed(SCallbackVars &sCallback) {
    SCallbackMeta &sCallbackMeta = sCallback.m_sCallbackMeta;
    float d;

    d = sCallbackMeta.m_DatWei.Get(sCallback.m_iOne, sCallback.m_iTwo);
    sCallbackMeta.m_DatDeltaeD.Get(sCallback.m_iOne, sCallback.m_iTwo) += d;
    sCallback.m_pDatOut->Get(sCallback.m_iOne, sCallback.m_iTwo) += d * d;
}

void callback_ybare(SCallbackVars &sCallback) {
    SCallbackMeta &sCallbackMeta = sCallback.m_sCallbackMeta;

    sCallbackMeta.m_DatYBare.Get(sCallback.m_iOne, sCallback.m_iTwo) +=
            sCallbackMeta.m_DatWei.Get(sCallback.m_iOne, sCallback.m_iTwo) *
            sCallback.m_dValue;
    sCallback.m_pMatCounts->Get(sCallback.m_iOne, sCallback.m_iTwo) += sCallbackMeta.m_DatWei.Get(sCallback.m_iOne,
                                                                                                  sCallback.m_iTwo);
}

void callback_qe(SCallbackVars &sCallback) {
    SCallbackMeta &sCallbackMeta = sCallback.m_sCallbackMeta;
    float d;

    d = sCallback.m_dValue - sCallbackMeta.m_DatYBare.Get(sCallback.m_iOne, sCallback.m_iTwo);
    sCallbackMeta.m_DatQe.Get(sCallback.m_iOne, sCallback.m_iTwo) +=
            sCallbackMeta.m_DatWei.Get(sCallback.m_iOne, sCallback.m_iTwo) * d * d;
}

void initialize_meta(SCallback &sCallback) {
    SCallbackMeta &sCallbackMeta = sCallback.m_sCallbackVars.m_sCallbackMeta;
    const vector <string> &vecstrGenes = sCallback.m_sCallbackVars.m_pDatOut->GetGeneNames();
    size_t i, j;
    CDat *pDatOut;
    CHalfMatrix<float> *pMatCounts;
    CHalfMatrix<float> MatCounts;
    float d;
    CDat DatTmp;

    sCallbackMeta.m_DatDeltaeD.Open(vecstrGenes);
    sCallbackMeta.m_DatDeltaeD.Clear(0);
    DatTmp.Open(vecstrGenes, false);
    DatTmp.Clear(0);
    pDatOut = sCallback.m_sCallbackVars.m_pDatOut;
    sCallback.m_sCallbackVars.m_pDatOut = &DatTmp;
    iterate_inputs(sCallback, callback_deltaed, callback_wei);
    sCallback.m_sCallbackVars.m_pDatOut = pDatOut;
    for (i = 0; i < sCallbackMeta.m_DatDeltaeD.GetGenes(); ++i)
        for (j = (i + 1); j < sCallbackMeta.m_DatDeltaeD.GetGenes(); ++j)
            if ((d = sCallbackMeta.m_DatDeltaeD.Get(i, j)))
                sCallbackMeta.m_DatDeltaeD.Get(i, j) -= DatTmp.Get(i, j) / d;

    sCallbackMeta.m_DatYBare.Open(vecstrGenes, false);
    sCallbackMeta.m_DatYBare.Clear(0);
    MatCounts.Initialize(vecstrGenes.size());
    MatCounts.Clear();
    pMatCounts = sCallback.m_sCallbackVars.m_pMatCounts;
    sCallback.m_sCallbackVars.m_pMatCounts = &MatCounts;
    iterate_inputs(sCallback, callback_ybare, callback_wei);
    sCallback.m_sCallbackVars.m_pMatCounts = pMatCounts;
    for (i = 0; i < sCallbackMeta.m_DatYBare.GetGenes(); ++i)
        for (j = (i + 1); j < sCallbackMeta.m_DatYBare.GetGenes(); ++j)
            if ((d = MatCounts.Get(i, j)))
                sCallbackMeta.m_DatYBare.Get(i, j) /= d;

    sCallbackMeta.m_DatQe.Open(vecstrGenes, false);
    sCallbackMeta.m_DatQe.Clear(0);
    iterate_inputs(sCallback, callback_qe, callback_wei);
}

int MainDATs(const gengetopt_args_info &sArgs) {
    CDataset Dataset;
    CDat DatOut, DatCur, DatQW;
    CHalfMatrix<float> MatCounts;
    size_t i, j;
    float d;
    vector <string> vecstrTerms;
    CPCL PCLWeights(false);
    CGenome Genome;
    CGenes GenesIn(Genome);
    vector < CGenes * > vecpTerms;
    EMethod eMethod;
    vector<float> vecdWis;
    int iRet;
    SCallback sCallback(sArgs, GenesIn, PCLWeights, vecstrTerms, vecpTerms);
    void (*pfnInitialize)(SCallbackVars &);

    if (!input_files.size())
        return 1;

    if (!Dataset.OpenGenes(input_files)) {
        cerr << "Couldn't open: " << input_files[0];
        for (i = 1; i < input_files.size(); ++i)
            cerr << ", " << input_files[i];
        cerr << endl;
        return 1;
    }
    if (sArgs.weights_arg && !PCLWeights.Open(sArgs.weights_arg, 0)) {
        cerr << "Could not open: " << sArgs.weights_arg << endl;
        return 1;
    }
    if (sArgs.genes_arg && !GenesIn.Open(sArgs.genes_arg)) {
        cerr << "Could not open: " << sArgs.genes_arg << endl;
        return 1;
    }
    if (sArgs.terms_arg && !CGenes::Open(sArgs.terms_arg, Genome, vecstrTerms, vecpTerms)) {
        cerr << "Could not open: " << sArgs.terms_arg << endl;
        return 1;
    }

    for (eMethod = EMethodBegin; eMethod < EMethodEnd; eMethod = (EMethod) (eMethod + 1))
        if (!strcmp(c_aszMethods[eMethod], sArgs.method_arg))
            break;
    if (eMethod >= EMethodEnd) {
        cmdline_parser_print_help();
        return 1;
    }

    DatOut.Open(vecstrTerms.empty() ? Dataset.GetGeneNames() : vecstrTerms, false,
                sArgs.memmap_flag ? sArgs.output_arg : NULL);
    switch (eMethod) {
        case EMethodMax:
            d = -FLT_MAX;
            break;

        case EMethodMin:
            d = FLT_MAX;
            break;

        case EMethodGMean:
            d = 1;
            break;

        case EMethodProd:
            d = 1;
            break;

        default:
            d = 0;
    }
    for (i = 0; i < DatOut.GetGenes(); ++i)
        for (j = (i + 1); j < DatOut.GetGenes(); ++j)
            DatOut.Set(i, j, d);
    if (fabs(d) < 2) {
        MatCounts.Initialize(DatOut.GetGenes());
        MatCounts.Clear();
    }

    sCallback.m_sCallbackVars.m_eMethod = eMethod;
    sCallback.m_sCallbackVars.m_pDatOut = &DatOut;
    sCallback.m_sCallbackVars.m_pMatCounts = &MatCounts;
    if (eMethod == EMethodMeta)
        initialize_meta(sCallback);

    switch (eMethod) {
        case EMethodMeta:
            pfnInitialize = callback_wei;
            break;

        case EMethodQMeta:
            pfnInitialize = callback_andersondarling;
            break;

        default:
            pfnInitialize = NULL;
    }
    if ((iRet = iterate_inputs(sCallback, callback_combine, pfnInitialize)))
        return iRet;
    for (i = 0; i < DatOut.GetGenes(); ++i)
        for (j = (i + 1); j < DatOut.GetGenes(); ++j)
            switch (eMethod) {
                case EMethodMean:
                case EMethodMeta:
                case EMethodQMeta:
                    DatOut.Set(i, j, (d = MatCounts.Get(i, j)) ? (DatOut.Get(i, j) / (sArgs.reweight_flag ? 1 : d)) :
                                     CMeta::GetNaN());
                    break;

                case EMethodGMean:
                    DatOut.Set(i, j, (d = MatCounts.Get(i, j)) ?
                                     (float) pow((double) DatOut.Get(i, j), 1.0 / (sArgs.reweight_flag ? 1 : d))
                                                               : CMeta::GetNaN());
                    break;

                case EMethodHMean:
                    DatOut.Set(i, j, (d = MatCounts.Get(i, j)) ? ((sArgs.reweight_flag ? 1 : d) / DatOut.Get(i, j)) :
                                     CMeta::GetNaN());
                    break;

                case EMethodMax:
                    if (DatOut.Get(i, j) == -FLT_MAX)
                        DatOut.Set(i, j, CMeta::GetNaN());
                    break;

                case EMethodMin:
                    if (DatOut.Get(i, j) == FLT_MAX)
                        DatOut.Set(i, j, CMeta::GetNaN());
                    break;

                case EMethodProd:
                    DatOut.Set(i, j, (d = MatCounts.Get(i, j)) ?
                                     (float) (double) DatOut.Get(i, j) : CMeta::GetNaN());
                default:
                    break;
            }

    if (sArgs.zscore_flag)
        DatOut.Normalize(CDat::ENormalizeZScore);
    if (!sArgs.memmap_flag)
        DatOut.Save(sArgs.output_arg);

    for (i = 0; i < vecpTerms.size(); ++i)
        delete vecpTerms[i];
    return 0;
}

static int MainDABs(const gengetopt_args_info &sArgs) {
    CDatasetCompact Dataset;
    size_t i;
    ofstream ofsm;

    if (input_files.empty())
        return 1;

    if (!Dataset.Open(input_files)) {
        cerr << "Couldn't open: " << input_files[0];
        for (i = 1; i < input_files.size(); ++i)
            cerr << ", " << input_files[i];
        cerr << endl;
        return 1;
    }

    if (sArgs.output_arg) {
        ofsm.open(sArgs.output_arg);
        Dataset.Save(ofsm, true);
        ofsm.close();
    } else {
        Dataset.Save(cout, false);
        cout.flush();
    }

    return 0;
}

struct SSimpleome {
    map <string, size_t> m_mapstriGenes;
    map <size_t, string> m_mapistrGenes;

    size_t Get(const string &strGene) {
        map<string, size_t>::const_iterator iterGene;
        size_t iRet;

        if ((iterGene = m_mapstriGenes.find(strGene)) != m_mapstriGenes.end())
            return iterGene->second;

        m_mapstriGenes[strGene] = iRet = m_mapstriGenes.size();
        m_mapistrGenes[iRet] = strGene;
        return iRet;
    }

    const string &Get(size_t iGene) const {

        return m_mapistrGenes.find(iGene)->second;
    }
};

struct SModule {
    float m_dSpecificity{};
    set <size_t> m_setiGenes;
    set<const SModule *> m_setpsChildren;

    static float Open(const char *szFile, SSimpleome &sSimpleome, vector<SModule *> &vecpsModules) {
        static const size_t c_iBuffer = 131072;
        ifstream ifsm;
        char acBuffer[c_iBuffer];
        vector <string> vecstrLine;
        float dRet;
        size_t i, iModule;
        SModule *psModule;

        ifsm.open(szFile);
        if (!ifsm.is_open()) {
            cerr << "Could not open: " << szFile << endl;
            return CMeta::GetNaN();
        }
        for (iModule = 0; !ifsm.eof(); ++iModule) {
            ifsm.getline(acBuffer, c_iBuffer - 1);
            if (!acBuffer[0])
                iModule--;
        }
        ifsm.close();
        if (!iModule) {
            cerr << "No modules found in: " << szFile << endl;
            return CMeta::GetNaN();
        }
        cerr << "Found " << --iModule << " modules in: " << szFile << endl;
        vecpsModules.resize(iModule);

        dRet = CMeta::GetNaN();
        ifsm.clear();
        ifsm.open(szFile);
        for (iModule = 0; !ifsm.eof(); ++iModule) {
            ifsm.getline(acBuffer, c_iBuffer - 1);
            acBuffer[c_iBuffer - 1] = 0;
            if (CMeta::IsNaN(dRet)) {
                dRet = (float) strtod(acBuffer, nullptr);
                iModule--;
                continue;
            }
            vecstrLine.clear();
            CMeta::Tokenize(acBuffer, vecstrLine);
            if (vecstrLine.size() < 2) {
                iModule--;
                continue;
            }
            vecpsModules[iModule] = psModule = new SModule();
            psModule->m_dSpecificity = (float) strtod(vecstrLine[0].c_str(), nullptr);
            for (i = 1; i < vecstrLine.size(); ++i)
                psModule->m_setiGenes.insert(sSimpleome.Get(vecstrLine[i]));
        }

        return dRet;
    }

    float Jaccard(const SModule &sModule) const {
        set<size_t>::const_iterator iterGene;
        size_t iUnion, iIntersection;

        iIntersection = 0;
        iUnion = sModule.m_setiGenes.size();
        for (iterGene = m_setiGenes.begin(); iterGene != m_setiGenes.end(); ++iterGene)
            if (sModule.m_setiGenes.find(*iterGene) != sModule.m_setiGenes.end())
                iIntersection++;
            else
                iUnion++;

        return ((float) iIntersection / iUnion);
    }

    size_t Intersection(const SModule &sModule) const {
        size_t iRet;
        set<size_t>::const_iterator iterGene;

        for (iRet = 0, iterGene = m_setiGenes.begin(); iterGene != m_setiGenes.end(); ++iterGene)
            if (sModule.m_setiGenes.find(*iterGene) != sModule.m_setiGenes.end())
                iRet++;

        return iRet;
    }

    void Merge(const SModule &sModule) {
        set<size_t>::const_iterator iterGene;

        for (iterGene = sModule.m_setiGenes.begin(); iterGene != sModule.m_setiGenes.end(); ++iterGene)
            m_setiGenes.insert(*iterGene);
        m_dSpecificity = (m_dSpecificity + sModule.m_dSpecificity) / 2;
    }

    bool IsChild(const SModule *psModule) const {

        return (m_setpsChildren.find(psModule) != m_setpsChildren.end());
    }

    void Save(ostream &ostm, float dCutoff, const SSimpleome &sSimpleome) const {
        set<const SModule *>::const_iterator iterChild;
        set<size_t>::const_iterator iterGene;

        ostm << this << '\t' << dCutoff << '\t' << m_dSpecificity << '\t';
        for (iterChild = m_setpsChildren.begin(); iterChild != m_setpsChildren.end(); ++iterChild) {
            if (iterChild != m_setpsChildren.begin())
                ostm << '|';
            ostm << *iterChild;
        }
        for (iterGene = m_setiGenes.begin(); iterGene != m_setiGenes.end(); ++iterGene)
            ostm << '\t' << sSimpleome.Get(*iterGene);
        ostm << endl;
    }
};

struct SSorterModules {
    const vector<float> &m_vecdModules;

    explicit SSorterModules(const vector<float> &vecdModules) : m_vecdModules(vecdModules) {}

    bool operator()(size_t iOne, size_t iTwo) const {

        return (m_vecdModules[iOne] > m_vecdModules[iTwo]);
    }
};

int MainModules(const gengetopt_args_info &sArgs) {
    vector<float> vecdModules;
    vector <vector<SModule *>> vecvecpsModules;
    vector <size_t> veciIndices;
    size_t i, j, iOuter, iCutoffOne, iCutoffTwo, iModuleOne, iModuleTwo;
    SSimpleome sSimpleome;
    bool fDone;
    float d;
    ofstream ofsm;
    ostream *postm;

    vecdModules.resize(input_files.size());
    vecvecpsModules.resize(input_files.size());
    veciIndices.resize(input_files.size());
    for (i = 0; i < vecvecpsModules.size(); ++i) {
        veciIndices[i] = i;
        if (CMeta::IsNaN(vecdModules[i] = SModule::Open(input_files[i].c_str(), sSimpleome,
                                                        vecvecpsModules[i])))
            return 1;
    }
    sort(veciIndices.begin(), veciIndices.end(), SSorterModules(vecdModules));

    for (iOuter = 0, fDone = false; !fDone; ++iOuter) {
        fDone = true;
        cerr << "Outer loop: " << iOuter << endl;
        for (iCutoffOne = 0; iCutoffOne < veciIndices.size(); ++iCutoffOne) {
            vector < SModule * > &vecpsModulesOne = vecvecpsModules[veciIndices[iCutoffOne]];

            cerr << "Merging cutoff: " << vecdModules[veciIndices[iCutoffOne]] << endl;
            for (iModuleOne = 0; iModuleOne < vecpsModulesOne.size(); ++iModuleOne) {
                SModule *psOne = vecpsModulesOne[iModuleOne];

                if (!psOne)
                    continue;
                for (iModuleTwo = (iModuleOne + 1); iModuleTwo < vecpsModulesOne.size(); ++iModuleTwo) {
                    SModule *psTwo = vecpsModulesOne[iModuleTwo];

                    if (!psTwo)
                        continue;
                    if ((d = psOne->Jaccard(*psTwo)) >= sArgs.jaccard_arg) {
                        cerr << "Merging @" << d << ' ' << iCutoffOne << ':' << iModuleOne << " (" <<
                             psOne->m_setiGenes.size() << ") " << iCutoffOne << ':' << iModuleTwo << " (" <<
                             psTwo->m_setiGenes.size() << ')' << endl;
                        psOne->Merge(*psTwo);
                        delete vecpsModulesOne[iModuleTwo];
                        vecpsModulesOne[iModuleTwo] = nullptr;
                        iModuleTwo--;
                        fDone = false;
                    }
                }
            }
            for (iCutoffTwo = (iCutoffOne + 1); iCutoffTwo < veciIndices.size(); ++iCutoffTwo) {
                vector < SModule * > &vecpsModulesTwo = vecvecpsModules[veciIndices[iCutoffTwo]];

                for (iModuleOne = 0; iModuleOne < vecpsModulesOne.size(); ++iModuleOne) {
                    SModule *psOne = vecpsModulesOne[iModuleOne];

                    if (!psOne)
                        continue;
                    for (iModuleTwo = 0; iModuleTwo < vecpsModulesTwo.size(); ++iModuleTwo) {
                        SModule *psTwo = vecpsModulesTwo[iModuleTwo];

                        if (!psTwo)
                            continue;
                        if ((d = psOne->Jaccard(*psTwo)) >= sArgs.jaccard_arg) {
                            cerr << "Merging @" << d << ' ' << iCutoffOne << ':' << iModuleOne << " (" <<
                                 psOne->m_setiGenes.size() << ") " << iCutoffTwo << ':' << iModuleTwo <<
                                 " (" << psTwo->m_setiGenes.size() << ')' << endl;
                            psOne->Merge(*psTwo);
                            delete vecpsModulesTwo[iModuleTwo];
                            vecpsModulesTwo[iModuleTwo] = nullptr;
                            iModuleTwo--;
                            fDone = false;
                        }
                    }
                }
            }
        }
    }

    for (iCutoffOne = 0; iCutoffOne < veciIndices.size(); ++iCutoffOne) {
        vector < SModule * > &vecpsModulesOne = vecvecpsModules[veciIndices[iCutoffOne]];

        for (iCutoffTwo = (iCutoffOne + 1); iCutoffTwo < veciIndices.size(); ++iCutoffTwo) {
            vector < SModule * > &vecpsModulesTwo = vecvecpsModules[veciIndices[iCutoffTwo]];

            for (iModuleOne = 0; iModuleOne < vecpsModulesOne.size(); ++iModuleOne) {
                SModule *psOne = vecpsModulesOne[iModuleOne];

                if (!psOne)
                    continue;
                for (iModuleTwo = 0; iModuleTwo < vecpsModulesTwo.size(); ++iModuleTwo) {
                    SModule *psTwo = vecpsModulesTwo[iModuleTwo];

                    if (!psTwo)
                        continue;
                    if (!psTwo->IsChild(psOne) && ((d = ((float) psOne->Intersection(*psTwo) /
                                                         psOne->m_setiGenes.size())) >= sArgs.intersection_arg)) {
                        cerr << vecdModules[veciIndices[iCutoffOne]] << ':' << iModuleOne <<
                             " child of " << vecdModules[veciIndices[iCutoffTwo]] << ':' << iModuleTwo <<
                             " (" << psOne->m_setiGenes.size() << ", " << psTwo->m_setiGenes.size() << ", " <<
                             d << ')' << endl;
                        psTwo->m_setpsChildren.insert(psOne);
                    }
                }
            }
        }
    }

    if (sArgs.output_arg) {
        ofsm.open(sArgs.output_arg);
        postm = &ofsm;
    } else
        postm = &cout;
    for (iCutoffOne = 0; iCutoffOne < veciIndices.size(); ++iCutoffOne) {
        vector < SModule * > &vecpsModulesOne = vecvecpsModules[veciIndices[iCutoffOne]];

        for (iModuleOne = 0; iModuleOne < vecpsModulesOne.size(); ++iModuleOne) {
            SModule *psOne = vecpsModulesOne[iModuleOne];

            if (psOne)
                psOne->Save(*postm, vecdModules[veciIndices[iCutoffOne]], sSimpleome);
        }
    }
    if (sArgs.output_arg)
        ofsm.close();

    for (i = 0; i < vecvecpsModules.size(); ++i)
        for (j = 0; j < vecvecpsModules[i].size(); ++j)
            if (vecvecpsModules[i][j])
                delete vecvecpsModules[i][j];

    return 0;
}

int MainRevDATs(const gengetopt_args_info &sArgs) {
    vector <string> vecstrTerms, vecstrGenes;
    vector < CGenes * > vecpTerms;
    CGenome Genome;
    CDat DatOut;
    CHalfMatrix <uint32_t> MatCounts;
    size_t i, j, k, iOne, iTwo, iArg, iA, iB;
    float d;

    if (!sArgs.terms_arg) {
        cerr << "Terms argument required" << endl;
        return 1;
    }
    if (!CGenes::Open(sArgs.terms_arg, Genome, vecstrTerms, vecpTerms)) {
        cerr << "Could not open: " << sArgs.terms_arg << endl;
        return 1;
    }

    {
        set <string> setstrGenes;

        for (i = 0; i < vecpTerms.size(); ++i)
            for (j = 0; j < vecpTerms[i]->GetGenes(); ++j)
                setstrGenes.insert(vecpTerms[i]->GetGene(j).GetName());
        vecstrGenes.resize(setstrGenes.size());
        copy(setstrGenes.begin(), setstrGenes.end(), vecstrGenes.begin());
    }

    DatOut.Open(vecstrGenes, false);
    for (i = 0; i < DatOut.GetGenes(); ++i)
        memset(DatOut.Get(i), 0, (DatOut.GetGenes() - i - 1) * sizeof(*DatOut.Get(i)));
    MatCounts.Initialize(DatOut.GetGenes());
    MatCounts.Clear();
    for (iArg = 0; iArg < input_files.size(); ++iArg) {
        CDat DatIn;
        vector <vector<size_t>> vecveciGenes;

        if (!DatIn.Open(input_files[iArg].c_str())) {
            cerr << "Could not open: " << input_files[iArg] << endl;
            return 1;
        }
        vecveciGenes.resize(DatIn.GetGenes());
        for (i = 0; i < DatIn.GetGenes(); ++i) {
            for (j = 0; j < vecstrTerms.size(); ++j)
                if (DatIn.GetGene(i) == vecstrTerms[j])
                    break;
            if (j < vecstrTerms.size()) {
                for (k = 0; k < vecpTerms[j]->GetGenes(); ++k)
                    if ((iOne = DatOut.GetGene(vecpTerms[j]->GetGene(k).GetName())) != -1)
                        vecveciGenes[i].push_back(iOne);
            } else
                cerr << "Unrecognized gene: " << DatIn.GetGene(i) << endl;
        }
        for (i = 0; i < DatIn.GetGenes(); ++i) {
            if (vecveciGenes[i].empty())
                continue;
            for (j = (i + 1); j < DatIn.GetGenes(); ++j) {
                if (CMeta::IsNaN(d = DatIn.Get(i, j)))
                    continue;
                for (iA = 0; iA < vecveciGenes[i].size(); ++iA) {
                    iOne = vecveciGenes[i][iA];
                    for (iB = 0; iB < vecveciGenes[j].size(); ++iB) {
                        iTwo = vecveciGenes[j][iB];
                        MatCounts.Get(iOne, iTwo)++;
                        DatOut.Get(iOne, iTwo) += d;
                    }
                }
            }
        }
    }

    for (i = 0; i < DatOut.GetGenes(); ++i)
        for (j = (i + 1); j < DatOut.GetGenes(); ++j)
            if ((k = MatCounts.Get(i, j)))
                DatOut.Get(i, j) /= k;
            else
                DatOut.Set(i, j, CMeta::GetNaN());
    DatOut.Save(sArgs.output_arg);

    for (i = 0; i < vecpTerms.size(); ++i)
        delete vecpTerms[i];
    return 0;
}
