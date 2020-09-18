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

static const char c_acDab[] = ".dab";
static const char c_acQDab[] = ".qdab";
static const char c_acDat[] = ".dat";
static const char c_acSVM[] = ".svm";
static const char c_acOUT[] = ".out";

enum EMethod {
    EMethodBegin = 0,
    EMethodMean = EMethodBegin,
    EMethodMax = EMethodMean + 1,
    EMethodQuant = EMethodMax + 1,
    EMethodMedian = EMethodQuant + 1,
    EMethodSelectMean = EMethodMedian + 1,
    EMethodEnd = EMethodSelectMean + 1
};

static const char *c_aszMethods[] = {
        "mean", "max", "quant", "median", "selectmean", nullptr
};


float Percentile(vector<float> &vecVals, float quartile) {
    //size_t iOne, iTwo, iSize;
    //float d, dFrac;
    size_t iSize;
    float value, value2;

    iSize = vecVals.size();
    if (iSize == 0)
        return CMeta::GetNaN();

    if (iSize == 1)
        return vecVals[0];

    std::sort(vecVals.begin(), vecVals.end());

    float index = quartile * (iSize + 1);
    float remainder = index - (size_t) index;

    index = (size_t) index - 1;

    if (remainder == 0) {
        return vecVals[(size_t) index];
    } else {
        if (remainder > 0.5) {
            value = vecVals[(size_t) index - 1];
            value2 = vecVals[(size_t) index];
        } else {
            value = vecVals[(size_t) index];
            value2 = vecVals[(size_t) index + 1];
        }
        float interpolationValue = (value2 - value) * remainder;

        cerr << value << " " << value2 << " " << remainder << endl;
        //float interpolationValue = (value - vecVals[((size_t)index-1)]) * remainder;
        //cerr << interpolationValue << endl;

        //return (vecVals[(size_t)index] + interpolationValue);
        return (value + interpolationValue);
    }
}

float Median(vector<float> &vecVals, float quartile) {
    size_t iSize, idx;

    iSize = vecVals.size();
    if (iSize == 0)
        return CMeta::GetNaN();

    if (iSize == 1)
        return vecVals[0];

    std::sort(vecVals.begin(), vecVals.end());

    if (quartile == 0.5) {
        idx = vecVals.size() / 2;
        if (vecVals.size() % 2 != 0)
            return vecVals[idx];
        else
            return ((vecVals[(idx - 1)] + vecVals[idx]) * 0.5);
    } else {
        idx = round(vecVals.size() * quartile);
        return vecVals[idx];
    }

}

float SelectMean(vector<float> &vecVals) {
    size_t iSize, idx, i, j;
    float sum;

    iSize = vecVals.size();
    if (iSize == 0)
        return CMeta::GetNaN();

    if (iSize == 1)
        return vecVals[0];

    std::sort(vecVals.begin(), vecVals.end());

    // return the mean of the top quartile
    idx = (vecVals.size() / 4) * 3;

    j = 0;
    sum = 0.0;
    for (i = idx; i < iSize; ++i) {
        ++j;
        sum += vecVals[i];
    }

    return sum / j;


    // DEBUG
    if (1 == 2) {
        idx = vecVals.size() / 2;

        j = 0;
        sum = 0.0;
        for (i = idx; i < iSize; ++i) {
            ++j;
            sum += vecVals[i];
        }

        return sum / j;
    }

    // DEBUG
    if (1 == 2) {

        idx = (vecVals.size() / 10);
        idx = iSize - idx;

        j = 0;
        sum = 0.0;
        for (i = idx; i < iSize; ++i) {
            ++j;
            sum += vecVals[i];
        }

        return sum / j;
    }

    // DEBUG
    if (1 == 2) {

        idx = (vecVals.size() / 4) * 3;
        return vecVals[idx];
    }

}

int main(int iArgs, char **aszArgs) {
    gengetopt_args_info sArgs{};
    int iRet = 0;
    size_t i, j, k, l;
    size_t iRuns, inputNums;
    float d, dout;
    DIR *dp;
    struct dirent *ep;
    CDat DatOut, DatTrack, DatCur;
    vector <size_t> veciGenesCur;
    // store all input file names
    vector <string> input_files;
    vector <string> vecstrInputs, vecstrDatasets;
    EMethod eMethod;
    ifstream ifsm;
    vector<float> vecWeights;
    map <string, size_t> mapstriDatasets;
    CGenome Genome;
    CGenes Genes(Genome);

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }
    CMeta Meta(sArgs.verbosity_arg);

    for (eMethod = EMethodBegin; eMethod < EMethodEnd; eMethod = (EMethod) (eMethod + 1))
        if (!strcmp(c_aszMethods[eMethod], sArgs.method_arg)) {
            cerr << "combine method: " << c_aszMethods[eMethod] << endl;
            break;
        }

    // now collect the data files from directory if given
    if (sArgs.directory_given) {
        dp = opendir(sArgs.directory_arg);
        i = 0;
        if (dp != nullptr) {
            while ((ep = readdir(dp))) {
                if (strstr(ep->d_name, c_acDab) != nullptr ||
                    strstr(ep->d_name, c_acQDab) != nullptr ||
                    strstr(ep->d_name, c_acDat) != nullptr
                        ) {
                    if (strstr(ep->d_name, c_acSVM) != nullptr)
                        continue;
                    if (strstr(ep->d_name, c_acOUT) != nullptr)
                        continue;
                    vecstrInputs.push_back((string) sArgs.directory_arg + "/" + ep->d_name);
                    vecstrDatasets.push_back(vecstrInputs[i]);
                    mapstriDatasets[ep->d_name] = i;
                    ++i;
                }
            }
            (void) closedir(dp);

            inputNums = vecstrDatasets.size();

            // sort by ASCI
            //std::sort( vecstrInputs.begin(), vecstrInputs.end() );

            cerr << "Number of datasets: " << inputNums << '\n';
        } else {
            cerr << "Couldn't open the directory: " << sArgs.directory_arg << '\n';
            return 1;
        }
    }

    // read dataset weight file
    if (sArgs.weight_given) {
        static const size_t c_iBuffer = 1024;
        char szBuffer[c_iBuffer];
        string strFile;
        vector <string> vecstrLine;
        map<string, size_t>::const_iterator iterDataset;
        size_t numDataset;
        vecWeights.resize(vecstrInputs.size());

        // set default weight
        for (i = 0; i < vecWeights.size(); i++) {
            vecWeights[i] = 0;
        }

        ifsm.clear();
        ifsm.open(sArgs.weight_arg);
        if (!ifsm.is_open()) {
            cerr << "Could not open: " << sArgs.weight_arg << endl;
            return 1;
        }
        while (!ifsm.eof()) {
            ifsm.getline(szBuffer, c_iBuffer - 1);
            szBuffer[c_iBuffer - 1] = 0;
            if (!szBuffer[0])
                continue;
            vecstrLine.clear();
            CMeta::Tokenize(szBuffer, vecstrLine);
            if (vecstrLine.size() != 2) {
                cerr << "Illegal weight line: " << szBuffer << endl;
                return 1;
            }

            if ((iterDataset = mapstriDatasets.find(vecstrLine[0])) == mapstriDatasets.end())
                cerr << "Dataset in weights but not database: " << vecstrLine[0] << endl;
            else
                vecWeights[iterDataset->second] = (float) strtol(vecstrLine[1].c_str(), nullptr, 10);
        }
        ifsm.close();

        // now only keep datasets with non-zero weighting
        vecstrDatasets.clear();
        numDataset = 0;
        for (i = 0; i < vecstrInputs.size(); ++i) {
            if (vecWeights[i] == 0)
                continue;
            vecstrDatasets.push_back(vecstrInputs[i]);
            vecWeights[numDataset] = vecWeights[i];
            numDataset++;
        }
        vecstrDatasets.resize(numDataset);
        vecWeights.resize(numDataset);

        cerr << "Total number of datasets combining with non-zero weights: " << numDataset << endl;
    }

    if (eMethod == EMethodMedian or eMethod == EMethodSelectMean) {
        vector<float> vecVals;
        float val;
        vector < CDat * > vecData;

        vecData.resize(vecstrDatasets.size());
        // now iterate dat/dab networks
        for (i = 0; i < vecstrDatasets.size(); ++i) {
            vecData[i] = new CDat();
            if (!vecData[i]->Open(vecstrDatasets[i].c_str())) {
                cerr << "Couldn't open: " << vecstrDatasets[i] << endl;
                return 1;
            }

            if (sArgs.rank_flag)
                vecData[i]->Rank();
            if (sArgs.zscore_flag)
                vecData[i]->Normalize(CDat::ENormalizeZScore);

            cerr << "open: " << vecstrDatasets[i] << endl;
        }

        // initialized the output dab and pair value vector
        DatOut.Open(vecData[0]->GetGeneNames());
        vecVals.resize(vecstrDatasets.size());

        // debug
        cerr << "num dataset: " << vecstrDatasets.size() << endl;

        for (i = 0; i < DatOut.GetGenes(); ++i) {
            for (j = i + 1; j < DatOut.GetGenes(); ++j) {
                // iterate over each dataset
                vecVals.clear();
                for (k = 0; k < vecstrDatasets.size(); ++k) {
                    if (CMeta::IsNaN(val = vecData[k]->Get(i, j)))
                        continue;

                    if (sArgs.weight_given) {
                        val *= vecWeights[k];
                    }
                    vecVals.push_back(val);
                }

                if (vecVals.empty())
                    continue;

                if (eMethod == EMethodMedian)
                    // find median or quantile
                    DatOut.Set(i, j, Median(vecVals, sArgs.quantile_arg));
                else DatOut.Set(i, j, SelectMean(vecVals));
            }
        }

        DatOut.Save(sArgs.output_arg);
        return 0;
    }

    /// IF combine method is Quantile
    /// Beaware that values are qunatized to allow full read in of the input datsets
    if (eMethod == EMethodQuant) {
        CDatasetCompact Dataset;
        vector<float> vecVals;
        size_t val;

        if (!Dataset.Open(vecstrDatasets)) {
            cerr << "Couldn't open input datasets for quantile combine." << endl;
            return 1;
        }

        // initialized the output dab and pair value vector
        DatOut.Open(Dataset.GetGeneNames());
        vecVals.resize(Dataset.GetExperiments());

        // debug
        cerr << "num dataset: " << Dataset.GetExperiments() << endl;

        for (i = 0; i < DatOut.GetGenes(); ++i) {
            for (j = i + 1; j < DatOut.GetGenes(); ++j) {
                // iterate over each dataset
                vecVals.clear();
                for (k = 0; k < Dataset.GetExperiments(); ++k) {
                    // debug
                    //cerr << "Orig val: " << val << endl;

                    if ((val = Dataset.GetDiscrete(i, j, k)) == -1)
                        continue;

                    // debug
                    //cerr << "Got val: " << val << endl;
                    if (sArgs.weight_given) {
                        val *= vecWeights[k];
                    }

                    vecVals.push_back(val);
                }

                if (vecVals.empty())
                    continue;

                // debug
                //cerr << "quant val: " << CStatistics::Percentile(vecVals.begin(), vecVals.end(), sArgs.quantile_arg) << endl;
                //cerr << "stats: " << vecVals.begin() << " " << vecVals.end() << " " << endl; //(vecVals.end()-vecVals.begin()) << endl;
                //cerr << "stats: " << (vecVals.end()-vecVals.begin()) << endl;
                DatOut.Set(i, j, CStatistics::Percentile(vecVals.begin(), vecVals.end(), sArgs.quantile_arg));

                //DatOut.Set(i, j, Percentile(vecVals, sArgs.quantile_arg));
                //for(size_t t = 0; t < vecVals.size(); t++ )
                //	cerr << vecVals[t] << ' ';
                //cerr << endl;
            }
        }

        DatOut.Save(sArgs.output_arg);
        return 0;
    }


    // now iterate dat/dab networks
    for (i = 0; i < vecstrDatasets.size(); ++i) {

        // open first network, we will just add edge weights to this CDat
        if (i == 0) {
            if (!DatCur.Open(vecstrDatasets[i].c_str())) {
                cerr << "Couldn't open: " << vecstrDatasets[i] << endl;
                return 1;
            }

            if (sArgs.rank_flag)
                DatCur.Rank();
            if (sArgs.zscore_flag)
                DatCur.Normalize(CDat::ENormalizeZScore);

            DatOut.Open(DatCur);
            DatTrack.Open(DatCur);

            // this Dat is used to track various values (count, max)
            for (j = 0; j < DatTrack.GetGenes(); ++j)
                for (k = (j + 1); k < DatTrack.GetGenes(); ++k) {
                    if (eMethod == EMethodMean)
                        DatTrack.Set(j, k, 0.0);
                    if (CMeta::IsNaN(d = DatCur.Get(j, k)))
                        continue;

                    if (sArgs.weight_given) {
                        d *= vecWeights[i];
                        // store weighted value (numerator)
                        DatOut.Set(j, k, d);
                    }

                    switch (eMethod) {
                        case EMethodMax:
                            DatOut.Set(j, k, d);
                            break;
                        case EMethodMean:
                            // keep track of denominator
                            DatTrack.Set(j, k, 1.0);
                        default:
                            break;
                    }
                }

            cerr << "opened: " << vecstrDatasets[i] << endl;
            continue;
        }

        if (!DatCur.Open(vecstrDatasets[i].c_str())) {
            cerr << "Couldn't open: " << vecstrDatasets[i] << endl;
            return 1;
        }
        cerr << "opened: " << vecstrDatasets[i] << endl;

        if (sArgs.rank_flag)
            DatCur.Rank();
        if (sArgs.zscore_flag)
            DatCur.Normalize(CDat::ENormalizeZScore);

        if (sArgs.map_flag) {
            // Get gene index match
            cerr << "inside map flag" << endl;
            veciGenesCur.clear();
            veciGenesCur.resize(DatOut.GetGenes());
            for (l = 0; l < DatOut.GetGenes(); l++) {
                veciGenesCur[l] = DatCur.GetGene(DatOut.GetGene(l));
                if (veciGenesCur[l] == -1) {
                    cerr << "ERROR: missing gene: " << DatOut.GetGene(l) << ", " << vecstrDatasets[i] << endl;
                    return 1;
                }
            }
        }

        // now add edges to Dat
        for (j = 0; j < DatOut.GetGenes(); ++j)
            for (k = (j + 1); k < DatOut.GetGenes(); ++k) {

                if (sArgs.map_flag) {
                    if (CMeta::IsNaN(d = DatCur.Get(veciGenesCur[j], veciGenesCur[k]))) {
                        continue;
                    }
                } else {
                    if (CMeta::IsNaN(d = DatCur.Get(j, k))) {
                        continue;
                    }
                }

                if (sArgs.weight_given)
                    d *= vecWeights[i];

                switch (eMethod) {
                    case EMethodMax:
                        if (CMeta::IsNaN((dout = DatOut.Get(j, k))) || d > dout)
                            DatOut.Set(j, k, d);
                        break;

                    case EMethodMean:
                        DatTrack.Set(j, k, DatTrack.Get(j, k) + 1);

                        if (CMeta::IsNaN((dout = DatOut.Get(j, k))))
                            DatOut.Set(j, k, d);
                        else
                            DatOut.Set(j, k, DatOut.Get(j, k) + d);
                        break;

                    default:
                        break;
                }
            }
    }

    switch (eMethod) {
        case EMethodMean:
            // now convert sum to mean
            for (j = 0; j < DatOut.GetGenes(); ++j)
                for (k = (j + 1); k < DatOut.GetGenes(); ++k) {
                    if (CMeta::IsNaN(d = DatOut.Get(j, k)))
                        continue;
                    DatOut.Set(j, k, d / DatTrack.Get(j, k));
                }

        default:
            break;
    }

    // Filter dat
    if (sArgs.genes_given) {
        ifsm.clear();
        ifsm.open(sArgs.genes_arg);
        if (!Genes.Open(ifsm)) {
            cerr << "Could not open: " << sArgs.genes_arg << endl;
            return 1;
        }
        ifsm.close();

        DatOut.FilterGenes(Genes, CDat::EFilterInclude);
    }

    if (sArgs.genee_given) {
        cerr << "Filter only edges with genes in: " << sArgs.genee_arg << endl;
        if (!DatOut.FilterGenes(sArgs.genee_arg, CDat::EFilterEdge)) {
            cerr << "ERROR can't open file: " << sArgs.genee_arg << endl;
            return 1;
        }
    }

    DatOut.Save(sArgs.output_arg);
    return iRet;
}
