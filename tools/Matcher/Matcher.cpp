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

static const char c_acDAB[] = ".dab";

bool read_values(const CDat &Dat, size_t iMin, size_t iMax, vector<float> &vecdValues) {
    size_t i, j, iCount;
    float d;

    for (iCount = i = 0; i < Dat.GetGenes(); ++i)
        for (j = (i + 1); j < Dat.GetGenes(); ++j) {
            if (CMeta::IsNaN(d = Dat.Get(i, j)))
                continue;
            iCount++;
            if (vecdValues.size() < iMax)
                vecdValues.push_back(d);
            else if (((float) rand() / RAND_MAX) < ((float) iMax / iCount))
                vecdValues[rand() % vecdValues.size()] = d;
        }
    if (!(i = vecdValues.size()))
        return false;
    while (vecdValues.size() < iMin)
        vecdValues.push_back(vecdValues[rand() % i]);
    sort(vecdValues.begin(), vecdValues.end());

    return true;
}

int main(int iArgs, char **aszArgs) {
    gengetopt_args_info sArgs;
    IMeasure *pMeasure;
    CMeasurePearson Pearson;
    CMeasureEuclidean Euclidean;
    CMeasureKendallsTau KendallsTau;
    CMeasureKolmogorovSmirnov KolmSmir;
    CMeasureHypergeometric Hypergeom;
    CMeasureQuickPearson PearQuick;
    CMeasureInnerProduct InnerProd;
    CMeasureBinaryInnerProduct BinInnerProd;
    size_t i, iTwo;
    string strFile;

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }
    CMeta Meta(sArgs.verbosity_arg, sArgs.random_arg);

    CMeasureSigmoid EuclideanSig(&Euclidean, false, 1.0f / sArgs.inputs_num);
    IMeasure *apMeasures[] = {&Pearson, &EuclideanSig, &KendallsTau,
                              &KolmSmir, &Hypergeom, &PearQuick, &InnerProd, &BinInnerProd, NULL};

    pMeasure = NULL;
    if (sArgs.distance_arg)
        for (i = 0; apMeasures[i]; ++i)
            if (!strcmp(apMeasures[i]->GetName(), sArgs.distance_arg)) {
                pMeasure = apMeasures[i];
                break;
            }
    if (!pMeasure) {
        cerr << "Could not recognize: " << sArgs.distance_arg << endl;
        return 1;
    }

    if (sArgs.table_flag) {
        for (i = 0; i < sArgs.inputs_num; ++i)
            cout << '\t' << sArgs.inputs[i];
        cout << endl;
    }
    if (!sArgs.inputs_num)
        return 0;
    FOR_EACH_DIRECTORY_FILE((string) sArgs.input_arg, strFile)
    CDat DatOne;
    vector<float> vecdOne;

    if (!CMeta::IsExtension(strFile, c_acDAB))
        continue;

    strFile = (string) sArgs.input_arg + '/' + strFile;
    if (!DatOne.Open(strFile.c_str(), !!sArgs.memmap_flag)) {
        cerr << "Could not open: " << strFile << endl;
        return 1;
    }
    if (!read_values(DatOne, sArgs.size_min_arg, sArgs.size_max_arg, vecdOne))
        continue;

    for (iTwo = 0; iTwo < sArgs.inputs_num; ++iTwo) {
        CDat DatTwo;
        vector<float> vecdTwo;

        if (!DatTwo.Open(sArgs.inputs[iTwo], !!sArgs.memmap_flag)) {
            cerr << "Could not open: " << sArgs.inputs[iTwo] << endl;
            return 1;
        }
        if (sArgs.table_flag)
            cout << '\t';
        if (!read_values(DatTwo, sArgs.size_min_arg, sArgs.size_max_arg, vecdTwo))
            continue;
        if (!sArgs.table_flag)
            cout << strFile << '\t' << sArgs.inputs[iTwo] << '\t';
        cout << pMeasure->Measure(&vecdOne.front(), vecdOne.size(), &vecdTwo.front(), vecdTwo.size(),
                                  IMeasure::EMapNone);
        if (!sArgs.table_flag)
            cout << endl;
    }
    if (sArgs.table_flag)
        cout << endl;
}

return 0; }
