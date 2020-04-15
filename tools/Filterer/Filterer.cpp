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

int main(int iArgs, char **aszArgs) {
    gengetopt_args_info sArgs;
    CDat Dat;
    vector <string> vecstrTokens;
    size_t iArg, i, j;
    float d, dMin, dMax;
    CHalfMatrix<char> MatStatus;
    bool fDefaultExclude;
    char c;

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }
    CMeta Meta(sArgs.verbosity_arg);

    if (!Dat.Open(sArgs.input_arg)) {
        cerr << "Could not open: " << sArgs.input_arg << endl;
        return 1;
    }

    MatStatus.Initialize(Dat.GetGenes());
    MatStatus.Clear();
    fDefaultExclude = false;
    for (iArg = 0; iArg < sArgs.inputs_num; ++iArg) {
        vecstrTokens.clear();
        CMeta::Tokenize(sArgs.inputs[iArg] + 1, vecstrTokens, "=", true);
        d = (float) atof(vecstrTokens[0].c_str());
        if (vecstrTokens.size() == 1) {
            if ((sArgs.inputs[iArg][1]) == '=') {
                dMin = -FLT_MAX;
                dMax = d;
            } else {
                dMin = d;
                dMax = FLT_MAX;
            }
        } else {
            dMin = d;
            dMax = (float) atof(vecstrTokens[1].c_str());
        }
        switch (sArgs.inputs[iArg][0]) {
            case 'i':
                fDefaultExclude = true;
                for (i = 0; i < Dat.GetGenes(); ++i)
                    for (j = (i + 1); j < Dat.GetGenes(); ++j)
                        if (!CMeta::IsNaN(d = Dat.Get(i, j)) && (d <= dMax) && (d >= dMin) &&
                            (MatStatus.Get(i, j) == 0))
                            MatStatus.Set(i, j, 1);
                break;

            case 'x':
                for (i = 0; i < Dat.GetGenes(); ++i)
                    for (j = (i + 1); j < Dat.GetGenes(); ++j)
                        if (!CMeta::IsNaN(d = Dat.Get(i, j)) && (d <= dMax) && (d >= dMin))
                            MatStatus.Set(i, j, -1);
                break;

            default:
                cerr << "Unrecognized command: " << sArgs.inputs[iArg] << endl;
                return 1;
        }
    }

    for (i = 0; i < Dat.GetGenes(); ++i)
        for (j = (i + 1); j < Dat.GetGenes(); ++j) {
            c = MatStatus.Get(i, j);
            if ((c < 0) || (fDefaultExclude && (c < 1)))
                Dat.Set(i, j, CMeta::GetNaN());
        }
    Dat.Save(sArgs.output_arg);

    return 0;
}
